// L-systems: strings rewritten in parallel by rules, then drawn by a turtle.
//
// The strings and the counts of each symbol are exact (the counts by a
// matrix power, so they are cheap even when the string is far too long to
// write out). The turtle's path is exact when the angle is a multiple of 30°:
// every coordinate is then a + b·√3 with a and b rational, and the points are
// returned in that form as well as in floats.
import { Rational } from "@ratmath/core";
import { refuse, ZERO, ONE, Q } from "./tree.js";

/** Numbers a + b√3 with rational a, b. */
export class Q3 {
  constructor(a = ZERO, b = ZERO) {
    this.a = a;
    this.b = b;
  }
  add(o) {
    return new Q3(this.a.add(o.a), this.b.add(o.b));
  }
  mul(o) {
    // (a + b√3)(c + d√3) = ac + 3bd + (ad + bc)√3
    return new Q3(this.a.multiply(o.a).add(Q(3).multiply(this.b).multiply(o.b)), this.a.multiply(o.b).add(this.b.multiply(o.a)));
  }
  equals(o) {
    return this.a.equals(o.a) && this.b.equals(o.b);
  }
  toNumber() {
    return this.a.toNumber() + this.b.toNumber() * Math.sqrt(3);
  }
  toString() {
    const f = (q) => (q.denominator === 1n ? `${q.numerator}` : `${q.numerator}/${q.denominator}`);
    if (this.b.equals(ZERO)) return f(this.a);
    const s = this.b.equals(ONE) ? "√3" : this.b.equals(ONE.negate()) ? "-√3" : `${f(this.b)}√3`;
    if (this.a.equals(ZERO)) return s;
    return `${f(this.a)} ${this.b.numerator < 0n ? "-" : "+"} ${s.replace(/^-/, "")}`;
  }
}
const half = Q(1, 2);
/** cos and sin of k·30°, exactly. */
function turn30(k) {
  const m = ((k % 12) + 12) % 12;
  const C = [new Q3(ONE), new Q3(ZERO, half), new Q3(half), new Q3(ZERO), new Q3(half.negate()), new Q3(ZERO, half.negate()), new Q3(ONE.negate()), new Q3(ZERO, half.negate()), new Q3(half.negate()), new Q3(ZERO), new Q3(half), new Q3(ZERO, half)];
  return { cos: C[m], sin: C[(m + 9) % 12] };
}

export class LSystem {
  constructor({ axiom, rules, angle = 90, name = "L-system", draw = "FG", move = "f" }) {
    if (typeof axiom !== "string" || !axiom.length) refuse("bad-argument", "give an axiom string");
    this.axiom = axiom;
    this.rules = Object.freeze({ ...rules });
    for (const [k, v] of Object.entries(this.rules)) if (k.length !== 1 || typeof v !== "string") refuse("bad-argument", "rules map single symbols to strings");
    this.angle = angle;
    this.name = name;
    this.draw = draw;
    this.move = move;
    Object.freeze(this);
  }
  /** The string after n rewritings (refused past 5,000,000 symbols: use counts). */
  generate(n) {
    let s = this.axiom;
    for (let i = 0; i < n; i++) {
      let out = "";
      for (const ch of s) {
        out += this.rules[ch] ?? ch;
        if (out.length > 5_000_000) refuse("too-large", `the string passes 5,000,000 symbols at step ${i + 1}; counts(n) gives its make-up exactly`);
      }
      s = out;
    }
    return s;
  }
  /**
   * How many of each symbol the n-th string has, exactly, without writing it
   * out: the counts go forward by the matrix M[a][b] = (number of b in the
   * rule for a), raised to the n-th power by repeated squaring.
   */
  counts(n) {
    const syms = [...new Set([...this.axiom, ...Object.keys(this.rules), ...Object.values(this.rules).join("")])];
    const idx = new Map(syms.map((s, i) => [s, i]));
    const k = syms.length;
    const M = syms.map((a) => {
      const row = new Array(k).fill(0n);
      for (const ch of this.rules[a] ?? a) row[idx.get(ch)] += 1n;
      return row;
    });
    const mul = (A, B) => A.map((r) => B[0].map((_, j) => r.reduce((s, x, t) => s + x * B[t][j], 0n)));
    let P = syms.map((_, i) => syms.map((__, j) => (i === j ? 1n : 0n)));
    let base = M;
    let e = BigInt(n);
    while (e > 0n) {
      if (e & 1n) P = mul(P, base);
      base = mul(base, base);
      e >>= 1n;
    }
    const v = new Array(k).fill(0n);
    for (const ch of this.axiom) v[idx.get(ch)] += 1n;
    const out = new Array(k).fill(0n);
    for (let i = 0; i < k; i++) for (let j = 0; j < k; j++) out[j] += v[i] * P[i][j];
    const counts = Object.fromEntries(syms.map((s, i) => [s, out[i]]));
    return Object.freeze({ counts, length: out.reduce((s, x) => s + x, 0n) });
  }
  /**
   * The turtle's path for the n-th string: "F" and "G" draw a unit step, "f"
   * moves without drawing, "+" and "-" turn by the angle, "[" and "]" save and
   * restore the state. Returns { segments, exact } where segments are float
   * [x1, y1, x2, y2] and, when the angle is a multiple of 30°, `points` holds the
   * exact end points as Q3 pairs (a + b√3).
   */
  turtle(n, { step = 1 } = {}) {
    const s = this.generate(n);
    const exact = Number.isInteger(this.angle) && this.angle % 30 === 0;
    const k30 = exact ? this.angle / 30 : 0;
    const segs = [];
    const exactPts = exact ? [] : null;
    let x = 0;
    let y = 0;
    let h = 0; // heading in units of the angle (exact) or degrees
    let X = new Q3();
    let Y = new Q3();
    const st = [];
    const stepQ = new Q3(new Rational(BigInt(step)));
    for (const ch of s) {
      if (this.draw.includes(ch) || this.move.includes(ch)) {
        let nx;
        let ny;
        if (exact) {
          const { cos, sin } = turn30(h * k30);
          const NX = X.add(cos.mul(stepQ));
          const NY = Y.add(sin.mul(stepQ));
          nx = NX.toNumber();
          ny = NY.toNumber();
          if (this.draw.includes(ch)) exactPts.push([X, Y, NX, NY]);
          X = NX;
          Y = NY;
        } else {
          const r = (h * this.angle * Math.PI) / 180;
          nx = x + step * Math.cos(r);
          ny = y + step * Math.sin(r);
        }
        if (this.draw.includes(ch)) segs.push([x, y, nx, ny]);
        x = nx;
        y = ny;
      } else if (ch === "+") h++;
      else if (ch === "-") h--;
      else if (ch === "|") h += Math.round(180 / this.angle);
      else if (ch === "[") st.push([x, y, h, X, Y]);
      else if (ch === "]") {
        if (!st.length) refuse("bad-argument", "a ] without a [");
        [x, y, h, X, Y] = st.pop();
      }
    }
    return Object.freeze({ kind: exact ? "exact" : "float", segments: segs, exactSegments: exactPts, end: exact ? [X, Y] : [x, y] });
  }
}
export const lsystem = (spec) => new LSystem(spec);

/** Classic L-systems. */
export const lsystems = Object.freeze({
  koch: () => new LSystem({ name: "Koch curve", axiom: "F", rules: { F: "F+F--F+F" }, angle: 60 }),
  snowflake: () => new LSystem({ name: "Koch snowflake", axiom: "F--F--F", rules: { F: "F+F--F+F" }, angle: 60 }),
  arrowhead: () => new LSystem({ name: "Sierpinski arrowhead", axiom: "A", rules: { A: "B-A-B", B: "A+B+A" }, angle: 60, draw: "AB" }),
  dragon: () => new LSystem({ name: "Heighway dragon", axiom: "FX", rules: { X: "X+YF+", Y: "-FX-Y" }, angle: 90 }),
  hilbert: () => new LSystem({ name: "Hilbert curve", axiom: "A", rules: { A: "+BF-AFA-FB+", B: "-AF+BFB+FA-" }, angle: 90 }),
  levy: () => new LSystem({ name: "Lévy C curve", axiom: "F", rules: { F: "+F--F+" }, angle: 45 }),
  plant: () => new LSystem({ name: "fractal plant", axiom: "X", rules: { X: "F+[[X]-X]-F[-FX]+X", F: "FF" }, angle: 25 }),
});
