// Fractals from exact pieces: iterated function systems of affine maps with
// rational coefficients, Cantor sets given by digits in a base, and their
// dimensions.
//
// Three kinds of answer, as everywhere in this package:
//   exact      the points of the n-th stage, membership of a rational number
//              in a Cantor set, the similarity dimension when it is rational
//   certified  the similarity dimension as an interval (the root of Moran's
//              equation), a contraction bound
//   estimated  points from the chaos game and box-counting dimensions, said
//              to be estimates
import { Rational, BaseSystem } from "@ratmath/core";
import { parse } from "@ratmath/analytic";
import { refuse, toQ, ZERO, ONE, Q, boxValue } from "./tree.js";

const lt = (a, b) => a.compareTo(b) < 0;
const isqrtBig = (n) => {
  if (n < 0n) return null;
  if (n < 2n) return n;
  let x = BigInt(Math.floor(Math.sqrt(Number(n))));
  while (x * x > n) x--;
  while ((x + 1n) * (x + 1n) <= n) x++;
  return x;
};
/** A rational ≥ sqrt(q): exact when the root is rational, else rounded up to 2^-40. */
const sqrtUp = (q) => {
  const r = exactSqrt(q);
  if (r) return r;
  const scale = 1n << 40n;
  return new Rational(isqrtBig((q.numerator * scale * scale) / q.denominator) + 1n, scale);
};
/** sqrt(q) when it is rational, else null. */
const exactSqrt = (q) => {
  const a = isqrtBig(q.numerator);
  const b = isqrtBig(q.denominator);
  return a !== null && b !== null && a * a === q.numerator && b * b === q.denominator ? new Rational(a, b) : null;
};
const fmt = (q) => (q.denominator === 1n ? `${q.numerator}` : `${q.numerator}/${q.denominator}`);

/** A small seeded generator (mulberry32), so a chaos game can be repeated. */
export function seeded(seed = 1) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * An affine map of the plane with rational coefficients:
 *   (x, y) → (a x + b y + e,  c x + d y + f)
 */
export class Affine {
  constructor({ a = 0, b = 0, c = 0, d = 0, e = 0, f = 0 } = {}) {
    this.a = toQ(a);
    this.b = toQ(b);
    this.c = toQ(c);
    this.d = toQ(d);
    this.e = toQ(e);
    this.f = toQ(f);
    Object.freeze(this);
  }
  /** The similarity x → r·R(x) + t, R a rotation by `turn` quarter turns (exact) or reflected. */
  static similarity(r, { e = 0, f = 0, rotate = 0, reflect = false } = {}) {
    const R = toQ(r);
    const q = ((rotate % 4) + 4) % 4;
    const [cs, sn] = [[1, 0], [0, 1], [-1, 0], [0, -1]][q];
    const C = R.multiply(Q(cs));
    const S = R.multiply(Q(sn));
    return reflect ? new Affine({ a: C, b: S, c: S, d: C.negate(), e, f }) : new Affine({ a: C, b: S.negate(), c: S, d: C, e, f });
  }
  apply([x, y]) {
    return [this.a.multiply(x).add(this.b.multiply(y)).add(this.e), this.c.multiply(x).add(this.d.multiply(y)).add(this.f)];
  }
  applyFloat([x, y]) {
    return [this.a.toNumber() * x + this.b.toNumber() * y + this.e.toNumber(), this.c.toNumber() * x + this.d.toNumber() * y + this.f.toNumber()];
  }
  get det() {
    return this.a.multiply(this.d).subtract(this.b.multiply(this.c));
  }
  /** The ratio r when the map is a similarity (it scales every distance by r), else null. r may be irrational: then { r2 } with r² exact. */
  similarityRatio() {
    const rot = this.a.equals(this.d) && this.b.equals(this.c.negate());
    const ref = this.a.equals(this.d.negate()) && this.b.equals(this.c);
    if (!rot && !ref) return null;
    const r2 = this.a.multiply(this.a).add(this.c.multiply(this.c));
    return { r2, r: exactSqrt(r2) };
  }
  /**
   * An exact rational upper bound on the factor by which the map can stretch a
   * distance: the largest singular value, sqrt((T + sqrt(T² − 4D²))/2) with
   * T = a² + b² + c² + d² and D the determinant, with each square root rounded up.
   */
  lipschitzBound() {
    const s = this.similarityRatio();
    if (s) return s.r ?? sqrtUp(s.r2);
    const T = [this.a, this.b, this.c, this.d].reduce((t, x) => t.add(x.multiply(x)), ZERO);
    const D = this.det;
    const disc = T.multiply(T).subtract(Q(4).multiply(D).multiply(D));
    return sqrtUp(T.add(sqrtUp(disc)).divide(Q(2)));
  }
  toString() {
    return `(x, y) → (${fmt(this.a)}x + ${fmt(this.b)}y + ${fmt(this.e)}, ${fmt(this.c)}x + ${fmt(this.d)}y + ${fmt(this.f)})`.replace(/\+ -/g, "− ");
  }
}

const affine = (m) => (m instanceof Affine ? m : new Affine(m));

/**
 * An iterated function system: a finite list of contracting affine maps of
 * the plane with rational coefficients. Its attractor is the one compact set
 * A with A = f₁(A) ∪ … ∪ f_N(A).
 */
export class IFS {
  constructor(maps, { name = "IFS", probabilities = null } = {}) {
    if (!Array.isArray(maps) || maps.length < 1) refuse("bad-argument", "an IFS needs at least one map");
    this.maps = Object.freeze(maps.map(affine));
    this.name = name;
    if (probabilities) {
      if (probabilities.length !== maps.length) refuse("bad-argument", "one probability for each map");
      const ps = probabilities.map(toQ);
      if (!ps.reduce((s, p) => s.add(p), ZERO).equals(ONE) || ps.some((p) => lt(p, ZERO))) refuse("bad-argument", "the probabilities must be at least 0 and add up to 1");
      this.probabilities = Object.freeze(ps);
    } else this.probabilities = Object.freeze(maps.map(() => Q(1, maps.length)));
    Object.freeze(this);
  }
  /** Each map's contraction bound (exact rationals); all must be below 1 for an attractor. */
  contraction() {
    const bounds = this.maps.map((m) => m.lipschitzBound());
    return { bounds, contracting: bounds.every((b) => lt(b, ONE)) };
  }
  /** One step: the images of the points under every map (exact). */
  step(points) {
    return this.maps.flatMap((m) => points.map((p) => m.apply(p)));
  }
  /**
   * The n-th stage from a starting set of points (default: the fixed point of the first map): exact
   * rational points. There are N^n of them times the start, so n is limited.
   */
  stage(n, start = null) {
    let pts = (start ?? [this.fixedPoint(0)]).map(([x, y]) => [toQ(x), toQ(y)]);
    if (pts.length * this.maps.length ** n > 2_000_000) refuse("too-large", `stage ${n} has more than 2,000,000 points; use the chaos game`);
    for (let i = 0; i < n; i++) pts = this.step(pts);
    return pts;
  }
  /** The fixed point of map i, exactly: (I − A)⁻¹ t. */
  fixedPoint(i = 0) {
    const m = this.maps[i];
    const a = ONE.subtract(m.a);
    const d = ONE.subtract(m.d);
    const b = m.b.negate();
    const c = m.c.negate();
    const det = a.multiply(d).subtract(b.multiply(c));
    if (det.equals(ZERO)) refuse("degenerate", `map ${i} has no single fixed point`);
    return [d.multiply(m.e).subtract(b.multiply(m.f)).divide(det), a.multiply(m.f).subtract(c.multiply(m.e)).divide(det)];
  }
  /**
   * The chaos game: start at a point, apply a map chosen at random with the
   * given probabilities, and repeat. Returns float points (an estimate of the
   * attractor, with the seed to repeat it).
   */
  chaosGame(n = 10000, { seed = 1, skip = 20 } = {}) {
    const rnd = seeded(seed);
    const cum = [];
    let s = 0;
    for (const p of this.probabilities) cum.push((s += p.toNumber()));
    let pt = this.fixedPoint(0).map((q) => q.toNumber());
    const out = new Float64Array(2 * n);
    for (let i = 0; i < n + skip; i++) {
      const u = rnd();
      let k = cum.findIndex((c) => u < c);
      if (k < 0) k = cum.length - 1;
      pt = this.maps[k].applyFloat(pt);
      if (i >= skip) {
        out[2 * (i - skip)] = pt[0];
        out[2 * (i - skip) + 1] = pt[1];
      }
    }
    return { kind: "estimated", seed, points: out, count: n };
  }
  /**
   * The similarity dimension: the s with Σ rᵢ^s = 1 (Moran's equation), for an
   * IFS of similarities. Exact when all ratios are equal and log N / log(1/r)
   * is rational; otherwise a certified interval. It equals the Hausdorff
   * dimension of the attractor when the pieces do not overlap too much (the
   * open set condition), which is not checked here.
   */
  similarityDimension({ bits = 40 } = {}) {
    const ratios = this.maps.map((m) => m.similarityRatio());
    if (ratios.some((r) => r === null)) refuse("not-similar", "the similarity dimension needs every map to be a similarity (a rotation or reflection times a scale)");
    return moran(ratios.map((r) => r.r2), { bits, squared: true });
  }
}

/**
 * Solve Σ rᵢ^s = 1 for s. `ratios` are rationals in (0, 1) (or their squares
 * with { squared: true }). Returns { exact, value, interval, formula }.
 */
export function moran(ratios, { bits = 40, squared = false } = {}) {
  const r2 = ratios.map((r) => (squared ? toQ(r) : toQ(r).multiply(toQ(r))));
  if (r2.some((q) => !lt(ZERO, q) || !lt(q, ONE))) refuse("bad-argument", "each ratio must be between 0 and 1");
  const N = r2.length;
  const same = r2.every((q) => q.equals(r2[0]));
  if (same) {
    // s = log N / log(1/r) = 2 log N / log(1/r²)
    const inv = ONE.divide(r2[0]);
    const rational = rationalLogRatio(BigInt(N), inv);
    const formula = `log(${N})/log(${fmtSqrt(inv)})`;
    const box = boxValue(parse(`2*log(${N})/log(${fmt(inv)})`), {}, bits + 16);
    if (rational) {
      const s = rational.multiply(Q(2));
      return Object.freeze({ kind: "exact", value: s, interval: [s, s], formula });
    }
    return Object.freeze({ kind: "certified", value: null, interval: [box.lo, box.hi], formula });
  }
  // different ratios: bisection on the decreasing function g(s) = Σ r²ᵢ^(s/2) − 1, each term enclosed
  const g = (s) => {
    let lo = ZERO;
    let hi = ZERO;
    for (const q of r2) {
      const b = boxValue(parse(`exp(${fmt(s)}/2*log(${fmt(q)}))`), {}, bits + 24);
      lo = lo.add(b.lo);
      hi = hi.add(b.hi);
    }
    return [lo.subtract(ONE), hi.subtract(ONE)];
  };
  let a = ZERO;
  let b = Q(1);
  while (lt(ZERO, g(b)[0])) b = b.multiply(Q(2));
  const eps = new Rational(1n, 1n << BigInt(bits));
  while (lt(eps, b.subtract(a))) {
    const m = a.add(b).divide(Q(2));
    const [lo, hi] = g(m);
    if (lt(ZERO, lo)) a = m;
    else if (lt(hi, ZERO)) b = m;
    else {
      // too close to decide at this precision: keep both ends honest
      return Object.freeze({ kind: "certified", value: null, interval: [a, b], formula: "the s with Σ rᵢ^s = 1" });
    }
  }
  return Object.freeze({ kind: "certified", value: null, interval: [a, b], formula: "the s with Σ rᵢ^s = 1" });
}
const fmtSqrt = (q) => {
  const r = exactSqrt(q);
  return r ? fmt(r) : `sqrt(${fmt(q)})`;
};

/** Prime factorisation of a positive bigint (small numbers only). */
function factorBig(n) {
  const out = new Map();
  let m = n;
  for (let p = 2n; p * p <= m; p++) {
    while (m % p === 0n) {
      out.set(p, (out.get(p) ?? 0n) + 1n);
      m /= p;
    }
    if (p > 1_000_000n) return null;
  }
  if (m > 1n) out.set(m, (out.get(m) ?? 0n) + 1n);
  return out;
}
/** log N / log x as a rational when it is one (x a rational > 1), else null. */
function rationalLogRatio(N, x) {
  const fN = factorBig(N);
  const fa = factorBig(x.numerator);
  const fb = factorBig(x.denominator);
  if (!fN || !fa || !fb) return null;
  const fx = new Map(fa);
  for (const [p, e] of fb) fx.set(p, (fx.get(p) ?? 0n) - e);
  const primes = new Set([...fN.keys(), ...fx.keys()]);
  // need N = x^t: exponents proportional, t = eN/ex for every prime
  let t = null;
  for (const p of primes) {
    const eN = fN.get(p) ?? 0n;
    const ex = fx.get(p) ?? 0n;
    if (ex === 0n) {
      if (eN !== 0n) return null;
      continue;
    }
    const r = new Rational(eN, ex);
    if (t === null) t = r;
    else if (!t.equals(r)) return null;
  }
  return t;
}

/**
 * The famous attractors, with exact maps.
 *   sierpinski()   three maps of ratio 1/2: dimension log 3/log 2
 *   carpet()       eight maps of ratio 1/3: log 8/log 3
 *   vicsek()       five maps of ratio 1/3: log 5/log 3
 *   cantorDust()   four maps of ratio 1/3: log 4/log 3
 *   barnsley()     Barnsley's fern: four affine maps with his coefficients, as exact decimals
 *   levyDragon(), heighway()  two maps of ratio 1/√2: dimension 2
 * (The Koch curve needs rotations by 60°, which have √3 in them; it is an
 * L-system in lsystem.js.)
 */
export const attractors = Object.freeze({
  sierpinski: () => new IFS([{ a: "1/2", d: "1/2" }, { a: "1/2", d: "1/2", e: "1/2" }, { a: "1/2", d: "1/2", e: "1/4", f: "1/2" }], { name: "Sierpinski triangle" }),
  carpet: () => {
    const maps = [];
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) if (i !== 1 || j !== 1) maps.push({ a: "1/3", d: "1/3", e: Q(i, 3), f: Q(j, 3) });
    return new IFS(maps, { name: "Sierpinski carpet" });
  },
  vicsek: () => new IFS([[0, 0], [2, 0], [1, 1], [0, 2], [2, 2]].map(([i, j]) => ({ a: "1/3", d: "1/3", e: Q(i, 3), f: Q(j, 3) })), { name: "Vicsek fractal" }),
  cantorDust: () => new IFS([[0, 0], [2, 0], [0, 2], [2, 2]].map(([i, j]) => ({ a: "1/3", d: "1/3", e: Q(i, 3), f: Q(j, 3) })), { name: "Cantor dust" }),
  heighway: () => new IFS([{ a: "1/2", b: "-1/2", c: "1/2", d: "1/2" }, { a: "-1/2", b: "-1/2", c: "1/2", d: "-1/2", e: 1 }], { name: "Heighway dragon" }),
  levyDragon: () => new IFS([{ a: "1/2", b: "-1/2", c: "1/2", d: "1/2" }, { a: "1/2", b: "1/2", c: "-1/2", d: "1/2", e: "1/2", f: "1/2" }], { name: "Lévy C curve" }),
  barnsley: () =>
    new IFS(
      [
        { a: 0, b: 0, c: 0, d: "0.16", e: 0, f: 0 },
        { a: "0.85", b: "0.04", c: "-0.04", d: "0.85", e: 0, f: "1.6" },
        { a: "0.2", b: "-0.26", c: "0.23", d: "0.22", e: 0, f: "1.6" },
        { a: "-0.15", b: "0.28", c: "0.26", d: "0.24", e: 0, f: "0.44" },
      ],
      { name: "Barnsley fern", probabilities: ["0.01", "0.85", "0.07", "0.07"] },
    ),
});

// ---- Cantor sets through digits -------------------------------------------------------

/**
 * The set of numbers in [0, 1] that have a base-b expansion using only the
 * given digits. cantor() is the middle-thirds set: base 3, digits 0 and 2.
 *
 *   contains(q)    exact, for a rational q: is there such an expansion?
 *   expansion(q)   q's base-b expansion through Core (repeating part marked with #)
 *   stage(n)       the b-adic intervals of the n-th stage, exactly
 *   measure(n)     the length of the n-th stage, (k/b)^n
 *   dimension()    log k / log b: exact when rational, certified otherwise
 */
export class DigitCantor {
  constructor({ base = 3, digits = [0, 2] } = {}) {
    if (!Number.isInteger(base) || base < 2) refuse("bad-argument", "the base is a whole number ≥ 2");
    const ds = [...new Set(digits)].sort((a, b) => a - b);
    if (!ds.length || ds.some((d) => !Number.isInteger(d) || d < 0 || d >= base)) refuse("bad-argument", `digits must be whole numbers from 0 to ${base - 1}`);
    this.base = base;
    this.digits = Object.freeze(ds);
    Object.freeze(this);
  }
  /**
   * Is q in the set? x is in it exactly when some allowed digit d has
   * b·x − d in [0, 1] and in the set again. For a rational x those numbers
   * share x's denominator, so there are finitely many of them, and the
   * question is whether the graph of these steps has an infinite path from
   * x: a search that always ends. Returns { inside, expansion } with an
   * expansion that uses only the allowed digits when inside.
   */
  contains(q) {
    const x0 = toQ(q);
    if (lt(x0, ZERO) || lt(ONE, x0)) return Object.freeze({ inside: false, reason: "outside [0, 1]" });
    const b = new Rational(BigInt(this.base));
    const D = x0.denominator;
    const key = (x) => x.numerator.toString();
    // states are x = n/D with 0 ≤ n ≤ D; find a cycle reachable from x0 using allowed steps
    const next = (x) => this.digits.map((d) => ({ d, y: b.multiply(x).subtract(new Rational(BigInt(d))) })).filter(({ y }) => !lt(y, ZERO) && !lt(ONE, y));
    const color = new Map(); // 1 on the stack, 2 done (no infinite path)
    const path = [];
    const dfs = (x) => {
      const k = key(x.multiply(new Rational(D)));
      const c = color.get(k);
      if (c === 1) return { cycleAt: k };
      if (c === 2) return null;
      color.set(k, 1);
      for (const { d, y } of next(x)) {
        path.push({ k, d });
        const r = dfs(y);
        if (r) return r;
        path.pop();
      }
      color.set(k, 2);
      return null;
    };
    if (D > 10n ** 7n) refuse("too-large", "the denominator is too large for the exact search");
    const found = dfs(x0);
    if (!found) return Object.freeze({ inside: false, reason: `no base-${this.base} expansion of ${fmt(x0)} uses only the digits ${this.digits.join(", ")}` });
    const start = path.findIndex((p) => p.k === found.cycleAt);
    const pre = path.slice(0, start).map((p) => digitChar(p.d)).join("");
    const rep = path.slice(start).map((p) => digitChar(p.d)).join("");
    return Object.freeze({ inside: true, expansion: `0.${pre}#${rep}`, reason: `0.${pre}(${rep}) repeating, in base ${this.base}` });
  }
  /** q's base-b expansion as Core writes it ("#" starts the repeating part). */
  expansion(q) {
    return toQ(q).toRepeatingBase(BaseSystem.fromBase(this.base));
  }
  /** The k^n intervals of the n-th stage, as [lo, hi] pairs of rationals. */
  stage(n) {
    const k = this.digits.length;
    if (k ** n > 1_000_000) refuse("too-large", "more than a million intervals");
    let iv = [[ZERO, ONE]];
    const b = new Rational(BigInt(this.base));
    for (let i = 0; i < n; i++) {
      const out = [];
      for (const [lo, hi] of iv) {
        const w = hi.subtract(lo).divide(b);
        for (const d of this.digits) {
          const a = lo.add(w.multiply(new Rational(BigInt(d))));
          out.push([a, a.add(w)]);
        }
      }
      iv = mergeTouching(out);
    }
    return iv;
  }
  /** The total length of the n-th stage: (k/b)^n exactly (touching intervals overlap only at ends). */
  measure(n) {
    return new Rational(BigInt(this.digits.length) ** BigInt(n), BigInt(this.base) ** BigInt(n));
  }
  /** log k / log b. */
  dimension({ bits = 40 } = {}) {
    const k = this.digits.length;
    if (k === 1) return Object.freeze({ kind: "exact", value: ZERO, interval: [ZERO, ZERO], formula: "0" });
    const r = rationalLogRatio(BigInt(k), new Rational(BigInt(this.base)));
    const formula = `log(${k})/log(${this.base})`;
    if (r) return Object.freeze({ kind: "exact", value: r, interval: [r, r], formula });
    const box = boxValue(parse(formula), {}, bits + 16);
    return Object.freeze({ kind: "certified", value: null, interval: [box.lo, box.hi], formula });
  }
  toString() {
    return `numbers in [0, 1] with base-${this.base} digits ${this.digits.join(", ")}`;
  }
}
const digitChar = (d) => "0123456789abcdefghijklmnopqrstuvwxyz"[d] ?? `(${d})`;
function mergeTouching(iv) {
  iv.sort((p, q) => p[0].compareTo(q[0]));
  const out = [];
  for (const [lo, hi] of iv) {
    const last = out[out.length - 1];
    if (last && !lt(last[1], lo)) {
      if (lt(last[1], hi)) last[1] = hi;
    } else out.push([lo, hi]);
  }
  return out;
}
export const cantor = (options) => new DigitCantor(options);

// ---- box counting -----------------------------------------------------------------------

/**
 * The box-counting dimension, estimated: count the boxes of side 2^-k that
 * the points meet, for each k, and fit a line to log(count) against k·log 2.
 * Points are float pairs ({ points } from chaosGame, or [[x, y], …]). An
 * estimate: it depends on the sample and the scales, and says so.
 */
export function boxCounting(points, { scales = [2, 3, 4, 5, 6, 7], bounds = null } = {}) {
  const P = points.points ?? Float64Array.from(points.flat());
  const n = P.length / 2;
  let [x0, y0, x1, y1] = bounds ?? [Infinity, Infinity, -Infinity, -Infinity];
  if (!bounds)
    for (let i = 0; i < n; i++) {
      x0 = Math.min(x0, P[2 * i]);
      x1 = Math.max(x1, P[2 * i]);
      y0 = Math.min(y0, P[2 * i + 1]);
      y1 = Math.max(y1, P[2 * i + 1]);
    }
  const side = Math.max(x1 - x0, y1 - y0) || 1;
  const rows = scales.map((k) => {
    const m = 2 ** k;
    const seen = new Set();
    for (let i = 0; i < n; i++) {
      const a = Math.min(m - 1, Math.floor(((P[2 * i] - x0) / side) * m));
      const b = Math.min(m - 1, Math.floor(((P[2 * i + 1] - y0) / side) * m));
      seen.add(a * m + b);
    }
    return { k, boxes: seen.size };
  });
  // least squares slope of log2(count) on k
  const xs = rows.map((r) => r.k);
  const ys = rows.map((r) => Math.log2(r.boxes));
  const mx = xs.reduce((s, x) => s + x, 0) / xs.length;
  const my = ys.reduce((s, y) => s + y, 0) / ys.length;
  const slope = xs.reduce((s, x, i) => s + (x - mx) * (ys[i] - my), 0) / xs.reduce((s, x) => s + (x - mx) ** 2, 0);
  return Object.freeze({ kind: "estimated", dimension: slope, rows, points: n });
}
