// Three ways to evaluate one formula tree: exactly in rationals (when the
// formula allows it), as an interval that is certain to hold the true value,
// and in doubles as an ordinary program would. The trees are the ones
// analytic's parser produces.
import { Rational, RationalInterval } from "@ratmath/core";
import { parse, format, differentiate, simplify, compileFloat, enclose, constants } from "@ratmath/analytic";

export class ChaosError extends Error {
  constructor(code, message) {
    super(message);
    this.name = "ChaosError";
    this.code = code;
  }
}
export const refuse = (code, message) => {
  throw new ChaosError(code, message);
};

export const ZERO = new Rational(0n);
export const ONE = new Rational(1n);
export const Q = (n, d = 1) => new Rational(BigInt(n), BigInt(d));
/** a Rational from a Rational, a whole number, a bigint or a string such as "7/2" or "0.35" */
export function toQ(x) {
  if (x instanceof Rational) return x;
  if (typeof x === "bigint") return new Rational(x);
  if (typeof x === "number") {
    if (Number.isInteger(x)) return new Rational(BigInt(x));
    refuse("inexact", `${x} is a floating-point number; write it as a string ("${x}") so its exact value is the one you typed`);
  }
  if (typeof x === "string") {
    const s = x.trim();
    const m = /^(-?)(\d*)\.(\d+)$/.exec(s);
    if (m) {
      const n = BigInt((m[2] || "0") + m[3]);
      const q = new Rational(n, 10n ** BigInt(m[3].length));
      return m[1] ? q.negate() : q;
    }
    return new Rational(s);
  }
  refuse("bad-number", `cannot read ${String(x)} as an exact number`);
}

// ---- intervals with exact endpoints, rounded outward ---------------------------

const floorDiv = (a, b) => (a >= 0n ? a / b : -((-a + b - 1n) / b));
const down = (q, bits) => {
  const s = 1n << BigInt(bits);
  return new Rational(floorDiv(q.numerator * s, q.denominator), s);
};
const up = (q, bits) => down(q.negate(), bits).negate();
const lt = (a, b) => a.compareTo(b) < 0;
const minQ = (...xs) => xs.reduce((a, b) => (lt(b, a) ? b : a));
const maxQ = (...xs) => xs.reduce((a, b) => (lt(a, b) ? b : a));

export class Box {
  constructor(lo, hi = lo) {
    this.lo = lo;
    this.hi = hi;
  }
  static of(x) {
    if (x instanceof Box) return x;
    if (x instanceof RationalInterval) return new Box(x.low, x.high);
    return new Box(toQ(x));
  }
  get width() {
    return this.hi.subtract(this.lo);
  }
  get mid() {
    return this.lo.add(this.hi).divide(Q(2));
  }
  get isPoint() {
    return this.lo.equals(this.hi);
  }
  /** the largest |x| in the interval */
  get mag() {
    return maxQ(this.lo.abs(), this.hi.abs());
  }
  /** the smallest |x| in the interval */
  get mig() {
    return this.containsZero() ? ZERO : minQ(this.lo.abs(), this.hi.abs());
  }
  containsZero() {
    return this.lo.numerator <= 0n && this.hi.numerator >= 0n;
  }
  contains(q) {
    return this.lo.compareTo(q) <= 0 && q.compareTo(this.hi) <= 0;
  }
  inside(o) {
    return lt(o.lo, this.lo) && lt(this.hi, o.hi);
  }
  round(bits) {
    return this.isPoint && this.lo.denominator.toString(2).length <= bits + 8 ? this : new Box(down(this.lo, bits), up(this.hi, bits));
  }
  add(o) {
    return new Box(this.lo.add(o.lo), this.hi.add(o.hi));
  }
  sub(o) {
    return new Box(this.lo.subtract(o.hi), this.hi.subtract(o.lo));
  }
  neg() {
    return new Box(this.hi.negate(), this.lo.negate());
  }
  mul(o) {
    const p = [this.lo.multiply(o.lo), this.lo.multiply(o.hi), this.hi.multiply(o.lo), this.hi.multiply(o.hi)];
    return new Box(minQ(...p), maxQ(...p));
  }
  div(o) {
    if (o.containsZero()) refuse("division", "division by an interval that contains zero");
    return this.mul(new Box(ONE.divide(o.hi), ONE.divide(o.lo)));
  }
  abs() {
    return new Box(this.mig, this.mag);
  }
  pow(n) {
    if (n === 0) return new Box(ONE);
    if (n < 0) return new Box(ONE).div(this.pow(-n));
    if (n % 2 === 0) {
      const a = this.abs();
      return new Box(a.lo.pow(n), a.hi.pow(n));
    }
    return new Box(this.lo.pow(n), this.hi.pow(n));
  }
  toRationalInterval() {
    return new RationalInterval(this.lo, this.hi);
  }
}

// ---- the three evaluators -------------------------------------------------------

const wholeNumber = (node) => {
  const neg = node.type === "op" && node.op === "neg";
  const n = neg ? node.args[0] : node;
  if (n.type !== "num" || n.value.denominator !== 1n) return null;
  const k = Number(n.value.numerator);
  return neg ? -k : k;
};

/** exact value in rationals, or null when the formula needs a function that has no rational value */
export function exactValue(ast, env) {
  const go = (n) => {
    switch (n.type) {
      case "num":
        return n.value;
      case "var": {
        const v = env[n.name];
        if (v === undefined) refuse("unbound", `no value for ${n.name}`);
        return v;
      }
      case "const":
        return null;
      case "op": {
        if (n.op === "^") {
          const k = wholeNumber(n.args[1]);
          const b = go(n.args[0]);
          if (k === null || b === null) return null;
          if (k < 0 && b.numerator === 0n) refuse("division", "division by zero");
          return b.pow(k);
        }
        const [a, b] = n.args.map(go);
        if (a === null || (n.op !== "neg" && b === null)) return null;
        switch (n.op) {
          case "neg": return a.negate();
          case "+": return a.add(b);
          case "-": return a.subtract(b);
          case "*": return a.multiply(b);
          case "/":
            if (b.numerator === 0n) refuse("division", "division by zero");
            return a.divide(b);
        }
        return null;
      }
      case "call": {
        const a = n.args.map(go);
        if (a.includes(null)) return null;
        if (n.name === "abs") return a[0].abs();
        if (n.name === "min") return minQ(...a);
        if (n.name === "max") return maxQ(...a);
        return null;
      }
    }
    return null;
  };
  return go(ast);
}

/** an interval certain to contain the value, with endpoints kept to about `bits` binary digits */
export function boxValue(ast, env, bits = 96) {
  const memo = new Map();
  const go = (n) => {
    if (n.type !== "op" && n.type !== "call") return one(n);
    let r = memo.get(n);
    if (r === undefined) memo.set(n, (r = one(n).round(bits)));
    return r;
  };
  const one = (n) => {
    switch (n.type) {
      case "num":
        return new Box(n.value);
      case "var": {
        const v = env[n.name];
        if (v === undefined) refuse("unbound", `no value for ${n.name}`);
        return v;
      }
      case "const":
        return Box.of(constants[n.name](bits + 8));
      case "op": {
        if (n.op === "^") {
          const k = wholeNumber(n.args[1]);
          if (k === null) refuse("unsupported", "powers here need a whole-number exponent");
          return go(n.args[0]).pow(k);
        }
        const [a, b] = n.args.map(go);
        switch (n.op) {
          case "neg": return a.neg();
          case "+": return a.add(b);
          case "-": return a.sub(b);
          case "*": return a.mul(b);
          case "/": return a.div(b);
        }
        break;
      }
      case "call": {
        const a = n.args.map(go);
        if (n.name === "abs") return a[0].abs();
        if (n.name === "min") return new Box(minQ(a[0].lo, a[1].lo), minQ(a[0].hi, a[1].hi));
        if (n.name === "max") return new Box(maxQ(a[0].lo, a[1].lo), maxQ(a[0].hi, a[1].hi));
        const f = enclose[n.name === "ln" ? "log" : n.name];
        if (!f || a.length !== 1) refuse("unsupported", `${n.name} cannot be enclosed here`);
        return Box.of(f(a[0].toRationalInterval(), bits + 8));
      }
    }
    refuse("unsupported", "this formula cannot be enclosed");
  };
  return go(ast);
}

/** log|x| for an interval that stays away from zero */
export function logAbs(box, bits = 96) {
  const a = box.abs();
  if (a.lo.numerator === 0n) refuse("log-zero", "the derivative may vanish on the orbit, so the logarithm has no finite enclosure");
  return Box.of(enclose.log(a.toRationalInterval(), bits + 8));
}

/** replace a variable by a tree */
export function substitute(ast, name, replacement) {
  switch (ast.type) {
    case "var":
      return ast.name === name ? replacement : ast;
    case "op":
      return { type: "op", op: ast.op, args: ast.args.map((a) => substitute(a, name, replacement)) };
    case "call":
      return { type: "call", name: ast.name, args: ast.args.map((a) => substitute(a, name, replacement)) };
    default:
      return ast;
  }
}

// ---- polynomials with rational coefficients, constant term first ----------------

export const P = {
  trim(p) {
    const q = p.slice();
    while (q.length > 1 && q[q.length - 1].numerator === 0n) q.pop();
    return q;
  },
  deg: (p) => (p.length === 1 && p[0].numerator === 0n ? -1 : p.length - 1),
  add(a, b) {
    const r = [];
    for (let i = 0; i < Math.max(a.length, b.length); i++) r.push((a[i] ?? ZERO).add(b[i] ?? ZERO));
    return P.trim(r);
  },
  scale: (a, q) => P.trim(a.map((c) => c.multiply(q))),
  sub: (a, b) => P.add(a, P.scale(b, Q(-1))),
  mul(a, b) {
    const r = new Array(a.length + b.length - 1).fill(ZERO);
    for (let i = 0; i < a.length; i++) {
      if (a[i].numerator === 0n) continue;
      for (let j = 0; j < b.length; j++) r[i + j] = r[i + j].add(a[i].multiply(b[j]));
    }
    return P.trim(r);
  },
  pow(a, k) {
    let r = [ONE];
    for (let i = 0; i < k; i++) r = P.mul(r, a);
    return r;
  },
  divmod(a, b) {
    let r = P.trim(a);
    const db = P.deg(b);
    if (db < 0) refuse("division", "division by the zero polynomial");
    if (P.deg(r) < db) return { quo: [ZERO], rem: r };
    const quo = new Array(P.deg(r) - db + 1).fill(ZERO);
    while (P.deg(r) >= db) {
      const k = P.deg(r) - db;
      const c = r[r.length - 1].divide(b[b.length - 1]);
      quo[k] = c;
      const next = r.slice();
      for (let j = 0; j <= db; j++) next[j + k] = next[j + k].subtract(c.multiply(b[j]));
      next.pop();
      r = P.trim(next.length ? next : [ZERO]);
    }
    return { quo: P.trim(quo), rem: r };
  },
  gcd(a, b) {
    let [x, y] = [P.trim(a), P.trim(b)];
    while (P.deg(y) >= 0) [x, y] = [y, P.divmod(x, y).rem];
    return P.deg(x) < 0 ? x : P.scale(x, ONE.divide(x[x.length - 1]));
  },
  derivative: (a) => (a.length === 1 ? [ZERO] : a.slice(1).map((c, i) => c.multiply(Q(i + 1)))),
  compose(a, b) {
    let r = [ZERO];
    for (let i = a.length - 1; i >= 0; i--) r = P.add(P.mul(r, b), [a[i]]);
    return r;
  },
  evaluate(a, x) {
    let r = ZERO;
    for (let i = a.length - 1; i >= 0; i--) r = r.multiply(x).add(a[i]);
    return r;
  },
  evaluateBox(a, x, bits) {
    let r = new Box(ZERO);
    for (let i = a.length - 1; i >= 0; i--) r = r.mul(x).add(new Box(a[i])).round(bits);
    return r;
  },
  /** integer coefficients with the same roots */
  integers(a) {
    let L = 1n;
    const g = (x, y) => (y ? g(y, x % y) : x);
    for (const c of a) L = (L * c.denominator) / g(L, c.denominator);
    return a.map((c) => c.numerator * (L / c.denominator));
  },
  text(a, name = "x") {
    const parts = [];
    for (let k = a.length - 1; k >= 0; k--) {
      const c = a[k];
      if (c.numerator === 0n) continue;
      const mag = c.abs();
      const x = k === 0 ? "" : k === 1 ? name : `${name}^${k}`;
      const body = k === 0 ? mag.toString() : mag.equals(ONE) ? x : `${mag.denominator === 1n ? mag : `(${mag})`}*${x}`;
      parts.push(parts.length ? `${c.numerator < 0n ? "-" : "+"} ${body}` : `${c.numerator < 0n ? "-" : ""}${body}`);
    }
    return parts.length ? parts.join(" ") : "0";
  },
};

/** the formula as a polynomial in `name` (other names must have rational values in env), or null */
export function polynomialOf(ast, name, env = {}) {
  const go = (n) => {
    switch (n.type) {
      case "num":
        return [n.value];
      case "var":
        if (n.name === name) return [ZERO, ONE];
        return env[n.name] instanceof Rational ? [env[n.name]] : null;
      case "op": {
        if (n.op === "^") {
          const k = wholeNumber(n.args[1]);
          const b = go(n.args[0]);
          if (k === null || k < 0 || b === null || k * Math.max(P.deg(b), 0) > 4096) return null;
          return P.pow(b, k);
        }
        const [a, b] = n.args.map(go);
        if (!a || (n.op !== "neg" && !b)) return null;
        switch (n.op) {
          case "neg": return P.scale(a, Q(-1));
          case "+": return P.add(a, b);
          case "-": return P.sub(a, b);
          case "*": return P.mul(a, b);
          case "/": return P.deg(b) === 0 ? P.scale(a, ONE.divide(b[0])) : null;
        }
        return null;
      }
      default:
        return null;
    }
  };
  return go(ast);
}

export { parse, format, differentiate, simplify, compileFloat };
