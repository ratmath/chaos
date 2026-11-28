// A map of the line, iterated. An orbit comes in three kinds: exact (every
// term a rational number, the true orbit), enclosed (every term an interval
// that is certain to hold the true value) and floating point (what an
// ordinary program computes). Putting them side by side shows the step at
// which the doubles stop describing the orbit they claim to.
import { Rational, RationalInterval } from "@ratmath/core";
import { roots as algebraicRoots, format as algebraicFormat } from "@ratmath/algebraic";
import { zeros } from "@ratmath/analytic";

/** About log2 of a nonnegative rational's size: the numerator's bits less the denominator's. */
const bitSize = (q) => q.numerator.toString(2).length - q.denominator.toString(2).length;
import { ChaosError, refuse, ZERO, ONE, Q, toQ, Box, P, parse, format, differentiate, simplify, compileFloat, exactValue, boxValue, logAbs, substitute, polynomialOf } from "./tree.js";

const bitLength = (q) => q.numerator.toString(2).length + q.denominator.toString(2).length;
const num = (q) => ({ type: "num", value: q });
const divisors = (n) => Array.from({ length: n }, (_, i) => i + 1).filter((d) => n % d === 0);

export class IteratedMap {
  /**
   * @param formula the map as a formula in `variable`, e.g. "r*x*(1 - x)"
   * @param options { variable = "x", params: { r: "4" }, critical } (critical: the turning point, for maps such as the tent whose derivative jumps there)
   */
  constructor(formula, options = {}) {
    this.variable = options.variable ?? "x";
    this.params = Object.fromEntries(Object.entries(options.params ?? {}).map(([k, v]) => [k, toQ(v)]));
    this.text = typeof formula === "string" ? formula : format(formula);
    const raw = typeof formula === "string" ? parse(formula) : formula;
    // parameters are written into the tree, so every evaluator sees the same exact map
    let tree = raw;
    for (const [k, v] of Object.entries(this.params)) tree = substitute(tree, k, num(v));
    this.tree = tree;
    this.family = raw;
    this.derivativeTree = simplify(differentiate(tree, this.variable));
    this._float = compileFloat(tree);
    this._dfloat = compileFloat(this.derivativeTree);
    this.polynomial = polynomialOf(tree, this.variable);
    this.name = options.name ?? null;
    this.critical = options.critical === undefined ? null : toQ(options.critical);
  }

  /** the same family with other parameter values */
  with(params) {
    return new IteratedMap(this.family, { variable: this.variable, params: { ...this.params, ...params }, name: this.name, critical: this.critical ?? undefined });
  }

  toString() {
    const ps = Object.entries(this.params).map(([k, v]) => `${k} = ${v}`).join(", ");
    return `${this.variable} → ${this.text}${ps ? ` (${ps})` : ""}`;
  }

  /** f(x), exactly. Refuses when the formula has no rational value there (use enclose). */
  at(x) {
    const v = exactValue(this.tree, { [this.variable]: toQ(x) });
    if (v === null) refuse("not-exact", `${this.text} has no exact rational value; use enclose(x) or enclosedOrbit`);
    return v;
  }
  /** f(x) in doubles */
  float(x) {
    return this._float({ [this.variable]: x });
  }
  /** an interval certain to contain f(X) for every x in X */
  enclose(x, bits = 96) {
    return boxValue(this.tree, { [this.variable]: Box.of(x) }, bits).toRationalInterval();
  }
  /** f′(x): exact when the formula allows, otherwise an enclosure */
  derivativeAt(x, bits = 96) {
    if (x instanceof RationalInterval) return boxValue(this.derivativeTree, { [this.variable]: Box.of(x) }, bits).toRationalInterval();
    const v = exactValue(this.derivativeTree, { [this.variable]: toQ(x) });
    return v ?? boxValue(this.derivativeTree, { [this.variable]: Box.of(x) }, bits).toRationalInterval();
  }

  /**
   * The true orbit x₀, f(x₀), …, fⁿ(x₀) in exact rationals. The numbers can
   * grow fast (for a quadratic map the digits double every step), so the
   * orbit stops with an error once a term passes `maxBits` binary digits.
   */
  orbit(x0, n, { maxBits = 200000 } = {}) {
    const out = [toQ(x0)];
    for (let i = 0; i < n; i++) {
      const next = this.at(out[i]);
      if (bitLength(next) > maxBits) refuse("too-large", `the exact orbit needs more than ${maxBits} binary digits by step ${i + 1}; use enclosedOrbit for long orbits`);
      out.push(next);
    }
    return out;
  }

  /** the orbit as an ordinary program computes it */
  floatOrbit(x0, n) {
    const out = [typeof x0 === "number" ? x0 : toQ(x0).toNumber()];
    for (let i = 0; i < n; i++) out.push(this.float(out[i]));
    return out;
  }

  /**
   * The true orbit as intervals. Working digits are added until the last
   * term is known to `digits` decimal places (or `maxBits` is reached), so
   * the answer is the true orbit to that accuracy however chaotic the map.
   * { terms: RationalInterval[], bits, width: the widest last term }.
   */
  enclosedOrbit(x0, n, { digits = 12, bits, maxBits = 1 << 17 } = {}) {
    const start = Box.of(x0);
    const want = Q(1).divide(Q(10).pow(digits));
    let work = bits ?? Math.max(64, Math.ceil(digits * 3.33) + 32);
    for (;;) {
      const terms = [start];
      let ok = true;
      try {
        for (let i = 0; i < n && ok; i++) {
          const next = boxValue(this.tree, { [this.variable]: terms[i] }, work);
          terms.push(next);
          // an interval this wide says nothing any more: stop and start again with more digits
          if (next.width.compareTo(ONE) > 0) ok = false;
          else if (bitSize(next.mag) > 4000) refuse("escapes", `the orbit leaves every bounded range (it passes 2^4000 by step ${i + 1})`);
        }
      } catch (e) {
        if (!(e instanceof ChaosError) || e.code !== "division") throw e;
        ok = false;
      }
      const width = ok ? terms[n].width : null;
      if ((ok && width.compareTo(want) <= 0) || bits !== undefined || work >= maxBits) {
        if (!ok) refuse("not-enclosed", `the orbit could not be enclosed with ${work} binary digits: an interval grew past width 1 or met a division by zero`);
        return { terms: terms.map((b) => b.toRationalInterval()), bits: work, width };
      }
      work *= 2;
    }
  }

  /**
   * The float orbit beside the true one: rows { n, float, exact, error,
   * digits } where `exact` is the true value (an enclosure), `error` the
   * distance from the double to it and `digits` how many decimal digits of
   * the double are right. `parted` is the first step with no correct digit.
   */
  compare(x0, n, { digits = 20 } = {}) {
    const floats = this.floatOrbit(x0, n);
    const truth = this.enclosedOrbit(x0, n, { digits }).terms;
    let parted = null;
    const rows = floats.map((f, i) => {
      const mid = truth[i].low.add(truth[i].high).divide(Q(2)).toNumber();
      const error = Math.abs(f - mid);
      const scale = Math.max(Math.abs(mid), 1e-300);
      const good = error === 0 ? 17 : Math.max(0, Math.min(17, -Math.log10(error / scale)));
      if (parted === null && good < 1) parted = i;
      return { n: i, float: f, exact: truth[i], error, digits: good };
    });
    return { rows, parted };
  }

  /** the points of a cobweb diagram for plotting: (x₀, 0), (x₀, x₁), (x₁, x₁), (x₁, x₂), … in doubles */
  cobweb(x0, n) {
    const xs = this.floatOrbit(x0, n);
    const points = [[xs[0], 0]];
    for (let i = 0; i < n; i++) points.push([xs[i], xs[i + 1]], [xs[i + 1], xs[i + 1]]);
    return points;
  }

  /** fⁿ as a polynomial (constant term first); needs a polynomial map */
  iterate(n) {
    if (!this.polynomial) refuse("not-polynomial", "this needs a polynomial map with rational coefficients");
    let p = [ZERO, ONE];
    for (let i = 0; i < n; i++) p = P.compose(this.polynomial, p);
    return p;
  }

  /**
   * Every real point of least period n, exactly. For a polynomial map these
   * are algebraic numbers: the real roots of fⁿ(x) − x that are not roots
   * of f^d(x) − x for a smaller d dividing n. Each is returned with its
   * cycle's multiplier (fⁿ)′(x) as an enclosure and the stability that
   * follows; a multiplier of exactly 0, 1 or −1 is recognised exactly.
   */
  periodicPoints(n = 1, { bits = 96 } = {}) {
    if (!Number.isInteger(n) || n < 1) refuse("bad-period", "the period is a positive whole number");
    if (!this.polynomial) refuse("not-polynomial", "exact periodic points need a polynomial map with rational coefficients; periodicPointsIn(n, [a, b]) finds certified ones for any formula");
    if (P.deg(this.polynomial) ** n > 4096) refuse("too-large", "fⁿ(x) − x has too high a degree to solve exactly");
    const fn = this.iterate(n);
    let dyn = P.sub(fn, [ZERO, ONE]);
    if (P.deg(dyn) < 1) refuse("degenerate", "every point has this period");
    for (const d of divisors(n)) {
      if (d === n) continue;
      const lower = P.sub(this.iterate(d), [ZERO, ONE]);
      for (;;) {
        const g = P.gcd(dyn, lower);
        if (P.deg(g) < 1) break;
        dyn = P.divmod(dyn, g).quo;
      }
    }
    if (P.deg(dyn) < 1) return [];
    const dfn = P.derivative(fn);
    const found = algebraicRoots(P.integers(dyn)).map(({ value }) => value);
    return found.map((value) => {
      const multiplier = this._multiplierAt(value, dfn, bits);
      return { value, text: algebraicFormat(value), period: n, ...multiplier };
    });
  }
  fixedPoints(options) {
    return this.periodicPoints(1, options);
  }

  /**
   * The linear test says nothing when the multiplier is ±1; the next terms of
   * g = fⁿ at the point decide. For g′ = −1 the Schwarzian
   * Sg = g‴/g′ − (3/2)(g″/g′)² = −g‴ − (3/2)g″² decides: negative attracts,
   * positive repels (g∘g − p = −2·Sg·u³/3 + …). For g′ = +1, g″ ≠ 0 attracts on
   * one side only ("semistable"); with g″ = 0, the sign of g‴ decides.
   * → { nonlinear: "attracting" | "repelling" | "semistable" | "undecided", nonlinearReason }
   */
  _beyondLinear(value, dfn, minimal, target, bits) {
    const d2 = P.derivative(dfn);
    const d3 = P.derivative(d2);
    const signAt = (q) => {
      if (P.deg(P.divmod(q, minimal).rem) < 0) return 0; // q vanishes at the root exactly
      for (let b = bits; b <= 4 * bits; b *= 2) {
        const [lo, hi] = value.interval(b);
        const e = P.evaluateBox(q, new Box(lo, hi), b + 32);
        if (e.lo.numerator > 0n) return 1;
        if (e.hi.numerator < 0n) return -1;
      }
      return null;
    };
    if (target.equals(ONE)) {
      const s2 = signAt(d2);
      if (s2 === null) return { nonlinear: "undecided", nonlinearReason: "the sign of (fⁿ)″ at the point was not settled" };
      if (s2 !== 0) return { nonlinear: "semistable", nonlinearReason: `multiplier +1 and (fⁿ)″ ${s2 > 0 ? ">" : "<"} 0: points on one side are attracted, on the other repelled` };
      const s3 = signAt(d3);
      if (s3 === null || s3 === 0) return { nonlinear: "undecided", nonlinearReason: "multiplier +1, (fⁿ)″ = 0 and the sign of (fⁿ)‴ was not settled (or it is 0)" };
      return { nonlinear: s3 < 0 ? "attracting" : "repelling", nonlinearReason: `multiplier +1, (fⁿ)″ = 0 and (fⁿ)‴ ${s3 < 0 ? "<" : ">"} 0` };
    }
    // Sg = −g‴ − (3/2)g″² as a polynomial in the point
    const S = P.sub(P.scale(d3, Q(-1)), P.scale(P.mul(d2, d2), Q(3, 2)));
    const sS = signAt(S);
    if (sS === null || sS === 0) return { nonlinear: "undecided", nonlinearReason: "multiplier −1 and the Schwarzian derivative at the point was not settled (or it is 0)" };
    return { nonlinear: sS < 0 ? "attracting" : "repelling", nonlinearReason: `multiplier −1 and Schwarzian derivative ${sS < 0 ? "negative" : "positive"} at the point, so the second iterate ${sS < 0 ? "attracts" : "repels"} (slowly, at cubic order)` };
  }

  _multiplierAt(value, dfn, bits) {
    const minimal = value.poly.map((c) => new Rational(c));
    const exactly = (target) => P.deg(P.gcd(minimal, P.sub(dfn, [target]))) >= 1 && P.deg(P.divmod(P.sub(dfn, [target]), minimal).rem) < 0;
    const [lo, hi] = value.interval(bits);
    const enc = P.evaluateBox(dfn, new Box(lo, hi), bits + 32);
    const multiplier = enc.toRationalInterval();
    // the multiplier is a polynomial in the root, so it equals a rational exactly when that polynomial minus the rational vanishes at the root
    for (const [target, stability] of [[ZERO, "superattracting"], [ONE, "neutral"], [Q(-1), "neutral"]]) {
      if (enc.contains(target) && exactly(target)) {
        const out = { multiplier: new RationalInterval(target, target), stability, exactMultiplier: target };
        if (stability === "neutral") Object.assign(out, this._beyondLinear(value, dfn, minimal, target, bits));
        return out;
      }
    }
    const mag = enc.abs();
    const stability = mag.hi.compareTo(ONE) < 0 ? "attracting" : mag.lo.compareTo(ONE) > 0 ? "repelling" : "undecided";
    return { multiplier, stability, exactMultiplier: null };
  }

  /**
   * The cycles of least period n: the periodic points grouped into orbits.
   * [{ period, points: [text], values, multiplier, stability, lyapunov }]
   * with lyapunov = log|multiplier|/n enclosed (−Infinity for a
   * superattracting cycle).
   */
  cycles(n = 1, { bits = 96 } = {}) {
    const points = this.periodicPoints(n, { bits });
    const left = points.map((p, i) => i);
    const out = [];
    const image = (i) => {
      for (let b = bits; b <= bits * 64; b *= 2) {
        const [lo, hi] = points[i].value.interval(b);
        const y = P.evaluateBox(this.polynomial, new Box(lo, hi), b + 32);
        const hits = points.map((p, j) => j).filter((j) => {
          const [l, h] = points[j].value.interval(b);
          return l.compareTo(y.hi) <= 0 && y.lo.compareTo(h) <= 0;
        });
        if (hits.length === 1) return hits[0];
      }
      refuse("undecided", "could not tell the points of a cycle apart");
    };
    while (left.length) {
      const members = [left.shift()];
      for (let j = image(members[0]); j !== members[0]; j = image(j)) {
        members.push(j);
        left.splice(left.indexOf(j), 1);
      }
      const first = points[members[0]];
      let lyapunov = null;
      if (first.stability === "superattracting") lyapunov = -Infinity;
      else {
        try {
          const l = logAbs(Box.of(first.multiplier), bits);
          lyapunov = new Box(l.lo.divide(Q(n)), l.hi.divide(Q(n))).toRationalInterval();
        } catch {
          lyapunov = null;
        }
      }
      out.push({ period: n, points: members.map((i) => points[i].text), values: members.map((i) => points[i].value), multiplier: first.multiplier, stability: first.stability, lyapunov });
    }
    return out;
  }

  /**
   * Points of period n inside [a, b] for any formula, as certified
   * enclosures: the zeros of fⁿ(x) − x there, each proven to hold exactly
   * one. Points whose least period is a proper divisor of n are included and
   * marked with it when that can be decided.
   */
  periodicPointsIn(n, [a, b], { digits = 20 } = {}) {
    const X = { type: "var", name: this.variable };
    let comp = X;
    for (let i = 0; i < n; i++) comp = substitute(this.tree, this.variable, comp);
    const eq = { type: "op", op: "-", args: [comp, X] };
    const found = zeros(format(eq), [toQ(a), toQ(b)], { variable: this.variable, tolerance: 10 ** -Math.min(digits, 30) }).exact();
    const bits = Math.ceil(digits * 3.33) + 64;
    return found.map((z) => {
      const start = new Box(z.lo, z.hi);
      let x = start;
      let m = new Box(ONE);
      let least = n;
      for (let i = 0; i < n; i++) {
        m = m.mul(boxValue(this.derivativeTree, { [this.variable]: x }, bits)).round(bits);
        x = boxValue(this.tree, { [this.variable]: x }, bits);
        // the orbit left the enclosure of the starting point: it has not come back yet
        const back = x.lo.compareTo(start.hi) <= 0 && start.lo.compareTo(x.hi) <= 0;
        if (back && i + 1 < least && n % (i + 1) === 0) least = i + 1;
      }
      const mag = m.abs();
      return {
        value: new RationalInterval(z.lo, z.hi),
        period: n,
        leastPeriod: least === n ? n : `${least} or ${n}`,
        multiplier: m.toRationalInterval(),
        stability: mag.hi.compareTo(ONE) < 0 ? "attracting" : mag.lo.compareTo(ONE) > 0 ? "repelling" : "undecided",
      };
    });
  }

  /** the critical points in [a, b] (zeros of f′), as rationals when exact and enclosures otherwise */
  criticalPoints([a, b] = [0, 1]) {
    if (this.critical) return [this.critical];
    const dp = this.polynomial && P.derivative(this.polynomial);
    if (dp && P.deg(dp) === 1) {
      const c = dp[0].negate().divide(dp[1]);
      return c.compareTo(toQ(a)) >= 0 && c.compareTo(toQ(b)) <= 0 ? [c] : [];
    }
    return zeros(format(this.derivativeTree), [toQ(a), toQ(b)], { variable: this.variable }).exact().map((z) => (z.lo.equals(z.hi) ? z.lo : new RationalInterval(z.lo, z.hi)));
  }
}

/** A map of the line from a formula: map("r*x*(1 - x)", { params: { r: 4 } }). */
export const map = (formula, options) => new IteratedMap(formula, options);
/** The logistic map x → r·x·(1 − x). */
export const logistic = (r = 4) => new IteratedMap("r*x*(1 - x)", { params: { r }, name: "logistic" });
/** The tent map x → μ·(1/2 − |x − 1/2|), which is μ·x up to 1/2 and μ·(1 − x) after. */
export const tent = (mu = 2) => new IteratedMap("mu*(1/2 - abs(x - 1/2))", { params: { mu }, name: "tent", critical: "1/2" });
/** The quadratic family x → x² + c. */
export const quadratic = (c = -2) => new IteratedMap("x^2 + c", { params: { c }, name: "quadratic" });

/**
 * The doubling map x → 2x mod 1 on [0, 1), done exactly. A rational p/q is
 * eventually periodic: { orbit, preperiod, period } with the orbit listed up
 * to the first repeat.
 */
export function doubling(x0, { base = 2, limit = 100000 } = {}) {
  const b = toQ(base);
  if (b.denominator !== 1n || b.numerator < 2n) refuse("bad-base", "the multiplier is a whole number, at least 2");
  const step = (x) => {
    const y = x.multiply(b);
    const whole = y.numerator / y.denominator;
    return y.subtract(new Rational(whole));
  };
  let x = toQ(x0);
  if (x.numerator < 0n || x.compareTo(ONE) >= 0) refuse("bad-point", "the doubling map lives on [0, 1)");
  const seen = new Map();
  const orbit = [];
  while (!seen.has(x.toString())) {
    if (orbit.length > limit) refuse("too-long", `no repeat within ${limit} steps`);
    seen.set(x.toString(), orbit.length);
    orbit.push(x);
    x = step(x);
  }
  const preperiod = seen.get(x.toString());
  // what doubles do: each step drops one binary digit, so the orbit reaches 0 and stays
  let f = toQ(x0).toNumber();
  let floatDies = null;
  for (let i = 1; i <= 2000 && floatDies === null; i++) {
    f = (Number(b.numerator) * f) % 1;
    if (f === 0) floatDies = i;
  }
  return { orbit, preperiod, period: orbit.length - preperiod, cycle: orbit.slice(preperiod), floatReachesZeroAt: orbit.slice(preperiod).some((q) => q.numerator === 0n) ? null : floatDies };
}
