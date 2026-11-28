// Lyapunov exponents, sensitive dependence, bifurcations and symbolic
// dynamics for maps of the line.
//
// Three kinds of answer, and each says which it is. Exact: a rational or an
// algebraic number. Certified: an interval proven to contain the true value
// (a finite-time Lyapunov exponent along the true orbit, a parameter at
// which a cycle is born). Estimated: a bifurcation diagram, or the limit a
// sequence of certified numbers appears to approach.
import { Rational, RationalInterval } from "@ratmath/core";
import { zeros } from "@ratmath/analytic";
import { IteratedMap } from "./map.js";
import { refuse, ZERO, ONE, Q, toQ, Box, format, differentiate, simplify, compileFloat, exactValue, boxValue, logAbs, substitute } from "./tree.js";

const asMap = (m) => (m instanceof IteratedMap ? m : refuse("bad-map", "give a map made by map(), logistic(), tent() or quadratic()"));
const out = (b) => b.toRationalInterval();

/**
 * The Lyapunov exponent of the orbit of x₀ over n steps,
 * (1/n)·Σ log|f′(xᵢ)|, as a certified enclosure: the average is taken
 * along the true orbit, not a floating-point one. `transient` steps are
 * discarded first. The exponent of the map (the limit as n grows) is what
 * this estimates; for a periodic cycle use cycles(), which gives it exactly.
 * { value, steps, transient, bits, kind, float }.
 */
export function lyapunov(m, x0, { steps = 200, transient = 0, digits = 10 } = {}) {
  asMap(m);
  const total = steps + transient;
  const { terms, bits } = m.enclosedOrbit(x0, total, { digits: digits + 4 });
  // Σ log|f′(xᵢ)| is the log of a product, so the derivatives are multiplied in short runs and one logarithm is taken per run
  const logBits = Math.ceil(digits * 3.33) + 48;
  let sum = new Box(ZERO);
  let run = new Box(ONE);
  for (let i = transient; i < total; i++) {
    run = run.mul(boxValue(m.derivativeTree, { [m.variable]: Box.of(terms[i]) }, bits).abs());
    if ((i - transient) % 8 === 7 || i === total - 1) {
      sum = sum.add(logAbs(run, logBits)).round(logBits);
      run = new Box(ONE);
    }
  }
  const value = out(new Box(sum.lo.divide(Q(steps)), sum.hi.divide(Q(steps))));
  // the same average the way a program would compute it, for comparison
  let x = toQ(x0 instanceof RationalInterval ? x0.low : x0).toNumber();
  let fsum = 0;
  for (let i = 0; i < total; i++) {
    if (i >= transient) fsum += Math.log(Math.abs(m._dfloat({ [m.variable]: x })));
    x = m.float(x);
  }
  return { value, steps, transient, bits, kind: "certified finite-time average", float: fsum / steps };
}

/**
 * Sensitive dependence, measured on true orbits: two starting points a
 * distance `delta` apart, followed together. rows { n, separation } with
 * the separation enclosed; `doublings` is the average factor by which the
 * separation grows per step while it is still small, and `reaches` the
 * first step at which it passes `threshold`.
 */
export function sensitivity(m, x0, delta, { steps = 60, threshold = "1/10", digits = 12 } = {}) {
  asMap(m);
  const a = m.enclosedOrbit(x0, steps, { digits: digits + 6 }).terms;
  const b = m.enclosedOrbit(toQ(x0).add(toQ(delta)), steps, { digits: digits + 6 }).terms;
  const limit = toQ(threshold);
  let reaches = null;
  const rows = a.map((t, n) => {
    const sep = Box.of(b[n]).sub(Box.of(t)).abs();
    if (reaches === null && sep.lo.compareTo(limit) > 0) reaches = n;
    return { n, separation: out(sep) };
  });
  const upTo = reaches ?? steps;
  const d0 = toQ(delta).abs().toNumber();
  const dn = rows[upTo].separation.low.add(rows[upTo].separation.high).divide(Q(2)).toNumber();
  const rate = upTo > 0 && dn > 0 ? Math.log(dn / d0) / upTo : null;
  return { rows, reaches, rate, factor: rate === null ? null : Math.exp(rate) };
}

// ---- bifurcations ----------------------------------------------------------------

/**
 * The data of a bifurcation diagram: for each parameter value, the points
 * the orbit visits after a transient. This is a picture made from
 * floating-point orbits and is labelled as an estimate; the proven
 * statements are bifurcationPoint, superstable and window.
 */
export function bifurcationDiagram(family, { parameter, range, count = 200, start = 0.5, transient = 400, keep = 64 }) {
  asMap(family);
  const f = compileFloat(family.family);
  const others = Object.fromEntries(Object.entries(family.params).filter(([k]) => k !== parameter).map(([k, v]) => [k, v.toNumber()]));
  const [a, b] = range.map((x) => (typeof x === "number" ? x : toQ(x).toNumber()));
  const columns = [];
  for (let i = 0; i < count; i++) {
    const value = count === 1 ? a : a + ((b - a) * i) / (count - 1);
    const env = { ...others, [parameter]: value, [family.variable]: start };
    for (let k = 0; k < transient; k++) env[family.variable] = f(env);
    const seen = new Set();
    for (let k = 0; k < keep; k++) {
      env[family.variable] = f(env);
      seen.add(Math.round(env[family.variable] * 1e6) / 1e6);
    }
    columns.push({ parameter: value, points: [...seen].sort((p, q) => p - q) });
  }
  return { kind: "estimated", parameter, columns };
}

/** partial derivatives of the one-step map that the chain rule for (fⁿ)′, and its derivatives, needs */
function jets(family, parameter) {
  const x = family.variable;
  let tree = family.family;
  for (const [k, v] of Object.entries(family.params)) if (k !== parameter) tree = substitute(tree, k, { type: "num", value: v });
  const d = (t, name) => simplify(differentiate(t, name));
  const fx = d(tree, x);
  const trees = { f: tree, fx, fr: d(tree, parameter), fxx: d(fx, x), fxr: d(fx, parameter) };
  const floats = Object.fromEntries(Object.entries(trees).map(([k, t]) => [k, compileFloat(t)]));
  return { trees, floats, x };
}
/** g = fⁿ and its partials g_x, g_r, g_xx, g_xr at (x, r), with any number type */
function compose(n, evalOne, x0, ops) {
  let g = x0;
  let gx = ops.one;
  let gr = ops.zero;
  let gxx = ops.zero;
  let gxr = ops.zero;
  for (let i = 0; i < n; i++) {
    const v = evalOne(g);
    // order matters: the second derivatives use the old first derivatives
    const nxx = ops.add(ops.mul(v.fxx, ops.mul(gx, gx)), ops.mul(v.fx, gxx));
    const nxr = ops.add(ops.add(ops.mul(v.fxx, ops.mul(gx, gr)), ops.mul(v.fxr, gx)), ops.mul(v.fx, gxr));
    const nr = ops.add(ops.mul(v.fx, gr), v.fr);
    const nx = ops.mul(v.fx, gx);
    g = v.f;
    gx = nx;
    gr = nr;
    gxx = nxx;
    gxr = nxr;
  }
  return { g, gx, gr, gxx, gxr };
}
const FLOAT_OPS = { one: 1, zero: 0, add: (a, b) => a + b, mul: (a, b) => a * b };
const boxOps = (bits) => ({ one: new Box(ONE), zero: new Box(ZERO), add: (a, b) => a.add(b).round(bits), mul: (a, b) => a.mul(b).round(bits) });
const fromDouble = (x) => {
  // the exact rational value of a double
  if (!Number.isFinite(x)) refuse("no-convergence", "the search left the range of doubles");
  let d = 1n;
  let y = x;
  while (!Number.isInteger(y)) {
    y *= 2;
    d *= 2n;
  }
  return new Rational(BigInt(y), d);
};

/**
 * A parameter value at which a cycle of period n has multiplier exactly
 * `multiplier` (−1: the cycle is about to double its period; +1: a cycle is
 * born in a tangent bifurcation), proven to exist and to be the only one in
 * the returned box. Start the search near (parameter, x) = `near`.
 * { parameter, point, period, multiplier, kind } with enclosures.
 */
export function bifurcationPoint(family, { parameter, period = 1, multiplier = -1, near, digits = 20 }) {
  asMap(family);
  const { trees, floats, x } = jets(family, parameter);
  const s = Number(multiplier);
  if (s !== 1 && s !== -1) refuse("bad-multiplier", "the multiplier at a bifurcation is 1 or -1");
  const [r0, x0] = near.map((v) => (typeof v === "number" ? v : toQ(v).toNumber()));
  const fl = (r) => (g) => {
    const env = { [x]: g, [parameter]: r };
    return { f: floats.f(env), fx: floats.fx(env), fr: floats.fr(env), fxx: floats.fxx(env), fxr: floats.fxr(env) };
  };
  // Newton's method in doubles to find a good starting box
  let X = x0;
  let R = r0;
  let J = null;
  for (let it = 0; it < 60; it++) {
    const c = compose(period, fl(R), X, FLOAT_OPS);
    const F1 = c.g - X;
    const F2 = c.gx - s;
    J = [[c.gx - 1, c.gr], [c.gxx, c.gxr]];
    const det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    if (!Number.isFinite(det) || det === 0) refuse("no-convergence", "the search met a singular step; start nearer the bifurcation");
    const dx = (F1 * J[1][1] - F2 * J[0][1]) / det;
    const dr = (J[0][0] * F2 - J[1][0] * F1) / det;
    X -= dx;
    R -= dr;
    if (Math.abs(dx) + Math.abs(dr) < 1e-15 * (1 + Math.abs(X) + Math.abs(R))) break;
    if (it === 59) refuse("no-convergence", "Newton's method did not settle; start nearer the bifurcation");
  }
  // Krawczyk's test in exact interval arithmetic: K(B) inside B proves one solution in B
  const bits = Math.ceil(digits * 3.33) + 48 + 8 * period;
  const ops = boxOps(bits);
  const bx = (r) => (g) => {
    const env = { [x]: g, [parameter]: r };
    return Object.fromEntries(Object.entries(trees).map(([k, t]) => [k, boxValue(t, env, bits)]));
  };
  const det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
  const C = [[J[1][1] / det, -J[0][1] / det], [-J[1][0] / det, J[0][0] / det]].map((row) => row.map(fromDouble));
  let mx = fromDouble(X);
  let mr = fromDouble(R);
  let rad = Q(1, 1n << 20n);
  let proven = null;
  for (let round = 0; round < 60; round++) {
    const pm = compose(period, bx(new Box(mr)), new Box(mx), ops);
    const Fm = [pm.g.sub(new Box(mx)), pm.gx.sub(new Box(Q(s)))];
    const BX = new Box(mx.subtract(rad), mx.add(rad));
    const BR = new Box(mr.subtract(rad), mr.add(rad));
    const pb = compose(period, bx(BR), BX, ops);
    const JB = [[pb.gx.sub(new Box(ONE)), pb.gr], [pb.gxx, pb.gxr]];
    const cb = (q) => new Box(q);
    // the preconditioner: the inverse of the midpoint Jacobian, recomputed as the box shrinks
    const jm = JB.map((row) => row.map((e) => e.mid));
    const dm = jm[0][0].multiply(jm[1][1]).subtract(jm[0][1].multiply(jm[1][0]));
    if (proven && dm.numerator !== 0n) {
      const short = (q) => new Box(q).round(96).lo;
      C[0][0] = short(jm[1][1].divide(dm));
      C[0][1] = short(jm[0][1].divide(dm).negate());
      C[1][0] = short(jm[1][0].divide(dm).negate());
      C[1][1] = short(jm[0][0].divide(dm));
    }
    // I − C·J(B)
    const M = [0, 1].map((i) => [0, 1].map((j) => new Box(i === j ? ONE : ZERO).sub(cb(C[i][0]).mul(JB[0][j]).add(cb(C[i][1]).mul(JB[1][j])))));
    const sym = new Box(rad.negate(), rad);
    const K = [0, 1].map((i) => {
      const centre = new Box(i === 0 ? mx : mr).sub(cb(C[i][0]).mul(Fm[0]).add(cb(C[i][1]).mul(Fm[1])));
      return centre.add(M[i][0].mul(sym)).add(M[i][1].mul(sym)).round(bits);
    });
    if (K[0].inside(BX) && K[1].inside(BR)) {
      proven = K;
      const want = Q(1).divide(Q(10).pow(digits));
      if (K[1].width.compareTo(want) <= 0 || round > 40) break;
      // contract: recentre on K and shrink the box to just cover it
      mx = K[0].mid;
      mr = K[1].mid;
      const w = K[0].width.compareTo(K[1].width) > 0 ? K[0].width : K[1].width;
      // K holds the solution, so a box just covering K does too; the test contracts faster on a smaller box
      rad = w.multiply(Q(11, 20));
      continue;
    }
    if (proven) break;
    // not yet proven: take the exact Newton step when it is a small one (the doubles' solution is only
    // approximate), and try a smaller box, where the derivatives vary less
    const small = Q(1, 1024);
    if (K[0].mid.subtract(mx).abs().compareTo(small) < 0 && K[1].mid.subtract(mr).abs().compareTo(small) < 0) {
      mx = K[0].mid;
      mr = K[1].mid;
    }
    rad = rad.divide(Q(16));
    if (rad.compareTo(Q(1, 1n << 200n)) < 0) break;
  }
  if (!proven) refuse("not-proven", "a solution was found numerically but could not be proven; the bifurcation may be degenerate here");
  return {
    parameter: out(proven[1]),
    point: out(proven[0]),
    period,
    multiplier: s,
    kind: s === -1 ? "period doubling" : "tangent (a cycle is born)",
  };
}

/**
 * Parameters in `range` at which the critical point lies on a cycle of
 * period n (superstable cycles), each a certified enclosure. Those whose
 * critical orbit already closes after a proper divisor of n are left out
 * when that can be decided.
 */
export function superstable(family, { parameter, period, range, critical, digits = 20 }) {
  asMap(family);
  let tree = family.family;
  for (const [k, v] of Object.entries(family.params)) if (k !== parameter) tree = substitute(tree, k, { type: "num", value: v });
  const c = critical === undefined ? family.with({ [parameter]: range[0] }).criticalPoints([-1e6, 1e6])[0] : toQ(critical);
  if (!(c instanceof Rational)) refuse("no-critical-point", "give the critical point as { critical }");
  const start = { type: "num", value: c };
  const closes = (n) => {
    let comp = start;
    for (let i = 0; i < n; i++) comp = substitute(tree, family.variable, comp);
    return { type: "op", op: "-", args: [comp, start] };
  };
  const found = zeros(format(closes(period)), range.map(toQ), { variable: parameter, tolerance: 10 ** -Math.min(digits, 30) }).exact();
  const bits = Math.ceil(digits * 3.33) + 64;
  const keep = [];
  for (const z of found) {
    const r = new Box(z.lo, z.hi);
    let lower = false;
    for (let d = 1; d < period; d++) {
      if (period % d) continue;
      const v = boxValue(closes(d), { [parameter]: r }, bits);
      if (v.containsZero()) lower = true;
    }
    if (!lower) keep.push(new RationalInterval(z.lo, z.hi));
  }
  return keep;
}

/**
 * The period-doubling cascade of a family: the parameters at which the
 * cycle of period 1, 2, 4, … doubles, each proven by bifurcationPoint, and
 * the ratios of successive gaps, which approach Feigenbaum's constant
 * 4.669201609…. The parameters and the ratios are certified; that the
 * ratios converge to that constant is not proven here.
 */
export function periodDoubling(family, { parameter, levels = 4, start, digits = 16 }) {
  asMap(family);
  const points = [];
  const { floats, x } = jets(family, parameter);
  let here = null;
  let gap = null;
  for (let k = 0; k < levels; k++) {
    const period = 2 ** k;
    const guess = k === 0 ? start : settle(floats, x, parameter, here, period, gap);
    let found;
    try {
      found = bifurcationPoint(family, { parameter, period, multiplier: -1, near: guess, digits });
    } catch (e) {
      // the composed map's exact coefficients grow quickly with the period; some engines cap a BigInt at 2^20 bits
      if (e instanceof RangeError) refuse("budget", `the doubling of the period-${period} cycle needs numbers larger than this JavaScript engine can hold; levels up to ${k} were found (ask for levels: ${k})`);
      throw e;
    }
    points.push(found);
    const r = found.parameter.low.toNumber();
    gap = here === null ? null : r - here;
    here = r;
  }
  const ratios = [];
  for (let k = 2; k < points.length; k++) {
    const a = Box.of(points[k - 1].parameter).sub(Box.of(points[k - 2].parameter));
    const b = Box.of(points[k].parameter).sub(Box.of(points[k - 1].parameter));
    ratios.push(out(a.div(b)));
  }
  return { points, parameters: points.map((p) => p.parameter), ratios, feigenbaum: "4.669201609102990…", kind: "parameters and ratios certified; the limit of the ratios is not proven here" };
}
/** a float guess (parameter, x) for the doubling of the cycle of the given period beyond `here` */
function settle(floats, x, parameter, here, period, gap) {
  // walk the parameter forward, following the attracting cycle, until its multiplier nears −1
  const step = gap === null ? 2e-3 * Math.max(1, Math.abs(here)) : Math.abs(gap) / 4.669 / 100;
  const dir = gap === null ? 1 : Math.sign(gap);
  let xx = 0.5;
  let last = null;
  for (let i = 1; i <= 3000; i++) {
    const r = here + dir * step * i;
    for (let k = 0; k < 4000; k++) xx = floats.f({ [x]: xx, [parameter]: r });
    let m = 1;
    let y = xx;
    for (let k = 0; k < period; k++) {
      m *= floats.fx({ [x]: y, [parameter]: r });
      y = floats.f({ [x]: y, [parameter]: r });
    }
    if (!Number.isFinite(m)) break;
    const closed = Math.abs(y - xx) < 1e-7;
    if (closed && Math.abs(m) < 1) {
      last = [r, xx];
      if (m < -0.9) return last;
    } else if (last) return last;
  }
  return last ?? refuse("no-convergence", "could not follow the cycle to its next doubling; give fewer levels");
}

/**
 * A periodic window: the stretch of parameters in which an attracting
 * cycle of period n exists, from its birth in a tangent bifurcation
 * (multiplier +1) to its first period doubling (multiplier −1), with the
 * superstable parameter between them. All three are certified.
 * `near` is a rough (parameter, x) inside the window.
 */
export function periodicWindow(family, { parameter, period, near, digits = 16 }) {
  asMap(family);
  const born = bifurcationPoint(family, { parameter, period, multiplier: 1, near, digits });
  const doubles = bifurcationPoint(family, { parameter, period, multiplier: -1, near, digits });
  const lo = born.parameter.low.compareTo(doubles.parameter.low) < 0 ? born.parameter.low : doubles.parameter.low;
  const hi = born.parameter.high.compareTo(doubles.parameter.high) > 0 ? born.parameter.high : doubles.parameter.high;
  let middle = [];
  try {
    middle = superstable(family, { parameter, period, range: [lo, hi], digits });
  } catch {
    middle = [];
  }
  return { period, born: born.parameter, doubles: doubles.parameter, superstable: middle[0] ?? null, width: out(Box.of(doubles.parameter).sub(Box.of(born.parameter)).abs()) };
}

// ---- symbolic dynamics -------------------------------------------------------------

/**
 * The itinerary of x₀: L, C or R for each of the first n points of its
 * orbit according to which side of the critical point c it lies on. The
 * orbit is the true one; a point too close to c to place gets more working
 * digits until it is decided (or is C exactly, for exact orbits).
 */
export function itinerary(m, x0, n, { critical } = {}) {
  asMap(m);
  const c = critical === undefined ? m.criticalPoints([-1e6, 1e6])[0] : toQ(critical);
  if (!(c instanceof Rational)) refuse("no-critical-point", "give the turning point as { critical }");
  // exact while the numbers stay small, enclosures after that
  const symbols = [];
  let x = toQ(x0);
  let exact = true;
  let digits = 20;
  for (let i = 0; i < n; i++) {
    if (exact) {
      const s = x.compareTo(c);
      symbols.push(s < 0 ? "L" : s > 0 ? "R" : "C");
      const next = exactValue(m.tree, { [m.variable]: x });
      if (next === null || next.denominator.toString(2).length > 4000) exact = false;
      else x = next;
      continue;
    }
    for (;;) {
      const t = Box.of(m.enclosedOrbit(x0, i, { digits }).terms[i]);
      if (t.hi.compareTo(c) < 0) symbols.push("L");
      else if (t.lo.compareTo(c) > 0) symbols.push("R");
      else if (digits < 2000) {
        digits *= 2;
        continue;
      } else symbols.push("?");
      break;
    }
  }
  return symbols.join("");
}

/** The kneading sequence: the itinerary of the critical value f(c). */
export function kneading(m, n, { critical } = {}) {
  asMap(m);
  const c = critical === undefined ? m.criticalPoints([-1e6, 1e6])[0] : toQ(critical);
  if (!(c instanceof Rational)) refuse("no-critical-point", "give the turning point as { critical }");
  return itinerary(m, m.at(c), n, { critical: c });
}

/** where n sits in Sharkovskii's order: odd part and power of two */
const split = (n) => {
  if (!Number.isInteger(n) || n < 1) refuse("bad-period", "periods are positive whole numbers");
  let k = 0;
  let m = n;
  while (m % 2 === 0) {
    m /= 2;
    k++;
  }
  return { odd: m, two: k };
};
/**
 * Sharkovskii's order 3 ▷ 5 ▷ 7 ▷ … ▷ 2·3 ▷ 2·5 ▷ … ▷ 4·3 ▷ … ▷ 8 ▷ 4 ▷ 2 ▷ 1.
 * sharkovskii.compare(a, b) is negative when a comes first; a continuous map
 * of an interval with a point of period a has points of every period after a.
 */
export const sharkovskii = {
  compare(a, b) {
    const A = split(a);
    const B = split(b);
    if (a === b) return 0;
    const pa = A.odd === 1;
    const pb = B.odd === 1;
    if (pa && pb) return B.two - A.two; // powers of two come last, in decreasing order
    if (pa) return 1;
    if (pb) return -1;
    if (A.two !== B.two) return A.two - B.two;
    return A.odd - B.odd;
  },
  /** does a point of period a force a point of period b? */
  forces(a, b) {
    return sharkovskii.compare(a, b) <= 0;
  },
  /** the periods up to `limit` that a point of period n forces, in Sharkovskii's order */
  forced(n, limit = 32) {
    return sharkovskii.order(limit).filter((k) => sharkovskii.forces(n, k));
  },
  /** 1 … limit sorted into Sharkovskii's order */
  order(limit = 16) {
    return Array.from({ length: limit }, (_, i) => i + 1).sort(sharkovskii.compare);
  },
};
