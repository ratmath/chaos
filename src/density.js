// The Lyapunov exponent of a map from its invariant measure.
//
// For a map f with an ergodic invariant probability measure μ, almost every
// orbit has the same exponent, the space average
//
//     λ = ∫ log|f′(x)| dμ(x).
//
// A finite orbit average (lyapunov()) only approaches this. When the measure
// is known, the integral can be enclosed instead, and the enclosure is a
// statement about the limit itself.
//
// The measure is given by its distribution function F(x) = μ([a, x]); the
// density is F′. Working with F has a point: the density of the logistic map
// at r = 4, 1/(π·sqrt(x(1 − x))), is unbounded at both ends, but F is
// continuous, and the measure of a cell is just F(r) − F(l).
//
// On a cell [l, r] with midpoint m, width h and g = log|f′|:
//
//     ∫ g dF = g(m)·ΔF + g′(m)·M₁ + ½·g″(ξ)·∫ (x − m)² dF,
//     M₁ = ∫ (x − m) dF = (h/2)·(F(l) + F(r)) − ∫ F dx        (by parts),
//
// with ∫ (x − m)² dF in [0, h²/4]·ΔF, and ∫ F dx between the midpoint and
// trapezoid rules when F is convex or concave on the cell (between h·F(l)
// and h·F(r) otherwise). Every quantity is an interval.
//
// At a critical point c (f′(c) = 0) g has a logarithmic singularity. On a
// cell with c at one end, f′(x) = f″(ξ)·(x − c), so
//
//     ∫ g dF = log|f″(ξ)|·ΔF + ∫ log|x − c|·F′(x) dx,
//
// and the last integral lies between (sup F′)·J and (inf F′)·J with
// J = ∫₀^δ log u du = δ·log δ − δ < 0.
import { Rational, RationalInterval } from "@ratmath/core";
import { zeros } from "@ratmath/analytic";
import { IteratedMap } from "./map.js";
import { refuse, ZERO, ONE, Q, toQ, Box, parse, format, differentiate, simplify, exactValue, boxValue, logAbs } from "./tree.js";

const HALF = Q(1, 2);
const soft = (f) => {
  try {
    return f();
  } catch (e) {
    if (e?.name === "ChaosError" || e?.name === "AnalyticError" || e instanceof RangeError) return null;
    throw e;
  }
};
const scaleBox = (b, q) => (q.numerator >= 0n ? new Box(b.lo.multiply(q), b.hi.multiply(q)) : new Box(b.hi.multiply(q), b.lo.multiply(q)));
const hullQ = (p, q) => (p.compareTo(q) <= 0 ? new Box(p, q) : new Box(q, p));

/** abs(u) replaced by u or −u according to the sign of u at a point (for one side of a corner). */
function resolveAbs(t, name, at) {
  if (t.type === "call" && t.name === "abs") {
    const inner = resolveAbs(t.args[0], name, at);
    const v = exactValue(inner, { [name]: at });
    if (v === null || v.numerator === 0n) refuse("unsupported", "the sign inside an absolute value could not be settled on this piece");
    return v.numerator > 0n ? inner : { type: "op", op: "neg", args: [inner] };
  }
  return t.args ? { ...t, args: t.args.map((a) => resolveAbs(a, name, at)) } : t;
}

/**
 * Whether the tree is ≥ 0 wherever it is defined on the cell: by its shape
 * (products, quotients and sums of things that are ≥ 0, square roots, even
 * powers), or else by its enclosure. The shape test works where a
 * denominator touches 0 at an end of the cell and no enclosure exists.
 */
function nonnegOn(t, name, cell, bits) {
  const byBox = () => {
    const b = soft(() => boxValue(t, { [name]: cell }, bits));
    return b !== null && b.lo.numerator >= 0n;
  };
  switch (t.type) {
    case "num":
      return t.value.numerator >= 0n;
    case "const":
      return true; // pi and e
    case "op":
      if (t.op === "*" || t.op === "/") return (t.args.every((a) => nonnegOn(a, name, cell, bits)) || t.args.every((a) => nonnegOn({ type: "op", op: "neg", args: [a] }, name, cell, bits))) || byBox();
      if (t.op === "+") return t.args.every((a) => nonnegOn(a, name, cell, bits)) || byBox();
      if (t.op === "^" && t.args[1].type === "num" && t.args[1].value.denominator === 1n && t.args[1].value.numerator % 2n === 0n) return true;
      if (t.op === "neg") {
        const a = t.args[0];
        if (a.type === "num") return a.value.numerator <= 0n;
        if (a.type === "op" && a.op === "neg") return nonnegOn(a.args[0], name, cell, bits);
        return byBox();
      }
      return byBox();
    case "call":
      if (t.name === "sqrt" || t.name === "abs" || t.name === "exp") return true;
      return byBox();
  }
  return byBox();
}

const KNOWN = [
  {
    test: (m) => m.name === "logistic" && m.params.r?.equals(Q(4)),
    distribution: "2/pi*asin(sqrt(x))",
    density: "1/(pi*sqrt(x*(1 - x)))",
    interval: [ZERO, ONE],
    exponent: "log(2)",
    reason: "x = sin(pi*θ/2)^2 turns the map into the tent map θ → 1 − |1 − 2θ|, which keeps length; the image of length under that change of variable has this distribution",
  },
  {
    test: (m) => m.name === "tent" && m.params.mu?.equals(Q(2)),
    distribution: "x",
    density: "1",
    interval: [ZERO, ONE],
    exponent: "log(2)",
    reason: "each interval has two preimages of half its length, so length is kept",
  },
  {
    test: (m) => m.name === "quadratic" && m.params.c?.equals(Q(-2)),
    distribution: "1/2 + asin(x/2)/pi",
    density: "1/(pi*sqrt(4 - x^2))",
    interval: [Q(-2), Q(2)],
    exponent: "log(2)",
    reason: "x = 2cos(pi*θ) turns the map into angle doubling, which keeps length",
  },
];

/**
 * The invariant density of a map where it is known in closed form: the
 * logistic map at r = 4, the tent map at slope 2 and x² − 2.
 * → { density, distribution, interval, exponent, reason }.
 */
export function invariantDensity(m) {
  if (!(m instanceof IteratedMap)) refuse("bad-map", "give a map made by map(), logistic(), tent() or quadratic()");
  const k = KNOWN.find((e) => e.test(m));
  if (!k) refuse("unknown-density", "no closed form for the invariant density of this map is known here; pass { distribution, interval } to lyapunovExponent yourself");
  return { density: k.density, distribution: k.distribution, interval: k.interval, exponent: k.exponent, reason: k.reason };
}

/**
 * λ = ∫ log|f′| dμ for the measure with distribution function F, as a
 * certified enclosure.
 *
 *   lyapunovExponent(logistic(4))                       // an interval around log 2
 *   lyapunovExponent(map("4*x^3 - 3*x"), { distribution: "1 - acos(x)/pi", interval: [-1, 1] })   // log 3
 *
 * → { value, kind, density, distribution, interval, criticalPoints, cells,
 *     invariance: { proven, reason, checks, holds }, exponent, note }.
 *
 * The integral is enclosed for the measure given. That the measure is
 * invariant is proven only for the built-in maps (by the change of variable
 * named in `invariance.reason`); for a distribution you supply it is tested:
 * at several levels y the measure of {x : f(x) ≤ y} is computed from certified
 * preimages and compared with F(y).
 */
export function lyapunovExponent(m, { distribution, interval, digits = 5, checks = 5, maxCells = 20000 } = {}) {
  if (!(m instanceof IteratedMap)) refuse("bad-map", "give a map made by map(), logistic(), tent() or quadratic()");
  let known = null;
  if (distribution === undefined) {
    known = invariantDensity(m);
    distribution = known.distribution;
    interval = known.interval;
  }
  if (!interval) refuse("bad-argument", "give the interval [a, b] that carries the measure");
  const x = m.variable;
  const [a, b] = interval.map(toQ);
  if (a.compareTo(b) >= 0) refuse("bad-argument", "the interval needs a < b");
  const bits = Math.ceil(digits * 3.33) + 56;
  const F = typeof distribution === "string" ? parse(distribution) : distribution;
  const F1 = simplify(differentiate(F, x));
  const F2 = simplify(differentiate(F1, x));
  const Fmemo = new Map();
  const Fat = (q) => {
    const key = q.toString();
    let v = Fmemo.get(key);
    if (v === undefined) Fmemo.set(key, (v = boxValue(F, { [x]: new Box(q) }, bits).round(bits)));
    return v;
  };
  // F is a distribution function on [a, b]
  const [Fa, Fb] = [Fat(a), Fat(b)];
  const tiny = Q(1, 10 ** 9);
  if (Fa.mag.compareTo(tiny) > 0 || Fb.sub(new Box(ONE)).mag.compareTo(tiny) > 0) refuse("bad-distribution", `the distribution must run from 0 at ${a} to 1 at ${b}; it runs from about ${Fa.mid.toNumber()} to ${Fb.mid.toNumber()}`);

  // pieces between corners of the map; on each, the derivatives of f as formulas
  const corners = m.critical !== null && /abs\(/.test(format(m.tree)) && m.critical.compareTo(a) > 0 && m.critical.compareTo(b) < 0 ? [m.critical] : [];
  const edges = [a, ...corners, b];
  const pieces = [];
  const criticalPoints = [];
  for (let i = 0; i + 1 < edges.length; i++) {
    const [l, r] = [edges[i], edges[i + 1]];
    const tree = corners.length ? simplify(resolveAbs(m.tree, x, l.add(r).multiply(HALF))) : m.tree;
    const d1 = simplify(differentiate(tree, x));
    const d2 = simplify(differentiate(d1, x));
    const d3 = simplify(differentiate(d2, x));
    const g1 = { type: "op", op: "/", args: [d2, d1] };
    const g2 = { type: "op", op: "/", args: [{ type: "op", op: "-", args: [{ type: "op", op: "*", args: [d3, d1] }, { type: "op", op: "^", args: [d2, { type: "num", value: Q(2) }] }] }, { type: "op", op: "^", args: [d1, { type: "num", value: Q(2) }] }] };
    const cuts = [];
    if (d1.type !== "num") {
      for (const z of zeros(format(d1), [l, r], { variable: x }).exact()) {
        if (z.kind === "everywhere") refuse("degenerate", "the map is constant on an interval");
        if (z.kind !== "exact") refuse("unsupported", `a critical point near ${z.lo.toNumber()} is not a rational number; the bound at a critical point needs its exact position`);
        cuts.push(z.at);
        criticalPoints.push(z.at);
      }
    } else if (d1.value.numerator === 0n) refuse("degenerate", "the map is constant on an interval");
    pieces.push({ l, r, d1, d2, g1, g2, cuts: cuts.sort((p, q) => p.compareTo(q)) });
  }

  // one cell's share of the integral
  const share = (cell) => {
    const { l, r, piece, at } = cell; // at: the critical point at one end, or null
    const h = r.subtract(l);
    const dF = Fat(r).sub(Fat(l));
    const box = new Box(l, r);
    if (at !== null) {
      const second = soft(() => boxValue(piece.d2, { [x]: box }, bits));
      const rho = soft(() => boxValue(F1, { [x]: box }, bits));
      if (!second || !rho || second.containsZero() || rho.lo.numerator < 0n || h.compareTo(ONE) > 0) return null;
      const J = scaleBox(logAbs(new Box(h), bits).sub(new Box(ONE)), h); // δ·log δ − δ, negative
      return logAbs(second, bits).mul(dF).add(new Box(rho.hi.multiply(J.lo), rho.lo.multiply(J.hi))).round(bits);
    }
    if (!nonnegOn(F1, x, box, bits)) return null; // the measure must not be negative anywhere on the cell
    const mid = l.add(r).multiply(HALF);
    const pm = { [x]: new Box(mid) };
    const d1m = soft(() => boxValue(piece.d1, pm, bits));
    const curve = soft(() => boxValue(piece.g2, { [x]: box }, bits));
    if (!d1m || d1m.containsZero()) return null;
    if (!curve) {
      // first order: the range of g times the measure
      const d1box = soft(() => boxValue(piece.d1, { [x]: box }, bits));
      return d1box && !d1box.containsZero() ? logAbs(d1box, bits).mul(dF).round(bits) : null;
    }
    const gm = logAbs(d1m, bits);
    const slope = boxValue(piece.g1, pm, bits);
    // ∫ F dx on the cell
    const trap = scaleBox(Fat(l).add(Fat(r)), h.multiply(HALF));
    const midRule = scaleBox(Fat(mid), h);
    const bend = soft(() => boxValue(F2, { [x]: box }, bits));
    let intF;
    if (bend && (bend.lo.numerator >= 0n || bend.hi.numerator <= 0n)) intF = new Box(trap.lo.compareTo(midRule.lo) < 0 ? trap.lo : midRule.lo, trap.hi.compareTo(midRule.hi) > 0 ? trap.hi : midRule.hi);
    else intF = new Box(Fat(l).lo.multiply(h), Fat(r).hi.multiply(h));
    const M1 = trap.sub(intF);
    const quad = scaleBox(curve, HALF).mul(new Box(ZERO, h.multiply(h).multiply(Q(1, 4)))).mul(dF);
    return gm.mul(dF).add(slope.mul(M1)).add(quad).round(bits);
  };

  let cells = [];
  for (const piece of pieces) {
    const marks = [piece.l, ...piece.cuts, piece.r].filter((q, i, arr) => i === 0 || !q.equals(arr[i - 1]));
    for (let i = 0; i + 1 < marks.length; i++) {
      const [l, r] = [marks[i], marks[i + 1]];
      const isCut = (q) => piece.cuts.some((c) => c.equals(q));
      // start with a few cells, the ones next to a critical point marked
      const n = 8;
      for (let k = 0; k < n; k++) {
        const cl = l.add(r.subtract(l).multiply(Q(k, n)));
        const cr = l.add(r.subtract(l).multiply(Q(k + 1, n)));
        cells.push({ l: cl, r: cr, piece, at: k === 0 && isCut(l) ? l : k === n - 1 && isCut(r) ? r : null });
      }
    }
  }
  const split = (c) => {
    const mid = c.l.add(c.r).multiply(HALF);
    return [
      { l: c.l, r: mid, piece: c.piece, at: c.at !== null && c.at.equals(c.l) ? c.at : null },
      { l: mid, r: c.r, piece: c.piece, at: c.at !== null && c.at.equals(c.r) ? c.at : null },
    ];
  };
  const weigh = (c) => {
    c.value = share(c);
    c.depth = c.depth ?? 0;
    return c;
  };
  cells = cells.map(weigh);
  const tol = Q(1, 10 ** digits);
  let total = null;
  let stuck = null;
  for (let round = 0; round < 60; round++) {
    const open = cells.filter((c) => c.value === null);
    const widths = cells.map((c) => (c.value === null ? null : c.value.width));
    const sum = widths.reduce((s, w) => (w === null ? s : s.add(w)), ZERO);
    if (!open.length && sum.compareTo(tol) <= 0) break;
    if (cells.length >= maxCells) break;
    const bar = sum.divide(Q(2 * cells.length));
    const next = [];
    for (const c of cells) {
      const needs = c.value === null || c.value.width.compareTo(bar) > 0;
      if (!needs) next.push(c);
      else if (c.depth >= 60) {
        stuck = c;
        next.push(c);
      } else for (const half of split(c)) next.push(weigh({ ...half, depth: c.depth + 1 }));
    }
    cells = next;
    if (stuck) break;
  }
  const failed = cells.find((c) => c.value === null);
  if (failed) refuse("no-enclosure", `log|f′| could not be integrated near ${x} = ${failed.l.toNumber()}: the derivative of the map or of the distribution has no enclosure there`);
  total = cells.reduce((s, c) => s.add(c.value), new Box(ZERO)).round(bits);

  // is the measure invariant? μ({f ≤ y}) against F(y) at a few levels
  const fx = (q) => m.float(q.toNumber());
  const levels = [];
  for (let i = 1; i <= checks; i++) {
    const y = a.add(b.subtract(a).multiply(Q(2 * i - 1, 2 * checks + 1)));
    const one = soft(() => {
      const roots = zeros(`${format(m.tree)} - (${y})`, [a, b], { variable: x }).exact();
      if (roots.some((z) => z.kind === "everywhere")) return null;
      const marks = [new Box(a), ...roots.map((z) => (z.at ? new Box(z.at) : new Box(z.lo, z.hi))).sort((p, q) => p.lo.compareTo(q.lo)), new Box(b)];
      let measure = new Box(ZERO);
      for (let k = 0; k + 1 < marks.length; k++) {
        const [p, q] = [marks[k], marks[k + 1]];
        if (p.hi.compareTo(q.lo) >= 0) continue;
        if (fx(p.hi.add(q.lo).multiply(HALF)) <= y.toNumber()) measure = measure.add(boxValue(F, { [x]: q }, bits).sub(boxValue(F, { [x]: p }, bits)));
      }
      const want = Fat(y);
      const holds = measure.lo.compareTo(want.hi) <= 0 && want.lo.compareTo(measure.hi) <= 0 && measure.width.compareTo(Q(1, 1000)) < 0;
      return { y, measure: measure.toRationalInterval(), distribution: want.toRationalInterval(), holds };
    });
    if (one) levels.push(one);
  }
  const holds = levels.length > 0 && levels.every((c) => c.holds);
  if (levels.length && !holds) refuse("not-invariant", `the measure is not invariant under the map: at level ${levels.find((c) => !c.holds).y} the measure of {f ≤ y} differs from F(y)`);

  return {
    value: total.toRationalInterval(),
    kind: "certified integral of log|f′| against the invariant measure",
    distribution: format(F),
    density: known ? known.density : format(F1),
    interval: [a, b],
    criticalPoints,
    cells: cells.length,
    exponent: known ? known.exponent : null,
    invariance: {
      proven: known !== null,
      reason: known ? known.reason : `tested at ${levels.length} levels: the measure of {f ≤ y} from certified preimages agrees with F(y)`,
      checks: levels,
      holds,
    },
    note: "this is the exponent of almost every orbit when the measure is ergodic; a finite orbit average (lyapunov) approaches it",
  };
}
