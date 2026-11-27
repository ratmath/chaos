// Maps of the plane, iterated. As on the line, an orbit comes exact (rational
// arithmetic, the true orbit), enclosed (intervals certain to hold it) or in
// doubles, and the three are put side by side. Periodic points are proven to
// exist, one in each small box, with their multipliers and the stable and
// unstable directions; for polynomial maps the low periods are also found
// exactly as algebraic numbers. Two maps are exact throughout: automorphisms
// of the torus such as the cat map, where every rational point is periodic,
// and Smale's horseshoe in its piecewise linear form.
import { Rational, RationalInterval } from "@ratmath/core";
import { roots as algebraicRoots, format as algebraicFormat, add as aAdd, mul as aMul, inv as aInv, fromRational, equals as aEquals } from "@ratmath/algebraic";
import { ChaosError, refuse, ZERO, ONE, Q, toQ, Box, P, parse, format, differentiate, simplify, compileFloat, exactValue, boxValue, logAbs, substitute } from "./tree.js";
import * as F from "./fi.js";

const num = (q) => ({ type: "num", value: q });
const bitLength = (q) => q.numerator.toString(2).length + q.denominator.toString(2).length;
const properDivisors = (n) => Array.from({ length: n - 1 }, (_, i) => i + 1).filter((d) => n % d === 0);
const showEigen = (ev) => ({ re: F.toRI(ev.re), im: F.toRI(ev.im) });
const modulus = (ev) => F.sqrt(F.add(F.sqr(ev.re), F.sqr(ev.im)));

// ---- polynomials in several symbols (a term is a product of powers, possibly negative) ----
// Used to expand a formula with its sines, cosines and π treated as symbols:
// enough to see that a Jacobian determinant is constant, and to eliminate a
// variable from the equations of a periodic point.

const MP = {
  constant: (q) => (q.numerator === 0n ? new Map() : new Map([["", { c: q, e: {} }]])),
  symbol: (name) => new Map([[`${name}^1`, { c: ONE, e: { [name]: 1 } }]]),
  key: (e) => Object.keys(e).sort().map((k) => `${k}^${e[k]}`).join(" "),
  add(a, b, sign = 1) {
    const out = new Map(a);
    for (const [k, t] of b) {
      const c = (out.get(k)?.c ?? ZERO).add(sign === 1 ? t.c : t.c.negate());
      if (c.numerator === 0n) out.delete(k);
      else out.set(k, { c, e: t.e });
    }
    return out;
  },
  mul(a, b) {
    let out = new Map();
    for (const s of a.values()) {
      const part = new Map();
      for (const t of b.values()) {
        const e = { ...s.e };
        for (const [k, v] of Object.entries(t.e)) {
          e[k] = (e[k] ?? 0) + v;
          if (e[k] === 0) delete e[k];
        }
        const key = MP.key(e);
        const c = (part.get(key)?.c ?? ZERO).add(s.c.multiply(t.c));
        if (c.numerator === 0n) part.delete(key);
        else part.set(key, { c, e });
      }
      out = MP.add(out, part);
    }
    return out;
  },
  pow(a, k) {
    let out = MP.constant(ONE);
    for (let i = 0; i < k; i++) out = MP.mul(out, a);
    return out;
  },
  /** the constant when the polynomial is one, else null */
  constantOf(a) {
    if (a.size === 0) return ZERO;
    return a.size === 1 && a.has("") ? a.get("").c : null;
  },
  degree: (a, name) => Math.max(0, ...[...a.values()].map((t) => t.e[name] ?? 0)),
  /** coefficients in one symbol, constant term first, each a polynomial in the others */
  coefficients(a, name) {
    const out = [];
    for (const t of a.values()) {
      const k = t.e[name] ?? 0;
      if (k < 0) return null;
      const e = { ...t.e };
      delete e[name];
      out[k] = MP.add(out[k] ?? new Map(), new Map([[MP.key(e), { c: t.c, e }]]));
    }
    for (let i = 0; i < out.length; i++) out[i] ??= new Map();
    return out;
  },
};
/** expand a formula tree; variables listed in `vars` stay symbols, as do function calls and π */
function expand(ast, vars) {
  const go = (n) => {
    switch (n.type) {
      case "num":
        return MP.constant(n.value);
      case "var":
        return vars.includes(n.name) ? MP.symbol(n.name) : null;
      case "const":
        return MP.symbol(`#${n.name}`);
      case "call": {
        const inner = n.args.map((a) => go(a));
        if (inner.includes(null)) return null;
        return MP.symbol(`#${format(n)}`);
      }
      case "op": {
        if (n.op === "^") {
          const e = n.args[1];
          const inner = e.type === "op" && e.op === "neg" ? e.args[0] : e;
          const b = go(n.args[0]);
          if (!b || inner.type !== "num" || inner.value.denominator !== 1n || inner !== e) return null;
          const k = Number(inner.value.numerator);
          return k > 64 ? null : MP.pow(b, k);
        }
        const a = go(n.args[0]);
        if (!a) return null;
        if (n.op === "neg") return MP.mul(MP.constant(Q(-1)), a);
        const b = go(n.args[1]);
        if (!b) return null;
        switch (n.op) {
          case "+": return MP.add(a, b);
          case "-": return MP.add(a, b, -1);
          case "*": return MP.mul(a, b);
          case "/": {
            // only by a single term
            if (b.size !== 1) return null;
            const [t] = b.values();
            const e = Object.fromEntries(Object.entries(t.e).map(([k, v]) => [k, -v]));
            return MP.mul(a, new Map([[MP.key(e), { c: ONE.divide(t.c), e }]]));
          }
        }
        return null;
      }
      default:
        return null;
    }
  };
  return go(ast);
}

// ---- maps from formulas ----------------------------------------------------------------------

export class PlaneMap {
  /**
   * @param field the two components as formulas: ["1 - a*x^2 + y", "b*x"]
   * @param options { variables: ["x", "y"], params: { a: "1.4", b: "0.3" }, name }
   */
  constructor(field, options = {}) {
    if (!Array.isArray(field) || field.length !== 2) refuse("bad-map", "a map of the plane has two components");
    this.variables = options.variables ?? ["x", "y"];
    if (this.variables.length !== 2) refuse("bad-map", "a map of the plane has two variables");
    this.params = Object.fromEntries(Object.entries(options.params ?? {}).map(([k, v]) => [k, toQ(v)]));
    this.family = field.map((f) => (typeof f === "string" ? parse(f) : f));
    this.text = field.map((f) => (typeof f === "string" ? f : format(f)));
    this.trees = this.family.map((t) => {
      for (const [k, v] of Object.entries(this.params)) t = substitute(t, k, num(v));
      return t;
    });
    this.name = options.name ?? null;
    this._float = this.trees.map(compileFloat);
    this._J = this.trees.map((t) => this.variables.map((v) => simplify(differentiate(t, v))));
    this._tapes = null;
    const det = expand({ type: "op", op: "-", args: [{ type: "op", op: "*", args: [this._J[0][0], this._J[1][1]] }, { type: "op", op: "*", args: [this._J[0][1], this._J[1][0]] }] }, this.variables);
    /** the Jacobian determinant when it is the same everywhere (a Rational), else null */
    this.determinant = det ? MP.constantOf(det) : null;
    /** true when the map keeps areas: the determinant is 1 or −1 everywhere */
    this.areaPreserving = this.determinant !== null && this.determinant.abs().equals(ONE);
  }
  with(params) {
    return new PlaneMap(this.family, { variables: this.variables, params: { ...this.params, ...params }, name: this.name });
  }
  toString() {
    const ps = Object.entries(this.params).map(([k, v]) => `${k} = ${v}`).join(", ");
    return `(${this.variables.join(", ")}) → (${this.text.join(", ")})${ps ? ` (${ps})` : ""}`;
  }
  #env(p) {
    return { [this.variables[0]]: p[0], [this.variables[1]]: p[1] };
  }
  /** F(p), exactly. Refuses when a formula has no rational value there. */
  at(p) {
    const env = this.#env(p.map(toQ));
    const out = this.trees.map((t) => exactValue(t, env));
    if (out.includes(null)) refuse("not-exact", "this map has no exact rational values; use enclose(p) or enclosedOrbit");
    return out;
  }
  /** F(p) in doubles */
  float(p) {
    const env = this.#env(p);
    return this._float.map((f) => f(env));
  }
  /** intervals certain to contain F(p) for every p in the box */
  enclose(p, bits = 96) {
    const env = this.#env(p.map((x) => Box.of(x)));
    return this.trees.map((t) => boxValue(t, env, bits).toRationalInterval());
  }
  /** the Jacobian matrix as formulas */
  jacobian() {
    return this._J.map((row) => row.map(format));
  }
  /** the true orbit in exact rationals; stops with an error once a term passes maxBits binary digits */
  orbit(p0, n, { maxBits = 200000 } = {}) {
    const out = [p0.map(toQ)];
    for (let i = 0; i < n; i++) {
      const next = this.at(out[i]);
      if (bitLength(next[0]) + bitLength(next[1]) > maxBits) refuse("too-large", `the exact orbit needs more than ${maxBits} binary digits by step ${i + 1}; use enclosedOrbit for long orbits`);
      out.push(next);
    }
    return out;
  }
  /** the orbit as an ordinary program computes it */
  floatOrbit(p0, n) {
    const out = [p0.map((x) => (typeof x === "number" ? x : toQ(x).toNumber()))];
    for (let i = 0; i < n; i++) out.push(this.float(out[i]));
    return out;
  }
  /**
   * The true orbit as intervals. Working digits are added until the last
   * point is known to `digits` decimal places, however chaotic the map:
   * { points: [[RationalInterval, RationalInterval]], bits, width }.
   */
  enclosedOrbit(p0, n, { digits = 12, bits, maxBits = 1 << 17 } = {}) {
    const start = p0.map((x) => Box.of(x));
    const want = ONE.divide(Q(10).pow(digits));
    let work = bits ?? Math.max(64, Math.ceil(digits * 3.33) + 32);
    for (;;) {
      const points = [start];
      let ok = true;
      try {
        for (let i = 0; i < n && ok; i++) {
          const env = this.#env(points[i]);
          const next = this.trees.map((t) => boxValue(t, env, work));
          points.push(next);
          if (next.some((b) => b.width.compareTo(ONE) > 0)) ok = false;
          else if (next.some((b) => b.mag.numerator.toString(2).length > 4096)) refuse("escapes", `the orbit leaves every bounded range (it passes 2^4000 by step ${i + 1})`);
        }
      } catch (e) {
        if (!(e instanceof ChaosError) || e.code !== "division") throw e;
        ok = false;
      }
      const width = ok ? (points[n][0].width.compareTo(points[n][1].width) > 0 ? points[n][0].width : points[n][1].width) : null;
      if ((ok && width.compareTo(want) <= 0) || bits !== undefined || work >= maxBits) {
        if (!ok) refuse("not-enclosed", `the orbit could not be enclosed with ${work} binary digits`);
        return { points: points.map((p) => p.map((b) => b.toRationalInterval())), bits: work, width };
      }
      work *= 2;
    }
  }
  /**
   * The float orbit beside the true one: rows { n, float, exact, error,
   * digits }, and `parted`, the first step at which the doubles have no
   * correct digit left.
   */
  compare(p0, n, { digits = 20 } = {}) {
    const floats = this.floatOrbit(p0, n);
    const truth = this.enclosedOrbit(p0, n, { digits }).points;
    let parted = null;
    const rows = floats.map((f, i) => {
      const mids = truth[i].map((iv) => iv.low.add(iv.high).divide(Q(2)).toNumber());
      const error = Math.hypot(f[0] - mids[0], f[1] - mids[1]);
      const size = Math.max(Math.hypot(...mids), 1e-300);
      const good = error === 0 ? 17 : Math.max(0, Math.min(17, -Math.log10(error / size)));
      if (parted === null && good < 1) parted = i;
      return { n: i, float: f, exact: truth[i], error, digits: good };
    });
    return { rows, parted };
  }

  #tape() {
    return (this._tapes ??= { f: F.compile(this.trees, this.variables), jac: F.compile(this._J.flat(), this.variables) });
  }
  /** Fⁿ and its derivative over a box of doubles: { image, derivative } */
  #iterate(X, n) {
    const { f, jac } = this.#tape();
    let cur = X;
    let D = F.ivM(F.identity(2));
    for (let i = 0; i < n; i++) {
      const j = F.values(jac, cur);
      D = F.mmul([[j[0], j[1]], [j[2], j[3]]], D);
      cur = F.values(f, cur);
    }
    return { image: cur, derivative: D };
  }
  #describe(box, n) {
    const { derivative: M } = this.#iterate(box, n);
    const ev = F.eigenvalues(M);
    let stability = "undecided";
    let stable = null;
    let unstable = null;
    if (ev) {
      const mods = ev.map(modulus);
      const complex = ev[0].im[0] > 0 || ev[0].im[1] < 0;
      if (mods.every((m) => m[1] < 1)) stability = "attracting";
      else if (mods.every((m) => m[0] > 1)) stability = "repelling";
      else if (!complex && mods.every((m) => m[1] < 1 || m[0] > 1)) stability = "saddle";
      else if (complex && this.areaPreserving) stability = "elliptic";
      if (!complex) {
        // an eigenvector for λ is (b, λ − a) or (λ − d, c): take the better conditioned, scaled to have a 1
        const direction = (lam) => {
          const u = [M[0][1], F.sub(lam, M[0][0])];
          const w = [F.sub(lam, M[1][1]), M[1][0]];
          const size = (q) => Math.max(F.mag(q[0]), F.mag(q[1]));
          const q = size(u) >= size(w) ? u : w;
          const k = Math.abs(F.mid(q[0])) >= Math.abs(F.mid(q[1])) ? 0 : 1;
          if (F.hasZero(q[k])) return null;
          const other = F.div(q[1 - k], q[k]);
          return (k === 0 ? [F.pt(1), other] : [other, F.pt(1)]).map(F.toRI);
        };
        if (stability === "saddle") {
          const s = mods[0][1] < 1 ? 0 : 1;
          stable = direction(ev[s].re);
          unstable = direction(ev[1 - s].re);
        }
      }
    }
    return { multipliers: ev ? ev.map(showEigen) : null, stability, stable, unstable };
  }
  /**
   * Every point of period n in the region { x: [a, b], y: [c, d] }, each
   * proven to exist and to be the only one in its enclosure:
   * [{ point, period, leastPeriod, multipliers, stability, stable, unstable }].
   * The multipliers are the eigenvalues of the derivative of Fⁿ there; at a
   * saddle, `stable` and `unstable` are its two directions. `complete` on the
   * list says the rest of the region was proven to hold no other. Points
   * whose least period is smaller are left out unless `all` is set.
   */
  periodicPoints(n, region, { all = false, tol = 1e-9, maxBoxes = 40000 } = {}) {
    if (!Number.isInteger(n) || n < 1) refuse("bad-period", "the period is a positive whole number");
    const start = this.variables.map((v) => {
      if (!region?.[v]) refuse("bad-region", `the region needs a range for ${v}`);
      return F.ofExact(region[v]);
    });
    const equations = (k) => [
      (X) => this.#iterate(X, k).image.map((x, i) => F.sub(x, X[i])),
      (X) => this.#iterate(X, k).derivative.map((row, i) => row.map((a, j) => (i === j ? F.sub(a, F.pt(1)) : a))),
    ];
    const [G, DG] = equations(n);
    const list = F.solveAll(G, DG, start, { tol, maxBoxes });
    const out = [];
    for (const box of list) {
      // the least period: a smaller one is ruled out when the image under F^d misses the box
      let least = n;
      for (const d of properDivisors(n)) {
        const img = this.#iterate(box, d).image;
        if (img.some((x, i) => F.meet(x, box[i]) === null)) continue;
        const wide = box.map((x) => [x[0] - 1e-9 * (1 + F.mag(x)), x[1] + 1e-9 * (1 + F.mag(x))]);
        least = F.krawczyk(...equations(d), wide) ? d : null;
        break;
      }
      if (!all && least !== n) continue;
      out.push({ point: box.map(F.toRI), period: n, leastPeriod: least, ...this.#describe(box, n) });
    }
    out.sort((p, q) => p.point[0].low.compareTo(q.point[0].low) || p.point[1].low.compareTo(q.point[1].low));
    out.complete = list.complete;
    return out;
  }
  fixedPoints(region, options) {
    return this.periodicPoints(1, region, options);
  }
  /** the periodic points of least period n grouped into orbits: [{ period, points, multipliers, stability }] */
  cycles(n, region, options) {
    const points = this.periodicPoints(n, region, options);
    const boxes = points.map((p) => p.point.map((iv) => [F.ofExact(iv.low)[0], F.ofExact(iv.high)[1]]));
    const left = points.map((_, i) => i);
    const out = [];
    while (left.length) {
      const first = left.shift();
      const members = [first];
      let cur = boxes[first];
      for (let step = 1; step < n; step++) {
        cur = this.#iterate(cur, 1).image;
        const hits = left.filter((i) => boxes[i].every((x, k) => F.meet(x, cur[k]) !== null));
        if (hits.length !== 1) break;
        members.push(hits[0]);
        left.splice(left.indexOf(hits[0]), 1);
      }
      out.push({ period: n, points: members.map((i) => points[i].point), multipliers: points[first].multipliers, stability: points[first].stability, whole: members.length === n });
    }
    out.complete = points.complete;
    return out;
  }
  /**
   * The points of period n exactly, as algebraic numbers, for a polynomial
   * map in which one coordinate can be eliminated (one of the two equations
   * of a periodic point is of first degree in it). For the Hénon family that
   * covers periods 1 and 2. [{ x, y, text, period, multipliers, stability }].
   */
  exactPeriodicPoints(n = 1) {
    const [xn, yn] = this.variables;
    let comp = [MP.symbol(xn), MP.symbol(yn)];
    const polys = this.trees.map((t) => expand(t, this.variables));
    const usable = (p) => p && [...p.values()].every((t) => Object.entries(t.e).every(([k, v]) => (k === xn || k === yn) && v > 0));
    if (!polys.every(usable)) refuse("not-polynomial", "exact periodic points need a polynomial map with rational coefficients; periodicPoints(n, region) finds proven ones for any formula");
    const substituteBoth = (p, into) => {
      let out = new Map();
      for (const t of p.values()) {
        let term = MP.constant(t.c);
        term = MP.mul(term, MP.pow(into[0], t.e[xn] ?? 0));
        term = MP.mul(term, MP.pow(into[1], t.e[yn] ?? 0));
        out = MP.add(out, term);
      }
      return out;
    };
    for (let i = 0; i < n; i++) {
      comp = polys.map((p) => substituteBoth(p, comp));
      if (comp.some((p) => p.size > 4000)) refuse("too-large", "Fⁿ is too large to expand");
    }
    const eqs = [MP.add(comp[0], MP.symbol(xn), -1), MP.add(comp[1], MP.symbol(yn), -1)];
    // one equation of first degree in one variable, with a constant coefficient there
    let plan = null;
    for (const e of [1, 0]) {
      for (const v of [yn, xn]) {
        const cs = MP.coefficients(eqs[e], v);
        if (plan || !cs || cs.length !== 2) continue;
        const lead = MP.constantOf(cs[1]);
        if (lead === null || lead.numerator === 0n) continue;
        plan = { e, v, other: v === xn ? yn : xn, solved: MP.mul(MP.constant(ONE.divide(lead).negate()), cs[0]) };
      }
    }
    if (!plan) refuse("not-eliminable", `neither equation of a point of period ${n} is of first degree in a variable, so it cannot be solved exactly here; periodicPoints(n, region) finds proven enclosures`);
    // v = solved(other); put it in the remaining equation
    const into = plan.v === xn ? [plan.solved, MP.symbol(yn)] : [MP.symbol(xn), plan.solved];
    const rest = substituteBoth(eqs[1 - plan.e], into);
    const toUni = (p) => {
      const cs = MP.coefficients(p, plan.other) ?? [];
      const out = cs.map((c) => MP.constantOf(c));
      return out.length ? out : [ZERO];
    };
    const uni = P.trim(toUni(rest));
    if (P.deg(uni) < 1) refuse("degenerate", "every point on a curve has this period");
    if (P.deg(uni) > 64) refuse("too-large", "the eliminated equation has too high a degree to solve exactly");
    const solvedUni = toUni(plan.solved);
    const evalAlg = (coeffs, a) => coeffs.reduceRight((acc, c) => aAdd(aMul(acc, a), fromRational(c)), fromRational(ZERO));
    const evalMap = (p, ax, ay) => {
      let out = fromRational(ZERO);
      for (const t of p.values()) {
        let term = fromRational(t.c);
        for (let i = 0; i < (t.e[xn] ?? 0); i++) term = aMul(term, ax);
        for (let i = 0; i < (t.e[yn] ?? 0); i++) term = aMul(term, ay);
        out = aAdd(out, term);
      }
      return out;
    };
    const out = [];
    for (const { value } of algebraicRoots(P.integers(uni))) {
      const second = evalAlg(solvedUni, value);
      const [ax, ay] = plan.v === xn ? [second, value] : [value, second];
      // leave out the points of smaller period
      let cx = ax;
      let cy = ay;
      let least = n;
      for (let d = 1; d < n; d++) {
        [cx, cy] = [evalMap(polys[0], cx, cy), evalMap(polys[1], cx, cy)];
        if (n % d === 0 && aEquals(cx, ax) && aEquals(cy, ay)) {
          least = d;
          break;
        }
      }
      if (least !== n) continue;
      const box = [ax, ay].map((a) => {
        const r = a.rational;
        if (r instanceof Rational) return F.ofQ(r);
        const [lo, hi] = a.interval(80);
        return [F.ofQ(lo)[0], F.ofQ(hi)[1]];
      });
      out.push({ x: ax, y: ay, text: [algebraicFormat(ax), algebraicFormat(ay)], period: n, ...this.#describe(box, n) });
    }
    return out;
  }
  /**
   * Lyapunov exponents along the orbit of p0, estimated in doubles by
   * re-orthogonalising the tangent vectors at every step: { exponents, kind:
   * "estimated", sum }. When the Jacobian determinant is constant the sum of
   * the two exponents is log|det| exactly, and `sum` is an enclosure of it.
   */
  lyapunov(p0, n = 10000, { skip = 100 } = {}) {
    const { jac } = this.#tape();
    let p = p0.map((x) => (typeof x === "number" ? x : toQ(x).toNumber()));
    for (let i = 0; i < skip; i++) p = this.float(p);
    let u = [1, 0];
    let w = [0, 1];
    const sums = [0, 0];
    for (let i = 0; i < n; i++) {
      const j = F.values(jac, p.map(F.pt)).map(F.mid);
      const apply = (v) => [j[0] * v[0] + j[1] * v[1], j[2] * v[0] + j[3] * v[1]];
      u = apply(u);
      w = apply(w);
      const lu = Math.hypot(...u);
      u = u.map((x) => x / lu);
      const d = u[0] * w[0] + u[1] * w[1];
      w = [w[0] - d * u[0], w[1] - d * u[1]];
      const lw = Math.hypot(...w);
      w = w.map((x) => x / lw);
      sums[0] += Math.log(lu);
      sums[1] += Math.log(lw);
      p = this.float(p);
      if (!p.every(Number.isFinite)) refuse("escapes", "the orbit in doubles left every bounded range");
    }
    let sum = null;
    if (this.determinant !== null && this.determinant.numerator !== 0n) sum = logAbs(new Box(this.determinant)).toRationalInterval();
    return { exponents: sums.map((s) => s / n), kind: "estimated", steps: n, sum };
  }
}

/** A map of the plane from two formulas. */
export const planeMap = (field, options) => new PlaneMap(field, options);
/** The Hénon map (x, y) → (1 − a·x² + y, b·x). */
export const henon = (a = "1.4", b = "0.3") => new PlaneMap(["1 - a*x^2 + y", "b*x"], { params: { a, b }, name: "Hénon" });
/**
 * Chirikov's standard map, lifted to the plane: y → y + (k/2π)·sin(2πx),
 * then x → x + y. Positions are in turns; take them modulo 1 to draw it on
 * the torus.
 */
export const standardMap = (k = 1) => new PlaneMap(["x + y + k/(2*pi)*sin(2*pi*x)", "y + k/(2*pi)*sin(2*pi*x)"], { params: { k }, name: "standard map" });

// ---- automorphisms of the torus ------------------------------------------------------------

const mod1 = (q) => {
  const f = q.numerator % q.denominator;
  return new Rational(f < 0n ? f + q.denominator : f, q.denominator);
};

export class TorusMap {
  /** @param matrix a 2 × 2 matrix of whole numbers, acting on the torus R²/Z² */
  constructor(matrix = [[2, 1], [1, 1]], options = {}) {
    this.matrix = matrix.map((row) => row.map((x) => BigInt(x)));
    if (this.matrix.length !== 2 || this.matrix.some((r) => r.length !== 2)) refuse("bad-map", "a torus map is a 2 × 2 matrix of whole numbers");
    const [[a, b], [c, d]] = this.matrix;
    this.determinant = a * d - b * c;
    this.trace = a + d;
    this.name = options.name ?? null;
    /** true when the map is invertible on the torus */
    this.invertible = this.determinant === 1n || this.determinant === -1n;
    // λ² − tr·λ + det
    const disc = this.trace * this.trace - 4n * this.determinant;
    const onCircle = disc <= 0n ? this.determinant === 1n || disc === 0n : [1n, -1n].some((s) => 1n - s * this.trace + this.determinant === 0n);
    /** true when no eigenvalue has modulus 1 */
    this.hyperbolic = !onCircle;
  }
  toString() {
    const term = (k, v) => (k === 0n ? "" : `${k === 1n ? "" : k === -1n ? "-" : k}${v}`);
    const row = ([p, q]) => [term(p, "x"), term(q, "y")].filter(Boolean).join(" + ").replace(/\+ -/g, "- ") || "0";
    return `(x, y) → (${row(this.matrix[0])}, ${row(this.matrix[1])}) mod 1`;
  }
  /** the image of a point, exactly */
  at(p) {
    const [x, y] = p.map(toQ);
    const [[a, b], [c, d]] = this.matrix;
    return [mod1(x.multiply(new Rational(a)).add(y.multiply(new Rational(b)))), mod1(x.multiply(new Rational(c)).add(y.multiply(new Rational(d))))];
  }
  orbit(p0, n) {
    const out = [p0.map((x) => mod1(toQ(x)))];
    for (let i = 0; i < n; i++) out.push(this.at(out[i]));
    return out;
  }
  floatOrbit(p0, n) {
    const m = this.matrix.map((r) => r.map(Number));
    const frac = (x) => x - Math.floor(x);
    const out = [p0.map((x) => (typeof x === "number" ? x : toQ(x).toNumber()))];
    for (let i = 0; i < n; i++) {
      const [x, y] = out[i];
      out.push([frac(m[0][0] * x + m[0][1] * y), frac(m[1][0] * x + m[1][1] * y)]);
    }
    return out;
  }
  /**
   * A rational point is eventually periodic, and periodic when the map is
   * invertible: { period, preperiod, cycle }. The orbit lives on the grid of
   * points with the same denominator.
   */
  period(p0, { limit = 1000000 } = {}) {
    let p = p0.map((x) => mod1(toQ(x)));
    const seen = new Map();
    const orbit = [];
    for (;;) {
      const key = `${p[0]},${p[1]}`;
      if (seen.has(key)) {
        const preperiod = seen.get(key);
        return { period: orbit.length - preperiod, preperiod, cycle: orbit.slice(preperiod) };
      }
      if (orbit.length >= limit) refuse("too-long", `no repeat within ${limit} steps`);
      seen.set(key, orbit.length);
      orbit.push(p);
      p = this.at(p);
    }
  }
  /** the periods of all points i/q, j/q: a Map from period to how many points have it (invertible maps) */
  periods(q) {
    if (!Number.isInteger(q) || q < 1 || q > 2000) refuse("bad-denominator", "the denominator is a whole number from 1 to 2000");
    if (!this.invertible) refuse("not-invertible", "periods(q) needs a map of determinant ±1, where every rational point is periodic");
    const m = this.matrix.map((r) => r.map((x) => Number(((x % BigInt(q)) + BigInt(q)) % BigInt(q))));
    const seen = new Uint8Array(q * q);
    const out = new Map();
    for (let s = 0; s < q * q; s++) {
      if (seen[s]) continue;
      let len = 0;
      let cur = s;
      do {
        seen[cur] = 1;
        const x = Math.floor(cur / q);
        const y = cur % q;
        cur = ((m[0][0] * x + m[0][1] * y) % q) * q + ((m[1][0] * x + m[1][1] * y) % q);
        len++;
      } while (cur !== s);
      out.set(len, (out.get(len) ?? 0) + len);
    }
    return new Map([...out].sort((a, b) => a[0] - b[0]));
  }
  /** the number of fixed points of the n-th iterate, |det(Aⁿ − I)|, for a hyperbolic map */
  fixedPointCount(n) {
    let [[a, b], [c, d]] = [[1n, 0n], [0n, 1n]];
    const [[p, q], [r, s]] = this.matrix;
    for (let i = 0; i < n; i++) [a, b, c, d] = [a * p + b * r, a * q + b * s, c * p + d * r, c * q + d * s];
    const det = (a - 1n) * (d - 1n) - b * c;
    return det < 0n ? -det : det;
  }
  /**
   * The eigenvalues as exact algebraic numbers with the directions they
   * belong to: [{ value, text, direction, expanding }], direction being the
   * slope dy/dx as an algebraic number (null for a vertical direction).
   */
  eigen() {
    const [[a, b], [, d]] = this.matrix;
    const found = algebraicRoots([this.determinant, -this.trace, 1n]);
    if (found.length < 2 && !(found.length === 1 && found[0].multiplicity === 2)) refuse("not-real", "the eigenvalues are not real: the map is a rotation of finite order or has no invariant directions");
    return found.map(({ value }) => {
      // (A − λ)v = 0 with v = (b, λ − a), or (λ − d, c) when b = 0
      let slope = null;
      if (b !== 0n) slope = aMul(aAdd(value, fromRational(new Rational(-a))), fromRational(new Rational(1n, b)));
      else if (aEquals(value, fromRational(new Rational(d)))) slope = null;
      else slope = aMul(fromRational(new Rational(this.matrix[1][0])), aInv(aAdd(value, fromRational(new Rational(-d)))));
      const [lo, hi] = value.rational instanceof Rational ? [value.rational, value.rational] : value.interval(64);
      const absLo = lo.numerator >= 0n ? lo : hi.negate();
      return { value, text: algebraicFormat(value), direction: slope, directionText: slope ? algebraicFormat(slope) : "vertical", expanding: absLo.compareTo(ONE) > 0 };
    });
  }
  /**
   * The topological entropy, which is also the larger Lyapunov exponent: the
   * logarithm of the larger eigenvalue in size. { text, value } with value an
   * enclosure.
   */
  entropy(bits = 96) {
    const big = this.eigen().find((e) => e.expanding);
    if (!big) return { text: "0", value: new RationalInterval(ZERO, ZERO) };
    const [lo, hi] = big.value.rational instanceof Rational ? [big.value.rational, big.value.rational] : big.value.interval(bits);
    const value = logAbs(new Box(lo, hi), bits).toRationalInterval();
    const neg = hi.numerator < 0n;
    return { text: `log(${neg ? "-(" : ""}${big.text}${neg ? ")" : ""})`, value };
  }
  /**
   * What doubles do to a periodic orbit: the first step at which the float
   * orbit of a rational point is farther than `tolerance` from the true one.
   */
  floatPartsAt(p0, { limit = 2000, tolerance = 1e-3 } = {}) {
    const truth = this.orbit(p0, limit);
    const floats = this.floatOrbit(p0, limit);
    const gap = (a, b) => {
      const d = Math.abs(a - b);
      return Math.min(d, 1 - d);
    };
    for (let i = 0; i <= limit; i++) {
      if (Math.max(gap(floats[i][0], truth[i][0].toNumber()), gap(floats[i][1], truth[i][1].toNumber())) > tolerance) return i;
    }
    return null;
  }
}
/** Arnold's cat map (x, y) → (2x + y, x + y) mod 1. */
export const catMap = () => new TorusMap([[2, 1], [1, 1]], { name: "cat map" });
/** A torus map from a 2 × 2 matrix of whole numbers. */
export const torusMap = (matrix, options) => new TorusMap(matrix, options);

// ---- the horseshoe -----------------------------------------------------------------------------

const THIRD = Q(1, 3);
const TWO_THIRDS = Q(2, 3);
const THREE = Q(3);

/**
 * Smale's horseshoe in its piecewise linear form on the unit square: the
 * square is squeezed to a third of its width, stretched to three times its
 * height and folded back over itself, so that
 *   H(x, y) = (x/3, 3y)            for y ≤ 1/3   (symbol 0),
 *   H(x, y) = (1 − x/3, 3 − 3y)    for y ≥ 2/3   (symbol 1),
 * and points with 1/3 < y < 2/3 leave the square. The points that never
 * leave, forwards or backwards, form the product of two Cantor sets, and on
 * it H is the shift on two-sided sequences of 0s and 1s.
 */
export class Horseshoe {
  /** the image of a point of the square, or null when it leaves */
  at(p) {
    const [x, y] = p.map(toQ);
    if (y.compareTo(THIRD) <= 0) return [x.divide(THREE), y.multiply(THREE)];
    if (y.compareTo(TWO_THIRDS) >= 0) return [ONE.subtract(x.divide(THREE)), THREE.subtract(y.multiply(THREE))];
    return null;
  }
  /** the forward itinerary: the symbols of p, H(p), …, up to n of them; shorter when the point leaves */
  itinerary(p, n) {
    let cur = p.map(toQ);
    let out = "";
    for (let i = 0; i < n && cur; i++) {
      if (cur[1].numerator < 0n || cur[1].compareTo(ONE) > 0 || cur[0].numerator < 0n || cur[0].compareTo(ONE) > 0) break;
      if (cur[1].compareTo(THIRD) <= 0) out += "0";
      else if (cur[1].compareTo(TWO_THIRDS) >= 0) out += "1";
      else break;
      cur = this.at(cur);
    }
    return out;
  }
  /** the periodic point whose itinerary is the word repeated for ever, exactly */
  periodicPoint(word) {
    if (!/^[01]+$/.test(word)) refuse("bad-word", "a word is a string of 0s and 1s");
    // each branch is x → αx + β, y → γy + δ; compose along the word and solve for the fixed point
    let [ax, bx, ay, by] = [ONE, ZERO, ONE, ZERO];
    for (const s of word) {
      if (s === "0") [ax, bx, ay, by] = [ax.divide(THREE), bx.divide(THREE), ay.multiply(THREE), by.multiply(THREE)];
      else [ax, bx, ay, by] = [ax.divide(THREE).negate(), ONE.subtract(bx.divide(THREE)), ay.multiply(THREE).negate(), THREE.subtract(by.multiply(THREE))];
    }
    return [bx.divide(ONE.subtract(ax)), by.divide(ONE.subtract(ay))];
  }
  /** every point fixed by Hⁿ (there are 2ⁿ), with its word and least period */
  periodicPoints(n) {
    if (!Number.isInteger(n) || n < 1 || n > 16) refuse("bad-period", "the period is a whole number from 1 to 16");
    const out = [];
    for (let k = 0; k < 2 ** n; k++) {
      const word = k.toString(2).padStart(n, "0");
      let least = n;
      for (const d of properDivisors(n)) {
        if (word.slice(0, d).repeat(n / d) === word) {
          least = d;
          break;
        }
      }
      out.push({ word, point: this.periodicPoint(word), leastPeriod: least });
    }
    return out;
  }
  /** the number of orbits of least period n */
  orbitCount(n) {
    let total = 0;
    for (let d = 1; d <= n; d++) {
      if (n % d) continue;
      // Möbius function of n/d
      let m = n / d;
      let mu = 1;
      for (let p = 2; p * p <= m; p++) {
        if (m % p) continue;
        m /= p;
        if (m % p === 0) mu = 0;
        mu = -mu;
      }
      if (m > 1) mu = -mu;
      total += mu * 2 ** d;
    }
    return total / n;
  }
  /**
   * The points still in the square after n steps forwards form 2ⁿ horizontal
   * strips, and after n steps backwards 2ⁿ vertical strips: the exact
   * intervals of each kind, { horizontal: [[lo, hi]], vertical: [[lo, hi]] }.
   */
  strips(n) {
    if (!Number.isInteger(n) || n < 0 || n > 12) refuse("bad-depth", "the depth is a whole number from 0 to 12");
    let ivs = [[ZERO, ONE]];
    for (let i = 0; i < n; i++) ivs = ivs.flatMap(([a, b]) => [[a.divide(THREE), b.divide(THREE)], [ONE.subtract(b.divide(THREE)), ONE.subtract(a.divide(THREE))]]);
    ivs.sort((p, q) => p[0].compareTo(q[0]));
    return { horizontal: ivs, vertical: ivs };
  }
  /** the dimension of the invariant set: twice that of the middle-thirds Cantor set, 2·log 2/log 3 */
  get dimension() {
    return { text: "2*log(2)/log(3)", value: 2 * Math.log(2) / Math.log(3) };
  }
}
/** Smale's horseshoe on the unit square. */
export const horseshoe = () => new Horseshoe();
