// Flows: systems of ordinary differential equations x′ = f(x), solved with
// guaranteed enclosures fast enough to follow a set of starting points once
// round an attractor.
//
// The solver is a Taylor method in outward-rounded double intervals with
// Lohner's representation of sets. One step from the set c + B·r:
//   1. a rough box R with X + [0, h]·f(R) ⊆ R holds the solution on the whole
//      step (Picard–Lindelöf), and the same for the variational equation
//      V′ = Df(x)·V, V(0) = I;
//   2. the Taylor polynomial of the solution from the centre c, with the next
//      coefficient taken over R as the remainder, gives a tight enclosure z of
//      the centre's image;
//   3. the Taylor polynomial of V over the whole set gives an interval matrix
//      A that contains the derivative of the time-h map at every point of it,
//      so the image of the set lies in z + (A·B)·r (mean value theorem);
//   4. the new frame is the orthogonal factor of A·B, which keeps rotating
//      sets from being wrapped in ever larger boxes.
// On top of it: first returns to a section with their derivative, periodic
// orbits proven by Krawczyk's test on the return map, equilibria proven by the
// same test on the field, and bifurcations of equilibria bracketed in a
// parameter.
import { refuse, toQ, parse, format, differentiate, simplify, compileFloat, substitute } from "./tree.js";
import * as F from "./fi.js";

const num = (q) => ({ type: "num", value: q });
const isZeroNode = (t) => t.type === "num" && t.value.numerator === 0n;
const defaultNames = (n) => (n === 1 ? ["x"] : n === 2 ? ["x", "y"] : n === 3 ? ["x", "y", "z"] : Array.from({ length: n }, (_, i) => `x${i + 1}`));

/** Solution of a small linear system in doubles. */
function solveFloat(M, b) {
  const inv = F.inverseFloat(M);
  return inv ? inv.map((row) => row.reduce((s, x, j) => s + x * b[j], 0)) : null;
}
const modulus = (ev) => F.sqrt(F.add(F.sqr(ev.re), F.sqr(ev.im)));
const showEigen = (ev) => ({ re: F.toRI(ev.re), im: F.toRI(ev.im) });

export class Flow {
  /**
   * @param field the right-hand sides as formulas, one per variable: ["s*(y - x)", "x*(r - z) - y", "x*y - b*z"]
   * @param options { variables: ["x", "y", "z"], params: { s: 10, r: 28, b: "8/3" }, name }
   */
  constructor(field, options = {}) {
    if (!Array.isArray(field) || !field.length) refuse("bad-field", "a flow needs a list of formulas, one for each variable");
    this.n = field.length;
    this.variables = options.variables ?? defaultNames(this.n);
    if (this.variables.length !== this.n) refuse("bad-field", `${this.n} formulas need ${this.n} variables`);
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
  }

  /** the same family with other parameter values */
  with(params) {
    return new Flow(this.family, { variables: this.variables, params: { ...this.params, ...params }, name: this.name });
  }
  toString() {
    const ps = Object.entries(this.params).map(([k, v]) => `${k} = ${v}`).join(", ");
    return this.variables.map((v, i) => `${v}′ = ${this.text[i]}`).join(", ") + (ps ? ` (${ps})` : "");
  }
  /** the Jacobian matrix as formulas */
  jacobian() {
    return this._J.map((row) => row.map(format));
  }
  /** the field at a point, in doubles */
  velocity(p) {
    const env = Object.fromEntries(this.variables.map((v, i) => [v, p[i]]));
    return this._float.map((f) => f(env));
  }

  #tapes() {
    if (this._tapes) return this._tapes;
    const { n, variables: vars } = this;
    const vname = (i, j) => `__V${i}_${j}`;
    const vnames = [];
    const vfield = [];
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        vnames.push(vname(i, j));
        let sum = null;
        for (let k = 0; k < n; k++) {
          if (isZeroNode(this._J[i][k])) continue;
          const term = { type: "op", op: "*", args: [this._J[i][k], { type: "var", name: vname(k, j) }] };
          sum = sum ? { type: "op", op: "+", args: [sum, term] } : term;
        }
        vfield.push(sum ?? num(toQ(0)));
      }
    }
    this._tapes = {
      y: F.compile(this.trees, vars),
      ext: F.compile([...this.trees, ...vfield], [...vars, ...vnames]),
      jac: F.compile(this._J.flat(), vars),
    };
    return this._tapes;
  }
  #field(box) {
    return F.field(this.#tapes().y, box);
  }
  #jacobian(box) {
    const flat = F.values(this.#tapes().jac, box);
    return Array.from({ length: this.n }, (_, i) => flat.slice(i * this.n, (i + 1) * this.n));
  }
  #point(p, what = "a point") {
    const list = Array.isArray(p) ? p : this.variables.map((v) => p?.[v]);
    if (list.length !== this.n || list.some((x) => x === undefined)) refuse("bad-point", `${what} needs a value for each of ${this.variables.join(", ")}`);
    return list.map(F.ofExact);
  }
  #setOf(X) {
    const c = F.midV(X);
    return { c, B: F.identity(this.n), r: X.map((x, i) => F.sub(x, F.pt(c[i]))), X };
  }

  // ---- one validated step -------------------------------------------------------------------

  #rough(Z0, h) {
    const { ext } = this.#tapes();
    const hI = [0, h];
    const image = (Rv) => {
      const f = F.field(ext, Rv);
      return Z0.map((z, i) => F.hull(z, F.add(z, F.mul(hI, f[i]))));
    };
    try {
      // a few Picard sweeps give the size of the motion; then pad and test
      let N = Z0;
      for (let it = 0; it < 3; it++) N = image(N);
      if (!N.every(F.finite)) return null;
      const size = N.reduce((s, x, i) => Math.max(s, F.width(x) - F.width(Z0[i])), 0);
      let by = 0.05;
      for (let it = 0; it < 8; it++, by *= 2) {
        const R = N.map((x) => {
          const pad = by * (F.width(x) + size) + 1e-15 * (1 + F.mag(x));
          return [x[0] - pad, x[1] + pad];
        });
        const next = image(R);
        if (!next.every(F.finite)) return null;
        if (next.every((x, i) => F.within(x, R[i]))) return next;
        N = next.map((x, i) => F.hull(x, R[i]));
      }
    } catch {
      return null;
    }
    return null;
  }
  #step(set, h, K, pc) {
    const { n } = this;
    const { y: ty, ext } = this.#tapes();
    const Z0 = [...set.X, ...F.identity(n).flat().map(F.pt)];
    const R = this.#rough(Z0, h);
    if (!R) return null;
    try {
      const pr = F.series(ty, R.slice(0, n), K + 1).map((c) => c[K + 1]);
      const jc = F.series(ext, Z0, K);
      const jr = F.series(ext, R, K + 1).map((c) => c[K + 1]);
      return { h, K, R, pc, pr, jc, jr, set };
    } catch {
      return null;
    }
  }
  /** the step's enclosures at time s (an interval inside [0, h]) after its start */
  #dense(st, s) {
    const { n } = this;
    const { K, pc, pr, jc, jr, set } = st;
    const horner = (coef, rem) => {
      let acc = rem;
      for (let k = K; k >= 0; k--) acc = F.add(coef[k], F.mul(s, acc));
      return acc;
    };
    const z = pc.map((c, i) => horner(c, pr[i]));
    const direct = Array.from({ length: n }, (_, i) => horner(jc[i], jr[i]));
    const A = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => horner(jc[n + i * n + j], jr[n + i * n + j])));
    const M = F.mmul(A, F.ivM(set.B));
    const Mr = F.mvec(M, set.r);
    const X = z.map((x, i) => {
      const viaSet = F.add(x, Mr[i]);
      return F.meet(viaSet, direct[i]) ?? viaSet;
    });
    if (!X.every(F.finite)) refuse("no-enclosure", "the enclosure of the solution grew without bound");
    return { z, A, M, X };
  }
  #reframe(z, M, r, X) {
    const c = F.midV(z);
    const zc = z.map((x, i) => F.sub(x, F.pt(c[i])));
    const B = F.frame(F.midM(M), r.map(F.mag));
    const Bi = F.inverseEnclosed(B);
    if (!Bi) refuse("degenerate", "could not invert the frame of the set");
    const r1 = F.mvec(F.mmul(Bi, M), r).map((x, i) => F.add(x, F.mvec(Bi, zc)[i]));
    return { c, B, r: r1, X };
  }
  #advance(st, s) {
    const d = this.#dense(st, s);
    return { set: this.#reframe(d.z, d.M, st.set.r, d.X), A: d.A };
  }
  /** a step from the set of a length the solver chooses (at most hmax): the step data */
  #trial(set, hmax, opts) {
    const K = opts.order;
    const pc = F.series(this.#tapes().y, set.c.map(F.pt), K);
    let h = hmax;
    const tol = 1e-17 * (1 + Math.max(...set.c.map(Math.abs)));
    for (const k of [K - 1, K]) {
      const m = Math.max(...pc.map((c) => F.mag(c[k])));
      if (m > 0) h = Math.min(h, 0.8 * (tol / m) ** (1 / k));
    }
    for (; h > 1e-13; h /= 2) {
      const st = this.#step(set, h, K, pc);
      if (st) return st;
    }
    refuse("no-enclosure", "no enclosure of the solution could be found, even for a very short step: the set of starting points has grown too wide, or the solution blows up");
  }

  // ---- solving ---------------------------------------------------------------------------------

  /**
   * The solution through p0 (a point, or a box of points given as intervals)
   * up to time `until`, enclosed: { final, at(t), steps, complete, reached }.
   * `final` and `at(t)` are RationalIntervals, one per variable, that hold
   * the true solution for every starting point in the box.
   */
  solve(p0, until, options = {}) {
    const opts = { order: 12, maxStep: 0.25, maxSteps: 20000, ...options };
    const T = F.ofExact(until);
    if (!(T[0] > 0)) refuse("bad-time", "the end time is after 0");
    let set = this.#setOf(this.#point(p0));
    let t = F.pt(0);
    const steps = [];
    let hlast = opts.maxStep;
    let complete = true;
    let reason = null;
    const variables = this.variables;
    const dense = (st, s) => this.#dense(st, s).X;
    try {
      for (;;) {
        if (steps.length >= opts.maxSteps) {
          complete = false;
          reason = `stopped after ${opts.maxSteps} steps`;
          break;
        }
        const left = F.sub(T, t);
        const st = this.#trial(set, Math.min(opts.maxStep, 2 * hlast, left[1]), opts);
        hlast = st.h;
        steps.push({ t0: t, st });
        if (st.h >= left[1]) {
          set = this.#advance(st, [Math.max(0, left[0]), left[1]]).set;
          t = T;
          break;
        }
        set = this.#advance(st, F.pt(st.h)).set;
        t = F.add(t, F.pt(st.h));
      }
    } catch (e) {
      if (e.code !== "no-enclosure") throw e;
      complete = false;
      reason = e.message;
    }
    const final = set.X.map(F.toRI);
    return {
      variables,
      final,
      complete,
      reason,
      reached: F.toRI(t),
      steps: steps.length,
      /** the widest of the final enclosures, as a double */
      width: Math.max(...set.X.map(F.width)),
      at(time) {
        const tau = F.ofExact(time);
        const hit = steps.find(({ t0, st }) => tau[0] >= t0[1] && tau[1] - t0[0] <= st.h);
        if (!hit || tau[1] > t[1]) refuse("outside-domain", "that time is outside the solved span");
        const s = F.sub(tau, hit.t0);
        return dense(hit.st, [Math.max(0, s[0]), Math.min(hit.st.h, s[1])]).map(F.toRI);
      },
    };
  }

  /** the orbit as an ordinary program computes it (classical Runge–Kutta): { t, points } */
  floatOrbit(p0, until, { steps = 2000 } = {}) {
    let y = this.#point(p0).map(F.mid);
    const T = typeof until === "number" ? until : toQ(until).toNumber();
    const h = T / steps;
    const f = (p) => this.velocity(p);
    const ts = [0];
    const points = [y];
    for (let k = 0; k < steps; k++) {
      const k1 = f(y);
      const k2 = f(y.map((v, i) => v + (h / 2) * k1[i]));
      const k3 = f(y.map((v, i) => v + (h / 2) * k2[i]));
      const k4 = f(y.map((v, i) => v + h * k3[i]));
      y = y.map((v, i) => v + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
      ts.push((k + 1) * h);
      points.push(y);
    }
    return { t: ts, points };
  }

  // ---- equilibria ------------------------------------------------------------------------------

  /**
   * Every equilibrium in the box { x: [a, b], y: [c, d], … }, each proven to
   * exist and to be the only one in its enclosure, with the eigenvalues of
   * the linearisation there: [{ point, eigenvalues, kind, stability,
   * unstableDimension }]. `complete` on the list says the rest of the box was
   * proven to hold no other.
   */
  equilibria(region, options = {}) {
    const start = this.variables.map((v) => {
      if (!region?.[v]) refuse("bad-region", `the region needs a range for ${v}`);
      return F.ofExact(region[v]);
    });
    const list = F.solveAll((X) => this.#field(X), (X) => this.#jacobian(X), start, { tol: options.tol ?? 1e-9, maxBoxes: options.maxBoxes ?? 20000 });
    const out = list.map((box) => this.#linearisation(box));
    out.sort((a, b) => {
      for (let i = 0; i < this.n; i++) {
        const c = a.point[i].low.compareTo(b.point[i].low);
        if (c) return c;
      }
      return 0;
    });
    out.complete = list.complete;
    return out;
  }
  #linearisation(box) {
    const J = this.#jacobian(box);
    const { n } = this;
    const ev = F.eigenvalues(J);
    let stable = 0;
    let unstable = 0;
    let spiral = false;
    if (ev) {
      for (const e of ev) {
        if (e.re[1] < 0) stable++;
        else if (e.re[0] > 0) unstable++;
        if (e.im[0] > 0 || e.im[1] < 0) spiral = true;
      }
    } else if (n === 2) {
      const tr = F.add(J[0][0], J[1][1]);
      const det = F.sub(F.mul(J[0][0], J[1][1]), F.mul(J[0][1], J[1][0]));
      if (det[1] < 0) [stable, unstable] = [1, 1];
      else if (det[0] > 0 && tr[1] < 0) stable = 2;
      else if (det[0] > 0 && tr[0] > 0) unstable = 2;
    }
    const decided = stable + unstable === n;
    const kind = !decided ? "undecided" : unstable === 0 ? (spiral ? "spiral sink" : "sink") : stable === 0 ? (spiral ? "spiral source" : "source") : spiral ? "saddle-focus" : "saddle";
    return {
      point: box.map(F.toRI),
      eigenvalues: ev ? ev.map(showEigen) : null,
      kind,
      stability: !decided ? "undecided" : unstable === 0 ? "stable" : "unstable",
      unstableDimension: decided ? unstable : null,
      spiral,
    };
  }

  // ---- sections and return maps -------------------------------------------------------------

  #section(section) {
    const j = this.variables.indexOf(section?.variable);
    if (j < 0) refuse("bad-section", `a section is { variable, value, direction }: the plane where one of ${this.variables.join(", ")} takes a value`);
    const dir = section.direction ?? 1;
    const d = dir === 1 || dir === "up" || dir === "increasing" ? 1 : dir === -1 || dir === "down" || dir === "decreasing" ? -1 : 0;
    if (!d) refuse("bad-section", 'the direction is 1 ("up": the variable increasing through the value) or -1 ("down")');
    return { j, v: F.ofExact(section.value ?? 0), d, names: this.variables.filter((_, i) => i !== j) };
  }
  #onSection(sec, start) {
    const list = Array.isArray(start) ? start : sec.names.map((v) => start?.[v]);
    if (list.length !== sec.names.length || list.some((x) => x === undefined)) refuse("bad-point", `a point of the section needs a value for each of ${sec.names.join(", ")}`);
    const X = list.map((x) => (typeof x === "number" ? F.pt(x) : F.ofExact(x)));
    X.splice(sec.j, 0, sec.v);
    return X;
  }
  /** follow the set, which lies on the section, to its next crossing in the section's direction */
  #toSection(set0, D0, t0, sec, opts) {
    const { n } = this;
    const { j, v, d } = sec;
    const sigma = (X) => (d === 1 ? F.sub(X[j], v) : F.sub(v, X[j]));
    let set = set0;
    let D = D0;
    let t = t0;
    let state = "leaving";
    let hlast = opts.maxStep;
    for (let count = 0; count < opts.maxSteps; count++) {
      let st = this.#trial(set, Math.min(opts.maxStep, 2 * hlast), opts);
      let touches = F.hasZero(F.sub(st.R[j], v));
      let sgn = 0;
      while (touches) {
        sgn = F.signOf(this.#field(st.R.slice(0, n))[j]);
        if (sgn !== 0) break;
        // the field may be tangent to the section somewhere on this step: look closer
        if (st.h < 1e-9) refuse("tangent", "the orbit may touch the section without crossing it; choose another section");
        st = this.#trial(set, st.h / 4, opts);
        touches = F.hasZero(F.sub(st.R[j], v));
      }
      hlast = st.h;
      const target = touches && state === "away" && sgn === d && !(sigma(set.X)[0] > 0);
      if (!target) {
        if (touches && state === "leaving" && sgn !== d) refuse("bad-section", "the starting points do not cross the section in the direction asked for");
        if (!touches && state === "leaving") state = "away";
        const next = this.#advance(st, F.pt(st.h));
        set = next.set;
        D = F.mmul(next.A, D);
        t = F.add(t, F.pt(st.h));
        continue;
      }
      // the crossing: first go to a time when every point is still short of the section
      if (!(sigma(set.X)[1] < 0)) refuse("no-enclosure", "the set is too wide to follow through the section; start from a smaller box");
      let lo = 0;
      let hi = st.h;
      if (sigma(this.#dense(st, F.pt(hi)).X)[1] < 0) lo = hi;
      else {
        for (let it = 0; it < 40; it++) {
          const m = (lo + hi) / 2;
          if (sigma(this.#dense(st, F.pt(m)).X)[1] < 0) lo = m;
          else hi = m;
        }
      }
      if (lo > 0) {
        const next = this.#advance(st, F.pt(lo));
        set = next.set;
        D = F.mmul(next.A, D);
        t = F.add(t, F.pt(lo));
      }
      if (lo === st.h) continue;
      // then one step that takes every point across
      const gap = -sigma(set.X)[0];
      const fNow = this.#field(set.X)[j];
      if (F.signOf(fNow) !== d) refuse("tangent", "the orbit may touch the section without crossing it; choose another section");
      const speed = F.div(F.pt(1), fNow);
      let h2 = 1.5 * gap * F.mag(speed) + 1e-12;
      let st2 = null;
      const pc = F.series(this.#tapes().y, set.c.map(F.pt), opts.order);
      for (let tries = 0; tries < 12; tries++, h2 *= 2) {
        st2 = this.#step(set, h2, opts.order, pc);
        if (!st2) break;
        if (sigma(this.#dense(st2, F.pt(h2)).X)[0] > 0) break;
        st2 = null;
      }
      if (!st2) refuse("no-enclosure", "the set is too wide to follow through the section; start from a smaller box");
      lo = 0;
      hi = h2;
      for (let it = 0; it < 40; it++) {
        const m = (lo + hi) / 2;
        if (sigma(this.#dense(st2, F.pt(m)).X)[0] > 0) hi = m;
        else lo = m;
      }
      const over = this.#dense(st2, [0, hi]);
      const W = over.X.map((x, i) => F.meet(x, st2.R[i]) ?? x);
      const fW = this.#field(W);
      if (F.signOf(fW[j]) !== d) refuse("tangent", "the orbit may touch the section without crossing it; choose another section");
      const rho = fW.map((f) => F.div(f, fW[j]));
      // each point moves from where it is to the section along the field: x + (v − x_j)·ρ
      const cj = F.sub(v, F.pt(set.c[j]));
      const z = set.c.map((ci, i) => (i === j ? v : F.add(F.pt(ci), F.mul(cj, rho[i]))));
      const Bj = set.B[j];
      const M = set.B.map((row, i) => row.map((b, k) => (i === j ? [0, 0] : F.sub(F.pt(b), F.scale(rho[i], Bj[k])))));
      const Mr = F.mvec(M, set.r);
      const X = z.map((x, i) => (i === j ? v : F.add(x, Mr[i])));
      const time = F.add(t, F.div(F.sub(v, set.X[j]), fW[j]));
      const AD = F.mmul(over.A, D);
      const Dn = AD.map((row, i) => row.map((a, k) => (i === j ? [0, 0] : F.sub(a, F.mul(rho[i], AD[j][k])))));
      return { set: this.#reframe(z, M, set.r, X), D: Dn, t: time };
    }
    refuse("no-return", `no return to the section within ${opts.maxSteps} steps`);
  }
  #returns(sec, X0, opts) {
    let cur = { set: this.#setOf(X0), D: F.ivM(F.identity(this.n)), t: F.pt(0) };
    for (let k = 0; k < opts.returns; k++) cur = this.#toSection(cur.set, cur.D, cur.t, sec, opts);
    const keep = this.variables.map((_, i) => i).filter((i) => i !== sec.j);
    return {
      point: keep.map((i) => cur.set.X[i]),
      derivative: keep.map((i) => keep.map((k) => cur.D[i][k])),
      time: cur.t,
    };
  }
  /**
   * The first-return map of a section. `section` is { variable, value,
   * direction }: the plane where that variable takes that value, crossed with
   * the variable increasing (direction 1) or decreasing (−1). `start` gives
   * the other variables, as numbers or intervals. Returns enclosures that
   * hold for every starting point: { point, time, derivative }, the
   * derivative being the matrix of the return map in the section's
   * coordinates. `returns: k` gives the k-th return.
   */
  returnMap(section, start, options = {}) {
    const opts = { order: 12, maxStep: 0.25, maxSteps: 20000, returns: 1, ...options };
    const sec = this.#section(section);
    const out = this.#returns(sec, this.#onSection(sec, start), opts);
    return {
      variables: sec.names,
      point: out.point.map(F.toRI),
      time: F.toRI(out.time),
      derivative: out.derivative.map((row) => row.map(F.toRI)),
    };
  }
  /**
   * Look for a periodic orbit through the section near `guess` (the section's
   * other variables, as doubles or exact numbers) and prove it: Newton's
   * method on the return map finds a candidate, and Krawczyk's test on a
   * small box round it shows that the return map has exactly one fixed point
   * there. { proven, point, period, multipliers, stability, radius }.
   */
  periodicOrbit(section, guess, options = {}) {
    const opts = { order: 12, maxStep: 0.25, maxSteps: 20000, returns: 1, ...options };
    const sec = this.#section(section);
    const m = sec.names.length;
    const list = Array.isArray(guess) ? guess : sec.names.map((v) => guess?.[v]);
    if (list.length !== m || list.some((x) => x === undefined)) refuse("bad-point", `the guess needs a value for each of ${sec.names.join(", ")}`);
    let x = list.map((g) => (typeof g === "number" ? g : F.mid(F.ofExact(g))));
    const call = (X) => {
      const full = X.slice();
      full.splice(sec.j, 0, sec.v);
      return this.#returns(sec, full, opts);
    };
    const failed = (reason) => ({ proven: false, variables: sec.names, point: x, reason });
    // Newton's method in doubles, on enclosures of the return of single points
    let thin;
    try {
      for (let it = 0; it < 30; it++) {
        thin = call(x.map(F.pt));
        const g = thin.point.map((p, i) => F.mid(p) - x[i]);
        if (Math.max(...g.map(Math.abs)) < 1e-13 * (1 + Math.max(...x.map(Math.abs)))) break;
        const J = thin.derivative.map((row, i) => row.map((a, k) => F.mid(a) - (i === k ? 1 : 0)));
        const dx = solveFloat(J, g);
        if (!dx) return failed("the return map's derivative minus the identity could not be inverted");
        x = x.map((xi, i) => xi - dx[i]);
        if (!x.every(Number.isFinite)) return failed("Newton's method on the return map did not settle");
        thin = null;
      }
      if (!thin) return failed("Newton's method on the return map did not settle");
    } catch (e) {
      return failed(e.message);
    }
    const scaleOf = 1 + Math.max(...x.map(Math.abs));
    const radii = options.radius !== undefined ? [toQ(options.radius).toNumber()] : [1e-9 * scaleOf, 1e-7 * scaleOf, 1e-11 * scaleOf];
    for (const w of radii) {
      let wide;
      try {
        wide = call(x.map((xi) => [xi - w, xi + w]));
      } catch {
        continue;
      }
      const Cm = F.inverseFloat(wide.derivative.map((row, i) => row.map((a, k) => F.mid(a) - (i === k ? 1 : 0))));
      if (!Cm) continue;
      const Ci = F.ivM(Cm);
      // K = x − C·(P(x) − x) + (I − C·(DP(X) − I))·(X − x)
      const Gx = thin.point.map((p, i) => F.sub(p, F.pt(x[i])));
      const DG = wide.derivative.map((row, i) => row.map((a, k) => (i === k ? F.sub(a, F.pt(1)) : a)));
      const E = F.mmul(Ci, DG).map((row, i) => row.map((a, k) => F.sub(F.pt(i === k ? 1 : 0), a)));
      const CG = F.mvec(Ci, Gx);
      const Ed = F.mvec(E, x.map(() => [-w, w]));
      const Kb = x.map((xi, i) => F.add(F.sub(F.pt(xi), CG[i]), Ed[i]));
      if (!Kb.every((k, i) => F.inside(k, [x[i] - w, x[i] + w]))) continue;
      // proven: one more pass over the smaller box for the period and the multipliers
      let fin = wide;
      let box = Kb;
      try {
        fin = call(Kb);
      } catch {
        box = x.map((xi) => [xi - w, xi + w]);
      }
      const ev = F.eigenvalues(fin.derivative);
      let stability = "undecided";
      if (ev) {
        const mods = ev.map(modulus);
        if (mods.every((a) => a[1] < 1)) stability = "attracting";
        else if (mods.every((a) => a[0] > 1)) stability = "repelling";
        else if (mods.every((a) => a[1] < 1 || a[0] > 1)) stability = "saddle";
      }
      const full = box.slice();
      full.splice(sec.j, 0, sec.v);
      return {
        proven: true,
        variables: this.variables,
        point: full.map(F.toRI),
        period: F.toRI(fin.time),
        returns: opts.returns,
        multipliers: ev ? ev.map(showEigen) : null,
        stability,
        radius: w,
      };
    }
    return failed("a candidate was found but Krawczyk's test did not prove it; the orbit may be close to neutral");
  }

  // ---- bifurcations of equilibria ---------------------------------------------------------

  /**
   * Bifurcations of equilibria in the box as the parameter runs over
   * [from, to]: the range is sampled, the equilibria are proven and counted
   * at each sample, and wherever their number or stability differs between
   * two samples the bracket is halved until it is `tol` wide or the proofs
   * stop succeeding. [{ parameter: [lo, hi], kind, before, after }] with kind
   * "saddle-node", "pitchfork", "transcritical", "Hopf" or "stability change".
   */
  bifurcations(param, [from, to], region, options = {}) {
    if (!(param in this.params)) refuse("bad-parameter", `${param} is not a parameter of this flow`);
    const { samples = 24, tol = "1/1000000", maxBoxes = 20000 } = options;
    const a = toQ(from);
    const b = toQ(to);
    const width = toQ(tol);
    const memo = new Map();
    const look = (mu) => {
      const key = mu.toString();
      if (!memo.has(key)) {
        let eq = null;
        try {
          eq = this.with({ [param]: mu }).equilibria(region, { maxBoxes });
          if (!eq.complete || eq.some((e) => e.kind === "undecided")) eq = null;
        } catch {
          eq = null;
        }
        memo.set(key, eq);
      }
      return memo.get(key);
    };
    const centre = (e) => e.point.map((p) => p.low.add(p.high).toNumber() / 2);
    // pair the equilibria of two parameter values by position
    const pairs = (A, B) => {
      const all = [];
      A.forEach((e, i) => {
        const c = centre(e);
        B.forEach((f, k) => all.push([Math.hypot(...centre(f).map((x, q) => x - c[q])), i, k]));
      });
      all.sort((p, q) => p[0] - q[0]);
      const usedA = new Set();
      const usedB = new Set();
      const out = [];
      for (const [, i, k] of all) {
        if (usedA.has(i) || usedB.has(k)) continue;
        usedA.add(i);
        usedB.add(k);
        out.push([A[i], B[k]]);
      }
      return out;
    };
    const tag = (e) => `${e.unstableDimension}${e.spiral ? "s" : ""}`;
    const differ = (A, B) => A.length !== B.length || pairs(A, B).some(([p, q]) => p.unstableDimension !== q.unstableDimension);
    const classify = (A, B) => {
      const [small, big] = A.length <= B.length ? [A, B] : [B, A];
      const changed = pairs(small, big).filter(([p, q]) => p.unstableDimension !== q.unstableDimension);
      if (big.length - small.length === 2) return changed.length ? "pitchfork" : "saddle-node";
      if (big.length !== small.length) return "equilibria appear";
      if (changed.length && changed.every(([p, q]) => p.spiral && q.spiral && Math.abs(p.unstableDimension - q.unstableDimension) === 2)) return "Hopf";
      if (changed.length === 2 && changed.every(([p, q]) => Math.abs(p.unstableDimension - q.unstableDimension) === 1)) return "transcritical";
      return "stability change";
    };
    const found = [];
    const narrow = (lo, hi, A, B) => {
      if (hi.subtract(lo).compareTo(width) <= 0) return found.push({ lo, hi, A, B });
      const mid = lo.add(hi).divide(toQ(2));
      const Mq = look(mid);
      if (!Mq) return found.push({ lo, hi, A, B, stalled: true });
      const first = differ(A, Mq);
      const second = differ(Mq, B);
      if (first) narrow(lo, mid, A, Mq);
      if (second) narrow(mid, hi, Mq, B);
    };
    let prev = null;
    for (let i = 0; i <= samples; i++) {
      // samples off the round values, where bifurcations of textbook systems tend to sit
      const mu = a.add(b.subtract(a).multiply(toQ(BigInt(1000 * i + (i === 0 ? 0 : i === samples ? 0 : 137)))).divide(toQ(BigInt(1000 * samples))));
      const eq = look(mu);
      if (!eq) continue;
      if (prev && differ(prev.eq, eq)) narrow(prev.mu, mu, prev.eq, eq);
      prev = { mu, eq };
    }
    return found.map(({ lo, hi, A, B, stalled }) => ({
      parameter: [lo, hi],
      kind: classify(A, B),
      before: A.map((e) => e.kind),
      after: B.map((e) => e.kind),
      /** true when the bracket could not be halved to the tolerance because the proofs stopped succeeding near the bifurcation */
      stalled: Boolean(stalled),
      tags: [A.map(tag).join(","), B.map(tag).join(",")],
    }));
  }
}

/** A flow from formulas: flow(["y", "-x"], { variables: ["x", "y"] }). */
export const flow = (field, options) => new Flow(field, options);
/** The Lorenz system. */
export const lorenz = (s = 10, r = 28, b = "8/3") => new Flow(["s*(y - x)", "x*(r - z) - y", "x*y - b*z"], { variables: ["x", "y", "z"], params: { s, r, b }, name: "Lorenz" });
/** The Rössler system. */
export const rossler = (a = "0.2", b = "0.2", c = "5.7") => new Flow(["-y - z", "x + a*y", "b + z*(x - c)"], { variables: ["x", "y", "z"], params: { a, b, c }, name: "Rössler" });
/** Van der Pol's oscillator x″ − μ(1 − x²)x′ + x = 0, as x′ = y, y′ = μ(1 − x²)y − x. */
export const vanDerPol = (mu = 1) => new Flow(["y", "mu*(1 - x^2)*y - x"], { variables: ["x", "y"], params: { mu }, name: "Van der Pol" });
/** The pendulum θ″ + c·θ′ + sin θ = 0, as x′ = y, y′ = −sin x − c·y. */
export const pendulum = (damping = 0) => new Flow(["y", "-sin(x) - c*y"], { variables: ["x", "y"], params: { c: damping }, name: "pendulum" });
