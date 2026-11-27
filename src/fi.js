// Interval arithmetic in doubles, rounded outward, and what the plane maps
// and flows build on it: formulas compiled to a tape, Taylor coefficients of
// the solution of an ordinary differential equation, small interval linear
// algebra, and a Krawczyk test that proves a system of equations has exactly
// one solution in a box.
//
// An interval is a pair [lo, hi] of doubles. After a correctly rounded
// operation the error is at most half a unit in the last place, which is at
// most |x|·2^-53; each end is moved out by |x|·2^-52 plus the smallest
// subnormal, which covers it with room for the rounding of the move itself.
// Doubles are exact rationals, so every enclosure can be handed back as a
// RationalInterval without loss.
import { Rational, RationalInterval } from "@ratmath/core";
import { refuse, toQ } from "./tree.js";

const EPS = 2 ** -52;
const TINY = Number.MIN_VALUE;
export const dn = (x) => x - (Math.abs(x) * EPS + TINY);
export const up = (x) => x + (Math.abs(x) * EPS + TINY);

export const pt = (x) => [x, x];
const ZERO = [0, 0];
const isZero = (a) => a[0] === 0 && a[1] === 0;
export const add = (a, b) => (isZero(a) ? b : isZero(b) ? a : [dn(a[0] + b[0]), up(a[1] + b[1])]);
export const sub = (a, b) => (isZero(b) ? a : [dn(a[0] - b[1]), up(a[1] - b[0])]);
export const neg = (a) => [-a[1], -a[0]];
export function mul(a, b) {
  if (isZero(a) || isZero(b)) return ZERO;
  const p = a[0] * b[0];
  const q = a[0] * b[1];
  const r = a[1] * b[0];
  const s = a[1] * b[1];
  return [dn(Math.min(p, q, r, s)), up(Math.max(p, q, r, s))];
}
export function sqr(a) {
  if (isZero(a)) return ZERO;
  const lo = a[0] > 0 ? a[0] : a[1] < 0 ? -a[1] : 0;
  const hi = Math.max(Math.abs(a[0]), Math.abs(a[1]));
  return [lo === 0 ? 0 : dn(lo * lo), up(hi * hi)];
}
export function div(a, b) {
  if (b[0] <= 0 && b[1] >= 0) refuse("division", "division by an interval that contains zero");
  return mul(a, [dn(1 / b[1]), up(1 / b[0])]);
}
export const scale = (a, x) => mul(a, [x, x]);
export const hull = (a, b) => [Math.min(a[0], b[0]), Math.max(a[1], b[1])];
export const meet = (a, b) => {
  const lo = Math.max(a[0], b[0]);
  const hi = Math.min(a[1], b[1]);
  return lo <= hi ? [lo, hi] : null;
};
export const mid = (a) => (a[0] === a[1] ? a[0] : a[0] / 2 + a[1] / 2);
export const width = (a) => a[1] - a[0];
export const mag = (a) => Math.max(Math.abs(a[0]), Math.abs(a[1]));
export const hasZero = (a) => a[0] <= 0 && a[1] >= 0;
export const inside = (a, b) => a[0] > b[0] && a[1] < b[1];
export const within = (a, b) => a[0] >= b[0] && a[1] <= b[1];
export const finite = (a) => Number.isFinite(a[0]) && Number.isFinite(a[1]);
/** −1 or 1 when every number in the interval has that sign, else 0. */
export const signOf = (a) => (a[0] > 0 ? 1 : a[1] < 0 ? -1 : 0);
export function sqrt(a) {
  if (a[0] < 0) refuse("domain", "square root of an interval that reaches below zero");
  return [Math.max(0, dn(Math.sqrt(a[0]))), up(Math.sqrt(a[1]))];
}

// ---- doubles and rationals ------------------------------------------------------------------

const buf = new DataView(new ArrayBuffer(8));
/** The exact value of a double. */
export function qOfFloat(x) {
  if (!Number.isFinite(x)) refuse("not-finite", "an enclosure grew without bound");
  if (x === 0) return new Rational(0n);
  buf.setFloat64(0, x);
  const hi = buf.getUint32(0);
  const lo = buf.getUint32(4);
  const e = (hi >>> 20) & 0x7ff;
  let m = (BigInt(hi & 0xfffff) << 32n) | BigInt(lo);
  let shift;
  if (e === 0) shift = -1074;
  else {
    m |= 1n << 52n;
    shift = e - 1075;
  }
  if (hi >>> 31) m = -m;
  return shift >= 0 ? new Rational(m << BigInt(shift)) : new Rational(m, 1n << BigInt(-shift));
}
/** An interval of doubles around an exact number or a rational interval. */
export function ofExact(v) {
  if (Array.isArray(v) && v.length === 2 && typeof v[0] === "number" && typeof v[1] === "number") return v;
  if (v instanceof RationalInterval) return [ofQ(v.low)[0], ofQ(v.high)[1]];
  if (Array.isArray(v) && v.length === 2) return [ofQ(toQ(v[0]))[0], ofQ(toQ(v[1]))[1]];
  return ofQ(toQ(v));
}
export function ofQ(q) {
  const x = q.toNumber();
  if (qOfFloat(x).equals(q)) return [x, x];
  const pad = 4 * (Math.abs(x) * EPS + TINY);
  return [x - pad, x + pad];
}
export const toRI = (a) => new RationalInterval(qOfFloat(a[0]), qOfFloat(a[1]));

// ---- sine and cosine ---------------------------------------------------------------------------
// From their Taylor series with the remainder bounded, after taking off a
// multiple of π/2 known as an interval: no trust is placed in Math.sin.

const PI = [Math.PI, up(Math.PI)]; // π = 3.14159265358979323846…, and Math.PI is 3.14159265358979311…
const HALF_PI = [PI[0] / 2, PI[1] / 2];
function sinCosSmall(r) {
  // |r| ≤ 1: alternating series, 12 terms each, the remainder below |r|^25/25! < 7e-26
  const r2 = sqr(r);
  let s = pt(0);
  let c = pt(0);
  for (let j = 11; j >= 0; j--) {
    s = add(pt(1), neg(mul(r2, div(s, pt((2 * j + 2) * (2 * j + 3))))));
    c = add(pt(1), neg(mul(r2, div(c, pt((2 * j + 1) * (2 * j + 2))))));
  }
  s = mul(r, s);
  const pad = 1e-24;
  return [[Math.max(-1, s[0] - pad), Math.min(1, s[1] + pad)], [Math.max(-1, c[0] - pad), Math.min(1, c[1] + pad)]];
}
function sinCosPoint(x) {
  const k = Math.round(x / (Math.PI / 2));
  if (Math.abs(k) > 2 ** 30) return [[-1, 1], [-1, 1]];
  const r = sub(pt(x), scale(HALF_PI, k));
  const [s, c] = sinCosSmall(r);
  switch (((k % 4) + 4) % 4) {
    case 0: return [s, c];
    case 1: return [c, neg(s)];
    case 2: return [neg(s), neg(c)];
    default: return [neg(c), s];
  }
}
/** [sin, cos] of an interval. */
export function sinCos(a) {
  if (!(a[1] - a[0] < 6)) return [[-1, 1], [-1, 1]];
  const [s0, c0] = sinCosPoint(a[0]);
  if (a[0] === a[1]) return [s0, c0];
  const [s1, c1] = sinCosPoint(a[1]);
  const s = hull(s0, s1);
  const c = hull(c0, c1);
  // a turning point inside the interval takes the value to ±1; the test errs on the side of saying yes
  const reaches = (offset) => Math.ceil((a[0] - offset) / (2 * Math.PI) - 1e-9) <= Math.floor((a[1] - offset) / (2 * Math.PI) + 1e-9);
  if (reaches(Math.PI / 2)) s[1] = 1;
  if (reaches(-Math.PI / 2)) s[0] = -1;
  if (reaches(0)) c[1] = 1;
  if (reaches(Math.PI)) c[0] = -1;
  return [s, c];
}

// ---- formulas as tapes ----------------------------------------------------------------------------

const C = 0, V = 1, ADD = 2, SUB = 3, MUL = 4, DIV = 5, NEG = 6, SQR = 7, SIN = 8, COS = 9;

/**
 * Compile formula trees (analytic's parser, parameters already written in) to
 * a tape over the named variables: { ops, roots, n }.
 */
export function compile(asts, vars) {
  const ops = [];
  const seen = new Map();
  const push = (op) => ops.push(op) - 1;
  const constant = (iv) => push({ k: C, v: iv });
  const power = (a, n) => {
    if (n === 1) return a;
    const h = power(a, Math.floor(n / 2));
    const s = push({ k: SQR, a: h });
    return n % 2 ? push({ k: MUL, a: s, b: a }) : s;
  };
  const go = (node) => {
    if (seen.has(node)) return seen.get(node);
    let id;
    switch (node.type) {
      case "num":
        id = constant(ofQ(node.value));
        break;
      case "const":
        if (node.name !== "pi") refuse("unsupported", `the constant ${node.name} is not available here`);
        id = constant(PI);
        break;
      case "var": {
        const i = vars.indexOf(node.name);
        if (i < 0) refuse("unbound", `${node.name} has no value: name it as a variable or give it in params`);
        id = push({ k: V, i });
        break;
      }
      case "op": {
        if (node.op === "^") {
          const e = node.args[1];
          const inner = e.type === "op" && e.op === "neg" ? e.args[0] : e;
          if (inner.type !== "num" || inner.value.denominator !== 1n) refuse("unsupported", "powers here need a whole-number exponent");
          const n = Number(inner.value.numerator) * (inner === e ? 1 : -1);
          const base = go(node.args[0]);
          if (n === 0) id = constant(pt(1));
          else if (n > 0) id = power(base, n);
          else id = push({ k: DIV, a: constant(pt(1)), b: power(base, -n) });
          break;
        }
        const a = go(node.args[0]);
        if (node.op === "neg") {
          id = push({ k: NEG, a });
          break;
        }
        const b = go(node.args[1]);
        id = push({ k: { "+": ADD, "-": SUB, "*": MUL, "/": DIV }[node.op], a, b });
        break;
      }
      case "call": {
        if ((node.name !== "sin" && node.name !== "cos") || node.args.length !== 1) refuse("unsupported", `${node.name} is not available here: the fast solver knows + − × ÷, whole powers, sin and cos`);
        const a = go(node.args[0]);
        const key = `sc${a}`;
        if (!seen.has(key)) {
          const s = push({ k: SIN, a, twin: ops.length + 1 });
          push({ k: COS, a, twin: s });
          seen.set(key, s);
        }
        id = seen.get(key) + (node.name === "cos" ? 1 : 0);
        break;
      }
      default:
        refuse("unsupported", "this formula cannot be used here");
    }
    seen.set(node, id);
    return id;
  };
  const roots = asts.map(go);
  return { ops, roots, n: vars.length };
}

/**
 * Taylor coefficients in time of the solution of y′ = f(y) through y0 (a
 * vector of intervals): Y[i][k] for k = 0 … K, each an interval that holds
 * the coefficient for every starting point in y0.
 */
export function series(tape, y0, K) {
  return run(tape, y0, K, true).Y;
}
/** The values of the tape's formulas over a box (they need not be a vector field). */
export function values(tape, y) {
  const { T } = run(tape, y, 1, false);
  return tape.roots.map((r) => T[r][0]);
}
function run(tape, y0, K, feed) {
  const { ops, roots, n } = tape;
  const N = ops.length;
  const T = Array.from({ length: N }, () => new Array(K));
  const Y = y0.map((v) => [v]);
  for (let k = 0; k < K; k++) {
    for (let i = 0; i < N; i++) {
      const op = ops[i];
      let out;
      switch (op.k) {
        case C:
          out = k === 0 ? op.v : ZERO;
          break;
        case V:
          out = Y[op.i][k];
          break;
        case ADD:
          out = add(T[op.a][k], T[op.b][k]);
          break;
        case SUB:
          out = sub(T[op.a][k], T[op.b][k]);
          break;
        case NEG:
          out = neg(T[op.a][k]);
          break;
        case MUL: {
          const a = T[op.a];
          const b = T[op.b];
          out = ZERO;
          for (let j = 0; j <= k; j++) out = add(out, mul(a[j], b[k - j]));
          break;
        }
        case SQR: {
          const a = T[op.a];
          if (k === 0) out = sqr(a[0]);
          else {
            out = ZERO;
            for (let j = 0; 2 * j < k; j++) out = add(out, mul(a[j], a[k - j]));
            out = scale(out, 2);
            if (k % 2 === 0) out = add(out, sqr(a[k / 2]));
          }
          break;
        }
        case DIV: {
          const a = T[op.a];
          const b = T[op.b];
          const c = T[i];
          out = a[k];
          for (let j = 0; j < k; j++) out = sub(out, mul(c[j], b[k - j]));
          out = div(out, b[0]);
          break;
        }
        case SIN:
        case COS: {
          const u = T[op.a];
          if (k === 0) out = sinCos(u[0])[op.k === SIN ? 0 : 1];
          else {
            // s_k = (1/k) Σ j·u_j·c_(k−j),  c_k = −(1/k) Σ j·u_j·s_(k−j)
            const other = T[op.twin];
            out = ZERO;
            for (let j = 1; j <= k; j++) out = add(out, scale(mul(u[j], other[k - j]), j));
            out = div(out, pt(op.k === SIN ? k : -k));
          }
          break;
        }
      }
      T[i][k] = out;
    }
    if (feed) for (let i = 0; i < n; i++) Y[i].push(div(T[roots[i]][k], pt(k + 1)));
  }
  return { Y, T };
}
/** f(y) over a box. */
export const field = values;

// ---- small linear algebra --------------------------------------------------------------------------

export const midV = (v) => v.map(mid);
export const ivV = (v) => v.map(pt);
export const ivM = (M) => M.map((row) => row.map(pt));
export const midM = (M) => M.map((row) => row.map(mid));
export function mmul(A, B) {
  return A.map((row) =>
    B[0].map((_, j) => {
      let s = ZERO;
      for (let k = 0; k < row.length; k++) s = add(s, mul(row[k], B[k][j]));
      return s;
    }),
  );
}
export const mvec = (A, v) => A.map((row) => row.reduce((s, a, k) => add(s, mul(a, v[k])), ZERO));
export const identity = (n) => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));

/** The inverse of a matrix of doubles, in doubles (no guarantee), or null. */
export function inverseFloat(M) {
  const n = M.length;
  const A = M.map((row, i) => [...row, ...identity(n)[i]]);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let i = c + 1; i < n; i++) if (Math.abs(A[i][c]) > Math.abs(A[p][c])) p = i;
    if (!(Math.abs(A[p][c]) > 0)) return null;
    [A[c], A[p]] = [A[p], A[c]];
    const d = A[c][c];
    for (let j = 0; j < 2 * n; j++) A[c][j] /= d;
    for (let i = 0; i < n; i++) {
      if (i === c) continue;
      const f = A[i][c];
      if (f !== 0) for (let j = 0; j < 2 * n; j++) A[i][j] -= f * A[c][j];
    }
  }
  return A.map((row) => row.slice(n));
}
/** An interval matrix that contains the inverse of a matrix of doubles, or null. */
export function inverseEnclosed(M) {
  const n = M.length;
  const A = M.map((row, i) => [...row.map(pt), ...identity(n)[i].map(pt)]);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let i = c + 1; i < n; i++) if (Math.abs(mid(A[i][c])) > Math.abs(mid(A[p][c]))) p = i;
    if (hasZero(A[p][c])) return null;
    [A[c], A[p]] = [A[p], A[c]];
    const d = A[c][c];
    A[c] = A[c].map((x) => div(x, d));
    A[c][c] = pt(1);
    for (let i = 0; i < n; i++) {
      if (i === c) continue;
      const f = A[i][c];
      if (isZero(f)) continue;
      A[i] = A[i].map((x, j) => sub(x, mul(f, A[c][j])));
      A[i][c] = ZERO;
    }
  }
  return A.map((row) => row.slice(n));
}
/**
 * A matrix of doubles with orthonormal columns (up to rounding) spanning the
 * columns of M, longest weighted column first; missing directions are filled
 * in. Used as the frame of a set, so any invertible matrix would be valid:
 * orthogonality only keeps the enclosures tight.
 */
export function frame(M, weights) {
  const n = M.length;
  const cols = Array.from({ length: n }, (_, j) => M.map((row) => row[j]));
  const order = cols.map((c, j) => j).sort((a, b) => Math.hypot(...cols[b]) * weights[b] - Math.hypot(...cols[a]) * weights[a]);
  const basis = [];
  const take = (v, floor) => {
    const w = v.slice();
    for (let pass = 0; pass < 2; pass++) {
      for (const b of basis) {
        const d = w.reduce((s, x, i) => s + x * b[i], 0);
        for (let i = 0; i < n; i++) w[i] -= d * b[i];
      }
    }
    const len = Math.hypot(...w);
    if (!(len > floor)) return false;
    basis.push(w.map((x) => x / len));
    return true;
  };
  for (const j of order) {
    const len = Math.hypot(...cols[j]);
    if (len > 0) take(cols[j], 1e-9 * len);
  }
  for (let pass = 0; basis.length < n && pass < 2; pass++) {
    for (let i = 0; i < n && basis.length < n; i++) take(identity(n)[i], pass ? 1e-12 : 0.3);
  }
  if (basis.length < n) refuse("degenerate", "could not build a frame for the set");
  return Array.from({ length: n }, (_, i) => basis.map((b) => b[i]));
}

// ---- equations in a box ---------------------------------------------------------------------------

/**
 * Krawczyk's test for G(x) = 0 on the box X (G and its Jacobian DG take
 * interval vectors). Returns { box } when a box inside X contains exactly one
 * zero, which is also the only zero in X; { empty: true } when X is proven to
 * hold no zero; { shrunk } when all that is known is that every zero of X
 * lies in the smaller box; null when nothing was learnt.
 */
export function krawczykTest(G, DG, X) {
  let box = X;
  let proven = false;
  for (let iter = 0; iter < 40; iter++) {
    const m = midV(box);
    let J;
    let Gm;
    try {
      J = DG(box);
      Gm = G(ivV(m));
    } catch {
      break;
    }
    const Cm = inverseFloat(midM(J));
    if (!Cm || !Cm.every((row) => row.every(Number.isFinite))) break;
    const Ci = ivM(Cm);
    const n = m.length;
    // K = m − C·G(m) + (I − C·J(X))·(X − m)
    const CJ = mmul(Ci, J);
    const E = CJ.map((row, i) => row.map((x, j) => sub(pt(i === j ? 1 : 0), x)));
    const CG = mvec(Ci, Gm);
    const dX = box.map((x, i) => sub(x, pt(m[i])));
    const Ed = mvec(E, dX);
    const K = Array.from({ length: n }, (_, i) => add(sub(pt(m[i]), CG[i]), Ed[i]));
    if (!K.every(finite)) break;
    if (K.every((k, i) => inside(k, box[i]))) proven = true;
    const next = K.map((k, i) => meet(k, box[i]));
    if (next.some((x) => x === null)) return { empty: true };
    const shrink = next.reduce((s, x, i) => Math.max(s, width(box[i]) > 0 ? width(x) / width(box[i]) : 1), 0);
    box = next;
    if (shrink > 0.9) break;
  }
  if (proven) return { box: polish(G, DG, box) };
  return box === X ? null : { shrunk: box };
}
/** A box known to hold exactly one zero, made small: Newton's method in doubles, then the test again round its answer. */
function polish(G, DG, box) {
  const big = (b) => b.some((x) => width(x) > 1e-9 * (1 + mag(x)));
  if (!big(box)) return box;
  try {
    let x = midV(box);
    for (let it = 0; it < 60; it++) {
      const Cm = inverseFloat(midM(DG(ivV(x))));
      if (!Cm) return box;
      const g = midV(G(ivV(x)));
      const dx = Cm.map((row) => row.reduce((s, c, j) => s + c * g[j], 0));
      x = x.map((xi, i) => xi - dx[i]);
      if (!x.every((xi, i) => xi >= box[i][0] && xi <= box[i][1])) return box;
      if (Math.max(...dx.map(Math.abs)) < 1e-15 * (1 + Math.max(...x.map(Math.abs)))) break;
    }
    for (const rad of [1e-13, 1e-11, 1e-9, 1e-7, 1e-5]) {
      const small = x.map((xi, i) => {
        const r = rad * (1 + Math.abs(xi));
        return [Math.max(box[i][0], xi - r), Math.min(box[i][1], xi + r)];
      });
      const m = x;
      const J = DG(small);
      const Cm = inverseFloat(midM(J));
      if (!Cm) continue;
      const Ci = ivM(Cm);
      const E = mmul(Ci, J).map((row, i) => row.map((a, j) => sub(pt(i === j ? 1 : 0), a)));
      const CG = mvec(Ci, G(ivV(m)));
      const Ed = mvec(E, small.map((s, i) => sub(s, pt(m[i]))));
      const K = m.map((mi, i) => add(sub(pt(mi), CG[i]), Ed[i]));
      // K inside the small box: the zero is there (and it is the only one in the big box)
      if (K.every((k, i) => inside(k, small[i]))) return K;
    }
  } catch {
    return box;
  }
  return box;
}
/** The box of krawczykTest when the zero is proven, else null. */
export function krawczyk(G, DG, X) {
  return krawczykTest(G, DG, X)?.box ?? null;
}

/**
 * Every zero of G in the box, each proven to exist and to be alone in its
 * small box. `complete` says the rest of the box was proven empty.
 */
export function solveAll(G, DG, start, { tol = 1e-9, maxBoxes = 20000 } = {}) {
  const found = [];
  let undecided = 0;
  const stack = [start];
  let budget = maxBoxes;
  while (stack.length && budget-- > 0) {
    const box = stack.pop();
    if (found.some((f) => box.every((x, i) => within(x, f.region[i])))) continue;
    let val = null;
    try {
      val = G(box);
    } catch {
      val = null;
    }
    if (val && val.some((x) => !hasZero(x))) continue;
    const w = box.reduce((s, x) => Math.max(s, width(x)), 0);
    if (val) {
      const region = box.map((x) => {
        const r = width(x) / 8 + 1e-300;
        return [x[0] - r, x[1] + r];
      });
      const test = krawczykTest(G, DG, region);
      if (test?.empty) continue;
      if (test?.box) {
        const hit = test.box;
        if (!found.some((f) => f.box.every((x, i) => meet(x, hit[i]) !== null))) found.push({ box: hit, region });
        continue;
      }
      if (test?.shrunk) {
        // every zero of the box lies in the smaller one: go on with their overlap when that is a real gain
        const cut = test.shrunk.map((x, i) => meet(x, box[i]));
        if (cut.some((x) => x === null)) continue;
        const gain = cut.reduce((sum, x, i) => sum + (width(box[i]) > 0 ? width(x) / width(box[i]) : 1), 0) / cut.length;
        if (gain < 0.7) {
          stack.push(cut);
          continue;
        }
      }
    }
    if (w < tol) {
      undecided++;
      continue;
    }
    let k = 0;
    for (let i = 1; i < box.length; i++) if (width(box[i]) > width(box[k])) k = i;
    const m = mid(box[k]);
    const left = box.slice();
    const right = box.slice();
    left[k] = [box[k][0], m];
    right[k] = [m, box[k][1]];
    stack.push(left, right);
  }
  // keep the zeros that can lie in the box asked about
  const list = found.filter((f) => f.box.every((x, i) => meet(x, start[i]) !== null)).map((f) => f.box);
  list.complete = stack.length === 0 && undecided === 0;
  return list;
}

/**
 * Eigenvalues of a small interval matrix (order 1, 2 or 3) as enclosures:
 * [{ re, im }] with im = [0, 0] for a real one and a conjugate pair listed
 * with im ≥ 0 first; null when they cannot be separated.
 */
export function eigenvalues(M) {
  const n = M.length;
  const quadratic = (b, c) => {
    // λ² + bλ + c
    const disc = sub(sqr(b), scale(c, 4));
    const half = scale(neg(b), 0.5);
    if (disc[0] > 0) {
      const s = scale(sqrt(disc), 0.5);
      return [{ re: sub(half, s), im: ZERO }, { re: add(half, s), im: ZERO }];
    }
    if (disc[1] < 0) {
      const s = scale(sqrt(neg(disc)), 0.5);
      return [{ re: half, im: s }, { re: half, im: neg(s) }];
    }
    return null;
  };
  if (n === 1) return [{ re: M[0][0], im: ZERO }];
  if (n === 2) {
    const tr = add(M[0][0], M[1][1]);
    const det = sub(mul(M[0][0], M[1][1]), mul(M[0][1], M[1][0]));
    return quadratic(neg(tr), det);
  }
  if (n !== 3) return null;
  // λ³ + a2λ² + a1λ + a0
  const a2 = neg(add(add(M[0][0], M[1][1]), M[2][2]));
  const minor = (i, j) => sub(mul(M[i][i], M[j][j]), mul(M[i][j], M[j][i]));
  const a1 = add(add(minor(0, 1), minor(0, 2)), minor(1, 2));
  const det3 = add(
    sub(mul(M[0][0], sub(mul(M[1][1], M[2][2]), mul(M[1][2], M[2][1]))), mul(M[0][1], sub(mul(M[1][0], M[2][2]), mul(M[1][2], M[2][0])))),
    mul(M[0][2], sub(mul(M[1][0], M[2][1]), mul(M[1][1], M[2][0]))),
  );
  const a0 = neg(det3);
  const p = (x) => add(mul(add(mul(add(x, a2), x), a1), x), a0);
  const dp = (x) => add(mul(add(scale(x, 3), scale(a2, 2)), x), a1);
  const [b2, b1, b0] = [mid(a2), mid(a1), mid(a0)];
  const pf = (x) => ((x + b2) * x + b1) * x + b0;
  // real roots of the midpoint cubic, by bisection between its turning points
  const bound = 1 + Math.max(Math.abs(b2), Math.abs(b1), Math.abs(b0));
  const cuts = [-bound];
  const dd = 4 * b2 * b2 - 12 * b1;
  if (dd > 0) cuts.push((-2 * b2 - Math.sqrt(dd)) / 6, (-2 * b2 + Math.sqrt(dd)) / 6);
  cuts.push(bound);
  const guesses = [];
  for (let i = 0; i + 1 < cuts.length; i++) {
    let [lo, hi] = [cuts[i], cuts[i + 1]];
    if (!(pf(lo) * pf(hi) <= 0)) continue;
    for (let it = 0; it < 200 && hi - lo > 1e-15 * (1 + Math.abs(lo)); it++) {
      const m = (lo + hi) / 2;
      if (pf(lo) * pf(m) <= 0) hi = m;
      else lo = m;
    }
    guesses.push((lo + hi) / 2);
  }
  // prove one real root by interval Newton, the best separated first
  guesses.sort((x, y) => Math.abs(3 * y * y + 2 * b2 * y + b1) - Math.abs(3 * x * x + 2 * b2 * x + b1));
  for (const g of guesses) {
    for (const rad of [1e-10, 1e-7, 1e-4]) {
      const r = rad * (1 + Math.abs(g));
      let X = [g - r, g + r];
      let ok = false;
      try {
        for (let it = 0; it < 6; it++) {
          const d = dp(X);
          if (hasZero(d)) break;
          const Nw = sub(pt(g), div(p(pt(g)), d));
          if (inside(Nw, X)) ok = true;
          const next = meet(Nw, X);
          if (!next) {
            ok = false;
            break;
          }
          X = next;
        }
      } catch {
        ok = false;
      }
      if (!ok) continue;
      // deflate: λ² + (a2 + r)λ + (a1 + r(a2 + r))
      const b = add(a2, X);
      const rest = quadratic(b, add(a1, mul(X, b)));
      if (!rest) return null;
      return [{ re: X, im: ZERO }, ...rest];
    }
  }
  return null;
}
