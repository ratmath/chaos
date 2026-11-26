// Julia sets, the Mandelbrot set and Newton basins, with honest pixels.
//
// A pixel is a small rectangle of the plane, not a point. Each one is put in
// one of three classes, and the first two are proofs about every point of
// the rectangle:
//
//   outside    the orbit escapes: some iterate is certainly farther out than
//              the escape radius, for every point of the pixel at once
//   inside     the orbit is trapped: it enters a box B that p steps of the map
//              carry into B, for every point of the pixel (Mandelbrot: also
//              the exact tests for the main cardioid and the period-2 disc)
//   undecided  neither was shown within the iteration limit. The boundary,
//              where all the interest is, is always undecided at some zoom.
//
// Grids use interval arithmetic in floating point with every result moved
// outward by more than its rounding error, so the classes are certified while
// staying fast. A single exact point is classified with the complex package's
// rational rectangles instead (mandelbrotPoint, juliaPoint).
import { Rational } from "@ratmath/core";
import { Rect, Complex } from "@ratmath/complex";
import { refuse, toQ, ZERO, ONE, Q } from "./tree.js";

// ---- outward-rounded float intervals ------------------------------------------------------
// After a correctly rounded operation the error is at most half a unit in the last
// place, which is at most |x|·2^-53; moving each end out by |x|·2^-52 + the smallest
// subnormal covers it, with room to spare for the rounding of the move itself.
const EPS = 2 ** -52;
const TINY = Number.MIN_VALUE;
const dn = (x) => x - (Math.abs(x) * EPS + TINY);
const upx = (x) => x + (Math.abs(x) * EPS + TINY);
/** [lo, hi] products */
const mulLo = (a0, a1, b0, b1) => dn(Math.min(a0 * b0, a0 * b1, a1 * b0, a1 * b1));
const mulHi = (a0, a1, b0, b1) => upx(Math.max(a0 * b0, a0 * b1, a1 * b0, a1 * b1));
const sqLo = (a0, a1) => (a0 > 0 ? dn(a0 * a0) : a1 < 0 ? dn(a1 * a1) : 0);
const sqHi = (a0, a1) => upx(Math.max(a0 * a0, a1 * a1));
/** A float interval around a rational. */
const qLo = (q) => dn(q.toNumber());
const qHi = (q) => upx(q.toNumber());

/** Complex interval boxes as [x0, x1, y0, y1]. */
function csq([x0, x1, y0, y1]) {
  // (x + iy)² = x² − y² + 2xy·i
  return [dn(sqLo(x0, x1) - sqHi(y0, y1)), upx(sqHi(x0, x1) - sqLo(y0, y1)), dn(2 * mulLo(x0, x1, y0, y1)), upx(2 * mulHi(x0, x1, y0, y1))];
}
const cadd = (a, b) => [dn(a[0] + b[0]), upx(a[1] + b[1]), dn(a[2] + b[2]), upx(a[3] + b[3])];
const csub = (a, b) => [dn(a[0] - b[1]), upx(a[1] - b[0]), dn(a[2] - b[3]), upx(a[3] - b[2])];
export function cmul(a, b) {
  const rr0 = mulLo(a[0], a[1], b[0], b[1]);
  const rr1 = mulHi(a[0], a[1], b[0], b[1]);
  const ii0 = mulLo(a[2], a[3], b[2], b[3]);
  const ii1 = mulHi(a[2], a[3], b[2], b[3]);
  const ri0 = mulLo(a[0], a[1], b[2], b[3]);
  const ri1 = mulHi(a[0], a[1], b[2], b[3]);
  const ir0 = mulLo(a[2], a[3], b[0], b[1]);
  const ir1 = mulHi(a[2], a[3], b[0], b[1]);
  return [dn(rr0 - ii1), upx(rr1 - ii0), dn(ri0 + ir0), upx(ri1 + ir1)];
}
/** a/b, or null when b may be 0. */
export function cdiv(a, b) {
  const n0 = dn(sqLo(b[0], b[1]) + sqLo(b[2], b[3]));
  const n1 = upx(sqHi(b[0], b[1]) + sqHi(b[2], b[3]));
  if (!(n0 > 0)) return null;
  const num = cmul(a, [b[0], b[1], -b[3], -b[2]]);
  const inv = [dn(1 / n1), upx(1 / n0)];
  const scale = (lo, hi) => [mulLo(lo, hi, inv[0], inv[1]), mulHi(lo, hi, inv[0], inv[1])];
  const [r0, r1] = scale(num[0], num[1]);
  const [i0, i1] = scale(num[2], num[3]);
  return [r0, r1, i0, i1];
}
const normLo = (z) => dn(sqLo(z[0], z[1]) + sqLo(z[2], z[3]));
const normHi = (z) => upx(sqHi(z[0], z[1]) + sqHi(z[2], z[3]));
const inside = (a, b) => a[0] > b[0] && a[1] < b[1] && a[2] > b[2] && a[3] < b[3];
const finite = (z) => z.every(Number.isFinite);
const inflate = (z, f, extra) => {
  const cx = (z[0] + z[1]) / 2;
  const cy = (z[2] + z[3]) / 2;
  const r = Math.max(z[1] - z[0], z[3] - z[2]) * f + extra;
  return [dn(cx - r), upx(cx + r), dn(cy - r), upx(cy + r)];
};

// ---- grids -----------------------------------------------------------------------------------

/**
 * Classify a box; when that fails, split it into four and classify those, down
 * to `depth` levels. The box gets a class only when every piece has that same
 * class, so the answer is still about every point of the box.
 */
function bySubdivision(B, depth, single) {
  const k = single(B);
  if (k || depth <= 0) return k;
  const mx = (B[0] + B[1]) / 2;
  const my = (B[2] + B[3]) / 2;
  let first = -1;
  for (const Q4 of [[B[0], mx, B[2], my], [mx, B[1], B[2], my], [B[0], mx, my, B[3]], [mx, B[1], my, B[3]]]) {
    const r = bySubdivision(Q4, depth - 1, single);
    if (!r || (first >= 0 && r !== first)) return 0;
    first = r;
  }
  return first;
}

const UNDECIDED = 0;
const INSIDE = 1;
const OUTSIDE = 2;

/** A classified grid of pixels. kind[i] is 0 (undecided), 1 (inside), 2 (outside); for Newton basins, k + 1 for root k. */
export class PixelGrid {
  constructor(fields) {
    Object.assign(this, fields);
    Object.freeze(this);
  }
  at(i, j) {
    return this.kind[j * this.width + i];
  }
  /** Rows of text, top row first: one character per pixel. */
  toText(chars = this.labels ? "?123456789" : "?█ ") {
    const lines = [];
    for (let j = this.height - 1; j >= 0; j--) {
      let s = "";
      for (let i = 0; i < this.width; i++) s += chars[this.kind[j * this.width + i]] ?? "?";
      lines.push(s);
    }
    return lines.join("\n");
  }
  toString() {
    const c = this.counts;
    return `${this.name}: ${this.width}×${this.height} pixels, ${Object.entries(c).map(([k, v]) => `${v} ${k}`).join(", ")}`;
  }
}

function edges(range, n, what) {
  if (!Array.isArray(range) || range.length !== 2) refuse("bad-argument", `${what} must be [low, high]`);
  const lo = toQ(range[0]);
  const hi = toQ(range[1]);
  if (lo.compareTo(hi) >= 0) refuse("bad-argument", `${what}: low must be below high`);
  const step = hi.subtract(lo).divide(new Rational(BigInt(n)));
  return Array.from({ length: n + 1 }, (_, i) => lo.add(step.multiply(new Rational(BigInt(i)))));
}
const pixelBox = (xs, ys, i, j) => [qLo(xs[i]), qHi(xs[i + 1]), qLo(ys[j]), qHi(ys[j + 1])];
const checkSize = (w, h) => {
  if (!Number.isInteger(w) || !Number.isInteger(h) || w < 1 || h < 1) refuse("bad-argument", "width and height are whole numbers");
  if (w * h > 1_000_000) refuse("too-large", "at most a million pixels");
};

/** Is the whole box of c values in the main cardioid or the period-2 disc? (exact tests, in intervals) */
function cardioidOrBulb(C) {
  // cardioid: q(q + (x − 1/4)) < y²/4 with q = (x − 1/4)² + y²
  const x = [dn(C[0] - 0.25), upx(C[1] - 0.25)];
  const y2 = [sqLo(C[2], C[3]), sqHi(C[2], C[3])];
  const q = [dn(sqLo(x[0], x[1]) + y2[0]), upx(sqHi(x[0], x[1]) + y2[1])];
  const t = [dn(q[0] + x[0]), upx(q[1] + x[1])];
  const lhsHi = mulHi(q[0], q[1], t[0], t[1]);
  if (dn(lhsHi - y2[0] / 4) < 0 && upx(lhsHi - y2[0] / 4) < 0) return "the main cardioid";
  // the period-2 disc: (x + 1)² + y² < 1/16
  const b = upx(sqHi(dn(C[0] + 1), upx(C[1] + 1)) + y2[1]);
  if (b < 1 / 16) return "the period-2 disc";
  return null;
}

/** Can a box around the orbit be trapped by p steps of z → z² + c (c in the box C)? */
function trapped(Z, C, maxPeriod) {
  for (const f of [1.5, 3]) {
    const B = inflate(Z, f, 1e-12);
    let W = B;
    for (let p = 1; p <= maxPeriod; p++) {
      W = cadd(csq(W), C);
      if (!finite(W) || normHi(W) > 16) break;
      if (inside(W, B)) return p;
    }
  }
  return 0;
}

/**
 * The Mandelbrot set on a grid: c ranges over re × im, each pixel classified.
 *   mandelbrot({ re: ["-2", "1/2"], im: ["-5/4", "5/4"], width: 60, height: 30 })
 * Options: maxIter (default 200), maxPeriod (default 12, for the trapping test),
 * depth (default 2: a pixel that cannot be classified whole is split into 4,
 * then 16 pieces, and gets a class when every piece has it).
 */
export function mandelbrot({ re = ["-2", "1/2"], im = ["-5/4", "5/4"], width = 80, height = 40, maxIter = 200, maxPeriod = 12, depth = 2 } = {}) {
  checkSize(width, height);
  const xs = edges(re, width, "re");
  const ys = edges(im, height, "im");
  const kind = new Uint8Array(width * height);
  const escape = new Int32Array(width * height);
  const period = new Uint8Array(width * height);
  const checks = new Set([8, 16, 32, 64, 128, 256, 512, 1024, 2048]);
  let lastEscape = 0;
  let lastPeriod = 0;
  const single = (C) => {
    if (cardioidOrBulb(C)) return INSIDE;
    let Z = [0, 0, 0, 0];
    for (let n = 1; n <= maxIter; n++) {
      Z = cadd(csq(Z), C);
      if (!finite(Z)) return UNDECIDED;
      if (normLo(Z) > 4) {
        lastEscape = Math.max(lastEscape, n);
        return OUTSIDE;
      }
      if (checks.has(n)) {
        if (normHi(Z) >= 64) return UNDECIDED; // the box has spread too wide to trap
        const p = trapped(Z, C, maxPeriod);
        if (p) {
          lastPeriod = p;
          return INSIDE;
        }
      }
    }
    return UNDECIDED;
  };
  for (let j = 0; j < height; j++)
    for (let i = 0; i < width; i++) {
      const at = j * width + i;
      lastEscape = 0;
      lastPeriod = 0;
      kind[at] = bySubdivision(pixelBox(xs, ys, i, j), depth, single);
      if (kind[at] === OUTSIDE) escape[at] = lastEscape;
      if (kind[at] === INSIDE) period[at] = lastPeriod;
    }
  return finish("Mandelbrot set", { width, height, kind, escape, period, re: xs[0], reHigh: xs[width], im: ys[0], imHigh: ys[height], maxIter });
}

function finish(name, f) {
  let inside = 0;
  let outside = 0;
  let undecided = 0;
  for (const k of f.kind) k === INSIDE ? inside++ : k === OUTSIDE ? outside++ : undecided++;
  return new PixelGrid({ name, ...f, counts: { inside, outside, undecided } });
}

/** A float approximation of the attracting cycle of z → z² + c, or null. */
function attractingCycle(cr, ci, maxPeriod = 24) {
  let x = 0;
  let y = 0;
  for (let n = 0; n < 4000; n++) {
    const t = x * x - y * y + cr;
    y = 2 * x * y + ci;
    x = t;
    if (x * x + y * y > 1e6) return null;
  }
  // the period: the first p with z_p ≈ z_0
  const z0 = [x, y];
  let [a, b] = z0;
  for (let p = 1; p <= maxPeriod; p++) {
    const t = a * a - b * b + cr;
    b = 2 * a * b + ci;
    a = t;
    if (Math.hypot(a - z0[0], b - z0[1]) < 1e-9 * (1 + Math.hypot(...z0))) {
      const cycle = [];
      let [u, w] = z0;
      for (let k = 0; k < p; k++) {
        cycle.push([u, w]);
        const s = u * u - w * w + cr;
        w = 2 * u * w + ci;
        u = s;
      }
      return cycle;
    }
  }
  return null;
}

/**
 * The filled Julia set of z → z² + c for an exact c, on a grid of z.
 *   julia("-1", { re: ["-2", "2"], im: ["-1", "1"] })
 * Outside: an iterate is certainly beyond max(2, |c|). Inside: the orbit
 * enters a box around a point of the attracting cycle that p steps carry into
 * itself (proved for the box). Without an attracting cycle nothing is inside.
 */
export function julia(c, { re = ["-2", "2"], im = ["-3/2", "3/2"], width = 80, height = 40, maxIter = 300, depth = 2 } = {}) {
  checkSize(width, height);
  const cz = Complex.from(c);
  const C = [qLo(cz.re), qHi(cz.re), qLo(cz.im), qHi(cz.im)];
  const R2 = Math.max(4, normHi(C)) * (1 + 1e-12);
  // trapping boxes around the attracting cycle, each proved
  const traps = [];
  const cyc = attractingCycle(cz.re.toNumber(), cz.im.toNumber());
  if (cyc) {
    const p = cyc.length;
    for (const [u, w] of cyc) {
      for (const r of [1e-2, 1e-3, 1e-4, 1e-6]) {
        const B = [dn(u - r), upx(u + r), dn(w - r), upx(w + r)];
        let W = B;
        for (let k = 0; k < p; k++) W = cadd(csq(W), C);
        if (inside(W, B)) {
          traps.push(B);
          break;
        }
      }
    }
  }
  const xs = edges(re, width, "re");
  const ys = edges(im, height, "im");
  const kind = new Uint8Array(width * height);
  const escape = new Int32Array(width * height);
  let lastEscape = 0;
  const single = (Z0) => {
    let Z = Z0;
    for (let n = 0; n <= maxIter; n++) {
      if (n) Z = cadd(csq(Z), C);
      if (!finite(Z)) return UNDECIDED;
      if (normLo(Z) > R2) {
        lastEscape = Math.max(lastEscape, n);
        return OUTSIDE;
      }
      if (traps.some((B) => inside(Z, B))) return INSIDE;
      if (normHi(Z) > 64 * R2) return UNDECIDED;
    }
    return UNDECIDED;
  };
  for (let j = 0; j < height; j++)
    for (let i = 0; i < width; i++) {
      const at = j * width + i;
      lastEscape = 0;
      kind[at] = bySubdivision(pixelBox(xs, ys, i, j), depth, single);
      if (kind[at] === OUTSIDE) escape[at] = lastEscape;
    }
  return finish(`filled Julia set for c = ${cz}`, { width, height, kind, escape, c: cz, cycle: cyc, traps: traps.length, re: xs[0], reHigh: xs[width], im: ys[0], imHigh: ys[height], maxIter });
}

// ---- Newton basins --------------------------------------------------------------------------

/** Complex polynomial evaluation in boxes (Horner). Coefficients highest first, as boxes. */
const hornerBox = (cs, Z) => cs.slice(1).reduce((acc, c) => cadd(cmul(acc, Z), c), cs[0]);

/** Roots of a polynomial in floats (Durand–Kerner), coefficients highest first as [re, im]. */
function floatRoots(cs) {
  const n = cs.length - 1;
  const [ar, ai] = cs[0];
  const d = ar * ar + ai * ai;
  const mon = cs.map(([r, i]) => [(r * ar + i * ai) / d, (i * ar - r * ai) / d]);
  let z = Array.from({ length: n }, (_, k) => [0.4 + 0.9 * Math.cos((2 * Math.PI * k) / n + 0.3), 0.9 * Math.sin((2 * Math.PI * k) / n + 0.3)]);
  const ev = ([x, y]) => mon.slice(1).reduce(([pr, pi], [cr, ci]) => [pr * x - pi * y + cr, pr * y + pi * x + ci], [1, 0]);
  for (let it = 0; it < 500; it++) {
    let moved = 0;
    z = z.map((zi, i) => {
      let [num0, num1] = ev(zi);
      let den = [1, 0];
      z.forEach((zj, j) => {
        if (i !== j) {
          const a = zi[0] - zj[0];
          const b = zi[1] - zj[1];
          den = [den[0] * a - den[1] * b, den[0] * b + den[1] * a];
        }
      });
      const dd = den[0] * den[0] + den[1] * den[1] || 1e-300;
      const q = [(num0 * den[0] + num1 * den[1]) / dd, (num1 * den[0] - num0 * den[1]) / dd];
      moved = Math.max(moved, Math.hypot(...q));
      return [zi[0] - q[0], zi[1] - q[1]];
    });
    if (moved < 1e-15) break;
  }
  return z;
}

/**
 * Basins of attraction of Newton's method for a polynomial, given by its
 * coefficients (highest first; exact, complex ones written "1/2 + 3i").
 *   newtonBasins(["1", "0", "0", "-1"])   z³ − 1, three basins
 * Each simple root gets a box B with N(B) inside B for Newton's map N, and p′
 * not 0 on B: that proves B holds exactly one root, which attracts all of B.
 * A pixel is labelled k + 1 when its Newton orbit (in intervals) enters the
 * k-th box; 0 when that was not shown.
 */
export function newtonBasins(coefficients, { re = ["-2", "2"], im = ["-3/2", "3/2"], width = 80, height = 40, maxIter = 60, depth = 3 } = {}) {
  checkSize(width, height);
  const cz = coefficients.map((c) => Complex.from(c));
  while (cz.length && cz[0].re.equals(ZERO) && cz[0].im.equals(ZERO)) cz.shift();
  const deg = cz.length - 1;
  if (deg < 1) refuse("bad-argument", "the polynomial needs degree at least 1");
  const box = (z) => [qLo(z.re), qHi(z.re), qLo(z.im), qHi(z.im)];
  const P = cz.map(box);
  const dP = cz.slice(0, -1).map((z, k) => box(z.scale(new Rational(BigInt(deg - k)))));
  const d2P = dP.length > 1 ? cz.slice(0, -2).map((z, k) => box(z.scale(new Rational(BigInt((deg - k) * (deg - k - 1)))))) : [[0, 0, 0, 0]];
  // Newton's map N(z) = z − p/p′ at a point (a thin box), and its derivative N′ = p·p″/p′² over a box
  const Npoint = (Z) => {
    const q = cdiv(hornerBox(P, Z), hornerBox(dP, Z));
    return q ? csub(Z, q) : null;
  };
  const dN = (Z) => {
    const d = hornerBox(dP, Z);
    return cdiv(cmul(hornerBox(P, Z), hornerBox(d2P, Z)), cmul(d, d));
  };
  // the mean-value form: N(Z) ⊆ N(m) + N′(Z)·(Z − m), much tighter than evaluating N on Z directly
  const N = (Z) => {
    const mx = (Z[0] + Z[1]) / 2;
    const my = (Z[2] + Z[3]) / 2;
    const m = [mx, mx, my, my];
    const Nm = Npoint(m);
    const D = dN(Z);
    if (!Nm || !D) return null;
    return cadd(Nm, cmul(D, csub(Z, m)));
  };
  const roots = floatRoots(cz.map((z) => [z.re.toNumber(), z.im.toNumber()]));
  const traps = [];
  const unproved = [];
  for (const [u, w] of roots) {
    let ok = null;
    for (const r of [1e-3, 1e-5, 1e-7, 1e-9]) {
      const s = r * (1 + Math.hypot(u, w));
      const B = [dn(u - s), upx(u + s), dn(w - s), upx(w + s)];
      const NB = N(B);
      if (NB && inside(NB, B)) {
        ok = B;
        break;
      }
    }
    if (ok && !traps.some((T) => !(T[1] < ok[0] || ok[1] < T[0] || T[3] < ok[2] || ok[3] < T[2]))) traps.push(ok);
    else unproved.push([u, w]);
  }
  const xs = edges(re, width, "re");
  const ys = edges(im, height, "im");
  const kind = new Uint8Array(width * height);
  const single = (Z0) => {
    let Z = Z0;
    for (let n = 0; n <= maxIter && Z; n++) {
      const k = traps.findIndex((B) => inside(Z, B));
      if (k >= 0) return k + 1;
      Z = N(Z);
      if (Z && (!finite(Z) || Z[1] - Z[0] > 1e3 || Z[3] - Z[2] > 1e3)) Z = null;
    }
    return 0;
  };
  for (let j = 0; j < height; j++)
    for (let i = 0; i < width; i++) kind[j * width + i] = bySubdivision(pixelBox(xs, ys, i, j), depth, single);
  const counts = { undecided: 0 };
  traps.forEach((_, k) => (counts[`root ${k + 1}`] = 0));
  for (const k of kind) k ? counts[`root ${k}`]++ : counts.undecided++;
  return new PixelGrid({
    name: `Newton basins of a degree-${deg} polynomial`,
    width,
    height,
    kind,
    labels: true,
    roots: traps.map((B) => ({ re: [B[0], B[1]], im: [B[2], B[3]], approx: [(B[0] + B[1]) / 2, (B[2] + B[3]) / 2] })),
    unproved,
    counts,
    re: xs[0],
    reHigh: xs[width],
    im: ys[0],
    imHigh: ys[height],
  });
}

// ---- single exact points, with the complex package's rectangles ------------------------------

/**
 * Is the exact point c in the Mandelbrot set? { verdict: "inside" | "outside" |
 * "undecided", reason }. Rational rectangles (rounded outward to `bits` bits)
 * carry the orbit of 0; the cardioid and period-2 disc tests are exact.
 */
export function mandelbrotPoint(c, { maxIter = 500, bits = 128, maxPeriod = 16 } = {}) {
  const z = Complex.from(c);
  const x = z.re;
  const y = z.im;
  const quarter = Q(1, 4);
  const xm = x.subtract(quarter);
  const q = xm.multiply(xm).add(y.multiply(y));
  if (q.multiply(q.add(xm)).compareTo(y.multiply(y).multiply(quarter)) < 0) return Object.freeze({ verdict: "inside", reason: "in the main cardioid (exact test)" });
  const x1 = x.add(ONE);
  if (x1.multiply(x1).add(y.multiply(y)).compareTo(Q(1, 16)) < 0) return Object.freeze({ verdict: "inside", reason: "in the period-2 disc (exact test)" });
  const C = Rect.from(z);
  let Z = Rect.from(0);
  const four = Q(4);
  for (let n = 1; n <= maxIter; n++) {
    Z = Z.mul(Z).add(C).round(bits);
    const nr = Z.normRange();
    if (nr.lo.compareTo(four) > 0) return Object.freeze({ verdict: "outside", reason: `|z_${n}| > 2`, escape: n });
    if ((n & (n - 1)) === 0 && n >= 16 && nr.hi.compareTo(four) < 0) {
      const p = trappedRect(Z, C, bits, maxPeriod);
      if (p) return Object.freeze({ verdict: "inside", reason: `the orbit of 0 enters a box that ${p} step${p === 1 ? "" : "s"} carry into itself (an attracting cycle of period dividing ${p})`, period: p });
    }
  }
  return Object.freeze({ verdict: "undecided", reason: `neither escape nor a trap within ${maxIter} steps` });
}
/** p when a box around Z is carried into itself by p steps of z → z² + C (rational rectangles), else 0. */
function trappedRect(Z, C, bits, maxPeriod) {
  for (const grow of [Q(3, 2), Q(3)]) {
    const w = Z.width.compareTo(Z.height) > 0 ? Z.width : Z.height;
    const r = w.multiply(grow).add(new Rational(1n, 1n << BigInt(bits - 8)));
    const ctr = Z.center;
    const B = new Rect({ lo: ctr.re.subtract(r), hi: ctr.re.add(r) }, { lo: ctr.im.subtract(r), hi: ctr.im.add(r) });
    let W = B;
    for (let p = 1; p <= maxPeriod; p++) {
      W = W.mul(W).add(C).round(bits);
      if (W.normRange().hi.compareTo(Q(16)) > 0) break;
      if (strictlyInside(W, B)) return p;
    }
  }
  return 0;
}
const strictlyInside = (A, B) => A.re.lo.compareTo(B.re.lo) > 0 && A.re.hi.compareTo(B.re.hi) < 0 && A.im.lo.compareTo(B.im.lo) > 0 && A.im.hi.compareTo(B.im.hi) < 0;

/** Is the exact point z0 in the filled Julia set of z → z² + c? { verdict: "inside" | "outside" | "undecided", reason } */
export function juliaPoint(c, z0, { maxIter = 500, bits = 128, maxPeriod = 16 } = {}) {
  const C = Rect.from(Complex.from(c));
  let Z = Rect.from(Complex.from(z0));
  const cn = C.normRange().hi;
  const R2 = cn.compareTo(Q(4)) > 0 ? cn : Q(4);
  for (let n = 0; n <= maxIter; n++) {
    if (n) Z = Z.mul(Z).add(C).round(bits);
    if (Z.normRange().lo.compareTo(R2) > 0) return Object.freeze({ verdict: "outside", reason: `|z_${n}| > max(2, |c|), after which the orbit escapes`, escape: n });
    if (n >= 16 && (n & (n - 1)) === 0) {
      const p = trappedRect(Z, C, bits, maxPeriod);
      if (p) return Object.freeze({ verdict: "inside", reason: `the orbit enters a box that ${p} step${p === 1 ? "" : "s"} carry into itself, so it stays bounded`, period: p });
    }
  }
  return Object.freeze({ verdict: "undecided", reason: `neither escape nor a trap within ${maxIter} steps` });
}
