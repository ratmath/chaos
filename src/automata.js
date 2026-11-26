// Cellular automata: the discrete end of the subject. Everything here is
// exact, since cells are 0 or 1.
//
//   elementary(rule)   Wolfram's 256 rules of one dimension, three neighbours
//   life(cells)        Conway's Life, or any rule written "B3/S23", on an
//                      unbounded plane (only live cells are stored)
import { refuse } from "./tree.js";

export class Elementary {
  constructor(rule) {
    if (!Number.isInteger(rule) || rule < 0 || rule > 255) refuse("bad-argument", "an elementary rule is a whole number from 0 to 255");
    this.rule = rule;
    Object.freeze(this);
  }
  /** The new cell for each neighbourhood 111, 110, …, 000 (Wolfram's order). */
  get table() {
    return Object.fromEntries([7, 6, 5, 4, 3, 2, 1, 0].map((n) => [n.toString(2).padStart(3, "0"), (this.rule >> n) & 1]));
  }
  /** The rule seen in a mirror (left and right swapped). */
  get mirror() {
    let r = 0;
    for (let n = 0; n < 8; n++) {
      const m = ((n & 1) << 2) | (n & 2) | ((n >> 2) & 1);
      if ((this.rule >> n) & 1) r |= 1 << m;
    }
    return r;
  }
  /** The rule with 0 and 1 swapped. */
  get complement() {
    let r = 0;
    for (let n = 0; n < 8; n++) if (!((this.rule >> (7 - n)) & 1)) r |= 1 << n;
    return r;
  }
  /** The rules that behave the same up to mirror and complement (at most four). */
  get equivalents() {
    const m = new Elementary(this.mirror);
    return [...new Set([this.rule, m.rule, this.complement, m.complement])].sort((a, b) => a - b);
  }
  /** Additive: the new cell is a sum mod 2 of neighbours (rules 0, 60, 90, 102, 150, 170, 204, 240). */
  get additive() {
    if (this.rule & 1) return false;
    for (let a = 0; a < 8; a++) for (let b = 0; b < 8; b++) if ((((this.rule >> a) ^ (this.rule >> b)) & 1) !== ((this.rule >> (a ^ b)) & 1)) return false;
    return true;
  }
  /** One step of a row (Uint8Array of 0s and 1s). Boundary "periodic" or "zero". */
  step(row, { boundary = "periodic" } = {}) {
    const n = row.length;
    const out = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
      const l = i > 0 ? row[i - 1] : boundary === "periodic" ? row[n - 1] : 0;
      const r = i < n - 1 ? row[i + 1] : boundary === "periodic" ? row[0] : 0;
      out[i] = (this.rule >> ((l << 2) | (row[i] << 1) | r)) & 1;
    }
    return out;
  }
  /** The rows from an initial row (an array, a Uint8Array or a string of 0/1), `steps` steps on. */
  run(initial, steps, options) {
    let row = typeof initial === "string" ? Uint8Array.from(initial, (c) => (c === "1" ? 1 : 0)) : Uint8Array.from(initial);
    const rows = [row];
    for (let t = 0; t < steps; t++) rows.push((row = this.step(row, options)));
    return rows;
  }
  /** From one live cell in the middle of a row wide enough that the edges are never reached. */
  fromCell(steps) {
    const w = 2 * steps + 1;
    const row = new Uint8Array(w);
    row[steps] = 1;
    return this.run(row, steps, { boundary: "zero" });
  }
  /**
   * The eventual behaviour of a periodic row: the step at which a row first
   * repeats and the period. Exact, since there are finitely many rows.
   */
  cycle(initial, { limit = 100000 } = {}) {
    let row = typeof initial === "string" ? Uint8Array.from(initial, (c) => (c === "1" ? 1 : 0)) : Uint8Array.from(initial);
    const seen = new Map();
    for (let t = 0; t <= limit; t++) {
      const k = row.join("");
      if (seen.has(k)) return Object.freeze({ transient: seen.get(k), period: t - seen.get(k) });
      seen.set(k, t);
      row = this.step(row);
    }
    refuse("too-large", `no repeat within ${limit} steps`);
  }
  toString() {
    return `rule ${this.rule}`;
  }
}
export const elementary = (rule) => new Elementary(rule);
/** Rows as text, one character per cell. */
export const showRows = (rows, on = "█", off = "·") => rows.map((r) => Array.from(r, (c) => (c ? on : off)).join("")).join("\n");

// ---- Life and its relatives ----------------------------------------------------------------

const parseRule = (rule) => {
  const m = /^B([0-8]*)\/S([0-8]*)$/i.exec(rule.replace(/\s/g, ""));
  if (!m) refuse("bad-argument", `write the rule as "B3/S23"`);
  return { born: new Set([...m[1]].map(Number)), survive: new Set([...m[2]].map(Number)) };
};
const K = (x, y) => `${x},${y}`;
const unK = (k) => k.split(",").map(Number);

/** A pattern on the unbounded plane: the set of live cells. */
export class Life {
  constructor(cells, { rule = "B3/S23" } = {}) {
    this.rule = rule;
    this._r = parseRule(rule);
    this.cells = new Set();
    if (typeof cells === "string") {
      cells.split("\n").forEach((line, y) => [...line].forEach((c, x) => {
        if (c === "O" || c === "#" || c === "1" || c === "█" || c === "*") this.cells.add(K(x, y));
      }));
    } else for (const [x, y] of cells) this.cells.add(K(x, y));
  }
  get population() {
    return this.cells.size;
  }
  /** The next generation, as a new Life. */
  step() {
    const counts = new Map();
    for (const k of this.cells) {
      const [x, y] = unK(k);
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) if (dx || dy) {
        const n = K(x + dx, y + dy);
        counts.set(n, (counts.get(n) ?? 0) + 1);
      }
    }
    const next = [];
    for (const [k, c] of counts) if (this.cells.has(k) ? this._r.survive.has(c) : this._r.born.has(c)) next.push(unK(k));
    return new Life(next, { rule: this.rule });
  }
  run(n) {
    let g = this;
    for (let i = 0; i < n; i++) g = g.step();
    return g;
  }
  bounds() {
    if (!this.cells.size) return null;
    let [x0, y0, x1, y1] = [Infinity, Infinity, -Infinity, -Infinity];
    for (const k of this.cells) {
      const [x, y] = unK(k);
      x0 = Math.min(x0, x);
      y0 = Math.min(y0, y);
      x1 = Math.max(x1, x);
      y1 = Math.max(y1, y);
    }
    return { x0, y0, x1, y1 };
  }
  /** The shape moved so its corner is at the origin, as text (for comparing shapes). */
  normalized() {
    const b = this.bounds();
    if (!b) return { shape: "", dx: 0, dy: 0 };
    const pts = [...this.cells].map(unK).map(([x, y]) => [x - b.x0, y - b.y0]).sort((p, q) => p[1] - q[1] || p[0] - q[0]);
    return { shape: pts.map(([x, y]) => K(x, y)).join(";"), dx: b.x0, dy: b.y0 };
  }
  /**
   * What the pattern does, found by running it: "dies", "still life",
   * "oscillator" (with its period), "spaceship" (period and displacement per
   * period) or "unknown" within the limit. Exact: a repeated shape is
   * compared cell by cell.
   */
  classify({ limit = 1000 } = {}) {
    const seen = new Map();
    let g = this;
    for (let t = 0; t <= limit; t++) {
      if (!g.population) return Object.freeze({ kind: "dies", at: t });
      const { shape, dx, dy } = g.normalized();
      if (seen.has(shape)) {
        const s = seen.get(shape);
        const period = t - s.t;
        const mx = dx - s.dx;
        const my = dy - s.dy;
        if (mx || my) {
          const d = Math.max(Math.abs(mx), Math.abs(my));
          let [a, b] = [d, period];
          while (b) [a, b] = [b, a % b];
          const num = d / a;
          const den = period / a;
          const speed = `${num === 1 ? "" : num}c${den === 1 ? "" : `/${den}`}`;
          return Object.freeze({ kind: "spaceship", period, displacement: [mx, my], from: s.t, speed });
        }
        return Object.freeze({ kind: period === 1 ? "still life" : "oscillator", period, from: s.t });
      }
      seen.set(shape, { t, dx, dy });
      g = g.step();
    }
    return Object.freeze({ kind: "unknown", limit });
  }
  toString(on = "█", off = "·") {
    const b = this.bounds();
    if (!b) return "(empty)";
    const lines = [];
    for (let y = b.y0; y <= b.y1; y++) {
      let s = "";
      for (let x = b.x0; x <= b.x1; x++) s += this.cells.has(K(x, y)) ? on : off;
      lines.push(s);
    }
    return lines.join("\n");
  }
}
export const life = (cells, options) => new Life(cells, options);

/** Patterns from the Life lexicon. */
export const patterns = Object.freeze({
  block: "OO\nOO",
  blinker: "OOO",
  toad: ".OOO\nOOO.",
  beacon: "OO..\nOO..\n..OO\n..OO",
  glider: ".O.\n..O\nOOO",
  lwss: ".O..O\nO....\nO...O\nOOOO.",
  pulsar: "..OOO...OOO..\n.............\nO....O.O....O\nO....O.O....O\nO....O.O....O\n..OOO...OOO..\n.............\n..OOO...OOO..\nO....O.O....O\nO....O.O....O\nO....O.O....O\n.............\n..OOO...OOO..",
  rpentomino: ".OO\nOO.\n.O.",
  gosperGun: "........................O...........\n......................O.O...........\n............OO......OO............OO\n...........O...O....OO............OO\nOO........O.....O...OO..............\nOO........O...O.OO....O.O...........\n..........O.....O.......O...........\n...........O...O....................\n............OO......................",
});
