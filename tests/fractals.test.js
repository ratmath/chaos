// Phases 4 and 5: fractals, L-systems, cellular automata and escape-time pictures.
// Each exact or certified answer is checked against something independent:
// brute force, a float orbit, or a known closed form.
import { describe, test, expect } from "bun:test";
import { Rational } from "@ratmath/core";
import {
  Affine, IFS, attractors, moran, cantor, boxCounting, seeded,
  LSystem, lsystem, lsystems, Q3,
  Elementary, elementary, showRows, Life, life, patterns,
  mandelbrot, julia, newtonBasins, mandelbrotPoint, juliaPoint,
} from "../index.js";

const q = (n, d = 1) => new Rational(BigInt(n), BigInt(d));

describe("iterated function systems", () => {
  test("Sierpinski stages have 3^n exact points inside the triangle", () => {
    const S = attractors.sierpinski();
    const pts = S.stage(5);
    expect(pts.length).toBe(3 ** 5);
    for (const [x, y] of pts) {
      // inside the triangle (0,0), (1,0), (1/2,1)
      expect(y.compareTo(q(0)) >= 0).toBe(true);
      expect(y.compareTo(x.multiply(q(2))) <= 0).toBe(true);
      expect(y.compareTo(q(2).subtract(x.multiply(q(2)))) <= 0).toBe(true);
    }
  });
  test("fixed points are fixed", () => {
    for (const S of [attractors.sierpinski(), attractors.heighway(), attractors.barnsley()])
      S.maps.forEach((m, i) => {
        const p = S.fixedPoint(i);
        const [x, y] = m.apply(p);
        expect(x.equals(p[0]) && y.equals(p[1])).toBe(true);
      });
  });
  test("contraction bounds are below 1 and at least the true ratio", () => {
    for (const k of Object.keys(attractors)) {
      const S = attractors[k]();
      const { bounds, contracting } = S.contraction();
      expect(contracting).toBe(true);
      // a bound is at least |A v| / |v| for sample unit vectors (float check)
      S.maps.forEach((m, i) => {
        for (let t = 0; t < 32; t++) {
          const v = [Math.cos(t / 5), Math.sin(t / 5)];
          const w = [m.a.toNumber() * v[0] + m.b.toNumber() * v[1], m.c.toNumber() * v[0] + m.d.toNumber() * v[1]];
          expect(Math.hypot(...w)).toBeLessThanOrEqual(bounds[i].toNumber() + 1e-12);
        }
      });
    }
  });
  test("similarity dimensions: exact and certified", () => {
    const s = attractors.sierpinski().similarityDimension();
    expect(s.kind).toBe("certified");
    const v = Math.log(3) / Math.log(2);
    expect(s.interval[0].toNumber()).toBeLessThanOrEqual(v);
    expect(s.interval[1].toNumber()).toBeGreaterThanOrEqual(v);
    const h = attractors.heighway().similarityDimension();
    expect(h.kind).toBe("exact");
    expect(h.value.equals(q(2))).toBe(true);
    // Moran with unequal ratios: 1/2 and 1/4 give (1/2)^s + (1/4)^s = 1, s = log2 of the golden ratio
    const m = moran(["1/2", "1/4"]);
    const g = Math.log2((1 + Math.sqrt(5)) / 2);
    expect(m.interval[0].toNumber()).toBeLessThanOrEqual(g + 1e-12);
    expect(m.interval[1].toNumber()).toBeGreaterThanOrEqual(g - 1e-12);
    expect(m.interval[1].toNumber() - m.interval[0].toNumber()).toBeLessThan(1e-9);
  });
  test("chaos game stays on the attractor and box counting is near the dimension", () => {
    const S = attractors.sierpinski();
    const g = S.chaosGame(20000, { seed: 7 });
    expect(g.kind).toBe("estimated");
    expect(S.chaosGame(5, { seed: 7 }).points).toEqual(g.points.slice(0, 10));
    const bc = boxCounting(g, { scales: [2, 3, 4, 5, 6] });
    expect(bc.kind).toBe("estimated");
    expect(Math.abs(bc.dimension - Math.log2(3))).toBeLessThan(0.1);
  });
  test("seeded is repeatable and in [0, 1)", () => {
    const a = seeded(3);
    const b = seeded(3);
    for (let i = 0; i < 100; i++) {
      const x = a();
      expect(x).toBe(b());
      expect(x >= 0 && x < 1).toBe(true);
    }
  });
  test("bad input is refused", () => {
    expect(() => new IFS([])).toThrow();
    expect(() => new IFS([{ a: "1/2" }], { probabilities: ["1/2"] })).toThrow();
    expect(() => moran(["2"])).toThrow();
    expect(() => attractors.barnsley().similarityDimension()).toThrow();
  });
});

describe("Cantor sets", () => {
  // brute force: follow the base-3 digits of n/d (with the alternative of a final 1 = 0222…)
  const bruteInside = (n, d, steps = 200) => {
    // x in C iff some path of allowed digits survives; check by BFS over numerators (finite)
    let frontier = new Set([n]);
    for (let s = 0; s < steps && frontier.size; s++) {
      const next = new Set();
      for (const a of frontier) for (const dig of [0, 2]) {
        const b = 3 * a - dig * d;
        if (b >= 0 && b <= d) next.add(b);
      }
      frontier = next;
    }
    return frontier.size > 0;
  };
  test("membership agrees with brute force for every n/d with d ≤ 40", () => {
    const C = cantor();
    for (let d = 1; d <= 40; d++)
      for (let n = 0; n <= d; n++) expect(C.contains(q(n, d)).inside).toBe(bruteInside(n, d));
  });
  test("famous members and non-members", () => {
    const C = cantor();
    expect(C.contains("1/4")).toMatchObject({ inside: true, expansion: "0.#02" });
    expect(C.contains("3/4").inside).toBe(true);
    expect(C.contains("1/3").inside).toBe(true);
    expect(C.contains("1/2").inside).toBe(false);
    expect(C.contains("2").inside).toBe(false);
  });
  test("stages, measure and dimension", () => {
    const C = cantor();
    const st = C.stage(4);
    expect(st.length).toBe(16);
    const total = st.reduce((s, [a, b]) => s.add(b.subtract(a)), q(0));
    expect(total.equals(C.measure(4))).toBe(true);
    expect(C.measure(4).equals(q(16, 81))).toBe(true);
    const dim = C.dimension();
    expect(dim.kind).toBe("certified");
    const v = Math.log(2) / Math.log(3);
    expect(dim.interval[0].toNumber() <= v && v <= dim.interval[1].toNumber()).toBe(true);
    // base 4 with two digits: log 2/log 4 = 1/2 exactly
    expect(cantor({ base: 4, digits: [0, 3] }).dimension().value.equals(q(1, 2))).toBe(true);
  });
});

describe("L-systems", () => {
  test("counts by matrix power match the written string", () => {
    for (const k of ["koch", "arrowhead", "dragon", "hilbert", "plant"]) {
      const L = lsystems[k]();
      for (const n of [0, 1, 2, 3, 4]) {
        const s = L.generate(n);
        const { counts, length } = L.counts(n);
        expect(length).toBe(BigInt(s.length));
        for (const [sym, c] of Object.entries(counts)) expect(c).toBe(BigInt([...s].filter((x) => x === sym).length));
      }
    }
  });
  test("Koch counts in closed form for large n", () => {
    const { counts } = lsystems.koch().counts(100);
    expect(counts.F).toBe(4n ** 100n);
  });
  test("exact turtle points agree with floats, and closed curves close", () => {
    const S = lsystems.snowflake().turtle(3);
    expect(S.kind).toBe("exact");
    for (let i = 0; i < S.segments.length; i++) {
      const [x1, y1, x2, y2] = S.segments[i];
      const [X1, Y1, X2, Y2] = S.exactSegments[i];
      expect(Math.abs(X1.toNumber() - x1) + Math.abs(Y2.toNumber() - y2) + Math.abs(X2.toNumber() - x2) + Math.abs(Y1.toNumber() - y1)).toBeLessThan(1e-9);
    }
    const [X, Y] = S.end;
    expect(X.equals(new Q3()) && Y.equals(new Q3())).toBe(true);
    expect(S.segments.length).toBe(3 * 4 ** 3);
  });
  test("the Hilbert curve visits every cell of the grid once", () => {
    const H = lsystems.hilbert().turtle(3);
    const seen = new Set(["0,0"]);
    for (const [, , x, y] of H.segments) seen.add(`${Math.round(x)},${Math.round(y)}`);
    expect(seen.size).toBe(64);
    expect(H.segments.length).toBe(63);
  });
  test("Q3 arithmetic", () => {
    const r3 = new Q3(q(0), q(1));
    expect(r3.mul(r3).equals(new Q3(q(3)))).toBe(true);
    expect(new Q3(q(1, 2), q(-1)).toString()).toBe("1/2 - √3");
  });
  test("bad input is refused", () => {
    expect(() => lsystem({ axiom: "" })).toThrow();
    expect(() => lsystem({ axiom: "F", rules: { FF: "F" } })).toThrow();
    expect(() => lsystem({ axiom: "]", rules: {} }).turtle(0)).toThrow();
    expect(() => lsystems.koch().generate(20)).toThrow();
  });
});

describe("elementary cellular automata", () => {
  test("rule 90 from one cell is Pascal's triangle mod 2", () => {
    const rows = elementary(90).fromCell(31);
    for (let t = 0; t <= 31; t++)
      for (let i = 0; i < rows[t].length; i++) {
        const k = i - 31 + t; // offset from the left edge of the light cone, in steps of 2
        const expected = k >= 0 && k <= 2 * t && k % 2 === 0 && ((t & (k / 2)) === k / 2) ? 1 : 0;
        expect(rows[t][i]).toBe(expected);
      }
  });
  test("tables, mirrors and complements against brute force", () => {
    for (let r = 0; r < 256; r++) {
      const E = new Elementary(r);
      const M = new Elementary(E.mirror);
      const C = new Elementary(E.complement);
      for (let n = 0; n < 8; n++) {
        const [l, c, rr] = [(n >> 2) & 1, (n >> 1) & 1, n & 1];
        const f = (R, a, b, d) => (R.rule >> ((a << 2) | (b << 1) | d)) & 1;
        expect(f(M, rr, c, l)).toBe(f(E, l, c, rr));
        expect(f(C, 1 - l, 1 - c, 1 - rr)).toBe(1 - f(E, l, c, rr));
      }
    }
    expect(elementary(110).equivalents).toEqual([110, 124, 137, 193]);
    expect(elementary(30).table["100"]).toBe(1);
  });
  test("additive rules are exactly the eight linear ones", () => {
    const found = [];
    for (let r = 0; r < 256; r++) if (elementary(r).additive) found.push(r);
    expect(found).toEqual([0, 60, 90, 102, 150, 170, 204, 240]);
  });
  test("cycles of a periodic row", () => {
    // rule 184 (traffic) on a ring of alternating cells just shifts: period 2
    expect(elementary(184).cycle("10101010")).toEqual({ transient: 0, period: 2 });
    const c = elementary(30).cycle("0000000001");
    // check by stepping
    const rows = elementary(30).run("0000000001", c.transient + c.period);
    expect(rows[c.transient].join("")).toBe(rows[c.transient + c.period].join(""));
  });
  test("showRows", () => {
    expect(showRows(elementary(90).fromCell(1), "#", ".")).toBe(".#.\n#.#");
    expect(() => elementary(256)).toThrow();
  });
});

describe("Life", () => {
  // independent brute force on a bounded torus big enough not to wrap
  const brute = (cells, steps, N = 64) => {
    let g = Array.from({ length: N }, () => new Uint8Array(N));
    for (const [x, y] of cells) g[y + 20][x + 20] = 1;
    for (let s = 0; s < steps; s++) {
      const h = Array.from({ length: N }, () => new Uint8Array(N));
      for (let y = 1; y < N - 1; y++)
        for (let x = 1; x < N - 1; x++) {
          let c = 0;
          for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) if (dx || dy) c += g[y + dy][x + dx];
          h[y][x] = c === 3 || (c === 2 && g[y][x]) ? 1 : 0;
        }
      g = h;
    }
    const out = [];
    for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) if (g[y][x]) out.push(`${x - 20},${y - 20}`);
    return out.sort();
  };
  test("agrees with brute force", () => {
    for (const p of ["rpentomino", "glider", "lwss", "toad"]) {
      const L = life(patterns[p]);
      const start = [...L.cells].map((k) => k.split(",").map(Number));
      expect([...L.run(12).cells].sort()).toEqual(brute(start, 12));
    }
  });
  test("classify", () => {
    expect(life(patterns.block).classify()).toMatchObject({ kind: "still life", period: 1 });
    expect(life(patterns.blinker).classify()).toMatchObject({ kind: "oscillator", period: 2 });
    expect(life(patterns.pulsar).classify()).toMatchObject({ kind: "oscillator", period: 3 });
    expect(life(patterns.glider).classify()).toMatchObject({ kind: "spaceship", period: 4, speed: "c/4" });
    expect(life(patterns.lwss).classify()).toMatchObject({ kind: "spaceship", period: 4, speed: "c/2" });
    expect(life("O").classify()).toMatchObject({ kind: "dies", at: 1 });
  });
  test("other rules: HighLife B36/S23 replicator grows", () => {
    const L = life("OOO\nO..\nO..", { rule: "B36/S23" });
    expect(L.run(12).population).toBeGreaterThan(0);
    expect(() => life("O", { rule: "3/23" })).toThrow();
  });
});

describe("escape-time pictures", () => {
  const sample = (G, i, j, re0, re1, im0, im1) => {
    // the pixel's centre as a float point
    return [re0 + ((i + 0.5) * (re1 - re0)) / G.width, im0 + ((j + 0.5) * (im1 - im0)) / G.height];
  };
  test("Mandelbrot pixels agree with float orbits at their centres", () => {
    const G = mandelbrot({ width: 48, height: 24 });
    expect(G.counts.inside).toBeGreaterThan(100);
    expect(G.counts.outside).toBeGreaterThan(400);
    for (let j = 0; j < G.height; j++)
      for (let i = 0; i < G.width; i++) {
        const k = G.at(i, j);
        if (!k) continue;
        const [cr, ci] = sample(G, i, j, -2, 0.5, -1.25, 1.25);
        let [x, y] = [0, 0];
        let n = 0;
        for (; n < 2000 && x * x + y * y <= 4; n++) [x, y] = [x * x - y * y + cr, 2 * x * y + ci];
        if (k === 1) expect(x * x + y * y).toBeLessThanOrEqual(4);
        else expect(n).toBeLessThan(2000);
      }
    expect(G.toText().split("\n").length).toBe(24);
  });
  test("Julia pixels agree with float orbits", () => {
    const G = julia("-1", { width: 40, height: 20, re: ["-2", "2"], im: ["-1", "1"] });
    expect(G.counts.inside).toBeGreaterThan(20);
    for (let j = 0; j < G.height; j++)
      for (let i = 0; i < G.width; i++) {
        const k = G.at(i, j);
        if (!k) continue;
        let [x, y] = sample(G, i, j, -2, 2, -1, 1);
        let n = 0;
        for (; n < 2000 && x * x + y * y <= 4; n++) [x, y] = [x * x - y * y - 1, 2 * x * y];
        if (k === 1) expect(x * x + y * y).toBeLessThanOrEqual(4);
        else expect(n).toBeLessThan(2000);
      }
  });
  test("Newton basins: labelled pixels converge to their root", () => {
    const G = newtonBasins(["1", "0", "0", "-1"], { width: 30, height: 20 });
    expect(G.roots.length).toBe(3);
    let labelled = 0;
    for (let j = 0; j < G.height; j++)
      for (let i = 0; i < G.width; i++) {
        const k = G.at(i, j);
        if (!k) continue;
        labelled++;
        let [x, y] = sample(G, i, j, -2, 2, -1.5, 1.5);
        for (let n = 0; n < 200; n++) {
          // z − (z³ − 1)/(3z²)
          const [a, b] = [x * x - y * y, 2 * x * y]; // z²
          const [c, d] = [a * x - b * y - 1, a * y + b * x]; // z³ − 1
          const den = 3 * (a * a + b * b);
          [x, y] = [x - (c * a + d * b) / den, y - (d * a - c * b) / den];
        }
        const [u, w] = G.roots[k - 1].approx;
        expect(Math.hypot(x - u, y - w)).toBeLessThan(1e-6);
      }
    expect(labelled).toBeGreaterThan(200);
    // the roots are the cube roots of 1
    for (const r of G.roots) expect(Math.abs(Math.hypot(...r.approx) - 1)).toBeLessThan(1e-6);
  });
  test("single exact points", () => {
    expect(mandelbrotPoint("0").verdict).toBe("inside");
    expect(mandelbrotPoint("-1").verdict).toBe("inside");
    expect(mandelbrotPoint("-13/10")).toMatchObject({ verdict: "inside", period: 4 });
    expect(mandelbrotPoint("1").verdict).toBe("outside");
    expect(mandelbrotPoint("1/4", { maxIter: 200 }).verdict).toBe("undecided");
    expect(juliaPoint("-1", "0").verdict).toBe("inside");
    expect(juliaPoint("-1", "2").verdict).toBe("outside");
  });
  test("bad grids are refused", () => {
    expect(() => mandelbrot({ width: 0 })).toThrow();
    expect(() => mandelbrot({ re: ["1", "0"] })).toThrow();
    expect(() => newtonBasins(["0", "3"])).toThrow();
  });
});
