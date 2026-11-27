import { describe, expect, test } from "bun:test";
import { henon, standardMap, planeMap, catMap, torusMap, horseshoe, ChaosError } from "../index.js";

const mid = (iv) => (iv.low.toNumber() + iv.high.toNumber()) / 2;
const width = (iv) => iv.high.subtract(iv.low).toNumber();
const region = { x: [-2, 2], y: [-1, 1] };

describe("the Hénon map", () => {
  const H = henon();
  test("exact orbits, and the step where doubles part from them", () => {
    expect(H.orbit([0, 0], 3).map((p) => p.map(String))).toEqual([["0", "0"], ["1", "0"], ["-2/5", "3/10"], ["269/250", "-3/25"]]);
    expect(() => H.orbit([0.5, 0], 2)).toThrow(ChaosError);
    const exact = H.orbit([0, 0], 8);
    const enc = H.enclosedOrbit([0, 0], 8).points;
    enc.forEach((p, i) => p.forEach((iv, k) => expect(iv.low.compareTo(exact[i][k]) <= 0 && exact[i][k].compareTo(iv.high) <= 0).toBe(true)));
    const c = H.compare([0, 0], 150);
    expect(c.parted).toBeGreaterThan(60);
    expect(c.parted).toBeLessThan(130);
    expect(c.rows[10].digits).toBeGreaterThan(12);
  });
  test("the Jacobian determinant is −b everywhere", () => {
    expect(String(H.determinant)).toBe("-3/10");
    expect(H.areaPreserving).toBe(false);
    expect(H.jacobian()).toEqual([["-14/5*x", "1"], ["3/10", "0"]]);
    expect(henon("1.4", 1).areaPreserving).toBe(true);
  });
  test("periodic points are proven and counted", () => {
    // the known numbers of orbits of the classical map: 1, 1, 0, 1, 0, 2, 4 for periods 1 to 7 (plus the second fixed point)
    const counts = [1, 2, 3, 4, 5, 6, 7].map((n) => {
      const pts = H.periodicPoints(n, region);
      expect(pts.complete).toBe(true);
      return pts.length;
    });
    expect(counts).toEqual([2, 2, 0, 4, 0, 12, 28]);
    const [p, q] = H.fixedPoints(region);
    expect(mid(q.point[0])).toBeCloseTo((-0.7 + Math.sqrt(6.09)) / 2.8, 12);
    expect(width(q.point[0])).toBeLessThan(1e-12);
    expect(p.stability).toBe("saddle");
    expect(q.stability).toBe("saddle");
    expect(q.multipliers.map((e) => mid(e.re))[0]).toBeCloseTo(-1.9237388581534, 9);
    // the unstable direction is an eigenvector: J·u is a multiple of u
    const x = mid(q.point[0]);
    const u = q.unstable.map(mid);
    const Ju = [-2.8 * x * u[0] + u[1], 0.3 * u[0]];
    expect(Ju[0] * u[1] - Ju[1] * u[0]).toBeCloseTo(0, 9);
    expect(H.periodicPoints(2, region, { all: true }).length).toBe(4);
  });
  test("cycles group the points into orbits", () => {
    const cy = H.cycles(7, region);
    expect(cy.length).toBe(4);
    expect(cy.every((c) => c.whole && c.points.length === 7 && c.stability === "saddle")).toBe(true);
  });
  test("low periods are exact algebraic numbers", () => {
    const fixed = H.exactPeriodicPoints(1);
    expect(fixed.map((p) => p.text)).toEqual([["(−7 − √609)/28", "(−21 − 3√609)/280"], ["(−7 + √609)/28", "(−21 + 3√609)/280"]]);
    const two = H.exactPeriodicPoints(2);
    expect(two.map((p) => p.text[0])).toEqual(["(7 − √413)/28", "(7 + √413)/28"]);
    expect(two[0].stability).toBe("saddle");
    expect(() => H.exactPeriodicPoints(3)).toThrow(/first degree/);
    // they agree with the proven enclosures
    const enc = H.periodicPoints(2, region);
    two.forEach((p, i) => expect(p.x.toNumber()).toBeCloseTo(mid(enc[i].point[0]), 12));
  });
  test("Lyapunov exponents: an estimate, with the sum known exactly", () => {
    const l = H.lyapunov([0, 0], 20000);
    expect(l.kind).toBe("estimated");
    expect(l.exponents[0]).toBeGreaterThan(0.39);
    expect(l.exponents[0]).toBeLessThan(0.45);
    expect(mid(l.sum)).toBeCloseTo(Math.log(0.3), 12);
    expect(l.exponents[0] + l.exponents[1]).toBeCloseTo(Math.log(0.3), 6);
  });
});

describe("other maps from formulas", () => {
  test("the standard map keeps areas, with one saddle and one elliptic point", () => {
    const S = standardMap("0.5");
    expect(String(S.determinant)).toBe("1");
    expect(S.areaPreserving).toBe(true);
    const pts = S.fixedPoints({ x: ["-0.25", "0.75"], y: ["-0.5", "0.5"] });
    expect(pts.complete).toBe(true);
    expect(pts.map((p) => p.stability)).toEqual(["saddle", "elliptic"]);
    expect(mid(pts[1].point[0])).toBeCloseTo(0.5, 12);
    expect(() => S.at([0, 0])).toThrow(/enclose/);
    expect(width(S.enclose(["1/4", "1/10"])[0])).toBeLessThan(1e-20);
  });
  test("an attracting fixed point of a linear map", () => {
    const M = planeMap(["x/2 + y/4 + 1", "x/4 - y/3"]);
    const [p] = M.fixedPoints({ x: [-10, 10], y: [-10, 10] });
    expect(p.stability).toBe("attracting");
    expect(mid(p.point[0])).toBeCloseTo(64 / 29, 12);
    expect(String(M.determinant)).toBe("-11/48");
  });
});

describe("the cat map", () => {
  const cat = catMap();
  test("every rational point is periodic", () => {
    expect(String(cat)).toBe("(x, y) → (2x + y, x + y) mod 1");
    expect(cat.at(["1/5", "2/5"]).map(String)).toEqual(["4/5", "3/5"]);
    expect(cat.period(["1/5", "2/5"]).period).toBe(2);
    expect(cat.period(["3/101", "7/101"])).toMatchObject({ period: 25, preperiod: 0 });
    const per = cat.periods(12);
    expect([...per.values()].reduce((s, x) => s + x, 0)).toBe(144);
    expect(per.get(1)).toBe(1);
  });
  test("fixed points of iterates are counted by a determinant", () => {
    expect([1, 2, 3, 4, 5].map((n) => cat.fixedPointCount(n))).toEqual([1n, 5n, 16n, 45n, 121n]);
    // the points of period dividing 2 with denominator 5: there are 5 fixed points of the square in all
    expect(cat.hyperbolic).toBe(true);
    expect(torusMap([[0, 1], [-1, 0]]).hyperbolic).toBe(false);
  });
  test("eigenvalues, directions and entropy are exact", () => {
    const e = cat.eigen();
    expect(e.map((x) => x.text)).toEqual(["(3 − √5)/2", "(3 + √5)/2"]);
    expect(e.map((x) => x.directionText)).toEqual(["(−1 − √5)/2", "(−1 + √5)/2"]);
    expect(e.map((x) => x.expanding)).toEqual([false, true]);
    const h = cat.entropy();
    expect(h.text).toBe("log((3 + √5)/2)");
    expect(mid(h.value)).toBeCloseTo(Math.log((3 + Math.sqrt(5)) / 2), 14);
  });
  test("doubles lose a periodic orbit within a few dozen steps", () => {
    const k = cat.floatPartsAt(["3/101", "7/101"]);
    expect(k).toBeGreaterThan(20);
    expect(k).toBeLessThan(60);
  });
});

describe("the horseshoe", () => {
  const hs = horseshoe();
  test("periodic points from their words", () => {
    const p = hs.periodicPoint("011");
    expect(p.map(String)).toEqual(["9/13", "3/13"]);
    expect(hs.itinerary(p, 9)).toBe("011011011");
    let q = p;
    for (let i = 0; i < 3; i++) q = hs.at(q);
    expect(q.map(String)).toEqual(p.map(String));
    expect(hs.at(["1/2", "1/2"])).toBe(null);
    expect(hs.itinerary(["1/2", "1/4"], 5)).toBe("01111");
    expect(hs.itinerary(["1/2", "1/6"], 9)).toBe("0");
  });
  test("counts", () => {
    const pts = hs.periodicPoints(4);
    expect(pts.length).toBe(16);
    expect(pts.filter((x) => x.leastPeriod === 4).length).toBe(12);
    for (const x of pts) expect(hs.itinerary(x.point, 8)).toBe(x.word + x.word);
    expect([1, 2, 3, 4, 5, 6].map((n) => hs.orbitCount(n))).toEqual([2, 1, 2, 3, 6, 9]);
    expect(hs.strips(2).horizontal.map((s) => s.map(String))).toEqual([["0", "1/9"], ["2/9", "1/3"], ["2/3", "7/9"], ["8/9", "1"]]);
  });
});
