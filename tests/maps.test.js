import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import { map, logistic, tent, quadratic, doubling, ChaosError } from "../index.js";

const mid = (iv) => (iv.low.toNumber() + iv.high.toNumber()) / 2;

describe("orbits", () => {
  test("exact orbits are the true orbit", () => {
    expect(logistic(4).orbit("1/3", 3).map(String)).toEqual(["1/3", "8/9", "32/81", "6272/6561"]);
    expect(tent(2).orbit("1/7", 4).map(String)).toEqual(["1/7", "2/7", "4/7", "6/7", "2/7"]);
    expect(logistic("7/2").at("1/2").toString()).toBe("7/8");
  });
  test("a floating-point number is refused as an exact input", () => {
    expect(() => logistic(4).orbit(0.1, 3)).toThrow(ChaosError);
    expect(() => logistic(4).orbit("1/3", 60)).toThrow(/enclosedOrbit/);
  });
  test("enclosed orbits contain the exact orbit", () => {
    const L = logistic(4);
    const exact = L.orbit("1/3", 12);
    const enc = L.enclosedOrbit("1/3", 12, { digits: 15 });
    enc.terms.forEach((t, i) => {
      expect(t.low.compareTo(exact[i]) <= 0 && exact[i].compareTo(t.high) <= 0).toBe(true);
    });
    expect(enc.width.toNumber()).toBeLessThan(1e-15);
    const long = L.enclosedOrbit("1/10", 200, { digits: 10 });
    expect(long.terms[200].high.subtract(long.terms[200].low).toNumber()).toBeLessThan(1e-10);
  });
  test("the float orbit parts from the true one", () => {
    const c = logistic(4).compare("1/10", 80);
    expect(c.parted).toBeGreaterThan(40);
    expect(c.parted).toBeLessThan(60);
    expect(c.rows[10].digits).toBeGreaterThan(10);
    expect(c.rows[70].digits).toBeLessThan(1);
  });
  test("an orbit that escapes is reported", () => {
    expect(() => quadratic(1).enclosedOrbit(1, 60)).toThrow(/leaves every bounded range/);
  });
  test("the doubling map on rationals", () => {
    const d = doubling("1/10");
    expect(d.cycle.map(String)).toEqual(["1/5", "2/5", "4/5", "3/5"]);
    expect(d.preperiod).toBe(1);
    expect(d.period).toBe(4);
    expect(d.floatReachesZeroAt).toBeGreaterThan(50);
    expect(doubling("1/8").cycle.map(String)).toEqual(["0"]);
    expect(doubling("1/7", { base: 10 }).period).toBe(6);
  });
  test("cobweb points", () => {
    const pts = logistic(2).cobweb("1/4", 2);
    expect(pts).toEqual([[0.25, 0], [0.25, 0.375], [0.375, 0.375], [0.375, 0.46875], [0.46875, 0.46875]]);
  });
});

describe("periodic points", () => {
  test("fixed points and their stability", () => {
    const f = logistic(4).fixedPoints();
    expect(f.map((p) => p.text)).toEqual(["0", "3/4"]);
    expect(f.map((p) => p.stability)).toEqual(["repelling", "repelling"]);
    expect(mid(f[1].multiplier)).toBeCloseTo(-2, 12);
  });
  test("multipliers of exactly 0 and -1 are recognised", () => {
    expect(logistic(2).fixedPoints()[1].stability).toBe("superattracting");
    const n = logistic(3).fixedPoints()[1];
    expect(n.stability).toBe("neutral");
    expect(n.exactMultiplier.equals(new Rational(-1n))).toBe(true);
  });
  test("a neutral point is classified by its higher-order terms", () => {
    // r = 3: x = 2/3 + u gives the second iterate u − 18u³ − 27u⁴, which attracts
    expect(logistic(3).fixedPoints()[1].nonlinear).toBe("attracting");
    expect(logistic(1).fixedPoints()[0].nonlinear).toBe("semistable"); // transcritical: x − x² attracts only from the right
    expect(map("x - x^3").fixedPoints()[0].nonlinear).toBe("attracting");
    expect(map("x + x^3").fixedPoints()[0].nonlinear).toBe("repelling");
    expect(map("-x - x^3").fixedPoints()[0].nonlinear).toBe("repelling"); // second iterate u + 2u³ + …
    expect(map("-x + x^3").fixedPoints().find((p) => p.text === "0").nonlinear).toBe("attracting");
    expect(logistic(4).fixedPoints()[0].nonlinear).toBeUndefined();
  });
  test("least period 2 excludes the fixed points", () => {
    const p = logistic(4).periodicPoints(2);
    expect(p.map((x) => x.text)).toEqual(["(5 − √5)/8", "(5 + √5)/8"]);
    expect(logistic(2).periodicPoints(2)).toEqual([]);
  });
  test("cycles of period 3 of the full logistic map", () => {
    const c = logistic(4).cycles(3);
    expect(c.length).toBe(2);
    expect(c.every((k) => k.points.length === 3 && k.stability === "repelling")).toBe(true);
    expect(mid(c[0].lyapunov)).toBeCloseTo(Math.log(2), 10);
  });
  test("an attracting 2-cycle", () => {
    const [c] = logistic("3.2").cycles(2);
    expect(c.stability).toBe("attracting");
    expect(mid(c.multiplier)).toBeCloseTo(0.16, 12);
    expect(quadratic(-1).cycles(2)[0].lyapunov).toBe(-Infinity);
  });
  test("certified periodic points of a map that is not a polynomial", () => {
    const s = map("sin(pi*x)");
    expect(() => s.periodicPoints(1)).toThrow(/polynomial/);
    const pts = s.periodicPointsIn(1, [0, 1]);
    expect(pts.length).toBe(2);
    expect(mid(pts[1].value)).toBeCloseTo(0.7364844482, 9);
    expect(pts[1].stability).toBe("repelling");
  });
  test("with() changes a parameter", () => {
    const L = logistic(4).with({ r: 2 });
    expect(L.at("1/2").toString()).toBe("1/2");
    expect(L.criticalPoints().map(String)).toEqual(["1/2"]);
  });
});

describe("long exact orbits", () => {
  test("a bounded orbit with long numerators is not mistaken for an escaping one", () => {
    // 1/5 under 4x(1 − x): denominators 5^(2^n) pass 2^4096 by step 12, though every term is in [0, 1]
    const r = logistic(4).enclosedOrbit("1/5", 14, { bits: 6000 });
    expect(r.terms.length).toBe(15);
    for (const t of r.terms) expect(t.high.toNumber()).toBeLessThanOrEqual(1);
    expect(() => logistic(5).enclosedOrbit("1/2", 40)).toThrow(/leaves every bounded range/);
  });
});
