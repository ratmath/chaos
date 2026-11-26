import { describe, expect, test } from "bun:test";
import { logistic, tent, lyapunov, sensitivity, bifurcationDiagram, bifurcationPoint, superstable, periodDoubling, periodicWindow, itinerary, kneading, sharkovskii } from "../index.js";

const mid = (iv) => (iv.low.toNumber() + iv.high.toNumber()) / 2;
const width = (iv) => iv.high.subtract(iv.low).toNumber();

describe("Lyapunov exponents and sensitivity", () => {
  test("the tent map has exponent log 2 on every orbit that misses 1/2", () => {
    const l = lyapunov(tent(2), "1/7", { steps: 50 });
    expect(mid(l.value)).toBeCloseTo(Math.log(2), 10);
  });
  test("an attracting cycle gives the cycle's exponent", () => {
    const l = lyapunov(logistic("3.2"), "1/10", { steps: 100, transient: 100 });
    expect(mid(l.value)).toBeCloseTo(Math.log(0.16) / 2, 8);
    expect(width(l.value)).toBeLessThan(1e-9);
  });
  test("a chaotic orbit: certified along the true orbit", () => {
    const l = lyapunov(logistic(4), "1/10", { steps: 200 });
    expect(width(l.value)).toBeLessThan(1e-9);
    expect(mid(l.value)).toBeGreaterThan(0.6);
    expect(mid(l.value)).toBeLessThan(0.75);
    expect(l.kind).toMatch(/finite-time/);
  });
  test("separation grows by about 2 a step", () => {
    const s = sensitivity(logistic(4), "1/10", "1/1000000000000", { steps: 50 });
    expect(s.reaches).toBeGreaterThan(30);
    expect(s.reaches).toBeLessThan(45);
    expect(s.factor).toBeGreaterThan(1.7);
    expect(s.factor).toBeLessThan(2.3);
  });
});

describe("bifurcations", () => {
  const family = logistic();
  test("the first two period doublings of the logistic family", () => {
    const b1 = bifurcationPoint(family, { parameter: "r", period: 1, near: [2.9, 0.6] });
    expect(mid(b1.parameter)).toBeCloseTo(3, 14);
    expect(mid(b1.point)).toBeCloseTo(2 / 3, 14);
    expect(width(b1.parameter)).toBeLessThan(1e-18);
    const b2 = bifurcationPoint(family, { parameter: "r", period: 2, near: [3.44, 0.44] });
    expect(mid(b2.parameter)).toBeCloseTo(1 + Math.sqrt(6), 13);
    expect(b2.kind).toBe("period doubling");
  });
  test("the period-3 window", () => {
    const w = periodicWindow(family, { parameter: "r", period: 3, near: [3.835, 0.5] });
    expect(mid(w.born)).toBeCloseTo(1 + Math.sqrt(8), 12);
    expect(mid(w.doubles)).toBeCloseTo(3.8414990075, 9);
    expect(mid(w.superstable)).toBeCloseTo(3.8318740552, 9);
  });
  test("superstable parameters", () => {
    const s = superstable(family, { parameter: "r", period: 2, range: [3, "7/2"] });
    expect(s.length).toBe(1);
    expect(mid(s[0])).toBeCloseTo(1 + Math.sqrt(5), 12);
  });
  test("the cascade and Feigenbaum's ratios", () => {
    const d = periodDoubling(family, { parameter: "r", levels: 5, start: [2.9, 0.6] });
    expect(d.parameters.map((p) => mid(p).toFixed(6))).toEqual(["3.000000", "3.449490", "3.544090", "3.564407", "3.568759"]);
    expect(mid(d.ratios[2])).toBeGreaterThan(4.6);
    expect(mid(d.ratios[2])).toBeLessThan(4.7);
  });
  test("a bifurcation diagram is labelled as an estimate", () => {
    const d = bifurcationDiagram(family, { parameter: "r", range: [2.8, 3.5], count: 3 });
    expect(d.kind).toBe("estimated");
    expect(d.columns.map((c) => c.points.length)).toEqual([1, 2, 4]);
  });
});

describe("symbolic dynamics", () => {
  test("itineraries on the true orbit", () => {
    expect(itinerary(tent(2), "1/7", 6)).toBe("LLRRLR");
    expect(itinerary(logistic(4), "1/2", 3)).toBe("CRL");
    expect(itinerary(logistic(4), "1/10", 30)).toBe("LLRLRRRLLRLRRLLLLLLRRRRLLRLLRR");
  });
  test("kneading sequences", () => {
    expect(kneading(logistic(4), 6)).toBe("RLLLLL");
    expect(kneading(logistic("3.83"), 9)).toBe("RLRRLRRLR");
  });
  test("Sharkovskii's order", () => {
    expect(sharkovskii.order(12)).toEqual([3, 5, 7, 9, 11, 6, 10, 12, 8, 4, 2, 1]);
    expect(sharkovskii.forces(3, 1000)).toBe(true);
    expect(sharkovskii.forces(4, 3)).toBe(false);
    expect(sharkovskii.forced(6, 12)).toEqual([6, 10, 12, 8, 4, 2, 1]);
    expect(() => sharkovskii.compare(0, 1)).toThrow();
  });
});
