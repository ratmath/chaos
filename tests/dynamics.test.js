import { describe, expect, test, setDefaultTimeout } from "bun:test";
import { logistic, tent, quadratic, map, lyapunovExponent, invariantDensity, lyapunov, sensitivity, bifurcationDiagram, bifurcationPoint, superstable, periodDoubling, periodicWindow, itinerary, kneading, sharkovskii } from "../index.js";

// certified searches: seconds each here, more on a slow machine
setDefaultTimeout(120000);

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

describe("the Lyapunov exponent from the invariant density", () => {
  const has = (iv, x) => iv.low.toNumber() <= x && x <= iv.high.toNumber();
  const width = (iv) => iv.high.subtract(iv.low).toNumber();
  test("logistic map at r = 4: an enclosure of log 2", () => {
    const e = lyapunovExponent(logistic(4), { digits: 4 });
    expect(has(e.value, Math.LN2)).toBe(true);
    expect(width(e.value)).toBeLessThan(2e-4);
    expect(e.criticalPoints.map(String)).toEqual(["1/2"]);
    expect(e.invariance.proven).toBe(true);
    expect(e.invariance.holds).toBe(true);
    expect(e.exponent).toBe("log(2)");
    // the finite orbit average is near it but is not it
    const avg = lyapunov(logistic(4), "1/5").value;
    expect(Math.abs(avg.low.toNumber() - Math.LN2)).toBeLessThan(5e-3);
  });
  test("tent map and x^2 - 2", () => {
    const t = lyapunovExponent(tent(2));
    expect(has(t.value, Math.LN2)).toBe(true);
    expect(width(t.value)).toBeLessThan(1e-12);
    const q = lyapunovExponent(quadratic(-2), { digits: 3 });
    expect(has(q.value, Math.LN2)).toBe(true);
    expect(q.criticalPoints.map(String)).toEqual(["0"]);
  });
  test("a density supplied by the caller is tested for invariance: Chebyshev T3 gives log 3", () => {
    const e = lyapunovExponent(map("4*x^3 - 3*x"), { distribution: "1 - acos(x)/pi", interval: [-1, 1], digits: 3 });
    expect(has(e.value, Math.log(3))).toBe(true);
    expect(e.invariance.proven).toBe(false);
    expect(e.invariance.checks.length).toBe(5);
    expect(e.invariance.holds).toBe(true);
  });
  test("refusals: a measure that is not invariant, a map with no known density", () => {
    expect(() => lyapunovExponent(logistic(4), { distribution: "x", interval: [0, 1], digits: 2 })).toThrow(/not invariant/);
    expect(() => lyapunovExponent(logistic("7/2"))).toThrow(/no closed form/);
    expect(() => lyapunovExponent(logistic(4), { distribution: "x/2", interval: [0, 1] })).toThrow(/from 0/);
    expect(invariantDensity(logistic(4)).density).toBe("1/(pi*sqrt(x*(1 - x)))");
  });
  test("known densities are recognised by the formula, not the name", () => {
    // the identity under the name "logistic": its exponent is 0, not log 2
    expect(() => invariantDensity(map("x", { name: "logistic", params: { r: 4 } }))).toThrow(/no closed form/);
    expect(() => invariantDensity(map("x^2", { name: "tent", params: { mu: 2 } }))).toThrow(/no closed form/);
    expect(invariantDensity(map("4*x*(1 - x)")).exponent).toBe("log(2)");
    expect(invariantDensity(map("x^2 - 2", { name: "logistic" })).density).toBe("1/(pi*sqrt(4 - x^2))");
    expect(invariantDensity(map("2*(1/2 - abs(x - 1/2))")).density).toBe("1");
    const t = lyapunovExponent(map("4*t - 4*t^2", { variable: "t" }), { digits: 2 });
    expect(t.invariance.proven).toBe(true);
  });
});
