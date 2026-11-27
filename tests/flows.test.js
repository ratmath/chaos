import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import { flow, lorenz, rossler, vanDerPol, pendulum } from "../index.js";

const mid = (iv) => (iv.low.toNumber() + iv.high.toNumber()) / 2;
const width = (iv) => iv.high.subtract(iv.low).toNumber();
const holds = (iv, x) => iv.low.toNumber() <= x && x <= iv.high.toNumber();

describe("solving", () => {
  test("a linear system against its exact solution", () => {
    // x′ = y, y′ = −x from (1, 0): (cos t, −sin t)
    const s = flow(["y", "-x"]).solve([1, 0], 2);
    expect(s.complete).toBe(true);
    expect(holds(s.final[0], Math.cos(2))).toBe(true);
    expect(holds(s.final[1], -Math.sin(2))).toBe(true);
    expect(s.width).toBeLessThan(1e-12);
    const half = s.at("1/2");
    expect(holds(half[0], Math.cos(0.5))).toBe(true);
    expect(width(half[0])).toBeLessThan(1e-12);
    expect(() => s.at(3)).toThrow(/outside/);
  });
  test("exponential decay from a box of starting points", () => {
    const s = flow(["-x"]).solve([[1, 2]], 1);
    expect(s.final[0].low.toNumber()).toBeLessThanOrEqual(Math.exp(-1));
    expect(s.final[0].high.toNumber()).toBeGreaterThanOrEqual(2 * Math.exp(-1));
    expect(width(s.final[0])).toBeLessThan(Math.exp(-1) + 1e-6);
  });
  test("the Lorenz system: the enclosure holds a careful float solution", () => {
    const L = lorenz();
    const s = L.solve([1, 1, 1], 3);
    expect(s.complete).toBe(true);
    expect(s.width).toBeLessThan(1e-9);
    const rk = L.floatOrbit([1, 1, 1], 3, { steps: 30000 }).points.at(-1);
    s.final.forEach((iv, i) => expect(Math.abs(mid(iv) - rk[i])).toBeLessThan(1e-7));
    expect(String(L)).toContain("x′ = s*(y - x)");
    expect(L.jacobian()[0]).toEqual(["-10", "10", "0"]);
  });
  test("a solution that blows up is reported, not invented", () => {
    const s = flow(["x^2"]).solve([1], 2);
    expect(s.complete).toBe(false);
    expect(s.reached.high.toNumber()).toBeLessThan(1.01);
  });
});

describe("equilibria", () => {
  test("the Lorenz equilibria and their kinds", () => {
    const eq = lorenz().equilibria({ x: [-30, 30], y: [-30, 30], z: [-10, 60] });
    expect(eq.complete).toBe(true);
    expect(eq.map((e) => e.kind)).toEqual(["saddle-focus", "saddle", "saddle-focus"]);
    expect(mid(eq[2].point[0])).toBeCloseTo(Math.sqrt(72), 12);
    expect(mid(eq[2].point[2])).toBeCloseTo(27, 12);
    // the origin's eigenvalues: −8/3 and the roots of λ² + 11λ − 270
    const ev = eq[1].eigenvalues.map((e) => mid(e.re)).sort((a, b) => a - b);
    expect(ev[0]).toBeCloseTo((-11 - Math.sqrt(1201)) / 2, 9);
    expect(ev[1]).toBeCloseTo(-8 / 3, 9);
    expect(ev[2]).toBeCloseTo((-11 + Math.sqrt(1201)) / 2, 9);
    expect(eq[1].unstableDimension).toBe(1);
  });
  test("the damped pendulum", () => {
    const eq = pendulum("1/2").equilibria({ x: [-4, 4], y: [-1, 1] });
    expect(eq.map((e) => e.kind)).toEqual(["saddle", "spiral sink", "saddle"]);
    expect(holds(eq[2].point[0], Math.PI)).toBe(true);
    expect(width(eq[2].point[0])).toBeLessThan(1e-12);
  });
});

describe("return maps and periodic orbits", () => {
  test("Van der Pol's limit cycle is proven, with its period", () => {
    const po = vanDerPol(1).periodicOrbit({ variable: "y", value: 0, direction: -1 }, { x: 2 });
    expect(po.proven).toBe(true);
    expect(po.stability).toBe("attracting");
    expect(mid(po.point[0])).toBeCloseTo(2.00861986087484, 10);
    expect(holds(po.period, 6.663286859323130)).toBe(true);
    expect(width(po.period)).toBeLessThan(1e-9);
    expect(po.multipliers.length).toBe(1);
    expect(mid(po.multipliers[0].re)).toBeCloseTo(0.00085969, 6);
  });
  test("the return map itself", () => {
    const r = vanDerPol(1).returnMap({ variable: "y", value: 0, direction: -1 }, { x: 2 });
    expect(r.variables).toEqual(["x"]);
    expect(mid(r.point[0])).toBeCloseTo(2.0086122778, 8);
    expect(mid(r.derivative[0][0])).toBeCloseTo(0.00089996559, 8);
    // a box of starting points is carried too
    const box = vanDerPol(1).returnMap({ variable: "y", value: 0, direction: -1 }, { x: ["1.9999", "2.0001"] });
    expect(holds(box.point[0], mid(r.point[0]))).toBe(true);
    expect(width(box.point[0])).toBeLessThan(1e-4);
    expect(() => vanDerPol(1).returnMap({ variable: "y", value: 0, direction: 1 }, { x: 2 })).toThrow(/direction/);
  });
  test("a saddle periodic orbit of the Lorenz attractor", () => {
    const po = lorenz().periodicOrbit({ variable: "z", value: 27 }, { x: -13.76, y: -19.58 }, { returns: 2 });
    expect(po.proven).toBe(true);
    expect(po.stability).toBe("saddle");
    // Viswanath's value for the shortest orbit
    expect(holds(po.period, 1.5586522107161747)).toBe(true);
    expect(width(po.period)).toBeLessThan(1e-8);
    expect(mid(po.point[0])).toBeCloseTo(-13.763610682134, 8);
    const big = Math.max(...po.multipliers.map((e) => Math.abs(mid(e.re))));
    expect(big).toBeCloseTo(4.7129, 3);
  });
  test("the Rössler system: an attracting cycle and a saddle one", () => {
    const sec = { variable: "y", value: 0, direction: -1 };
    const a = rossler("0.2", "0.2", "2.5").periodicOrbit(sec, { x: -3.7, z: 0.03 });
    expect(a.proven).toBe(true);
    expect(a.stability).toBe("attracting");
    expect(mid(a.period)).toBeCloseTo(5.7489911833, 8);
    const b = rossler().periodicOrbit(sec, { x: -7.2, z: 0.02 });
    expect(b.proven).toBe(true);
    expect(b.stability).toBe("saddle");
    expect(mid(b.period)).toBeCloseTo(5.8810884556, 7);
  });
  test("no orbit, no proof", () => {
    const po = flow(["y", "-x - y"]).periodicOrbit({ variable: "y", value: 0, direction: -1 }, { x: 1 });
    expect(po.proven).toBe(false);
    expect(typeof po.reason).toBe("string");
  });
});

describe("bifurcations of equilibria", () => {
  const show = (b) => b.map((x) => x.kind);
  const inside = (b, v) => b.parameter[0].toNumber() <= v && v <= b.parameter[1].toNumber();
  test("the three normal forms in one variable", () => {
    const fold = flow(["mu - x^2"], { params: { mu: 0 } }).bifurcations("mu", [-1, 1], { x: [-3, 3] });
    expect(show(fold)).toEqual(["saddle-node"]);
    expect(inside(fold[0], 0)).toBe(true);
    expect(fold[0].before).toEqual([]);
    expect(fold[0].after).toEqual(["source", "sink"]);
    expect(show(flow(["mu*x - x^3"], { params: { mu: 0 } }).bifurcations("mu", [-1, 1], { x: [-3, 3] }))).toEqual(["pitchfork"]);
    expect(show(flow(["mu*x - x^2"], { params: { mu: 0 } }).bifurcations("mu", [-1, 1], { x: [-3, 3] }))).toEqual(["transcritical"]);
  });
  test("Hopf bifurcations", () => {
    const v = vanDerPol(0).bifurcations("mu", [-1, 1], { x: [-3, 3], y: [-3, 3] });
    expect(show(v)).toEqual(["Hopf"]);
    expect(inside(v[0], 0)).toBe(true);
    expect(v[0].parameter[1].subtract(v[0].parameter[0]).compareTo(new Rational(1n, 1000000n)) <= 0).toBe(true);
  });
  test("the Lorenz system: the pitchfork at r = 1 and the Hopf bifurcation at 470/19", () => {
    const b = lorenz().bifurcations("r", ["1/2", 30], { x: [-30, 30], y: [-30, 30], z: [-10, 60] });
    expect(show(b)).toEqual(["pitchfork", "Hopf"]);
    expect(inside(b[0], 1)).toBe(true);
    expect(inside(b[1], 470 / 19)).toBe(true);
    expect(b[1].before).toEqual(["spiral sink", "saddle", "spiral sink"]);
    expect(b[1].after).toEqual(["saddle-focus", "saddle", "saddle-focus"]);
  });
});
