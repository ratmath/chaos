// Type declarations for @ratmath/chaos.
import type { Rational, RationalInterval } from "@ratmath/core";

/** An exact number: a whole number, a bigint, a Rational, or a string such as "1/3" or "3.2". A JavaScript number that is not whole is refused. */
export type Exact = Rational | bigint | number | string;

export interface PlannedPhase {
  phase: number;
  contents: string;
  /** Rough size in source lines, as text. */
  estimate: string;
}
/** The phases still to come. */
export const plan: readonly PlannedPhase[];

export class ChaosError extends Error {
  code: string;
}

export type Stability = "attracting" | "repelling" | "neutral" | "superattracting" | "undecided";
/** An algebraic real from @ratmath/algebraic. */
export interface AlgebraicValue {
  poly: bigint[];
  interval(bits?: number): [Rational, Rational];
  toNumber(): number;
  toString(digits?: number): string;
}
export interface PeriodicPoint {
  value: AlgebraicValue;
  text: string;
  period: number;
  multiplier: RationalInterval;
  stability: Stability;
  /** The multiplier when it is exactly 0, 1 or −1, else null. */
  exactMultiplier: Rational | null;
}
export interface Cycle {
  period: number;
  points: string[];
  values: AlgebraicValue[];
  multiplier: RationalInterval;
  stability: Stability;
  /** log|multiplier|/period, enclosed; −Infinity for a superattracting cycle. */
  lyapunov: RationalInterval | number | null;
}
export interface MapOptions {
  variable?: string;
  params?: Record<string, Exact>;
  /** The turning point, for maps whose derivative jumps there. */
  critical?: Exact;
  name?: string;
}
export class IteratedMap {
  constructor(formula: string, options?: MapOptions);
  readonly text: string;
  readonly variable: string;
  readonly params: Record<string, Rational>;
  /** Coefficients, constant term first, when the map is a polynomial; else null. */
  readonly polynomial: Rational[] | null;
  with(params: Record<string, Exact>): IteratedMap;
  at(x: Exact): Rational;
  float(x: number): number;
  enclose(x: Exact | RationalInterval, bits?: number): RationalInterval;
  derivativeAt(x: Exact | RationalInterval, bits?: number): Rational | RationalInterval;
  orbit(x0: Exact, n: number, options?: { maxBits?: number }): Rational[];
  floatOrbit(x0: Exact, n: number): number[];
  enclosedOrbit(x0: Exact | RationalInterval, n: number, options?: { digits?: number; bits?: number; maxBits?: number }): { terms: RationalInterval[]; bits: number; width: Rational };
  compare(x0: Exact, n: number, options?: { digits?: number }): { rows: Array<{ n: number; float: number; exact: RationalInterval; error: number; digits: number }>; parted: number | null };
  cobweb(x0: Exact, n: number): Array<[number, number]>;
  iterate(n: number): Rational[];
  periodicPoints(n?: number, options?: { bits?: number }): PeriodicPoint[];
  fixedPoints(options?: { bits?: number }): PeriodicPoint[];
  cycles(n?: number, options?: { bits?: number }): Cycle[];
  periodicPointsIn(n: number, interval: [Exact, Exact], options?: { digits?: number }): Array<{ value: RationalInterval; period: number; leastPeriod: number | string; multiplier: RationalInterval; stability: Stability }>;
  criticalPoints(interval?: [Exact, Exact]): Array<Rational | RationalInterval>;
  toString(): string;
}
export function map(formula: string, options?: MapOptions): IteratedMap;
export function logistic(r?: Exact): IteratedMap;
export function tent(mu?: Exact): IteratedMap;
export function quadratic(c?: Exact): IteratedMap;
export function doubling(x0: Exact, options?: { base?: Exact; limit?: number }): { orbit: Rational[]; preperiod: number; period: number; cycle: Rational[]; floatReachesZeroAt: number | null };

export function lyapunov(map: IteratedMap, x0: Exact, options?: { steps?: number; transient?: number; digits?: number }): { value: RationalInterval; steps: number; transient: number; bits: number; kind: string; float: number };
export function sensitivity(map: IteratedMap, x0: Exact, delta: Exact, options?: { steps?: number; threshold?: Exact; digits?: number }): { rows: Array<{ n: number; separation: RationalInterval }>; reaches: number | null; rate: number | null; factor: number | null };
export function bifurcationDiagram(family: IteratedMap, options: { parameter: string; range: [Exact, Exact]; count?: number; start?: number; transient?: number; keep?: number }): { kind: "estimated"; parameter: string; columns: Array<{ parameter: number; points: number[] }> };
export interface Bifurcation { parameter: RationalInterval; point: RationalInterval; period: number; multiplier: 1 | -1; kind: string }
export function bifurcationPoint(family: IteratedMap, options: { parameter: string; period?: number; multiplier?: 1 | -1; near: [Exact, Exact]; digits?: number }): Bifurcation;
export function superstable(family: IteratedMap, options: { parameter: string; period: number; range: [Exact, Exact]; critical?: Exact; digits?: number }): RationalInterval[];
export function periodDoubling(family: IteratedMap, options: { parameter: string; levels?: number; start: [Exact, Exact]; digits?: number }): { points: Bifurcation[]; parameters: RationalInterval[]; ratios: RationalInterval[]; feigenbaum: string; kind: string };
export function periodicWindow(family: IteratedMap, options: { parameter: string; period: number; near: [Exact, Exact]; digits?: number }): { period: number; born: RationalInterval; doubles: RationalInterval; superstable: RationalInterval | null; width: RationalInterval };
export function itinerary(map: IteratedMap, x0: Exact, n: number, options?: { critical?: Exact }): string;
export function kneading(map: IteratedMap, n: number, options?: { critical?: Exact }): string;
export const sharkovskii: {
  compare(a: number, b: number): number;
  forces(a: number, b: number): boolean;
  forced(n: number, limit?: number): number[];
  order(limit?: number): number[];
};

// ---- phases 4 and 5: fractals, L-systems, automata, escape-time pictures ----

export type Point = [Rational, Rational];
export interface Dimension { kind: "exact" | "certified"; value: Rational | null; interval: [Rational, Rational]; formula: string }
export function seeded(seed?: number): () => number;
export interface AffineSpec { a?: Exact; b?: Exact; c?: Exact; d?: Exact; e?: Exact; f?: Exact }
export class Affine {
  constructor(spec?: AffineSpec);
  readonly a: Rational; readonly b: Rational; readonly c: Rational; readonly d: Rational; readonly e: Rational; readonly f: Rational;
  /** r times a rotation by `rotate` quarter turns (optionally reflected), then shifted by (e, f). */
  static similarity(r: Exact, options?: { e?: Exact; f?: Exact; rotate?: number; reflect?: boolean }): Affine;
  readonly det: Rational;
  apply(p: [Exact, Exact]): Point;
  applyFloat(p: [number, number]): [number, number];
  similarityRatio(): { r2: Rational; r: Rational | null } | null;
  lipschitzBound(): Rational;
}
export class IFS {
  constructor(maps: Array<Affine | AffineSpec>, options?: { name?: string; probabilities?: Exact[] });
  readonly maps: readonly Affine[];
  readonly name: string;
  readonly probabilities: readonly Rational[];
  contraction(): { bounds: Rational[]; contracting: boolean };
  step(points: Point[]): Point[];
  stage(n: number, start?: Array<[Exact, Exact]> | null): Point[];
  fixedPoint(i?: number): Point;
  chaosGame(n?: number, options?: { seed?: number; skip?: number }): { kind: "estimated"; seed: number; points: Float64Array; count: number };
  similarityDimension(options?: { bits?: number }): Dimension;
}
export const attractors: Readonly<Record<"sierpinski" | "carpet" | "vicsek" | "cantorDust" | "heighway" | "levyDragon" | "barnsley", () => IFS>>;
export function moran(ratios: Exact[], options?: { bits?: number; squared?: boolean }): Dimension;
export class DigitCantor {
  constructor(options?: { base?: number; digits?: number[] });
  readonly base: number;
  readonly digits: readonly number[];
  contains(q: Exact): { inside: boolean; expansion?: string; reason: string };
  expansion(q: Exact): string;
  stage(n: number): Array<[Rational, Rational]>;
  measure(n: number): Rational;
  dimension(options?: { bits?: number }): Dimension;
}
export function cantor(options?: { base?: number; digits?: number[] }): DigitCantor;
export function boxCounting(points: { points: Float64Array } | Array<[number, number]>, options?: { scales?: number[]; bounds?: [number, number, number, number] | null }): { kind: "estimated"; dimension: number; rows: Array<{ k: number; boxes: number }>; points: number };

export class Q3 {
  constructor(a?: Rational, b?: Rational);
  readonly a: Rational;
  readonly b: Rational;
  add(o: Q3): Q3;
  mul(o: Q3): Q3;
  equals(o: Q3): boolean;
  toNumber(): number;
  toString(): string;
}
export interface LSystemSpec { axiom: string; rules?: Record<string, string>; angle?: number; name?: string; draw?: string; move?: string }
export class LSystem {
  constructor(spec: LSystemSpec);
  readonly axiom: string;
  readonly rules: Readonly<Record<string, string>>;
  readonly angle: number;
  readonly name: string;
  generate(n: number): string;
  counts(n: number): { counts: Record<string, bigint>; length: bigint };
  turtle(n: number, options?: { step?: number }): { kind: "exact" | "float"; segments: Array<[number, number, number, number]>; exactSegments: Array<[Q3, Q3, Q3, Q3]> | null; end: [Q3, Q3] | [number, number] };
}
export function lsystem(spec: LSystemSpec): LSystem;
export const lsystems: Readonly<Record<"koch" | "snowflake" | "arrowhead" | "dragon" | "hilbert" | "levy" | "plant", () => LSystem>>;

export class Elementary {
  constructor(rule: number);
  readonly rule: number;
  readonly table: Record<string, number>;
  readonly mirror: number;
  readonly complement: number;
  readonly equivalents: number[];
  readonly additive: boolean;
  step(row: Uint8Array, options?: { boundary?: "periodic" | "zero" }): Uint8Array;
  run(initial: string | ArrayLike<number>, steps: number, options?: { boundary?: "periodic" | "zero" }): Uint8Array[];
  fromCell(steps: number): Uint8Array[];
  cycle(initial: string | ArrayLike<number>, options?: { limit?: number }): { transient: number; period: number };
  toString(): string;
}
export function elementary(rule: number): Elementary;
export function showRows(rows: ArrayLike<number>[], on?: string, off?: string): string;
export type LifeClass =
  | { kind: "dies"; at: number }
  | { kind: "still life" | "oscillator"; period: number; from: number }
  | { kind: "spaceship"; period: number; displacement: [number, number]; from: number; speed: string }
  | { kind: "unknown"; limit: number };
export class Life {
  constructor(cells: string | Iterable<[number, number]>, options?: { rule?: string });
  readonly rule: string;
  readonly cells: Set<string>;
  readonly population: number;
  step(): Life;
  run(n: number): Life;
  bounds(): { x0: number; y0: number; x1: number; y1: number } | null;
  normalized(): { shape: string; dx: number; dy: number };
  classify(options?: { limit?: number }): LifeClass;
  toString(on?: string, off?: string): string;
}
export function life(cells: string | Iterable<[number, number]>, options?: { rule?: string }): Life;
export const patterns: Readonly<Record<"block" | "blinker" | "toad" | "beacon" | "glider" | "lwss" | "pulsar" | "rpentomino" | "gosperGun", string>>;

export type Range = [Exact, Exact];
export interface PictureOptions { re?: Range; im?: Range; width?: number; height?: number; maxIter?: number; depth?: number }
export class PixelGrid {
  readonly name: string;
  readonly width: number;
  readonly height: number;
  /** 0 undecided; 1 inside and 2 outside, or k for root k in Newton basins. */
  readonly kind: Uint8Array;
  readonly counts: Record<string, number>;
  readonly labels?: boolean;
  readonly roots?: Array<{ re: [number, number]; im: [number, number]; approx: [number, number] }>;
  readonly unproved?: Array<[number, number]>;
  readonly re: Rational;
  readonly reHigh: Rational;
  readonly im: Rational;
  readonly imHigh: Rational;
  at(i: number, j: number): number;
  toText(chars?: string): string;
  toString(): string;
}
export function mandelbrot(options?: PictureOptions & { maxPeriod?: number }): PixelGrid;
export function julia(c: Exact, options?: PictureOptions): PixelGrid;
export function newtonBasins(coefficients: Exact[], options?: PictureOptions): PixelGrid;
export interface PointVerdict { verdict: "inside" | "outside" | "undecided"; reason: string; period?: number; escape?: number }
export function mandelbrotPoint(c: Exact, options?: { maxIter?: number; bits?: number; maxPeriod?: number }): PointVerdict;
export function juliaPoint(c: Exact, z0: Exact, options?: { maxIter?: number; bits?: number; maxPeriod?: number }): PointVerdict;

// ---- maps of the plane

/** A complex number as two enclosures; im is [0, 0] for a real one. */
export interface EnclosedComplex { re: RationalInterval; im: RationalInterval }
export type PlaneStability = "attracting" | "repelling" | "saddle" | "elliptic" | "undecided";
/** A point given by exact numbers, or by intervals for a box of points. */
export type PlanePoint = [Exact | RationalInterval, Exact | RationalInterval];
export interface Linearisation {
  /** The eigenvalues of the derivative of Fⁿ at the point, or null when they could not be separated. */
  multipliers: EnclosedComplex[] | null;
  stability: PlaneStability;
  /** At a saddle: the contracting direction, scaled so that one component is 1. */
  stable: [RationalInterval, RationalInterval] | null;
  /** At a saddle: the expanding direction. */
  unstable: [RationalInterval, RationalInterval] | null;
}
export interface PlanePeriodicPoint extends Linearisation {
  point: [RationalInterval, RationalInterval];
  period: number;
  /** The least period; null when a smaller one could not be ruled out or proven. */
  leastPeriod: number | null;
}
export interface ExactPlanePoint extends Linearisation {
  x: AlgebraicValue;
  y: AlgebraicValue;
  text: [string, string];
  period: number;
}
export type PlaneRegion = Record<string, [Exact, Exact]>;
export interface PlaneMapOptions { variables?: [string, string]; params?: Record<string, Exact>; name?: string }
export class PlaneMap {
  constructor(field: [string, string], options?: PlaneMapOptions);
  readonly text: [string, string];
  readonly variables: [string, string];
  readonly params: Record<string, Rational>;
  /** The Jacobian determinant when it is the same everywhere, else null. */
  readonly determinant: Rational | null;
  readonly areaPreserving: boolean;
  with(params: Record<string, Exact>): PlaneMap;
  at(p: [Exact, Exact]): [Rational, Rational];
  float(p: [number, number]): [number, number];
  enclose(p: PlanePoint, bits?: number): [RationalInterval, RationalInterval];
  jacobian(): [[string, string], [string, string]];
  orbit(p0: [Exact, Exact], n: number, options?: { maxBits?: number }): Array<[Rational, Rational]>;
  floatOrbit(p0: [Exact, Exact], n: number): Array<[number, number]>;
  enclosedOrbit(p0: PlanePoint, n: number, options?: { digits?: number; bits?: number; maxBits?: number }): { points: Array<[RationalInterval, RationalInterval]>; bits: number; width: Rational };
  compare(p0: [Exact, Exact], n: number, options?: { digits?: number }): { rows: Array<{ n: number; float: [number, number]; exact: [RationalInterval, RationalInterval]; error: number; digits: number }>; parted: number | null };
  /** Every point of period n in the region, each proven; `complete` says none was missed. */
  periodicPoints(n: number, region: PlaneRegion, options?: { all?: boolean; tol?: number; maxBoxes?: number }): PlanePeriodicPoint[] & { complete: boolean };
  fixedPoints(region: PlaneRegion, options?: { tol?: number; maxBoxes?: number }): PlanePeriodicPoint[] & { complete: boolean };
  cycles(n: number, region: PlaneRegion, options?: { tol?: number; maxBoxes?: number }): Array<{ period: number; points: Array<[RationalInterval, RationalInterval]>; multipliers: EnclosedComplex[] | null; stability: PlaneStability; whole: boolean }> & { complete: boolean };
  /** Points of period n as exact algebraic numbers, when one coordinate can be eliminated. */
  exactPeriodicPoints(n?: number): ExactPlanePoint[];
  /** Lyapunov exponents estimated in doubles; `sum` is an enclosure of their sum when the determinant is constant. */
  lyapunov(p0: [Exact, Exact], n?: number, options?: { skip?: number }): { exponents: [number, number]; kind: "estimated"; steps: number; sum: RationalInterval | null };
  toString(): string;
}
export function planeMap(field: [string, string], options?: PlaneMapOptions): PlaneMap;
export function henon(a?: Exact, b?: Exact): PlaneMap;
export function standardMap(k?: Exact): PlaneMap;

export class TorusMap {
  constructor(matrix?: Array<Array<number | bigint>>, options?: { name?: string });
  readonly matrix: bigint[][];
  readonly determinant: bigint;
  readonly trace: bigint;
  readonly invertible: boolean;
  readonly hyperbolic: boolean;
  at(p: [Exact, Exact]): [Rational, Rational];
  orbit(p0: [Exact, Exact], n: number): Array<[Rational, Rational]>;
  floatOrbit(p0: [Exact, Exact], n: number): Array<[number, number]>;
  period(p0: [Exact, Exact], options?: { limit?: number }): { period: number; preperiod: number; cycle: Array<[Rational, Rational]> };
  /** For every period, how many points with denominator q have it. */
  periods(q: number): Map<number, number>;
  fixedPointCount(n: number): bigint;
  eigen(): Array<{ value: AlgebraicValue; text: string; direction: AlgebraicValue | null; directionText: string; expanding: boolean }>;
  entropy(bits?: number): { text: string; value: RationalInterval };
  floatPartsAt(p0: [Exact, Exact], options?: { limit?: number; tolerance?: number }): number | null;
  toString(): string;
}
export function catMap(): TorusMap;
export function torusMap(matrix: Array<Array<number | bigint>>, options?: { name?: string }): TorusMap;

export class Horseshoe {
  at(p: [Exact, Exact]): [Rational, Rational] | null;
  itinerary(p: [Exact, Exact], n: number): string;
  periodicPoint(word: string): [Rational, Rational];
  periodicPoints(n: number): Array<{ word: string; point: [Rational, Rational]; leastPeriod: number }>;
  orbitCount(n: number): number;
  strips(n: number): { horizontal: Array<[Rational, Rational]>; vertical: Array<[Rational, Rational]> };
  readonly dimension: { text: string; value: number };
}
export function horseshoe(): Horseshoe;

// ---- flows

/** A starting value: an exact number, or an interval for a box of starting points. */
export type FlowValue = Exact | RationalInterval | [Exact, Exact];
export type FlowPoint = FlowValue[] | Record<string, FlowValue>;
export interface FlowOptions { variables?: string[]; params?: Record<string, Exact>; name?: string }
export interface SolverOptions {
  /** Degree of the Taylor polynomials (default 12). */
  order?: number;
  maxStep?: number;
  maxSteps?: number;
}
export interface FlowSolution {
  variables: string[];
  /** Enclosures of the solution at the end time, one per variable. */
  final: RationalInterval[];
  /** False when the solver stopped early; `reason` says why and `reached` how far it got. */
  complete: boolean;
  reason: string | null;
  reached: RationalInterval;
  steps: number;
  /** The widest of the final enclosures, as a double. */
  width: number;
  at(time: Exact): RationalInterval[];
}
export interface Equilibrium {
  point: RationalInterval[];
  eigenvalues: EnclosedComplex[] | null;
  kind: "sink" | "spiral sink" | "source" | "spiral source" | "saddle" | "saddle-focus" | "undecided";
  stability: "stable" | "unstable" | "undecided";
  unstableDimension: number | null;
  spiral: boolean;
}
export interface Section {
  /** The section is the plane where this variable takes `value`. */
  variable: string;
  value?: Exact;
  /** 1 or "up": crossed with the variable increasing; −1 or "down": decreasing. */
  direction?: 1 | -1 | "up" | "down";
}
export type SectionPoint = Record<string, FlowValue | number> | Array<FlowValue | number>;
export interface ReturnMapResult {
  /** The section's coordinates: the variables other than the section's own. */
  variables: string[];
  point: RationalInterval[];
  time: RationalInterval;
  /** The derivative of the return map in the section's coordinates. */
  derivative: RationalInterval[][];
}
export interface ProvenOrbit {
  proven: true;
  variables: string[];
  /** A point of the orbit on the section, one enclosure per variable. */
  point: RationalInterval[];
  period: RationalInterval;
  returns: number;
  multipliers: EnclosedComplex[] | null;
  stability: "attracting" | "repelling" | "saddle" | "undecided";
  /** Half the side of the box in which the orbit is the only one. */
  radius: number;
}
export interface UnprovenOrbit { proven: false; variables: string[]; point: number[]; reason: string }
export interface EquilibriumBifurcation {
  /** The parameter value lies between these two. */
  parameter: [Rational, Rational];
  kind: "saddle-node" | "pitchfork" | "transcritical" | "Hopf" | "stability change" | "equilibria appear";
  /** The kinds of the equilibria just before and just after. */
  before: string[];
  after: string[];
  /** True when the bracket could not be narrowed to the tolerance. */
  stalled: boolean;
  tags: [string, string];
}
export class Flow {
  constructor(field: string[], options?: FlowOptions);
  readonly n: number;
  readonly text: string[];
  readonly variables: string[];
  readonly params: Record<string, Rational>;
  with(params: Record<string, Exact>): Flow;
  jacobian(): string[][];
  velocity(p: number[]): number[];
  /** The solution through a point or a box of points, enclosed. */
  solve(p0: FlowPoint, until: Exact, options?: SolverOptions): FlowSolution;
  floatOrbit(p0: FlowPoint, until: Exact | number, options?: { steps?: number }): { t: number[]; points: number[][] };
  equilibria(region: Record<string, [Exact, Exact]>, options?: { tol?: number; maxBoxes?: number }): Equilibrium[] & { complete: boolean };
  returnMap(section: Section, start: SectionPoint, options?: SolverOptions & { returns?: number }): ReturnMapResult;
  periodicOrbit(section: Section, guess: SectionPoint, options?: SolverOptions & { returns?: number; radius?: Exact }): ProvenOrbit | UnprovenOrbit;
  bifurcations(param: string, range: [Exact, Exact], region: Record<string, [Exact, Exact]>, options?: { samples?: number; tol?: Exact; maxBoxes?: number }): EquilibriumBifurcation[];
  toString(): string;
}
export function flow(field: string[], options?: FlowOptions): Flow;
/** The Lorenz system x′ = s(y − x), y′ = x(r − z) − y, z′ = xy − bz. */
export function lorenz(s?: Exact, r?: Exact, b?: Exact): Flow;
/** The Rössler system x′ = −y − z, y′ = x + ay, z′ = b + z(x − c). */
export function rossler(a?: Exact, b?: Exact, c?: Exact): Flow;
/** Van der Pol's oscillator x′ = y, y′ = μ(1 − x²)y − x. */
export function vanDerPol(mu?: Exact): Flow;
/** The pendulum x′ = y, y′ = −sin x − c·y. */
export function pendulum(damping?: Exact): Flow;
