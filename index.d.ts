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
