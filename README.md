# @ratmath/chaos

Dynamical systems and chaos for [RatMath](https://ratmath.com/), part of the
[RatMath Shelf](https://ratmath.com/shelf/).

**Status: set up, nothing built yet.** This repository holds the package
skeleton and the docs site. The plan is in the
[library roadmap](https://github.com/ratmath/ratmath/blob/main/docs/design/library-roadmap.md#10-a-dynamical-systems-package-ratmathchaos).

Dynamical systems are where floating point fails fastest: a chaotic map doubles every rounding error at each step. Here orbits are exact or enclosed, periodic points are algebraic numbers, and statements about stability and bifurcation are proved.

## Planned

| Phase | Contents |
|---|---|
| 1 | Maps of the line: orbits, periodic points, stability |
| 2 | Bifurcation diagrams, Lyapunov exponents, symbolic dynamics |
| 3 | Maps of the plane; flows, sections, periodic orbits, bifurcations |
| 4 | Iterated function systems, Cantor sets, dimensions, L-systems, cellular automata |
| 5 | Julia and Mandelbrot sets, Newton basins |
| 6 | Pictures |

It will depend on: Core, analytic (enclosures, the validated solver), algebraic (periodic points) and linalg (linearisation); complex for Julia and Mandelbrot sets.

## Layout

- `index.js`, `src/`: the library. For now it exports only `plan`, the list
  above, so that the tests and the docs build have something to load.
- `tests/`: `bun run test`.
- `docs/src/*.html`: the docs pages (overview, explore, reference, then one
  page per topic as topics are built). `bun run docs:build` writes `docs/`,
  typesets the mathematics and runs the code box at the end of every page.
  `docs/*.html` and `docs/assets/chaos.js` are build output.
- `docs/assets/runbox.js` and `scripts/docs-kit.js` are shared with the other
  RatMath docs sites and should be kept in step with them.

## Development

```sh
bun install
bun run test
bun run docs:build
```

Each topic that is built gets its own docs page with usage, the mathematics
explained, and a runnable code box at the end, and is then listed on the
[Shelf page](https://ratmath.com/shelf/).

MIT licensed.
