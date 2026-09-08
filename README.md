# aigg — Angular Integration Grid Generator

`aigg` is a research program for constructing and optimizing **spherical integration rules** (aka: *“angular integration grids”*).

It generates symmetry-adapted quadrature rules on the unit sphere, with particular emphasis on compact rules of high order and high numerical accuracy. `aigg` supports the polyhedral point groups \(T\), \(O\), \(I\), \(T_d\), \(T_h\), \(O_h\), and \(I_h\), and can perform grid optimizations using several floating-point precisions (by default, it builds `fp64`, … ,`fp512` binaries).

## Do I need aigg?

**Probably not!**

…if you only want to use the generated integration grids.

The normal way of using the grids is through the separately distributed [**aig libraries**](http://knizia.me/aig), which contain curated sets of precomputed and tested spherical integration rules. Using these libraries does not require compiling or running `aigg`.

* **If you want spherical integration grids for use in another program:** see the [aig libraries](http://knizia.me/aig).
* **If you want to generate, reproduce, refine, or investigate spherical integration rules themselves:** this repository contains the `aigg` optimizer.

`aigg` is published primarily to make the construction and verification of the integration rules reproducible, and as a research tool for developing new rules.

## What can aigg do?

Among other things, `aigg` can:

* construct and optimize spherical integration rules with fixed point-group symmetry;
* optimize grid points and weights for exact integration up to a target order \(L\);
* perform incremental upward or downward searches for the highest achievable order of a given grid ansatz;
* refine or independently verify existing integration rules;
* work with Cartesian or barycentric seed-point coordinates;
* use systematically generated or explicitly supplied initial guesses;
* perform calculations in several floating-point precisions, including high-precision optimization;
* export either symmetry-unique seed points or complete point sets for further processing.

The program is intended as specialist research software rather than as a general-purpose numerical integration library.

## Quick example

The `inputs/` directory contains commented example input files.

For example, `inputs/1-incremental-search-up.ini` shows the default usage of the program:
```text
{ point-group: icosahedral;
  degree: [4,step:+1,-1];
  max-step: 1e4*residual;
  initial-points: hex-grid{5;1;01c}
}
```
Explanation:

- instantiate the icosahedral symmetry point-group (`point-group: icosahedral`)
- …generate an initial guess for a chiral hexagonal tiling over the sphere (`initial-points: hex-grid{5;1;01c}`). `5;1` designates two integer parameters which identify the hexagonal chiral tiling; `01c` designates the order of vertices of the “*fundamental domain*” of the point group. (see [http://knizia.me/aigg](http://knizia.me/aigg) for details)
- …then, starting at $L=4$, incrementally optimize (`degree: [4,step:+1,…]`) quadrature root positions and weights (to ensure that all spherical harmonics $Y_{lm}$ with $0≤l≤L$ are integrated exactly)…
- …until the first $L$ is found which can no longer be integrated exactly with the given symmetry group and root orbits.
- …once that $L$ is found, step down by one step (`degree: […, -1]`), and re-optimize the grid.


It can be run, for example, with the 128-bit build:

```sh
./aigg_128 -c inputs/1-incremental-search-up.ini
```

This starts at order \(L=4\), repeatedly increases the target order, and uses each optimized grid as the initial guess for the next optimization. Once the ansatz can no longer realize the requested order exactly, the schedule backs off according to the specification in `degree`.

See `inputs/` for additional examples covering refinement of existing grids, different initial guesses, export, and visualization.

## Building aigg

### Requirements

The current build has been tested on Linux and under WSL using `g++`.

Required packages are:

* a C++17-capable GCC toolchain;
* [Eigen](https://eigen.tuxfamily.org/);
* [Boost](https://www.boost.org/), in particular Boost.Multiprecision;
* [SCons](https://scons.org/);
* OpenMP support.

On a system with the required development packages installed, the default build can be started with:

```sh
make -j "$(nproc)"
```

The default build creates several executable variants for different floating-point precisions, together with a 64-bit debug build. The individual executables are named according to their precision, for example:

```text
aigg_64
aigg_128
...
```

A simple build-and-run test is:

```sh
make -j "$(nproc)" &&
./aigg_128 -c inputs/1-incremental-search-up.ini
```

aigg uses `SCons` internally because the build generates several variants of the program with different scalar backends and floating-point precisions.

## Running aigg

You can control `aigg` through configuration files:
```sh
./aigg_128 -c inputs/1-incremental-search-up.ini
```
or, alternatively, pass the corresponding input arguments on the command line:
``` sh
./aigg_128 '{point-group: icosahedral; degree: [4,step:+1,-1]; max-step: 1e4*residual; initial-points: hex-grid{5;1;01c}}'
```

Run
```sh
./aigg_64
```
without arguments to display the main command-line options.

The commented files in `inputs/` are intended as the primary starting point for new calculations. Some specialist options are not separately documented; comments in `SearchOptions.cpp` provide the complete reference.

Typical workflows include:

* incremental searches for the maximum achievable order;
* refinement of an existing rule from Cartesian coordinates;
* refinement from barycentric coordinates in a fundamental domain;
* optimization from systematically generated initial grids;
* export of optimized seed points or complete grids.

## Precomputed grids: the AIG libraries

For almost all applications, particularly the use of angular grids in electronic-structure programs, the **precomputed AIG libraries** are more useful than the `aigg` optimizer itself:

[**http://knizia.me/aig**](http://knizia.me/aig)

The AIG libraries contain curated integration rules generated by `aigg`, which were tested, and selected, for stability and numerical performance. Interfaces are provided for C++98, Fortran (F90), and Python, including variants for basic `fp64` precision (“double precision”), and variants featuring full `fp128` precisions.

They allow applications to:

* enumerate the available grids;
* query properties such as point count \(N\), order \(L\), symmetry, and weight spread;
* instantiate the explicit grid points and weights.

### Weight convention

The tabulated AIG-library weights are normalized to

$$ \sum_g w_g = 1, $$

rather than \(4\pi\). Consequently, a full solid-angle integral is represented as

$$ \int_{S^2} f(\mathbf r)\,d\Omega \approx 4\pi \sum_g w_g f(\mathbf r_g). $$

This convention applies to the distributed AIG libraries; consult the full documentation for details and file-format conventions.

## Documentation

Full documentation, including background on integration rules, what `aigg` actually does — any why —, and info on point groups, symmetries, etc, is available at:

[**https://knizia.me/aigg/**](http://knizia.me/aigg)

It includes:

* mathematical background on spherical integration rules;
* the definition and interpretation of grid order;
* spherical harmonics and Cartesian-monomial target spaces;
* point-group-supported integration rules;
* fundamental domains, seed points, and symmetry orbits;
* optimization and incremental-search strategies;
* input and export formats;
* examples for generating and refining grids;
* documentation of the tabulated AIG libraries.

The long-form documentation is intended both as a user manual and as background material explaining the mathematical construction implemented by `aigg`.

## Citation

Citation information for the scientific article describing aigg and the associated integration rules will be added here when available.

If you use aigg or the distributed AIG grids in scientific work, a citation to the associated publication is appreciated.

## License

The **aigg program itself is distributed under GPLv2**.

The separately distributed **AIG libraries use a different, permissive license** intended to allow their incorporation into commercial, non-commercial, closed-source, and open-source software without placing the host program under the GPL.

Please see the respective license files for the exact terms. The license of `aigg` should not be assumed to apply to the separately distributed AIG libraries.
