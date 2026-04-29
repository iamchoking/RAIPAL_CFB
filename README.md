# RAIPAL CFB

This repository generates crossed four-bar linkage (CFB) solution data for RAIPAL.

The generated C++ header is published through a separate repository,
`railabatkaist/raipal_solution`, which is checked out here as a Git submodule at:

```text
SOLUTION/results/raipal_solution
```

That setup keeps this repository as the producer of the numerical solution while
allowing C++ projects to consume only the generated solution repository as their
own submodule.

## Repository Setup

Clone this repository with submodules:

```bash
git clone --recurse-submodules git@github.com:railabatkaist/RAIPAL_CFB.git
cd RAIPAL_CFB
```

If the repository was already cloned without submodules, initialize them with:

```bash
git submodule update --init --recursive
```

The solution submodule should then be available at:

```bash
ls SOLUTION/results/raipal_solution
```

## Generate The C++ Solution Header

The MATLAB scripts under `SOLUTION/` produce intermediate data files and finally
write the C++ header into the `raipal_solution` submodule.

Run the full pipeline from the `SOLUTION/` directory:

```bash
cd SOLUTION
matlab -batch "CFB_table"
matlab -batch "CFB_poly"
matlab -batch "CFB_format"
```

If `results/cfb_coeffs.mat` is already up to date and only the header needs to
be regenerated, run:

```bash
cd SOLUTION
matlab -batch "CFB_format"
```

The generated file is:

```text
SOLUTION/results/raipal_solution/cfbSolution.hpp
```

## Publish Updated Solutions

Because `raipal_solution` is a submodule, generated solution changes must be
committed in the submodule repository first, then the updated submodule pointer
must be committed in this repository.

Commit and push the generated header:

```bash
cd SOLUTION/results/raipal_solution
git status
git add cfbSolution.hpp README.md
git commit -m "Update CFB solution"
git push
```

Then commit the parent repository changes:

```bash
cd ../../..
git status
git add README.md SOLUTION/CFB_format.m SOLUTION/results/raipal_solution
git commit -m "Update CFB solution output"
git push
```

If this is the first time adding the submodule in a branch, include
`.gitmodules` in the parent commit:

```bash
git add .gitmodules SOLUTION/results/raipal_solution
```

## Consuming The Solution In CMake Projects

C++ repositories should add `raipal_solution` directly as their own submodule:

```bash
git submodule add git@github.com:railabatkaist/raipal_solution.git external/raipal_solution
git commit -m "Add RAIPAL solution submodule"
```

For existing clones:

```bash
git submodule update --init --recursive
```

Since the solution is currently header-only, a consuming CMake project can expose
it with an interface target:

```cmake
add_library(raipal_solution INTERFACE)
target_include_directories(raipal_solution INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/external/raipal_solution
)

target_link_libraries(your_target PRIVATE raipal_solution)
```

Then include the generated header from C++:

```cpp
#include <cfbSolution.hpp>
```

The generated symbols are in the `cfb` namespace.
