# SpinModel

A Julia implementation of Olivier Rivoire's spin model from https://doi.org/10.1103/PhysRevE.100.032411

A protein is modeled as a 2d spin-glass lattice of width `W` and length `L+1`. A sequence maps to local fields and couplings through a random lookup table `K`. Free energies are computed by transfer matrices, and evolution is simulated by Metropolis Monte Carlo.

## Requirements

- **Julia 1.11 or newer** (developed on Julia 1.12). The easiest install is [juliaup](https://github.com/JuliaLang/juliaup):
  ```bash
  curl -fsSL https://install.julialang.org | sh    # macOS / Linux
  winget install --name Julia --id 9NJNWW8PVKMN -e # Windows
  ```
- **matplotlib**, used for plotting through [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl). You usually don't need to do anything: PyPlot installs matplotlib into a private Conda Python the first time it loads. If PyCall is set up to use your system Python instead (the default on Linux), install it there with `pip install matplotlib`.

## Installation

SpinModel is not in the Julia General registry, so it installs from its GitHub URL. Pick one of the two options below.

### Option 1: add it as a package

Use this if you only want to use SpinModel from your own scripts. From the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/ericrouviere/SpinModel.git")
```

Or press `]` to open the package prompt and run:

```
pkg> add https://github.com/ericrouviere/SpinModel.git
```

This installs SpinModel and its dependencies into your active environment. The first `using SpinModel` precompiles the package, which takes a few minutes. Update to the latest version with `Pkg.update("SpinModel")`.

It is good practice to give each project its own environment instead of the global one. Activate it before adding:

```julia
using Pkg
Pkg.activate("path/to/my_project")   # creates Project.toml there if needed
Pkg.add(url="https://github.com/ericrouviere/SpinModel.git")
```

### Option 2: clone the repository

Use this if you want to read or change the source code.

```bash
git clone https://github.com/ericrouviere/SpinModel.git
cd SpinModel
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Pkg.instantiate()` installs the exact dependency versions recorded in `Manifest.toml`.

To use your clone from a different environment, add it with `Pkg.develop` instead of `Pkg.add`. The environment then tracks your local files, so edits to the source take effect:

```julia
using Pkg
Pkg.develop(path="path/to/SpinModel")
```

## Loading SpinModel

### In the REPL or a notebook

```julia
using SpinModel
```

When editing the source, load [Revise](https://github.com/timholy/Revise.jl) first so changes are picked up without restarting Julia:

```julia
using Revise
using SpinModel
```

### In a script

If you installed with Option 1, or used `Pkg.develop`, activate that environment and load the package:

```julia
using Pkg
Pkg.activate("path/to/my_project")
using SpinModel
```

If you cloned the repository (Option 2), activate the clone itself. Building the path from `@__DIR__`, the folder of the script, keeps the script working wherever it is run from:

```julia
using Pkg
Pkg.activate(joinpath(@__DIR__, "path/to/SpinModel"))   # path relative to this script
using SpinModel
```

Alternatively, leave `Pkg.activate` out of the script and choose the environment on the command line:

```bash
julia --project=path/to/SpinModel my_script.jl
```

Setting the `JULIA_PROJECT` environment variable to the same path has the same effect.

### With parallel workers

The ensemble functions (`evolvePopEnsemble`, `computeBindingDMSEnsemble`, ...) run on worker processes through `pmap`. Each worker must start in the same environment and load SpinModel:

```julia
using Distributed
addprocs(4; exeflags="--project=$(Base.active_project())")
@everywhere using SpinModel
```

## Quick start

This script builds a random model, computes binding free energies to two ligands, and evolves a sequence to bind both.

```julia
using Random
using SpinModel

Random.seed!(5)

# Model
W, L, q = 5, 10, 5              # lattice width, layers minus one, alphabet size
sh, sJ = 3, 3                   # std of the random fields and couplings
Q = Settings(W, L, q)
K = rand_table(W, L, q, sh, sJ) # sequence → couplings lookup table

# Solvent and two ligands, as fields on the middle spin of the last layer
site = [CartesianIndex(3, L + 1, 1)]
ligs = Ligands([
    Perturbation(site, [0.0]),   # 1: solvent (unbound)
    Perturbation(site, [1.0]),   # 2: ligand 1
    Perturbation(site, [-1.0]),  # 3: ligand 2
])

# Free energies of a random sequence
seq = randSeq(Q)
F = computeFreeEnergies(seq, K, Q, ligs)  # [F_solvent, F_1, F_2]
ΔF = F[2:3] .- F[1]                        # binding free energies; negative = tight binding

# Evolve the sequence to bind both ligands
assay = DoubleBinding()                    # fitness ϕ = min(-ΔF₁, -ΔF₂)
times, fitnesses, sequences = evolve(seq, K, Q, ligs, assay, 1000; ζ=100)

println("ΔF = ", ΔF)
println("fitness: ", fitnesses[1], " → ", fitnesses[end])
```

Other selective pressures are `Stability()`, `Binding()`, `Specificity()`, `Allostery()`, and `NegativeAllostery()`, passed in place of `DoubleBinding()`.
