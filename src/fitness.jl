# Selective-pressure scoring functions on raw ligand free energies.
# E_s is the solvent (unbound) free energy; the rest are bound states.
stability(E_s) = -E_s
specificity(E_s, E_r, E_w) = min(E_s-E_r, E_w-E_s)
binding(E_s, E_r) = E_s - E_r
doubleBinding(E_s, E_r, E_w) = min(E_s-E_r, E_s-E_w)
allostery(E00, E10, E01, E11) = (E10-E00) - (E11-E01)
negativeAllostery(E00, E10, E01, E11) = -allostery(E00, E10, E01, E11)
binding1(Es, E1, E2) = Es - E1
binding2(Es, E1, E2) = Es - E2


# One immutable singleton type per selective pressure. Selection dispatches on
# the assay instance instead of branching on a string, so adding a pressure is
# just a new struct plus a `computeFitness` method (no central if/elseif to edit).
abstract type Assay end
struct Stability         <: Assay end
struct Binding           <: Assay end
struct DoubleBinding     <: Assay end
struct Specificity       <: Assay end
struct Allostery         <: Assay end
struct NegativeAllostery <: Assay end
struct Binding1          <: Assay end
struct Binding2          <: Assay end

# Scalar fitness from precomputed ligand free energies, dispatched on the assay.
computeFitness(::Stability,         energies) = stability(energies...)
computeFitness(::Binding,           energies) = binding(energies...)
computeFitness(::DoubleBinding,     energies) = doubleBinding(energies...)
computeFitness(::Specificity,       energies) = specificity(energies...)
computeFitness(::Allostery,         energies) = allostery(energies...)
computeFitness(::NegativeAllostery, energies) = negativeAllostery(energies...)
computeFitness(::Binding1,          energies) = binding1(energies...)
computeFitness(::Binding2,          energies) = binding2(energies...)

"""
    computeFitness(seq, K, Q, ligs, assay::Assay; site_add, h_add)

Scalar fitness of `seq` under selective pressure `assay` (an `Assay` singleton,
e.g. `DoubleBinding()`). Computes the ligand free energies and dispatches on
`assay` to the matching scoring function.
"""
function computeFitness(seq::Sequence, K::Table, Q::Settings, ligs::Ligands,
                        assay::Assay; site_add=CartesianIndex(1,1,1), h_add=0)
    energies = computeFreeEnergies(seq, K, Q, ligs; site_add, h_add)
    return computeFitness(assay, energies)
end
