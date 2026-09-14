# Olivier Rivoire's code rewritten by Eric Rouviere in julia 1.7.


"""
True when every ligand perturbation is a field on the last layer, the case
`computeFreeEnergiesFast` (and the `TransferCache` mutant scan) can handle.
"""
function _isFieldsOnLastLayer(ligs::Ligands, L::Integer)
    allGood = true
    for i in 1:length(ligs)
        per = ligs.perturbs[i]
        for j in eachindex(per.sites)
            site = Tuple(per.sites[j])
            allGood = allGood && site[2]==(L+1) # test last layer
            allGood = allGood && site[3]==1 # test feild
        end
    end
    return allGood
end

"""
Calculates the Free energies of the model for each ligand bound
condition. Dispach the fast or general method depending on the ligand.
"""
function computeFreeEnergies(seq::Sequence,
                             K::Table,
                             Q::Settings,
                             ligs::Ligands;
                             site_add=CartesianIndex(1,1,1),
                             h_add=0)

    if _isFieldsOnLastLayer(ligs, Q.L)
        energies = computeFreeEnergiesFast(seq, K, Q, ligs; site_add, h_add)
    else
        energies = computeFreeEnergiesGeneral(seq, K, Q, ligs; site_add, h_add)
    end
    return energies
end


"""
calculates the Free energies of the model bound to eac
ligand.This is the general version that can work fo
all ligand cases.
"""
function computeFreeEnergiesGeneral(seq::Sequence,
                                    K::Table,
                                    Q::Settings,
                                    ligs::Ligands;
                                    site_add=CartesianIndex(1,1,1),
                                    h_add=0)
    W, L = Q.W, Q.L
    energies = zeros(length(ligs))
    J = seq2J(seq, K)
    J[site_add] += h_add # add extra field to site site_add
    for i in 1:length(ligs)
        perturb = ligs[i]
        J[perturb.sites] += perturb.fields # apply perturbation
        energies[i] =  computeFreeEnergy(J)
        J[perturb.sites] -= perturb.fields # remove perturbation
    end
    J[site_add] -= h_add
    return energies
end


"""
calculates the Free energies of the model bound to eac
ligand. This version can only be used when ligands apply
a field and only to the L+1 layer. 
"""
function computeFreeEnergiesFast(seq::Sequence,
                                    K::Table,
                                    Q::Settings,
                                    ligs::Ligands;
                                    site_add=CartesianIndex(1,1,1),
                                    h_add=0)

    # Ligands can only be feilds and can only be applied on the L+1 layer.
    W, L = Q.W, Q.L
    @assert all(x -> x == 1, last.(Tuple.(ligs.perturbs[1].sites));)
    @assert all(x -> x == L+1, getindex.(Tuple.(ligs.perturbs[1].sites), 2))

    energies = zeros(length(ligs))
    J = seq2J(seq, K)
    J[site_add] += h_add # add extra field to site site_add
    T = Matrix{Float64}(undef, 2^W,2^W) 
    C = Matrix{Float64}(undef, W, W) 
    C_states = Matrix{Float64}(undef, W, 2^W) 
    h_states = Vector{Float64}(undef, 2^W) 
    Z = ones(2^W)
    states = getSpinStates(W)
    
    # Compute the product of the transfer matrices of layer 1 to L
    for layer in 1:L
        updateTransferMatrix!(T, C, C_states, h_states, states, J, layer)
        mul!(h_states, T, Z) # use h_states as buffer
        Z .= h_states
    end

    # deal with the last layer for each different ligand
    for i in 1:length(ligs)
        perturb = ligs[i]
        J[perturb.sites] += perturb.fields # apply perturbation
        updateTransferMatrix!(T, C, C_states, h_states, states, J, L+1) # more writing than needed
        energies[i] = -log(sum(T * Z)) # allocates --------------------------------------
        J[perturb.sites] -= perturb.fields # remove perturbation
    end
    J[site_add] -= h_add # remove the extra field if applied.
    return energies
end

"""
Calculation of free energy by transfer matrices
"""
function computeFreeEnergy(J::Array{Float64,3})
    W = size(J, 1)
    L = size(J, 2) - 1
    T = Matrix{Float64}(undef, 2^W,2^W) # the transfer matrix.
    C = Matrix{Float64}(undef, W, W) #
    C_states = Matrix{Float64}(undef, W, 2^W) #
    h_states = Vector{Float64}(undef, 2^W) #
    Z = ones(2^W)
    states = getSpinStates(W)
    for layer in 1:L+1
        updateTransferMatrix!(T, C, C_states, h_states, states, J, layer)
        mul!(h_states, T, Z)
        Z .= h_states
    end
    return -log(sum(Z)) # taking sum instead of trace due to boundry conditions.
end


"""
Compute the transfer matrix at a given layer k
note: the botzmann factor purposefully omits the normal minus sign, exp( ...).
"""
function updateTransferMatrix!(T::Matrix{Y},
                               C::Matrix{Y},
                               C_states::Matrix{Y},
                               h_states::Vector{Y},
                               states::Matrix{Y},
                               J::Array{Y,3},
                               k::Int) where Y <: AbstractFloat
    W = size(J, 1)
    L = size(J, 2) - 1
    h = J[:,k,1] ###############################################
    mul!(h_states, states', h)

    if k < (L+1) # generic case
        updateCouplingMat!(C, J, k, W)
        mul!(C_states, C, states)
        mul!(T, states', C_states)
        @inbounds for i0 in 1:2^W, i1 in 1:2^W
            T[i1,i0] = exp( T[i1,i0] + h_states[i0] )
        end
    else # special case (only fields). ie k = L+1
        @inbounds for i0 in 1:2^W, i1 in 1:2^W
            i0 == i1 ? T[i1,i0] = exp(h_states[i0]) : T[i1,i0] = 0.0 
        end
    end
    return nothing
end


"""
Returns a W by W^2 matrix where the columns store
all spin configurations of layer. Each spin is -1 or +1.
"""
function getSpinStates(W::Int)
    states = zeros(W, 2^W)
    s = zeros(W)
    dig = zeros(Int, W)
    @inbounds for i in 1:2^W
        digits!(dig, i-1, base=2)
        states[:,i] = 2.0 .* (0.5 .- dig)
    end
    return states
end


"""
Fill the coupling matrix for a given layer.
"""
function updateCouplingMat!(C::Matrix{Y},
                            J::Array{Y,3},
                            k::Int,
                            W::Int) where Y <: AbstractFloat
    for j in 1:W, i in 1:W
        if i==j
            C[i,j] = J[i,k,2]
        elseif j-i==1
            C[i,j] = J[i,k,3]
        else
            C[i,j] = 0.0
        end
    end
    C[W,1] = J[W,k,3]
    return nothing
end


#################################################3
### Might be to put this in the analysis.jl #####
################################################

"""
Generate N_seqs random sequences and compute their binding energies.
m is the mean of the centroid
s is the standard deviation of the points about the centroid.
"""
function computeBindingRandSeqs(seqs::Vector{Vector{T}},
                                K::Table,
                                Q::Settings,
                                ligs::Ligands) where T <: Integer
    @assert length(ligs) == 3
    N_seqs = length(seqs)
    ΔF = zeros(N_seqs, length(ligs)-1)
    for i in 1:N_seqs
        energies = computeFreeEnergies(seqs[i], K, Q, ligs)
        ΔF[i,1] = energies[2] - energies[1]
        ΔF[i,2] = energies[3] - energies[1]
    end
    m, s = computeMeanStdOfPoints(ΔF)
    return ΔF, m, s
end

function computeBindingRandSeqs(N_seqs::Int, K::Table, Q::Settings, ligs::Ligands)
    seqs = map( x -> randSeq(Q), 1:N_seqs);
    return computeBindingRandSeqs(seqs, K, Q, ligs) 
end


"""
Compute the mean and variance of a 
Nx2 matrix where each row is a point in 2D space.
"""
function computeMeanStdOfPoints(points::Matrix{Float64})
    @assert size(points, 2) == 2 "Input must be a Nx2 matrix where each row is a point in 2D space."
    centroid = mean(points, dims=1)
    distances = sqrt.(sum((points .- centroid).^2, dims=2))
    return mean(centroid), std(distances)
end


####################################################
##### Scanning the mutants of one background sequence
####################################################

"""
Transfer matrices and partial products of one *background* sequence, so that the free
energies of its single mutants can be computed without rebuilding every layer.

A mutation at sequence position `p` sits in layer `k = cld(p, W)` and changes only the
fields of layer `k` and the couplings of layers `k-1` and `k`, so only the transfer
matrices of layers `k-1` and `k` differ from the background's. Building one transfer
matrix costs ~4^W exponentials and dominates a free-energy evaluation, while applying
a cached one is a single matrix-vector product; scanning mutants this way is several
times faster than recomputing each from scratch.

Fields:
- `Ts[m]` — background transfer matrix of layer `m`, for `m` in `1:L`
- `P[m]`  — the product of layers `1:m-1` applied to the all-ones vector, so `P[1]`
            is all ones and `P[m+1] = Ts[m] * P[m]`

Only valid when every ligand perturbation is a field on the last layer
(`_isFieldsOnLastLayer`), the same restriction as `computeFreeEnergiesFast`.
Use with `mutantFreeEnergies!`.
"""
struct TransferCache
    W::Int
    L::Int
    states::Matrix{Float64}
    J::Array{Float64,3}           # background couplings
    Ts::Vector{Matrix{Float64}}   # background transfer matrices, layers 1:L
    P::Vector{Vector{Float64}}    # partial products, P[m] = Ts[m-1]...Ts[1] * 1
    Jm::Array{Float64,3}          # scratch: the mutant's couplings
    T::Matrix{Float64}            # scratch: a rebuilt transfer matrix
    C::Matrix{Float64}
    C_states::Matrix{Float64}
    h_states::Vector{Float64}
    Z::Vector{Float64}
    Zt::Vector{Float64}
end

function TransferCache(seq::Sequence, K::Table, Q::Settings)
    W, L = Q.W, Q.L
    states = getSpinStates(W)
    J = seq2J(seq, K)
    T = Matrix{Float64}(undef, 2^W, 2^W)
    C = Matrix{Float64}(undef, W, W)
    C_states = Matrix{Float64}(undef, W, 2^W)
    h_states = Vector{Float64}(undef, 2^W)
    Ts = [Matrix{Float64}(undef, 2^W, 2^W) for _ in 1:L]
    P = [ones(2^W) for _ in 1:(L+1)]
    for m in 1:L
        updateTransferMatrix!(Ts[m], C, C_states, h_states, states, J, m)
        mul!(P[m+1], Ts[m], P[m])
    end
    return TransferCache(W, L, states, J, Ts, P, similar(J), T, C, C_states, h_states,
                         Vector{Float64}(undef, 2^W), Vector{Float64}(undef, 2^W))
end

"""
Free energies of the background sequence of `cache` with position `p` set to residue
`a`, written into `energies` (one entry per ligand). `seq` must be the background
sequence the cache was built from; it is mutated and restored, so it is left unchanged.

Equivalent to `computeFreeEnergies(mutated_seq, K, Q, ligs)` up to floating-point
round-off, but reuses the cached transfer matrices of the layers the mutation does not
touch. Requires ligand perturbations to be fields on the last layer.
"""
function mutantFreeEnergies!(energies::AbstractVector,
                             cache::TransferCache,
                             seq::Sequence,
                             K::Table,
                             ligs::Ligands,
                             p::Integer,
                             a::Integer)
    W, L = cache.W, cache.L
    type_wt = seq[p]
    seq[p] = a
    seq2J!(cache.Jm, seq, K)
    seq[p] = type_wt

    k = cld(p, W)             # layer holding position p
    kstart = max(1, k - 1)    # first layer whose transfer matrix the mutation changes
    Z, Zt = cache.Z, cache.Zt
    Z .= cache.P[kstart]
    for m in kstart:L
        if m == k || m == k - 1
            updateTransferMatrix!(cache.T, cache.C, cache.C_states, cache.h_states,
                                  cache.states, cache.Jm, m)
            mul!(Zt, cache.T, Z)
        else
            mul!(Zt, cache.Ts[m], Z)
        end
        Z, Zt = Zt, Z
    end

    # Last layer: its transfer matrix is diagonal, so the partition function is just
    # the weighted sum below rather than another matrix-vector product.
    for i in 1:length(ligs)
        per = ligs[i]
        for j in eachindex(per.sites)
            cache.Jm[per.sites[j]] += per.fields[j]
        end
        @views mul!(cache.h_states, cache.states', cache.Jm[:, L+1, 1])
        acc = 0.0
        @inbounds for j in eachindex(Z)
            acc += exp(cache.h_states[j]) * Z[j]
        end
        energies[i] = -log(acc)
        for j in eachindex(per.sites)
            cache.Jm[per.sites[j]] -= per.fields[j]
        end
    end
    return energies
end
