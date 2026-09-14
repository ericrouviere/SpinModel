
###################################################
##### Magnetization and participation #############
###################################################

function computeMagnetization(seq::Sequence, K::Table, Q::Settings, ligs::Ligands)
    # return the magnetization at each site for each ligand perturbation.
    W, L = Q.W, Q.L
    J = seq2J(seq, K)
    Z_right = Vector{Float64}(undef, 2^W)  
    Z_left  = Vector{Float64}(undef, 2^W)
    buff_vec = Vector{Float64}(undef, 2^W)
    Mmat = [zeros(W,L+1) for i in 1:length(ligs)]
    J_mat = Matrix{Float64}(undef, W, W) #
    J_mat_states = Matrix{Float64}(undef, W, 2^W) #
    h_states = Vector{Float64}(undef, 2^W) #
    T = [Matrix{Float64}(undef, 2^W, 2^W) for i in 1:(L+1)]
    states = getSpinStates(W)

    function mul_update!(B::AbstractArray, buff::AbstractArray, A::AbstractArray)
        # does the following opperation inplace: B <-- A*B
        mul!(buff, A, B)
        B .= buff
        return nothing
    end

    function computeMags!(j::Int)
        # compute all transfer matrices
        for k in 1:L+1 # each layer k
            updateTransferMatrix!(T[k], J_mat, J_mat_states, h_states, states, J, k)
        end

        # compute magnetizations
        for l = 1:(L+1) # for each layer
            Z_right .= 1.0
            for k in 1:(l-1)
                mul_update!(Z_right, buff_vec, T[k])
            end
            Z_left .= 1.0
            for k in (L+1):-1:l
                mul_update!(Z_left, buff_vec, T[k]')
            end
            Z = Z_left' * Z_right
            for i = 1:W
                @views buff_vec .= Z_left .* states[i,:] .* Z_right 
                Mmat[j][i,l] = sum(buff_vec) / Z
            end
        end
        return nothing
    end

    for j in 1:length(ligs) # for each perturbation
        perturb = ligs[j]
        J[perturb.sites] += perturb.fields # apply ligand perturbation
        computeMags!(j)
        J[perturb.sites] -= perturb.fields # remove ligan perturbation
    end
    return Mmat
end


function computeMagnetization(seq::Sequence, K::Table, Q::Settings)
    return computeMagnetization(seq, K, Q, ligs)
end

function computePartRatio(seq::Sequence, K::Table, Q::Settings, ligs::Ligands)
    # compute participation ratio upon binding each ligand and return the 
    # largest PR.
    M = computeMagnetization(seq, K, Q, ligs)
    return maximum([ mean(0.5*abs.(M[i] .- M[1])) for i in 2:length(M)])
end


######################################################
## Deep Mutational Scan Functions ######################
######################################################


"""
Binding energies of each ligand for every single mutant of `seq`, as an array of size
`(q, length(seq), length(ligs)-1)`; WT residues are left at zero.

Scans the mutants through a `TransferCache` when the ligands allow it, which rebuilds
only the one or two transfer matrices each mutation changes instead of all `L+1`.
`_computeBindingDMSDirect` is the same calculation with every free energy computed
from scratch; the two agree to floating-point round-off.
"""
function computeBindingDMS(seq::Sequence, K::Table, Q::Settings, ligs::Ligands)
    @assert length(ligs) > 1
    _isFieldsOnLastLayer(ligs, Q.L) || return _computeBindingDMSDirect(seq, K, Q, ligs)
    q = Q.q
    DMS = zeros(q, length(seq), length(ligs)-1)
    E_wt = computeFreeEnergies(seq, K, Q, ligs)
    ΔE_wt = (E_wt .- E_wt[1])[2:end]
    cache = TransferCache(seq, K, Q)
    E = Vector{Float64}(undef, length(ligs))
    for j in eachindex(seq)
        type_wt = seq[j]
        for type_mut in 1:q
            type_mut == type_wt && continue
            mutantFreeEnergies!(E, cache, seq, K, ligs, j, type_mut)
            @views DMS[type_mut, j, :] .= (E .- E[1])[2:end] .- ΔE_wt
        end
    end
    return DMS
end


function _computeBindingDMSDirect(seq::Sequence, K::Table, Q::Settings, ligs::Ligands)
    # return the binding energies of the each ligand for
    # each single mutant.
    # This is the fall back, slower method to use when ligand
    # feilds are not on last layer. 
    @assert length(ligs) > 1
    q = Q.q
    DMS = zeros(q, length(seq), length(ligs)-1)
    E_wt = computeFreeEnergies(seq, K, Q, ligs)
    ΔE_wt = (E_wt .- E_wt[1])[2:end]
    for j in eachindex(seq)
        type_wt = seq[j] # save WT aa
        for type_mut in 1:q
            type_mut == type_wt && continue
            seq[j] = type_mut # mutate seq
            E = computeFreeEnergies(seq, K, Q, ligs)
            DMS[type_mut, j, :]  .= (E .- E[1])[2:end] .- ΔE_wt
        end
        seq[j] = type_wt # revert back to WT
    end
    return DMS
end


"""
1. count the number of mutatant that change free energy by at least x.
2. count the number of positions that have at least 1 mutation that
changes fitness by at least x.
3. Measure the distance between the active site and the fartest mutation
that changes energies by at least x.
W is the number sites in a single layer, ie width of lattice.

`W` and `actSiteLayer` are required: sequence position `j` sits in layer
`ceil(j/W)`, so the returned distance is meaningless if `W` does not match the
lattice the scan came from. They used to default to the `W=5`, `L=10` lattice,
which silently returned wrong distances for every other geometry. Prefer the
`ligs` method below, which derives both from the data.
"""
function analyseBindingDMS(dms::AbstractArray{T,3}, # dims of (N_types, N_sites, N_ligs)
        thresh::Number;
        W::Integer,
        actSiteLayer::Integer) where T <: AbstractFloat

    N_types, N_sites, N_ligs = size(dms)
    # A scan covers every site of the lattice, so W must divide its length.
    @assert N_sites % W == 0 "W=$(W) does not divide the $(N_sites) scanned sites"
    N_sensMuts = N_sensPos = maxDist = 0
    acc = 0.0
    for j in axes(dms,2)
        pos = false
        dist = actSiteLayer - Int(ceil(j/W))
        for i in axes(dms,1)
            Δϕ = norm(dms[i,j,:])
            Δϕ2 = Δϕ^2
            acc += Δϕ2
            if Δϕ > thresh
                N_sensMuts += 1
                pos = true
                dist > maxDist && (maxDist = dist)
            end
        end
        pos && (N_sensPos +=1) 
    end
    N_muts = (N_types-1)*N_sites
    evolvability = sqrt(acc/N_muts)
    return N_sensMuts, N_sensPos, evolvability, maxDist
end





"""
`analyseBindingDMS` with the lattice geometry taken from the ligands and the scan
itself, so no caller has to supply it.

The perturbed layer gives `actSiteLayer`, and since a scan covers all `W*(L+1)`
sites and the ligands sit on the last layer (`actSiteLayer = L+1`, as
`_isFieldsOnLastLayer` requires of the fast free-energy path), the width follows
as `N_sites / actSiteLayer`. Pass `W` explicitly for the unusual case of a
perturbation that is *not* on the last layer, where that division does not hold.
"""
function analyseBindingDMS(dms::AbstractArray{T,3}, # dims of (N_types, N_sites, N_ligs)
                           ligs::Ligands,
                           thresh::Number;
                           W::Union{Integer,Nothing}=nothing) where T <: AbstractFloat
    # this should extract the layer of all perturbations of all ligands.
    bindingSiteLayers =  getindex.(vcat([Tuple.(ligs.perturbs[i].sites) for i in 1:length(ligs)]...), 2)
    @assert length(unique(bindingSiteLayers)) == 1 # make sure all site are on the same layer
    actSiteLayer = bindingSiteLayers[1]
    N_sites = size(dms, 2)
    if W === nothing
        @assert N_sites % actSiteLayer == 0 "cannot infer W: $(N_sites) sites is not a " *
            "multiple of the $(actSiteLayer) layers; pass W explicitly"
        W = N_sites ÷ actSiteLayer
    end
    return analyseBindingDMS(dms, thresh; W, actSiteLayer)
end



"""
Single and double mutant binding energy changes of `seq`, returned as `(dms, dms2)`
with dims `(q, L, N_ligs)` and `(q, L, q, L, N_ligs)`.

This is the most expensive routine in the package — one call is ~100 deep mutational
scans — so it scans each single-mutant background through its own `TransferCache`
when the ligands allow it. `_computeBindingDoubleDMSDirect` is the same calculation
without the cache; the two agree to floating-point round-off.
"""
function computeBindingDoubleDMS(seq::Sequence, K::Table, Q::Settings, ligs::Ligands)
    @assert length(ligs) > 1
    _isFieldsOnLastLayer(ligs, Q.L) ||
        return _computeBindingDoubleDMSDirect(seq, K, Q, ligs)
    q = Q.q
    L = length(seq)
    N_ligs = length(ligs) - 1
    dms = computeBindingDMS(seq, K, Q, ligs)
    dms2 = zeros(q, L, q, L, N_ligs)
    E_wt = computeFreeEnergies(seq, K, Q, ligs)
    ΔE_wt = (E_wt .- E_wt[1])[2:end]

    # fill diagonal with single mutant effects
    for i in 1:L, a in 1:q, b in 1:q
        dms2[a, i, b, i, :] .= dms[b, i, :]
    end

    E = Vector{Float64}(undef, length(ligs))
    for i in 1:L
        type_wt_i = seq[i]
        for a in 1:q
            a == type_wt_i && continue
            seq[i] = a
            # one cache per single-mutant background, reused by every second mutation
            cache = TransferCache(seq, K, Q)
            for j in (i+1):L
                type_wt_j = seq[j]
                for b in 1:q
                    b == type_wt_j && continue
                    mutantFreeEnergies!(E, cache, seq, K, ligs, j, b)
                    for l in 1:N_ligs
                        ΔΔE = (E[l+1] - E[1]) - ΔE_wt[l]
                        dms2[a, i, b, j, l] = ΔΔE
                        dms2[b, j, a, i, l] = ΔΔE  # symmetry
                    end
                end
            end
            seq[i] = type_wt_i
        end
    end
    return dms, dms2
end

function _computeBindingDoubleDMSDirect(seq::Sequence, K::Table, Q::Settings, ligs::Ligands)
    # return (dms, dms2): single and double mutant binding energy changes.
    # dms  dims: (q, L, N_ligs)
    # dms2 dims: (q, L, q, L, N_ligs)
    # exploits symmetry: only computes (i<j) pairs, mirrors to (j,i).
    # diagonal (i==j): DMS2[a,i,b,i,:] = dms[b,i,:] (single mutant effect of b at i)
    @assert length(ligs) > 1
    q = Q.q
    L = length(seq)
    N_ligs = length(ligs) - 1
    dms = _computeBindingDMSDirect(seq, K, Q, ligs)
    dms2 = zeros(q, L, q, L, N_ligs)
    E_wt = computeFreeEnergies(seq, K, Q, ligs)
    ΔE_wt = (E_wt .- E_wt[1])[2:end]

    # fill diagonal with single mutant effects
    for i in 1:L, a in 1:q, b in 1:q
        dms2[a, i, b, i, :] .= dms[b, i, :]
    end

    for i in 1:L
        type_wt_i = seq[i]
        for a in 1:q
            a == type_wt_i && continue
            seq[i] = a
            for j in (i+1):L
                type_wt_j = seq[j]
                for b in 1:q
                    b == type_wt_j && continue
                    seq[j] = b
                    E = computeFreeEnergies(seq, K, Q, ligs)
                    ΔΔE = (E .- E[1])[2:end] .- ΔE_wt
                    dms2[a, i, b, j, :] .= ΔΔE
                    dms2[b, j, a, i, :] .= ΔΔE  # symmetry
                    seq[j] = type_wt_j
                end
            end
            seq[i] = type_wt_i
        end
    end
    return dms, dms2
end


####################################################
##### functions for Coarse-graining sequence space #
####################################################

function computeSingleMutantEvolvabilities(seq::Sequence, K::Table, Q::Settings,
                                           ligs::Ligands,
                                           dms::AbstractArray, dms2::AbstractArray;
                                           ΔF_0::Real=1)
    # For each single mutant, compute its ΔΔF12 (= ΔF[1]-ΔF[2]) and its evolvability.
    # output: ΔΔF12s (q, L) with NaN for out-of-range and WT amino acids;
    #         evolvabilities (q, L) with zero for WT amino acids.
    # ΔF_0 is the scale of the binding energies, max|ℓ_k - ℓ_0| over the ligands;
    # a mutant with max|ΔF| > ΔF_0 is out of range.
    q = Q.q
    L = length(seq)
    N_ligs = length(ligs) - 1

    E_wt = computeFreeEnergies(seq, K, Q, ligs)
    ΔF_wt = (E_wt .- E_wt[1])[2:end]  # (N_ligs,)

    ΔΔF12s       = fill(NaN, q, L)
    evolvabilities = zeros(q, L)

    for i in 1:L
        type_wt_i = seq[i]
        for a in 1:q
            a == type_wt_i && continue

            ΔF_mut = ΔF_wt .+ dms[a, i, :]  # (N_ligs,)
            if maximum(abs, ΔF_mut) <= ΔF_0
                ΔΔF12s[a, i] = ΔF_mut[1] - ΔF_mut[2]
            end

            # DMS from this single mutant background: effect of each second mutation
            sm_dms = dms2[a, i, :, :, :] .- reshape(dms[a, i, :], 1, 1, N_ligs)  # (q, L, N_ligs)

            _, _, evolvabilities[a, i], _ = analyseBindingDMS(sm_dms, 0.0;
                                                              W=Q.W, actSiteLayer=Q.L+1)
        end
    end

    return ΔΔF12s, evolvabilities
end

function classifyMutantState(ΔΔF::Real, ev::Number,
                             ΔΔF12_hi_thresh::Real, ΔΔF12_lo_thresh::Real, EV_thresh::Number)
    # Returns state 1-6: 1-3 = low EV, 4-6 = high EV; within each, 1/4=right, 2/5=mid, 3/6=left
    # The two ΔΔF12 thresholds are the boundaries of the three binding-space regions.
    bin_region = ΔΔF > ΔΔF12_hi_thresh ? 1 : ΔΔF < ΔΔF12_lo_thresh ? 3 : 2
    ev_region  = ev >= EV_thresh ? 2 : 1
    return bin_region + (ev_region - 1) * 3
end

function countMutantStates(seq::Sequence, ΔΔF12s::AbstractMatrix{<:Real},
                            evolvabilities::AbstractMatrix,
                            ΔΔF12_hi_thresh::Real, ΔΔF12_lo_thresh::Real, EV_thresh::Number)
    # Count single mutants in each of 6 states; ΔΔF12s/evolvabilities are (q, L).
    # NaN entries in ΔΔF12s are treated as out-of-range and skipped.
    state_n = zeros(Int, 6)
    for i in eachindex(seq)
        for a in axes(ΔΔF12s, 1)
            a == seq[i]        && continue  # skip WT
            isnan(ΔΔF12s[a, i]) && continue  # skip out-of-range
            state_n[classifyMutantState(ΔΔF12s[a, i], evolvabilities[a, i],
                                        ΔΔF12_hi_thresh, ΔΔF12_lo_thresh, EV_thresh)] += 1
        end
    end
    return state_n
end

function computeAllMutantStates(seqs::Vector{Sequence}, K::Table, Q::Settings,
                                 ligs::Ligands,
                                 ΔΔF12_hi_thresh::Real, ΔΔF12_lo_thresh::Real,
                                 EV_thresh::Number; ΔF_0::Real=1)
    # Returns a 6×6 matrix A where A[i,j] is the mean fraction of single mutants
    # that take a sequence from state i to state j, averaged over seqs.
    # ΔF_0 is the scale of the binding energies; max|ΔF| > ΔF_0 is out of range.

    N_states = 6
    function analyze(seq, K, Q, ligs, ΔΔF12_hi_thresh, ΔΔF12_lo_thresh, EV_thresh, ΔF_0)
        dms, dms2 = computeBindingDoubleDMS(seq, K, Q, ligs)

        E_wt  = computeFreeEnergies(seq, K, Q, ligs)
        ΔF_wt = (E_wt .- E_wt[1])[2:end]
        wt_ΔΔF12 = maximum(abs, ΔF_wt) <= ΔF_0 ? ΔF_wt[1] - ΔF_wt[2] : NaN

        _, _, wt_ev, _ = analyseBindingDMS(dms, 0.0; W=Q.W, actSiteLayer=Q.L+1)
        wt_state = isnan(wt_ΔΔF12) ? 0 :
                    classifyMutantState(wt_ΔΔF12, wt_ev,
                                        ΔΔF12_hi_thresh, ΔΔF12_lo_thresh, EV_thresh)

        mut_ΔΔF12s, mut_ev = computeSingleMutantEvolvabilities(seq, K, Q, ligs, dms, dms2; ΔF_0)
        state_n = countMutantStates(seq, mut_ΔΔF12s, mut_ev,
                                    ΔΔF12_hi_thresh, ΔΔF12_lo_thresh, EV_thresh)

        return (wt_state=wt_state, state_n=state_n)
    end

    results = pmap(seq -> analyze(seq, K, Q, ligs,
                                  ΔΔF12_hi_thresh, ΔΔF12_lo_thresh, EV_thresh, ΔF_0), seqs)

    A = zeros(N_states, N_states)
    x = zeros(N_states)
    for r in results
        r.wt_state == 0 && continue
        A[r.wt_state, :] .+= r.state_n
        x[r.wt_state] += 1
    end
    N_single_muts = (Q.q - 1) * length(seqs[1])
    A ./= N_single_muts
    A ./= x
    return A, x
end



####################################################
##### Coarse-graining a recorded population trajectory
####################################################

"""
Binding energies and evolvability of a single sequence.
Returns `(ΔF, ev)` with `ΔF = F_bound - F_solvent` for each non-solvent ligand and `ev`
the evolvability from a deep mutational scan (the scan itself is discarded).
`analyseBindingDMS` computes the evolvability as `sqrt(mean‖dms‖²)`, which does not
depend on its threshold argument, so `0.0` is passed here (as in
`computeSingleMutantEvolvabilities`).
"""
function computeSeqPhenotype(seq::Sequence, K::Table, Q::Settings, ligs::Ligands)
    E = computeFreeEnergies(seq, K, Q, ligs)
    ΔF = (E .- E[1])[2:end]
    dms = computeBindingDMS(seq, K, Q, ligs)
    _, _, ev, _ = analyseBindingDMS(dms, 0.0; W=Q.W, actSiteLayer=Q.L+1)
    return ΔF, ev
end

"""
Index the distinct sequences of one or more recorded population trajectories.
Returns `(uniqueSeqs, maps)` where `maps[r][j,t]` is the index into `uniqueSeqs` of the
sequence `popTrajs[r][:,j,t]`. A deep mutational scan costs far more than hashing a
sequence, so scanning only the distinct sequences is a large saving whenever the
population is clonal or the mutation rate is low.
"""
function _indexUniqueSeqs(popTrajs::Vector{<:AbstractArray{<:Integer,3}})
    inds = Dict{Sequence,Int}()
    uniqueSeqs = Sequence[]
    maps = [Matrix{Int}(undef, size(pt,2), size(pt,3)) for pt in popTrajs]
    for r in eachindex(popTrajs)
        pt = popTrajs[r]
        for t in axes(pt,3), j in axes(pt,2)
            seq = Sequence(pt[:,j,t])
            maps[r][j,t] = get!(inds, seq) do
                push!(uniqueSeqs, seq)
                length(uniqueSeqs)
            end
        end
    end
    return uniqueSeqs, maps
end

"""
Binding energies and evolvabilities of every sequence in the recorded population
trajectories of an ensemble (output of `evolvePop_saveSeqsEnsemble`).
Returns `(ΔFs, EVs)`, one entry per replicate:
- `ΔFs[r]::Array{Float64,3}` of size `(N_ligs, P, N_gen)`
- `EVs[r]::Matrix{Float64}`  of size `(P, N_gen)`
Distinct sequences are scanned once and the results scattered back, so repeated
sequences (clones, and unmutated lineages across generations) cost nothing extra.
"""
function computeTrajPhenotypesEnsemble(popTrajs::Vector{<:AbstractArray{<:Integer,3}},
                                       K::Table, Q::Settings, ligs::Ligands;
                                       pool::AbstractWorkerPool=default_worker_pool())
    uniqueSeqs, maps = _indexUniqueSeqs(popTrajs)
    out = pmap(x -> computeSeqPhenotype(x, K, Q, ligs), pool, uniqueSeqs)

    N_ligs = length(ligs) - 1
    ΔFs = Vector{Array{Float64,3}}(undef, length(popTrajs))
    EVs = Vector{Matrix{Float64}}(undef, length(popTrajs))
    for r in eachindex(popTrajs)
        P, N_gen = size(maps[r])
        ΔF = Array{Float64,3}(undef, N_ligs, P, N_gen)
        EV = Matrix{Float64}(undef, P, N_gen)
        for t in 1:N_gen, j in 1:P
            ΔF[:,j,t] .= out[maps[r][j,t]][1]
            EV[j,t] = out[maps[r][j,t]][2]
        end
        ΔFs[r], EVs[r] = ΔF, EV
    end
    return ΔFs, EVs
end

"""
Single-trajectory version of `computeTrajPhenotypesEnsemble`.
"""
function computeTrajPhenotypes(popTraj::AbstractArray{<:Integer,3},
                               K::Table, Q::Settings, ligs::Ligands;
                               pool::AbstractWorkerPool=default_worker_pool())
    ΔFs, EVs = computeTrajPhenotypesEnsemble([popTraj], K, Q, ligs; pool)
    return ΔFs[1], EVs[1]
end

"""
Assign each sequence of a recorded trajectory to one of the 6 coarse-grained states
from its binding energies and evolvability (see `classifyMutantState`).
`ΔF` has size `(2, P, N_gen)` and `EV` size `(P, N_gen)`; returns a `(P, N_gen)` matrix
of states. Sequences out of range (`max|ΔF| > ΔF_0`, the scale of the binding energies)
are given state `0`, the convention used in `computeAllMutantStates`.
"""
function classifyTrajStates(ΔF::AbstractArray{<:Real,3}, EV::AbstractMatrix,
                            ΔΔF12_hi_thresh::Real, ΔΔF12_lo_thresh::Real,
                            EV_thresh::Number; ΔF_0::Real=1)
    @assert size(ΔF,1) == 2 # ΔΔF12 needs exactly two binding environments
    @assert size(ΔF)[2:3] == size(EV)
    states = zeros(Int, size(EV))
    for t in axes(EV,2), j in axes(EV,1)
        @views f = ΔF[:,j,t]
        maximum(abs, f) > ΔF_0 && continue # out of range: leave state 0
        states[j,t] = classifyMutantState(f[1] - f[2], EV[j,t],
                                          ΔΔF12_hi_thresh, ΔΔF12_lo_thresh, EV_thresh)
    end
    return states
end

function classifyTrajStates(ΔFs::Vector{<:AbstractArray{<:Real,3}},
                            EVs::Vector{<:AbstractMatrix},
                            ΔΔF12_hi_thresh::Real, ΔΔF12_lo_thresh::Real,
                            EV_thresh::Number; ΔF_0::Real=1)
    return map((ΔF, EV) -> classifyTrajStates(ΔF, EV, ΔΔF12_hi_thresh, ΔΔF12_lo_thresh,
                                              EV_thresh; ΔF_0),
               ΔFs, EVs)
end

"""
Fraction of the population in each coarse-grained state at each generation.
`states` is a `(P, N_gen)` matrix from `classifyTrajStates`; returns an
`(N_states, N_gen)` matrix. Columns sum to at most 1: the deficit is the fraction of
out-of-range (state 0) sequences.
"""
function computeStateOccupancy(states::AbstractMatrix{<:Integer}; N_states::Int=6)
    P, N_gen = size(states)
    occ = zeros(N_states, N_gen)
    for t in 1:N_gen, j in 1:P
        s = states[j,t]
        s == 0 && continue
        occ[s,t] += 1
    end
    return occ ./ P
end

function computeStateOccupancy(states::Vector{<:AbstractMatrix{<:Integer}}; N_states::Int=6)
    return map(x -> computeStateOccupancy(x; N_states), states)
end

"""
Average a per-generation quantity over the complete environmental periods it contains.
`X` has size `(N_rows, N_gen)` (e.g. the occupancy from `computeStateOccupancy`) and
`k_record` holds the environment index of each generation. The phase is taken from
`k_record` (the first generation at which the environment cycle restarts) rather than
assumed, so this is correct whatever `N % τ` is.
Returns an `(N_rows, N_env*τ)` matrix whose columns run over one period, starting at the
first generation of environment 1.
"""
function foldPeriods(X::AbstractMatrix, k_record::AbstractVector{<:Integer},
                     τ::Integer, N_env::Integer)
    period = N_env * τ
    if k_record[1] == 1 # `evolvePop_saveSeqs` already aligns the window
        start = 1
    else
        i = findfirst(i -> k_record[i] == 1 && k_record[i-1] == N_env, 2:length(k_record))
        i === nothing && error("k_record contains no start of an environmental period")
        start = i + 1 # findfirst indexed into 2:length(k_record)
    end
    N_periods = (length(k_record) - start + 1) ÷ period
    N_periods == 0 && error("k_record contains no complete environmental period")
    folded = zeros(size(X,1), period)
    for p in 1:N_periods
        @views folded .+= X[:, (start + (p-1)*period):(start + p*period - 1)]
    end
    return folded ./ N_periods
end

###############################################
## Allosteric Surface Scan functions #############
###############################################
function perturbationScan(seq::Sequence,
                          K::Table,
                          Q::Settings,
                          ligs::Ligands;
                          h_mag::Number=1,
                          fullScan=true) # 
    # compute the change in binding energies upon applying a feild at each site
    # set fullScan to true to perturbe all sites.
    # set fullScan to false to perturb just the sites on the allo surface.

    @assert length(ligs) == 3 # this is temperatory
    W, L = Q.W, Q.L
    fullScan ? P=L+1 : P=1 # P indicates the number of columns of sites to scan.

    # unperturbed system
    energies_0 = computeFreeEnergies(seq, K, Q, ligs)
    ΔF_0 = (energies_0 .- energies_0[1])[2:end] # take differences to get binding energies
    
    # perturbed systems
    energies_pert = zeros(length(ligs), W, P, 2) # 2 because two perturbs: +h and -h.
    for i in 1:W, j in 1:P
        energies_pert[:,i,j,1] = computeFreeEnergies(seq, K, Q, ligs; site_add=CartesianIndex(i,j,1), h_add=h_mag)
        energies_pert[:,i,j,2] = computeFreeEnergies(seq, K, Q, ligs; site_add=CartesianIndex(i,j,1), h_add=-h_mag)
    end
    ΔF_pert = (energies_pert .- reshape(energies_pert[1,:,:,:], 1,  W, P, 2) )[2:end,:,:,:]
    ΔΔF = reshape(ΔF_0, :,1,1,1) .- ΔF_pert 
    return ΔΔF
end

function scanAllostericSurface(seq::Sequence,
                               K::Table,
                               Q::Settings,
                               ligs::Ligands;
                               h_mag::Number=1)
    # compute the effect of applying an allosteric feild on the 
    # binding energy at the active site.
    return perturbationScan(seq, K, Q, ligs; h_mag, fullScan=false)[:,:,1,:]
end

#function getAlloHotSpot(seq::Sequence,
#                        K::Table,
#                        Q::Settings,
#                        sites::Sites,
#                        fields::Fields;
#                        h_perturb::Number=1)
#    # return the most sensitive site on the allosteric surface.
#    ΔG = scanAllostericSurface(seq, K, Q, sites, fields; h_perturb)
#    ΔG_max, ind = findmax(ΔG)
#    site = ind[1]
#    ind[2] == 1 ? field = h_perturb : field = -h_perturb
#    return ΔG_max, site, field
#end

function computeAllostery(seq::Sequence,
                          K::Table,
                          Q::Settings,
                          ligs::Ligands;
                          h_mag::Number=1)
    # Scans the allosteric surface and returns the largest change
    # in binding energy from an allosteric perturbation.
    ΔΔF = scanAllostericSurface(seq, K, Q, ligs; h_mag)
    return maximum(abs, ΔΔF)
end


#############################################################
## Functions to bin sequences in Binding space ##############
#############################################################

function binBindingSpace(bindingEnergies::AbstractMatrix, # has dims (popSize, N_ligs) ??
                         N_bins::Int;
                         N_collect=200,
                         useLims=false,
                         limits=[-2cos(pi/4), 2cos(pi/4)])
    # return vector of vector of indices. Each subvector
    # contains the indices of sequences that belong to the same
    # function bin (ie have the same binding energies to both ligs).
    # N_bins: number of bins to divide function space.
    # N_collect: the max number of indices to keep in each bin.
    # bindingEnergies: Matrix with rows containing the binding
    # energies to lig 1 (col 1) and lig 2 (col 2) of each sequence.
    
    # Rotate binding space 45 degrees
    θ = pi/4
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    rotBE = R * bindingEnergies
    coors = rotBE[1,:]
    m, M = extrema(coors)
    if useLims
        bins = LinRange(limits[1], limits[2], N_bins+1)
    else
        bins = LinRange(m, M, N_bins+1)
    end
    inds = Vector{Int}[]
    for i in 1:N_bins
        tmp_inds = findall( x -> bins[i]<=x<=bins[i+1], coors )
        if length(tmp_inds) > N_collect
            tmp_inds = tmp_inds[1:N_collect]
        end
        push!(inds, tmp_inds)
    end
    return inds
end


"""
Return a vector that stores the bin numbers of each sequence. 
default limits only applies to binding space in [-1, 1] × [-1, 1], Not General Fix this!!
"""

function getBins(bindingEnergies::AbstractMatrix; # has dims (N_ligs, popSize)
                 N_bins=50,
                 θ = pi/4,
                 limits=[-2cos(θ), 2cos(θ)])
    if maximum(abs, bindingEnergies) > 1
        error("Current version of `getBins` only handles binding space in [-1, 1] × [-1, 1]")
    end

    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    rotBE = R * bindingEnergies
    coors = rotBE[1,:]
    m, M = extrema(coors)
    bins = LinRange(limits[1], limits[2], N_bins+1)
    binvec = zeros(Int, size(bindingEnergies,2))
    for i in eachindex(coors)
        binvec[i] = findfirst(y -> y>=coors[i], bins)-1
    end
    return binvec
end

function getBins(bindingEnergies::AbstractArray{T, 3}; # has dims (N_ligs, popSize, N_reps)
                 N_bins=50,
                 θ = pi/4,
                 limits=[-2cos(θ), 2cos(θ)]) where T <: Number
    return mapslices(getBins, bindingEnergies, dims=(1,2))[:,1,:] # size (popSize, numReps)
end

function getBinnedSeqs(K::Table, Q::Settings, ligs::Ligands; N_randseq=100_000,
        N_bins=50, N_collect=200, useLims=false)
    # return N_collect random sequences for each bin in binding space.
    @assert length(ligs)==3 
    # gen rand seqs and measure binding phenotypes.
    randseqs = map( x -> randSeq(Q), 1:N_randseq)
    energies = pmap( x -> computeFreeEnergies(x, K, Q, ligs), randseqs)
    e1 = getindex.(energies, 1)
    e2 = getindex.(energies, 2)
    e3 = getindex.(energies, 3)
    bindEnergies = [e2 - e1 e3 - e1]'
    
    # bin binding space and get indices of seqs in each bin.
    binnedSeqInds = binBindingSpace(bindEnergies, N_bins; N_collect, useLims);
    binnedSeqs = [ randseqs[binnedSeqInds[i]] for i in 1:N_bins ]
    # return 3-Array with dims (seqLen, N_collect, N_bins) rather than Vector{Vector{Sequence}}
    return cat(map(x -> hcat(x...), binnedSeqs)..., dims=3)
end

"""
Returns the counts of each bin number
"""
function makeHistogram(bins::AbstractVector{T}; N_bins::Int) where T<:Integer
    # takes a vector of bins numbers for each sequence.
    # e.g. if bin[1] = 33, that means the first element is in the 33rd bin.
    # Returns a vector with the counts for that bin. 
    # e.g. if bins=[1,3,2,1,1,1] and N_bins=5
    # then counts = [4,1,1,0,0]
    counts = zeros(Int, N_bins)
    for b in bins
        counts[b] += 1
    end
    return counts
end

function makeHistogram(bins::AbstractMatrix{T}; N_bins::Int) where T<:Integer
    # Apply bincount(::AbstractVector{T}) on each col of bins.
    P, N_points = size(bins)
    counts = Matrix{Int}(undef, N_bins, N_points)
    for i in 1:N_points
        @views counts[:, i] = makeHistogram(bins[:,i]; N_bins)
    end
    return counts # has size (N_bins, N_points)
end


"""
Takes an array of binding energies, and returns the histgram along
the binding curve for each entry in the dim=3.
"""
function binding2Histogram(bindingEnergies::AbstractArray{T, 3}; # has dims  (popSize, numLigs, numReps)
                           N_bins=50,
                           θ = pi/4,
                           limits=[-2cos(θ), 2cos(θ)]) where T <: Number
    bins = getBins(bindingEnergies; N_bins, θ, limits)
    return makeHistogram(bins; N_bins)
end




