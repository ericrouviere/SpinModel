function evolvePopEnsemble(K::Table, Q::Settings, ligs::Ligands,
                           evo::EvoParams; N_reps=10,
                           k_init_vec=rand(1:2, N_reps), N_seqs2Keep::Int=evo.P)

    # Evolve the replicate trajectories in parallel over the worker pool.
    # Falls back to serial execution when no workers have been added.
    output = pmap( y -> evolvePop(K, Q, ligs, evo;
                                  k_init=y, N_seqs2Keep), k_init_vec)
    return output # returns Vector{Vector{Sequence}}
end

"""
Seeded, reproducible version of `evolvePopEnsemble` for a fixed table `K`.
Each replicate `i` is evolved after seeding its worker's RNG with `seeds[i]`,
so the ensemble is exactly reproducible: save `seeds` and `k_init_vec`
(see `run.jl`) and re-running yields identical populations regardless of how
the replicates are scheduled across workers.
"""
function evolvePopEnsemble(K::Table, Q::Settings, ligs::Ligands,
                           evo::EvoParams, seeds::AbstractVector{<:Integer};
                           k_init_vec=rand(1:(length(ligs)-1), length(seeds)),
                           N_seqs2Keep::Int=evo.P,
                           pool::AbstractWorkerPool=default_worker_pool())
    @assert length(seeds) == length(k_init_vec)
    # A CachingPool ships K/Q/ligs/evo to each worker once for the whole ensemble.
    output = pmap( (s, y) -> evolvePop(K, Q, ligs, evo, s;
                                       k_init=y, N_seqs2Keep), pool, seeds, k_init_vec)
    return output # returns Vector{Vector{Sequence}}
end

function evolvePopEnsemble(Q::Settings, ligs::Ligands,
                           evo::EvoParams, seeds::Vector, sh, sJ;
                           k_init_vec=rand(1:2, length(seeds)), N_seqs2Keep::Int=evo.P)

    # each element of populations is matrix of sequences,
    # The first column is the inital population.
    # the second column is the final population.
    output = map( (x,y) -> evolvePop(Q, ligs, evo, x, sh, sJ;
                                  k_init=y, N_seqs2Keep), seeds, k_init_vec)
    N_reps = length(seeds)
    tables = Vector{Table}(undef, N_reps)
    populations = Vector{ Vector{Sequence} }(undef, N_reps)
    for i in eachindex(output)
        tables[i], populations[i] = output[i]
    end
    return tables, populations
end

function evolvePop_saveBindingEnsemble(K::Table, Q::Settings, ligs::Ligands,
                                      evo::EvoParams; N_reps=10,
                                      k_init_vec=rand(1:2, N_reps))

    # each element of tables is a Table
    output = map( y -> evolvePop_saveBinding(K, Q, ligs, evo;
                                  k_init=y), k_init_vec)
    return output, k_init_vec 
end

"""
Seeded, reproducible ensemble of `evolvePopStatic` for a fixed table `K`.
Like the seeded `evolvePopEnsemble` but for static (non-fluctuating) selection,
so there is no `k_init_vec`. Each replicate `i` is evolved after seeding its
worker's RNG with `seeds[i]`, so the ensemble is exactly reproducible regardless
of how replicates are scheduled across workers. Returns `Vector{Vector{Sequence}}`
(same layout as `evolvePopEnsemble`, so the analysis ensembles apply unchanged).
"""
function evolvePopStaticEnsemble(K::Table, Q::Settings, ligs::Ligands,
                                 evo::EvoParams, assay::Assay,
                                 seeds::AbstractVector{<:Integer};
                                 N_seqs2Keep::Int=evo.P,
                                 pool::AbstractWorkerPool=default_worker_pool())
    # A CachingPool ships K/Q/ligs/evo to each worker once for the whole ensemble.
    output = pmap( s -> evolvePopStatic(K, Q, ligs, evo, assay, s;
                                        N_seqs2Keep), pool, seeds)
    return output # returns Vector{Vector{Sequence}}
end

"""
Seeded, reproducible ensemble of `evolvePop_saveSeqs` for a fixed table `K`.
Each replicate `i` is evolved after seeding its worker's RNG with `seeds[i]`, so the
recorded trajectories are exactly reproducible regardless of how replicates are
scheduled across workers.
Returns `(popTrajs, k_records)` with one `PopArray` of size `(seqLen, P, N_record)` and
one environment-index vector per replicate.
"""
function evolvePop_saveSeqsEnsemble(K::Table, Q::Settings, ligs::Ligands,
                                    evo::EvoParams, seeds::AbstractVector{<:Integer};
                                    k_init_vec=fill(1, length(seeds)),
                                    N_samplePeriods::Int=4,
                                    pool::AbstractWorkerPool=default_worker_pool())
    @assert length(seeds) == length(k_init_vec)
    # A CachingPool ships K/Q/ligs/evo to each worker once for the whole ensemble.
    output = pmap( (s, y) -> evolvePop_saveSeqs(K, Q, ligs, evo, s;
                                                k_init=y, N_samplePeriods), pool, seeds, k_init_vec)
    popTrajs = [o[1] for o in output]
    k_records = [o[2] for o in output]
    return popTrajs, k_records
end

## Ensemble Analysis functions ################################################

function computePartRatioEnsemble(sequences::Array{Sequence},
                                  K::Table,
                                  Q::Settings,
                                  ligs::Ligands)
    return pmap( x -> computePartRatio(x, K, Q, ligs), sequences)
end

function computePartRatioEnsemble(sequences::Array{T, 3},
                                  K::Table,
                                  Q::Settings,
                                  ligs::Ligands;
                                  pool::AbstractWorkerPool=default_worker_pool()) where T <: Integer
    # the i,j sequence in "sequences" is stored as "sequences[:, i, j]"
    # Pass a CachingPool to ship K/Q/ligs to each worker once instead of per task.
    seqsMat = [sequences[:, i, j] for i in axes(sequences,2), j in axes(sequences,3)]
    return pmap( x -> computePartRatio(x, K, Q, ligs), pool, seqsMat)
end

function computeBindingDMSEnsemble(sequences::Array{T,3},
                                   K::Table,
                                   Q::Settings,
                                   ligs::Ligands;
                                   pool::AbstractWorkerPool=default_worker_pool()) where T <: Integer
    # the i,j sequence in "sequences" is stored as "sequences[:, i, j]"
    # Pass a CachingPool to ship K/Q/ligs to each worker once instead of per task.
    seqsMat = [sequences[:, i, j] for i in axes(sequences,2), j in axes(sequences,3)]
    out = pmap( x -> computeBindingDMS(x, K, Q, ligs), pool, seqsMat);
    dms = zeros(size(out[1])..., size(seqsMat)...)
    for i in axes(seqsMat, 1), j in axes(seqsMat,2)
        dms[:,:,:,i,j] = out[i,j]
    end
    return dms
end

function computeBindingDMSEnsemble(sequences::Vector{Sequence},
                                   K::Table,
                                   Q::Settings,
                                   ligs::Ligands;
                                   pool::AbstractWorkerPool=default_worker_pool())
    # Pass a CachingPool to ship K/Q/ligs to each worker once instead of per task:
    # a scan is only a few ms, so per-task serialisation is a visible cost.
    out = pmap( x -> computeBindingDMS(x, K, Q, ligs), pool, sequences);
    # Splatting thousands of arrays into `cat` costs more than the scans themselves,
    # so fill a preallocated array instead.
    dms = Array{Float64,4}(undef, size(out[1])..., length(out))
    for i in eachindex(out)
        dms[:,:,:,i] = out[i]
    end
    return dms
end


function computeAllosteryEnsemble(sequences::Array{Sequence},
                                  K::Table,
                                  Q::Settings,
                                  ligs::Ligands;
                                  h_mag::Number=1)
    return pmap( x -> computeAllostery(x, K, Q, ligs; h_mag), sequences)
end

function computeAllosteryEnsemble(sequences::Array{T, 3},
                                  K::Table,
                                  Q::Settings,
                                  ligs::Ligands;
                                  h_mag::Number=1,
                                  pool::AbstractWorkerPool=default_worker_pool()) where T <: Integer
    # Pass a CachingPool to ship K/Q/ligs to each worker once instead of per task.
    seqsMat = [sequences[:, i, j] for i in axes(sequences,2), j in axes(sequences,3)]
    return pmap( x -> computeAllostery(x, K, Q, ligs; h_mag), pool, seqsMat)
end

## Populations-level analysis (single flat pmap over all replicates) ###########
# These methods take the full populations (`Vector{Vector{Sequence}}`) and run
# one `pmap` over every sequence at once. Compared with broadcasting the
# per-replicate methods, this (a) ships K/Q/ligs to each worker only once when a
# CachingPool is supplied, and (b) forms a single balanced task pool instead of
# N_reps serial rounds. The per-replicate `Vector{Vector{...}}` layout is
# reconstructed so callers/output files are unchanged.

"""
Apply `f` to every sequence across all replicates with a single `pmap`, then
split the flat result back into one vector per replicate (preserving order).
"""
function _mapEnsembleFlat(f, populations::Vector{Vector{Sequence}},
                          pool::AbstractWorkerPool)
    lens = length.(populations)
    flat = reduce(vcat, populations)
    res = pmap(f, pool, flat)
    out = Vector{Vector{eltype(res)}}(undef, length(populations))
    stop = 0
    for i in eachindex(populations)
        out[i] = res[(stop + 1):(stop + lens[i])]
        stop += lens[i]
    end
    return out
end

function computePartRatioEnsemble(populations::Vector{Vector{Sequence}},
                                  K::Table, Q::Settings, ligs::Ligands;
                                  pool::AbstractWorkerPool=default_worker_pool())
    return _mapEnsembleFlat(x -> computePartRatio(x, K, Q, ligs), populations, pool)
end

function computeAllosteryEnsemble(populations::Vector{Vector{Sequence}},
                                  K::Table, Q::Settings, ligs::Ligands;
                                  h_mag::Number=1,
                                  pool::AbstractWorkerPool=default_worker_pool())
    return _mapEnsembleFlat(x -> computeAllostery(x, K, Q, ligs; h_mag), populations, pool)
end

function computeBindingDMSEnsemble(populations::Vector{Vector{Sequence}},
                                   K::Table, Q::Settings, ligs::Ligands;
                                   pool::AbstractWorkerPool=default_worker_pool())
    perRep = _mapEnsembleFlat(x -> computeBindingDMS(x, K, Q, ligs), populations, pool)
    # match the per-replicate DMS layout: (dms..., N_seqs2Keep) for each replicate
    return [cat(rep..., dims=4) for rep in perRep]
end

"""
Binding energies `ΔF` of every kept sequence to each non-solvent environment.
Returns one `N_seqs2Keep × (length(ligs)-1)` matrix per replicate, column `j`
being `F_bound(j) - F_solvent`.
"""
function computeBindEnergiesEnsemble(populations::Vector{Vector{Sequence}},
                                     K::Table, Q::Settings, ligs::Ligands;
                                     pool::AbstractWorkerPool=default_worker_pool())
    f = x -> (e = computeFreeEnergies(x, K, Q, ligs); e[2:end] .- e[1])
    perRep = _mapEnsembleFlat(f, populations, pool)
    return [permutedims(reduce(hcat, rep)) for rep in perRep]
end

function perturbationScanEnsemble(sequences::Array{Sequence},
                                  K::Table,
                                  Q::Settings,
                                  ligs::Ligands;
                                  h_mag::Number=1)
    return map( x -> perturbationScan(x, K, Q, ligs; h_mag), sequences)
end

function perturbationScanEnsemble(sequences::Array{T, 3},
                                  K::Table,
                                  Q::Settings,
                                  ligs::Ligands;
                                  h_mag::Number=1) where T <: Integer
    seqsMat = [sequences[:, i, j] for i in axes(sequences,2), j in axes(sequences,3)]
    return pmap( x -> perturbationScan(x, K, Q, ligs; h_mag), seqsMat)
end


"""
`analyseBindingDMS` over an ensemble. `W` and `actSiteLayer` are required for the
same reason as there — see its docstring — and the `ligs` methods below derive
them for you.
"""
function analyseBindingDMSEnsemble(dms::AbstractArray{<:Number, 5},
                                   thresh::Number;
                                   W::Integer, actSiteLayer::Integer)
    N_types, N_nodes, N_ligands, popSize, N_reps = size(dms)
    N_sensMutsArray = Matrix{Int64}(undef, popSize, N_reps)
    N_sensPosArray = similar(N_sensMutsArray)
    maxDistArray = similar(N_sensMutsArray)
    evolvabilityArray = Matrix{Float64}(undef, popSize, N_reps)
    for i in 1:popSize, j in 1:N_reps
        @views dm = dms[:,:,:,i,j]
        N_sensMuts, N_sensPos, evolvability, maxDist = analyseBindingDMS(dm, thresh; W, actSiteLayer)
        N_sensMutsArray[i,j] = N_sensMuts
        N_sensPosArray[i,j] = N_sensPos
        evolvabilityArray[i,j] = evolvability
        maxDistArray[i,j] = maxDist
    end
    return N_sensMutsArray, N_sensPosArray, evolvabilityArray,  maxDistArray
end


function analyseBindingDMSEnsemble(dms::AbstractArray{<:Number, 4},
                                   thresh::Number;
                                   W::Integer, actSiteLayer::Integer)
    return analyseBindingDMSEnsemble(reshape(dms, size(dms)..., 1),  thresh; W, actSiteLayer)
end

function analyseBindingDMSEnsemble(dms::AbstractArray{<:Number, 5},
                                   ligs::Ligands,
                                   thresh::Number;
                                   W::Union{Integer,Nothing}=nothing)
    N_types, N_nodes, N_ligands, popSize, N_reps = size(dms)
    N_sensMutsArray = Matrix{Int64}(undef, popSize, N_reps)
    N_sensPosArray = similar(N_sensMutsArray)
    maxDistArray = similar(N_sensMutsArray)
    evolvabilityArray = Matrix{Float64}(undef, popSize, N_reps)
    for i in 1:popSize, j in 1:N_reps
        @views dm = dms[:,:,:,i,j]
        N_sensMuts, N_sensPos, evolvability, maxDist = analyseBindingDMS(dm, ligs, thresh; W)
        N_sensMutsArray[i,j] = N_sensMuts
        N_sensPosArray[i,j] = N_sensPos
        evolvabilityArray[i,j] = evolvability
        maxDistArray[i,j] = maxDist
    end
    return N_sensMutsArray, N_sensPosArray, evolvabilityArray,  maxDistArray
end

function analyseBindingDMSEnsemble(dms::AbstractArray{<:Number, 4},
                                   ligs::Ligands,
                                   thresh::Number;
                                   W::Union{Integer,Nothing}=nothing)
    return analyseBindingDMSEnsemble(reshape(dms, size(dms)..., 1), ligs, thresh; W)
end

function computeOverlapsEnsemble(sequences::Array{Sequence},
                         K::Table,
                         Q::Settings;
                         ligs::Ligands)
    out = map( x -> computeOverlaps(x, K, Q, ligs), sequences)
    OC = zeros(length(sequences[1]),size(sequences)...)
    vC = similar(OC)
    for i in CartesianIndices(sequences)
        OC[:,i] .= out[i][1] 
        vC[:,i] .= out[i][2] 
    end
    return OC, vC
end


##############################

function buildAndAnalyze(tables::Vector{Table},
                         populations::Vector{Vector{Sequence}},
                         Q::Settings,
                         ligs::Ligands,
                         assay::Assay)
    # build model from sequence and table then assay it for some phenotype.
    N_reps = length(populations)
    popSize = length(populations[1])
    spec = zeros(popSize, N_reps)
    for j in 1:N_reps
        seqs = populations[j]
        K = tables[j]
        for i in 1:popSize
            seq = seqs[i]
            spec[i,j] = computeFitness(seq, K, Q, ligs, assay)
        end
    end
    return spec
end

function buildAndAnalyze(K::Table,
                         populations::Vector{Vector{Sequence}},
                         Q::Settings,
                         ligs::Ligands,
                         assay::Assay)
    tables = fill(K, length(populations))
    return buildAndAnalyze(tables, populations, Q, ligs, assay)
end



## Helper functions ###############

function convertPopulations2PopArray(populations::Vector{Vector{Sequence}})::PopArray
    N_reps = length(populations)
    popSize = size(populations[1],1)
    seqLength = length(populations[1][1])
    popArray = Array{Int8,3}(undef, seqLength, popSize, N_reps)
    for l in 1:N_reps
        popArray[:,:,l] = hcat(populations[l]...)
    end
    return popArray
end

function convertPopArray2Populations(popArray::PopArray)::Vector{Vector{Sequence}}
    return [ Sequence[ eachcol(popArray[:,:,i])... ]  for i in 1:size(popArray, 3) ]
    return [eachcol( popArray[:,:,i])...]
end

function convertTables2TableArray(tables::Vector{Table})::TableArray
    return cat(tables..., dims=5)
end

function convertTableArray2Tables(tableArray::TableArray)::Vector{Table}
    Table[eachslice(tableArray, dims=5)...]
end
