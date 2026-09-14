



"""
    evolve(seq, K, Q, ligs, assay, N; ζ, ΔT=1)

Single-sequence Metropolis Monte Carlo under the selective pressure `assay`.
Each of the `N` steps mutates one random position of the sequence to a different
residue and accepts the mutation with probability `min(1, exp(ζ*(ϕ_mut - ϕ)))`, where
`ϕ = computeFitness(seq, K, Q, ligs, assay)` and `ζ` is the evolutionary inverse
temperature (larger `ζ`, stronger selection).

Records the step, fitness and sequence every `ΔT` steps and at the last step, and
returns `(times, fitnesses, sequences)`. `seq` itself is not modified.
"""
function evolve(seq::Sequence, K::Table, Q::Settings, ligs::Ligands, assay::Assay,
                N::Integer; ζ::Real, ΔT::Integer=1)
    seq = copy(seq)
    ϕ = computeFitness(seq, K, Q, ligs, assay)
    times, fitnesses, sequences = Int[], Float64[], Sequence[]
    for t in 0:N
        # recording:
        if t % ΔT == 0 || t == N
            push!(times, t)
            push!(fitnesses, ϕ)
            push!(sequences, copy(seq))
        end
        t == N && break

        # mutation:
        i, amut = rand(eachindex(seq)), rand(1:(Q.q - 1))
        awt = seq[i]
        seq[i] = ((awt - 1 + amut) % Q.q) + 1
        ϕ_mut = computeFitness(seq, K, Q, ligs, assay)

        # selection:
        if rand() < exp(ζ * (ϕ_mut - ϕ))
            ϕ = ϕ_mut
        else
            seq[i] = awt
        end
    end
    return times, fitnesses, sequences
end


"""
Evolve a popluation of sequences under a fluctuating selection.
Saves numSeqs2Keep from the last half of the trajectory.
"""
function evolvePop(seqs::Vector{Sequence},
                   K::Table,
                   Q::Settings,
                   ligs::Ligands,
                   evo::EvoParams;
                   k_init::Int=1, # Initial enviroment
                   N_seqs2Keep::Int=Int(evo.N/2) )

    N, μ, τ, α = evo.N, evo.μ, evo.τ, evo.α
    N_env = length(ligs) - 1
    @assert iseven(N) 
    @assert 0 <= μ <= 1 # constrain mutation rate to be a probability
    @assert 4τ <= N # Env periods can't be longer that 1/4 the total evo traj.
    @assert 0 < k_init <= N_env # Env index has to be positive
    @assert 2N_seqs2Keep < N # cant keep more seqs than half the number of generations
   
    Q = deepcopy(Q) # detach from outside function
    P = length(seqs) # population
    k = k_init # environment index
    fits = zeros(length(seqs), N_env)
    mutatedRecord = BitVector(undef, P); mutatedRecord .= 1 # keep track of who was mutated
    fitsTemp = similar(fits) # allocate memory
    seqsTemp = deepcopy(seqs) # allocate memory
    seqBucket = Sequence[] # list to save sequences

    # Get list of generations to sample a sequence from.
    gen2StartSampling = Int(N - (floor((N/2)/(2τ))*2τ)) + 1
    gens2Sample = sort((gen2StartSampling:N)[randperm(length(gen2StartSampling:N))[1:N_seqs2Keep]])

    for t in 1:N # generations
        # compute Fitness
        updateAllFitnesses!(fits, seqs, Q, K, ligs, mutatedRecord)

        # apply selection
        applySelection!(seqs, fits, α, k, fitsTemp, seqsTemp)

        # mutation.
        mutatedRecord = mutateAtRate!.(seqs, [μ], [Q]) # mutation all sequences every generation.

        # switch environments
        t % τ == 0 && (k = k % N_env + 1) # increment the env cyclically

        # Sample sequence
        t in gens2Sample && push!(seqBucket, deepcopy(seqs[rand(1:P)]))
    end
    return seqBucket
end

function evolvePop(K::Table, Q::Settings, ligs::Ligands,
                   evo::EvoParams; k_init=1, N_seqs2Keep::Int=evo.P)
    seqs = map(x -> randSeq(Q), 1:evo.P)
    return evolvePop(seqs, K, Q, ligs, evo; k_init, N_seqs2Keep)
end

"""
Seeded variant of `evolvePop` for a fixed table `K`.
Seeds the global RNG with `seed` before evolving, so the entire trajectory
(initial population, mutations, selection, and which sequences are sampled)
is exactly reproducible. Re-running with the same `K`, `evo`, `seed`, and
`k_init` reproduces the returned population bit-for-bit.
"""
function evolvePop(K::Table, Q::Settings, ligs::Ligands,
                   evo::EvoParams, seed::Integer;
                   k_init=1, N_seqs2Keep::Int=evo.P)
    Random.seed!(seed)
    return evolvePop(K, Q, ligs, evo; k_init, N_seqs2Keep)
end

function evolvePop(Q::Settings, ligs::Ligands,
                   evo::EvoParams, seed, sh, sJ;
                   k_init=1, N_seqs2Keep::Int=evo.P)
    Random.seed!(seed)
    K = rand_table(Q.W, Q.L, Q.q, sh, sJ)
    seqs  = evolvePop(K, Q, ligs, evo; k_init, N_seqs2Keep)
    return K, seqs
end

"""
A version of evolvePop but that saves the binding energies of
each sequence every generation.
"""
function evolvePop_saveBinding(seqs::Vector{Sequence},
                               K::Table,
                               Q::Settings,
                               ligs::Ligands,
                               evo::EvoParams;
                               k_init::Int=1 # Initial enviroment
                               )

    N, μ, τ, α = evo.N, evo.μ, evo.τ, evo.α
    N_env = length(ligs) - 1
    @assert iseven(N) 
    @assert 0 <= μ <= 1 # constrain mutation rate to be a probability
    @assert 4τ <= N # Env periods can't be longer that 1/4 the total evo traj.
    @assert 0 < k_init <= N_env # Env index has to be positive
   
    Q = deepcopy(Q) # detach from outside function
    P = length(seqs) # population
    k = k_init # environment index
    fits = zeros(length(seqs), N_env)
    mutatedRecord = BitVector(undef, P); mutatedRecord .= 1 # keep track of who was mutated
    fitsTemp = similar(fits) # allocate memory
    seqsTemp = deepcopy(seqs) # allocate memory
    binding_data = Array{Float64, 3}(undef, P, N_env, N)

    for t in 1:N # generations
        # compute Fitness
        updateAllFitnesses!(fits, seqs, Q, K, ligs, mutatedRecord)

        # save binding data
        binding_data[:,:,t] = fits

        # apply selection
        applySelection!(seqs, fits, α, k, fitsTemp, seqsTemp)

        # mutation.
        mutatedRecord = mutateAtRate!.(seqs, [μ], [Q]) # mutation all sequences every generation.

        # switch environments
        t % τ == 0 && (k = k % N_env + 1) # increment the env cyclically
    end
    return binding_data
end

function evolvePop_saveBinding(K::Table, Q::Settings, ligs::Ligands,
                   evo::EvoParams; k_init=1)
    seqs = map(x -> randSeq(Q), 1:evo.P)
    return evolvePop_saveBinding(seqs, K, Q, ligs, evo; k_init)
end


"""
Like `evolvePop` but records the *entire* population every generation over the last
`N_samplePeriods` environmental periods (one period is `N_env*τ` generations), so that
the within-period dynamics can be resolved.

Returns `(popTraj, k_record)`:
- `popTraj::PopArray`  of size `(seqLen, P, N_record)` with `N_record = N_samplePeriods*N_env*τ`
- `k_record::Vector{Int}` of length `N_record`, the environment index.

Recording happens at the end of the generation loop (after selection, mutation and the
environment switch), which is where `evolvePop` samples its sequences: `popTraj[:,:,i]` is
the population entering a generation and `k_record[i]` is the environment it will be
selected in.

The window is aligned to the environmental cycle: it holds exactly `N_samplePeriods`
complete periods and `k_record` starts at environment 1, so `foldPeriods` can average
over all of them. The last few generations of the run (fewer than one period) are
simulated but not recorded.
"""
function evolvePop_saveSeqs(seqs::Vector{Sequence},
                            K::Table,
                            Q::Settings,
                            ligs::Ligands,
                            evo::EvoParams;
                            k_init::Int=1, # Initial enviroment
                            N_samplePeriods::Int=4)

    N, μ, τ, α = evo.N, evo.μ, evo.τ, evo.α
    N_env = length(ligs) - 1
    N_record = N_samplePeriods * N_env * τ
    @assert iseven(N)
    @assert 0 <= μ <= 1 # constrain mutation rate to be a probability
    @assert 4τ <= N # Env periods can't be longer that 1/4 the total evo traj.
    @assert 0 < k_init <= N_env # Env index has to be positive
    @assert N_record <= N # can't record more generations than the trajectory has

    Q = deepcopy(Q) # detach from outside function
    P = length(seqs) # population
    k = k_init # environment index
    fits = zeros(length(seqs), N_env)
    mutatedRecord = BitVector(undef, P); mutatedRecord .= 1 # keep track of who was mutated
    fitsTemp = similar(fits) # allocate memory
    seqsTemp = deepcopy(seqs) # allocate memory

    # Latest generation at which recording can start and still fit N_record
    # generations in, walked back to the start of an environmental cycle: the
    # environment switches when t % τ == 0, and after t÷τ switches it is
    # ((k_init-1 + t÷τ) % N_env) + 1, so recording from a t with t % τ == 0 and
    # that expression equal to 1 makes the window exactly N_samplePeriods
    # complete periods beginning in environment 1.
    t_start = τ * ((N - N_record + 1) ÷ τ)
    while (k_init - 1 + t_start ÷ τ) % N_env != 0
        t_start -= τ
    end
    @assert t_start >= 1 "trajectory too short to record $N_samplePeriods aligned periods; increase N"
    t_stop = t_start + N_record - 1 # last generation that gets recorded
    popTraj = PopArray(undef, length(seqs[1]), P, N_record)
    k_record = Vector{Int}(undef, N_record)

    for t in 1:N # generations
        # compute Fitness
        updateAllFitnesses!(fits, seqs, Q, K, ligs, mutatedRecord)

        # apply selection
        applySelection!(seqs, fits, α, k, fitsTemp, seqsTemp)

        # mutation.
        mutatedRecord = mutateAtRate!.(seqs, [μ], [Q]) # mutation all sequences every generation.

        # switch environments
        t % τ == 0 && (k = k % N_env + 1) # increment the env cyclically

        # record the whole population
        if t_start <= t <= t_stop
            i = t - t_start + 1
            for j in 1:P
                @views popTraj[:,j,i] .= seqs[j]
            end
            k_record[i] = k
        end
    end
    return popTraj, k_record
end

function evolvePop_saveSeqs(K::Table, Q::Settings, ligs::Ligands,
                            evo::EvoParams; k_init=1, N_samplePeriods::Int=4)
    seqs = map(x -> randSeq(Q), 1:evo.P)
    return evolvePop_saveSeqs(seqs, K, Q, ligs, evo; k_init, N_samplePeriods)
end

"""
Seeded variant of `evolvePop_saveSeqs` for a fixed table `K`.
Seeds the global RNG with `seed` before evolving, so the recorded trajectory is exactly
reproducible (see the seeded `evolvePop`).
"""
function evolvePop_saveSeqs(K::Table, Q::Settings, ligs::Ligands,
                            evo::EvoParams, seed::Integer;
                            k_init=1, N_samplePeriods::Int=4)
    Random.seed!(seed)
    return evolvePop_saveSeqs(K, Q, ligs, evo; k_init, N_samplePeriods)
end


function updateAllFitnesses!(fits::Matrix, seqs::Vector{Sequence},
                             Q::Settings, K::Table, ligs::Ligands)
    mr = BitVector(undef, length(seqs)); mr .= 1
    return updateAllFitnesses!(fits, seqs, Q, K, ligs, mr)
end

"""
Compute fitness for all sequences for all enviroments.
omit those sequences that have not been mutated, (0 in mutatedRecord)
"""
function updateAllFitnesses!(fits::Matrix, seqs::Vector{Sequence}, Q::Settings,
                             K::Table, ligs::Ligands, mutatedRecord::BitVector)
    for i in eachindex(seqs)
        if mutatedRecord[i] # Bool value
            energies = computeFreeEnergies(seqs[i], K, Q, ligs)
            fits[i,:] .=  energies[1] .- energies[2:end]
        end
    end
    return nothing
end




"""
define exponential binding energy to fitness map.
α defines selection, larger α, stronger selection.
Subtract maxFits from fits to make the largest weight
equal to one to correct floating point error. wont effect selection.
"""
function fits2weights(fits::AbstractVector, α::Number)
    maxFits = maximum(fits)
    w = exp.( α .* (fits .- maxFits) )
    return w
end



"""
apply probablity selection.
1. fitnesses in env k pass through fits2weights to get a weight
2. sequences are drawn for survival according to weight
"""
function applySelection!(seqs::Vector{Sequence},
                         fits::Matrix,
                         α::Number,
                         k::Int,
                         fitsTemp::Matrix,
                         seqsTemp::Vector{Sequence})
    @views f = fits[:,k]
    P = length(f)
    w = fits2weights( f, α)
    winners = wsample(1:P, w, P)
    
    # make backup of fits and seqs
    fitsTemp .= fits
    for i in 1:P
        seqsTemp[i] .= seqs[i]
    end

    # reproduce winners
    for i in 1:P
        @views fits[i,:] = fitsTemp[winners[i],:]
        @views seqs[i] .= seqsTemp[winners[i]]
    end
    return nothing
end


"""
Evolve a population of sequences under a static (non-fluctuating) selection.
Unlike `evolvePop`, there is no environment index and no `τ` cycling: every
generation each sequence is scored by a single scalar fitness, `computeFitness`
dispatched on the `assay` singleton (e.g. `DoubleBinding()` selects on the
generalist score `min(-ΔF1, -ΔF2)`). Saves N_seqs2Keep sequences from the last
half of the run.
"""
function evolvePopStatic(seqs::Vector{Sequence},
                         K::Table,
                         Q::Settings,
                         ligs::Ligands,
                         evo::EvoParams,
                         assay::Assay;
                         N_seqs2Keep::Int=Int(evo.N/2) )

    N, μ, α = evo.N, evo.μ, evo.α
    @assert iseven(N)
    @assert 0 <= μ <= 1 # constrain mutation rate to be a probability
    @assert 2N_seqs2Keep < N # cant keep more seqs than half the number of generations

    Q = deepcopy(Q) # detach from outside function
    P = length(seqs) # population
    fits = zeros(P)
    mutatedRecord = BitVector(undef, P); mutatedRecord .= 1 # keep track of who was mutated
    fitsTemp = similar(fits) # allocate memory
    seqsTemp = deepcopy(seqs) # allocate memory
    seqBucket = Sequence[] # list to save sequences

    # Get list of generations to sample a sequence from (last half of the run).
    gen2StartSampling = Int(N/2) + 1
    gens2Sample = sort((gen2StartSampling:N)[randperm(length(gen2StartSampling:N))[1:N_seqs2Keep]])

    for t in 1:N # generations
        # compute Fitness
        updateAllFitnessesStatic!(fits, seqs, Q, K, ligs, assay, mutatedRecord)

        # apply selection
        applySelectionStatic!(seqs, fits, α, fitsTemp, seqsTemp)

        # mutation.
        mutatedRecord = mutateAtRate!.(seqs, [μ], [Q]) # mutation all sequences every generation.

        # Sample sequence
        t in gens2Sample && push!(seqBucket, deepcopy(seqs[rand(1:P)]))
    end
    return seqBucket
end

function evolvePopStatic(K::Table, Q::Settings, ligs::Ligands,
                         evo::EvoParams, assay::Assay; N_seqs2Keep::Int=evo.P)
    seqs = map(x -> randSeq(Q), 1:evo.P)
    return evolvePopStatic(seqs, K, Q, ligs, evo, assay; N_seqs2Keep)
end

"""
Seeded variant of `evolvePopStatic` for a fixed table `K`.
Seeds the global RNG with `seed` before evolving, so the entire trajectory
(initial population, mutations, selection, and which sequences are sampled)
is exactly reproducible.
"""
function evolvePopStatic(K::Table, Q::Settings, ligs::Ligands,
                         evo::EvoParams, assay::Assay, seed::Integer;
                         N_seqs2Keep::Int=evo.P)
    Random.seed!(seed)
    return evolvePopStatic(K, Q, ligs, evo, assay; N_seqs2Keep)
end


function updateAllFitnessesStatic!(fits::Vector, seqs::Vector{Sequence},
                                   Q::Settings, K::Table, ligs::Ligands, assay::Assay)
    mr = BitVector(undef, length(seqs)); mr .= 1
    return updateAllFitnessesStatic!(fits, seqs, Q, K, ligs, assay, mr)
end

"""
Compute the scalar fitness (`computeFitness` on `assay`) for all sequences.
Omit those sequences that have not been mutated (0 in mutatedRecord).
"""
function updateAllFitnessesStatic!(fits::Vector, seqs::Vector{Sequence}, Q::Settings,
                                   K::Table, ligs::Ligands, assay::Assay,
                                   mutatedRecord::BitVector)
    for i in eachindex(seqs)
        if mutatedRecord[i] # Bool value
            fits[i] = computeFitness(seqs[i], K, Q, ligs, assay)
        end
    end
    return nothing
end

"""
apply probability selection on a scalar fitness vector (static, no environment).
1. fitnesses pass through fits2weights to get a weight
2. sequences are drawn for survival according to weight
"""
function applySelectionStatic!(seqs::Vector{Sequence},
                               fits::Vector,
                               α::Number,
                               fitsTemp::Vector,
                               seqsTemp::Vector{Sequence})
    P = length(fits)
    w = fits2weights(fits, α)
    winners = wsample(1:P, w, P)

    # make backup of fits and seqs
    fitsTemp .= fits
    for i in 1:P
        seqsTemp[i] .= seqs[i]
    end

    # reproduce winners
    for i in 1:P
        fits[i] = fitsTemp[winners[i]]
        @views seqs[i] .= seqsTemp[winners[i]]
    end
    return nothing
end


"""
Returns the number of generations for the evolutionary trajectory.
Must have atleast totalMuts mutations per sequence.
and at least minNumEpoch's per simulation.
"""
function pickNumGens(μ::Number, τ::Integer, seqLen::Integer,
                     minNumMuts::Integer, minNumEpochs::Integer, 
                     minN::Integer)
    N_fromMuts = 2*Int(ceil( minNumMuts / (seqLen* μ*2)))
    N_fromEpochs =2*Int(round(τ * minNumEpochs/2))
    N_fromMinN = 2*Int(round(minN/2))
    return maximum([N_fromMuts, N_fromEpochs, N_fromMinN])
end

