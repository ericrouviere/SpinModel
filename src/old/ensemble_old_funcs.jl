# Commented-out code moved from src/ensemble.jl.

#function evolveEnsemble(tables::Vector{Table}, Q::Settings)
#    N_reps = length(tables)
#    output = pmap(x->evolve(x, Q, ΔT=0), tables)
#    sequences = Array{Sequence,2}(undef, N_reps, 2)
#    fitnesses = Array{Float64, 2}(undef, N_reps, 2)
#    for i in eachindex(output)
#        tims, fits, seqs = output[i]
#        sequences[i,:] = seqs
#        fitnesses[i,:] = fits
#    end
#    return sequences, fitnesses
#end
#
#function evolve2Ensemble(K::Table, Q::Settings, N::Int, N_batch::Int,
#                         N_rounds::Int, seeds::Vector{UInt}; ϕ_goal=Int,
#                         seq0=0)
#    # evolve many sequences, starting from random sequences
#    # all with the same interaction table for generating big MSAs.
#    # N = Max number of MCMC iteractions during evolve.
#    # N_batch = number of sequences evolved with single pmap call.
#    # N_rounds =  number of pmap calls.
#    # Use a sinlge starting sequence through seq0.
#
#    for i in 1:N_rounds
#        @printf("%s%03d\n","running batch ", i)
#        # evolve one batch of sequences
#        seeds_batch = seeds[ (i-1)*N_batch+1 : i*N_batch]
#        if typeof(seq0) == Sequence # start from seq0
#            output = pmap(seed -> evolve2(seq0, K, Q, N, seed; ϕ_goal), seeds_batch)
#        elseif seq0 == 0 # start from rand seq
#            output = pmap(seed -> evolve2(K, Q, N, seed; ϕ_goal), seeds_batch)
#        else
#            error("You did not pass a proper sequence seq0.")
#        end
#        # reformat output 
#        seq1, fit1 = output[1]
#        seqLen = length(seq1)
#        sequences = Matrix{eltype(seq1)}(undef, seqLen, N_batch)
#        fitnesses = Vector{typeof(fit1)}(undef, N_batch)
#
#        for i in eachindex(output)
#            seq, fit = output[i]
#            sequences[:,i] .= seq
#            fitnesses[i] = fit
#        end
#        label = @sprintf("%03d", i)
#        save("data_"*label*".jld2", "fitnesses", fitnesses, "sequences", sequences)
#    end
#end

#function computeEffectiveFieldEnsemble(sequences::Array{Sequence},
#                                       K::Table,
#                                       site,
#                                       Q::Settings)
#    return map( x -> computeEffectiveField(x, K, site, Q), sequences)
#end

#function analyseBindingDMSEnsemble(dms::AbstractArray{<:Number, 4},
#                                     x::Number;
#                                     W=5,
#                                     actSiteLayer=11)
#    numTypes, numNodes, numLigands, popSize = size(dms)
#    numSensMutsArray = Vector{Int64}(undef, popSize)
#    numSensPosArray = similar(numSensMutsArray)
#    maxDistArray = similar(numSensMutsArray)
#    for i in 1:popSize
#        @views dm = dms[:,:,:,i]
#        numSensMuts, numSensPos, maxDist = countMutSensitivity(dm, x; W, actSiteLayer)
#        numSensMutsArray[i] = numSensMuts
#        numSensPosArray[i] = numSensPos
#        maxDistArray[i] = maxDist
#    end
#    return numSensMutsArray, numSensPosArray, maxDistArray
#end

#function computeEvolvabilityEnsemble(dms::AbstractArray{T, 5}) where T <: Number
#    numTypes, numNodes, numLigands, popSize, numReps = size(dms)
#    evolvability = Matrix{Float64}(undef, popSize, numReps)
#    for i in 1:popSize, j in 1:numReps
#        @views dm = dms[:,:,:,i,j]
#        evolvability[i,j] = computeEvolvability(dm)
#    end
#    return evolvability
#end

# From inside computeOverlapsEnsemble (src/ensemble.jl), which now returns only (OC, vC):
#    vC = similar(OC) #, OJ, vJ = similar(OC), similar(OC), similar(OC)
        #OJ[:,i] .= out[i][3] 
        #vJ[:,i] .= out[i][4]
#    return OC, vC#, OJ, vJ 
