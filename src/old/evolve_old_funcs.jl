#function evolve(seq::Sequence, K::Table, Q::Settings, ligs::Ligands, N::Int; ΔT=1)
#    # Metropolis Monte Carlo in the space of couplings
#    # starting from sequence seq
#    # Q : settings object with all details of system.   
#    seq = copy(seq)
#    q, W, L, ζ = Q.q, Q.W, Q.L, Q.ζ
#    Lseq = W*(L+1) 
#    # compute initial fitness
#    ϕ = computeFitness(seq, K, Q, ligs)
#    times, fitnesses, sequences = Int[], Float64[], Sequence[]
#    t = 0    
#    ΔT == 0 && (ΔT = N-1)   
#
#    while t < N
#        # mutation:
#        imut, amut = rand(1:Lseq), rand(1:(q-1))
#        awt = seq[imut]
#        seq[imut] = ((awt-1 + amut) % q) +1
#        ϕ_mut = computeFitness(seq, K, Q, ligs)
#        # selection:
#        if rand() < exp(ζ*(ϕ_mut-ϕ))
#            ϕ = ϕ_mut         
#        else
#            seq[imut] = awt
#        end
#        # recording:
#        if mod(t, ΔT) == 0
#            push!(times, t)
#            push!(fitnesses, ϕ)
#            push!(sequences, deepcopy(seq))
#        end 
#        t += 1
#    end
#    return times, fitnesses, sequences
#end
#
#function evolve(K::Table, Q::Settings, ligs::Ligands, N::Int; seed=rand(UInt), ΔT=1)
#    seq = randSeq(Q) # start with random sequence.
#    return evolve(seq, K, Q, ligs, N; ΔT)
#end


#function evolve2(seq::Sequence, K::Table, Q::Settings, ligs::Ligands, N::Int; ϕ_goal=Inf)
#    # Metropolis Monte Carlo in the space of couplings
#    # starting from sequence seq and returns final seq and fitness
#    # Does not keep intermediate sequences or fitnesses.
#    # Q : settings object with all details of system.   
#    # This version terminates if ϕ>=ϕ_goal or t>=N.
#    q, W, L, ζ = Q.q, Q.W, Q.L, Q.ζ
#    Lseq = W*(L+1) 
#    # compute initial fitness
#    ϕ = computeFitness(seq, K, Q, ligs)
#    times, fitnesses, sequences = Int[], Float64[], Sequence[]
#    t = 0    
#    finished(ϕ, ϕ_goal, t, N) = (ϕ>=ϕ_goal || t>=N) ? (return true) : (return false)
#
#    while !finished(ϕ, ϕ_goal, t, N)
#        # mutation:
#        imut, amut = rand(1:Lseq), rand(1:(q-1))
#        awt = seq[imut]
#        seq[imut] = ((awt-1 + amut) % q) +1
#        ϕ_mut= computeFitness(seq, K, Q, ligs)
#        # selection:
#        if rand() < exp(ζ*(ϕ_mut-ϕ))
#            ϕ = ϕ_mut         
#        else
#            seq[imut] = awt
#        end
#        t+=1
#    end
#    return seq, ϕ
#end
#
#function evolve2(K::Table, Q::Settings, ligs::Ligands, N::Int, seed::UInt; ϕ_goal=Inf)
#    # set seed, gen rand sequence, evolve.
#    Random.seed!(seed)
#    seq = randSeq(Q) # start from random sequence
#    return evolve2(seq, K, Q, ligs, N; ϕ_goal)
#end
#function evolve2(seq::Sequence, K::Table, Q::Settings, ligs::Ligands, N::Int, seed::UInt; ϕ_goal=Inf)
#    # set seed, evolve from starting seq.
#    Random.seed!(seed)
#    return evolve2(seq, K, Q, ligs, N; ϕ_goal)
#end



#fits2weights(x, a, x0) = 1/(1 + exp(-a*(x-x0)))
#"""
#define sigmoid binding to fitness map.
#"""
#function fits2weights(fits::AbstractVector, a::Number, x0::Number)
#    mx  = maximum( x -> a*(x-x0), fits)
#    if mx < -100
#        w = exp.( a .* ( fits .- x0 ) .- mx)
#    else
#        w = 1 ./ ( 1 .+  exp.( -a .* (fits .- x0)) )
#    end
#    return w
#end
