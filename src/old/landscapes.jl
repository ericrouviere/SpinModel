# OLD CODE ###

#function computeMagnetization(seq::Sequence, K::Table, sites::Sites, fields::Fields, Q::Settings)
#    # return the magnetization at each site for each ligand perturbation.
#    W, L = Q.W, Q.L
#    J, β = seq2J(seq, K), Q.β
#    conf_left, conf_right = ones(2^W), ones(2^W)
#    Mmat = [zeros(W,L+1) for i in eachindex(sites)]
#    H = copy(J[1])
#    J_mat = Matrix{Float64}(undef, W, W) #
#    J_mat_states = Matrix{Float64}(undef, W, 2^W) #
#    h_states = Vector{Float64}(undef, 2^W) #
#    Mat = zeros(2^W,2^W)
#    states = getSpinStates(W)
#    
#    for j in eachindex(sites) # for each perturbation
#        site, field = sites[j], fields[j]
#        #h00 = J[1][site] # save old feild
#        J[1][site] += field # bind ligand by adding fields
#        for k = 1:(L+1) # for each layer
#            Z_right = conf_right
#            for layer in 1:(k-1)
#                MatLayer!(Mat, J_mat, J_mat_states, h_states, states, J, layer, β)
#                Z_right = Mat * Z_right
#            end
#            Z_left = conf_left'
#            for layer in (L+1):-1:k
#                MatLayer!(Mat, J_mat, J_mat_states, h_states, states, J, layer, β)
#                Z_left = Z_left * Mat
#            end
#            Z_den = Z_left*Z_right
#            for i = 1:W
#                seqq = zeros(W)
#                seqq[i] = 1
#                Z_num = Z_left * (2*diagm(seq2conf(seqq)) - I) * Z_right
#                Mmat[j][i,k] = Z_num / Z_den
#            end
#        end
#        J[1] .= H # unbind ligand
#    end
#    return Mmat
#end



#function computeMagnetization(seq::Sequence, K::Table, sites::Sites, fields::Fields, Q::Settings)
#    # return the magnetization at each site for each ligand perturbation.
#    W, L = Q.W, Q.L
#    J, β = seq2J(seq, K), Q.β
#    Z_right = Vector{Float64}(undef, 2^W)  
#    Z_left  = Vector{Float64}(undef, 2^W)
#    buff_vec = Vector{Float64}(undef, 2^W)
#    Mmat = [zeros(W,L+1) for i in eachindex(sites)]
#    J_mat = Matrix{Float64}(undef, W, W) #
#    J_mat_states = Matrix{Float64}(undef, W, 2^W) #
#    h_states = Vector{Float64}(undef, 2^W) #
#    T = [Matrix{Float64}(undef, 2^W, 2^W) for i in 1:(L+1)]
#    states = getSpinStates(W)
#
#    function mul_update!(B::AbstractArray, buff::AbstractArray, A::AbstractArray)
#        # does the following opperation inplace: B <-- A*B
#        mul!(buff, A, B)
#        B .= buff
#        return nothing
#    end
#
#    function computeMags!(j::Int)
#        # compute all transfer matrices
#        for k in 1:L+1 # each layer k
#            MatLayer!(T[k], J_mat, J_mat_states, h_states, states, J, k, β)
#        end
#
#        # compute magnetizations
#        for l = 1:(L+1) # for each layer
#            Z_right .= 1.0
#            for k in 1:(l-1)
#                mul_update!(Z_right, buff_vec, T[k])
#            end
#            Z_left .= 1.0
#            for k in (L+1):-1:l
#                mul_update!(Z_left, buff_vec, T[k]')
#            end
#            Z = Z_left' * Z_right
#            for i = 1:W
#                @views buff_vec .= Z_left .* states[i,:] .* Z_right 
#                Mmat[j][i,l] = sum(buff_vec) / Z
#            end
#        end
#        return nothing
#    end
#
#    for j in eachindex(sites) # for each perturbation
#        site, field = sites[j], fields[j]
#        h0 = J[1][site] # save field at site
#        J[1][site] += field # bind ligand by adding fields
#        computeMags!(j)
#        J[1][site] = h0 # put back the original field at site
#    end
#    return Mmat
#end





#function computeFreeEnergies(seq::Sequence,
#                         K::Table,
#                         Q::Settings,
#                         sites::Sites,
#                         fields::Fields;
#                         site_add=CartesianIndex(1,1),
#                         h_add=0)
#    # calculates the Free energies of the model 
#    # for each binding condition.
#    @assert length(sites) == length(fields)
#    W, L, β = Q.W, Q.L, Q.β
#    energies = zeros(length(sites)) 
#    J = seq2J(seq, K)
#    J[1][site_add] += h_add # add extra field to site site_add
#    H = copy(J[1]) # copy the fields
#    for i in eachindex(sites) 
#        site = sites[i]
#        field = fields[i]
#        J[1][site] += field  # binding ligand at site
#        energies[i] =  computeFreeEnergy(J, W, β)
#        J[1] .= H # undo mutation
#    end
#    J[1][site_add] -= h_add
#    return energies
#end


#function computeFreeEnergy(J::Vector{Matrix{Float64}}, W::Int, β::Number)
#    return computeFreeEnergy(J, ones(2^W), ones(2^W), β)
#end

#function MatLayer!(Mat::Matrix{Float64}, J_mat, J_mat_states, h_states,
#                    states::Matrix{Float64}, J::Vector{Matrix{Float64}}, k::Int, β::Float64)
#    # Compute the transfer matrix at a given layer k
#    # note: the botzmann factor purposefully omits the normal minus sign, exp(β ...).
#    W, L = size(J[2])
#    h_vec = β * J[1][:,k] ###############################################
#    #@show sum(h_vec)
#    mul!(h_states, states', h_vec)
#
#    if k < (L+1) # generic case
#        Coupling_mat!(J_mat, J, k, β, W) ###################################
#        mul!(J_mat_states, J_mat, states)
#        mul!(Mat, states', J_mat_states)
#        @inbounds for i0 in 1:2^W, i1 in 1:2^W
#            Mat[i1,i0] = exp( Mat[i1,i0] + h_states[i0] )
#        end
#    else # special case (only fields). ie k = L+1
#        @inbounds for i0 = 1:2^W, i1=1:2^W
#            i0 == i1 ? Mat[i0,i1] = exp(h_states[i0]) : Mat[i0,i1] = 0.0 
#        end
#    end
#    return nothing
#end

#function Coupling_mat!(Cmat, J, k, β, W)
#    # Matrix of coupling for a given layer
#    #M = Matrix{Float64}(undef, W, W)
#    for j in 1:W, i in 1:W
#        if i==j
#            Cmat[i,j] = β*J[2][i,k]
#        elseif j-i==1
#            Cmat[i,j] = β*J[3][i,k]
#        else
#            Cmat[i,j] = 0.0
#        end
#    end
#    Cmat[W,1] = β*J[3][W,k]
#    return nothing
#end



#function computeFreeEnergies(seq::Sequence,
#                         K::Table,
#                         Q::Settings;
#                         site_add=CartesianIndex(1,1),
#                         h_add=0)
#    return computeFreeEnergies(seq, K, Q, Q.sites, Q.fields; site_add, h_add)
#end

#function computeFreeEnergy(J::Vector{Matrix{Float64}},
#                           conf_left::Vector{Float64},
#                           conf_right::Vector{Float64}, β::Number)
#    # Calculation of free energy by transfer matrices
#    W, L = size(J[2])
#    Mat = Matrix{Float64}(undef, 2^W,2^W)
#    J_mat = Matrix{Float64}(undef, W, W) #
#    J_mat_states = Matrix{Float64}(undef, W, 2^W) #
#    h_states = Vector{Float64}(undef, 2^W) #
#
#    states = getSpinStates(W)
#    Z = conf_right
#    for layer in 1:L+1
#        MatLayer!(Mat, J_mat, J_mat_states, h_states, states, J, layer, β)
#        Z = Mat * Z
#    end
#    #return -log(conf_left'*Z)[1,1]/β
#    return -log(conf_left'*Z)/β
#end

#function computeFreeEnergy(J::Vector{Matrix{Float64}}, W::Int, β::Number)
#    # Calculation of free energy by transfer matrices
#    W, L = size(J[2])
#    Mat = Matrix{Float64}(undef, 2^W,2^W)
#    J_mat = Matrix{Float64}(undef, W, W) #
#    J_mat_states = Matrix{Float64}(undef, W, 2^W) #
#    h_states = Vector{Float64}(undef, 2^W) #
#    Z = ones(2^W)
#    tmp = similar(Z)
#    states = getSpinStates(W)
#    for layer in 1:L+1
#        MatLayer!(Mat, J_mat, J_mat_states, h_states, states, J, layer, β)
#        mul!(tmp, Mat, Z)
#        Z .= tmp
#    end
#    return -log(sum(Z))/β # taking sum instead of trace due to boundry conditions.
#end





#function computeFreeEnergy_old(J::Vector{Matrix{Float64}},
#                           conf_left::Vector{Float64},
#                           conf_right::Vector{Float64}, β::Number)
#    # Calculation of free energy by transfer matrices
#    Z = conf_right
#    for layer in 1:(size(J[3],2)+1)
#        Z = MatLayer(J, layer, β) * Z
#    end
#    return -log(conf_left'*Z)[1,1]/β
#end
#
#
#function MatLayer(J::Vector{Matrix{Float64}}, k::Int, β::Number)
#    # Computing the transfer matrix at a given layer k
#    # the last layer, k=L+1, is a special case (only fields)
#    W, L = size(J[2])
#    Mat = zeros(2^W,2^W)
#    h_vec = β * J[1][:,k]
#    if k < (L+1) # generic case
#        J_mat = β * Coupling_mat(J[2][:,k], J[3][:,k])
#         for i0 = 0:2^W-1
#            sigma0 = 2 .* (.5 .- digits(i0, base=2, pad = W))
#            h_sigma0 = dot(h_vec, sigma0)
#            J_mat_sigma0 = J_mat*sigma0
#            for i1 = 0:2^W-1
#                sigma1 = 2 .* ( .5 .- digits(i1, base = 2, pad = W))
#                Mat[i1+1,i0+1] = exp(dot(sigma1, J_mat_sigma0) + h_sigma0) # indices shifted by 1
#            end
#        end
#    else # special case (only fields)
#        for i0 = 0:2^W-1
#            sigma0 = 2 .* (.5 .- digits(i0, base=2, pad = W))
#            h_sigma0 = dot(h_vec, sigma0)
#            Mat[i0+1,i0+1] = exp(h_sigma0)
#        end
#    end
#    return Mat
#end



#function MatLayer!(Mat::Matrix{Float64}, J::Vector{Matrix{Float64}}, k::Int, β::Number)
#    # Computing the transfer matrix at a given layer k
#    # the last layer, k=L+1, is a special case (only fields)
#
#    W, L = size(J[2])
#    dig = zeros(Int, W)
#    sigma0 = zeros(W)
#    sigma1 = zeros(W)
#    h_vec = β * J[1][:,k]
#
#    if k < (L+1) # generic case
#        J_mat = β * Coupling_mat(J[2][:,k], J[3][:,k])
#        for i0 = 0:2^W-1
#            i2σ!(sigma0, dig, i0) 
#            h_sigma0 = fastdot(h_vec, sigma0)
#            J_mat_sigma0 = J_mat*sigma0
#            for i1 = 0:2^W-1
#                i2σ!(sigma1, dig, i1) 
#                Mat[i1+1,i0+1] = exp(fastdot(sigma1, J_mat_sigma0) + h_sigma0) # indices shifted by 1
#            end
#        end
#    else # special case (only fields). ie k = L+1
#        for i0 = 0:2^W-1
#            Mat[:,i0+1] .= 0
#            i2σ!(sigma0, dig, i0)
#            Mat[i0+1,i0+1] = exp(fastdot(h_vec, sigma0))
#        end
#    end
#    return nothing
#end
#
#
#
#function i2σ!(σ::Vector{Float64}, dig::Vector{Int}, i::Int)
#    digits!(dig, i, base=2)
#    σ .= 2 .* (.5 .- dig)
#    return nothing
#end

#function computeFreeEnergy(J::Vector{Matrix{Float64}},
#                           conf_left::Vector{Float64},
#                           conf_right::Vector{Float64}, β::Number)
#    # Calculation of free energy by transfer matrices
#    W, L = size(J[2])
#    Mat = zeros(2^W,2^W)
#    Z = conf_right
#    for layer in 1:L+1
#        MatLayer!(Mat, J, layer, β)
#        Z = Mat * Z
#    end
#    return -log(conf_left'*Z)[1,1]/β
#end

#function fastdot(a::AbstractVector, b::AbstractVector)
#    # A faster dot product.
#    s = zero(eltype(a))
#    @inbounds for i in eachindex(a)
#        s += a[i] * b[i]
#    end
#    return s
#end


#function Coupling_mat(Jleft, Jright)
#    # Matrix of coupling for a given layer
#    W = size(Jleft)[1]
#    M = Array(Bidiagonal(Jleft, Jright[1:(W-1)], :U)) # upper bidiagonal
#    M[W,1] = Jright[W]
#    return M
#end


# OLD CODE ###

#function computeFreeEnergy_old(J::Vector{Matrix{Float64}},
#                           conf_left::Vector{Float64},
#                           conf_right::Vector{Float64}, β::Number)
#    # Calculation of free energy by transfer matrices
#    Z = conf_right
#    for layer in 1:(size(J[3],2)+1)
#        Z = MatLayer(J, layer, β) * Z
#    end
#    return -log(conf_left'*Z)[1,1]/β
#end
#
#
#function MatLayer(J::Vector{Matrix{Float64}}, k::Int, β::Number)
#    # Computing the transfer matrix at a given layer k
#    # the last layer, k=L+1, is a special case (only fields)
#    W, L = size(J[2])
#    Mat = zeros(2^W,2^W)
#    h_vec = β * J[1][:,k]
#    if k < (L+1) # generic case
#        J_mat = β * Coupling_mat(J[2][:,k], J[3][:,k])
#         for i0 = 0:2^W-1
#            sigma0 = 2 .* (.5 .- digits(i0, base=2, pad = W))
#            h_sigma0 = dot(h_vec, sigma0)
#            J_mat_sigma0 = J_mat*sigma0
#            for i1 = 0:2^W-1
#                sigma1 = 2 .* ( .5 .- digits(i1, base = 2, pad = W))
#                Mat[i1+1,i0+1] = exp(dot(sigma1, J_mat_sigma0) + h_sigma0) # indices shifted by 1
#            end
#        end
#    else # special case (only fields)
#        for i0 = 0:2^W-1
#            sigma0 = 2 .* (.5 .- digits(i0, base=2, pad = W))
#            h_sigma0 = dot(h_vec, sigma0)
#            Mat[i0+1,i0+1] = exp(h_sigma0)
#        end
#    end
#    return Mat
#end



#function MatLayer!(Mat::Matrix{Float64}, J::Vector{Matrix{Float64}}, k::Int, β::Number)
#    # Computing the transfer matrix at a given layer k
#    # the last layer, k=L+1, is a special case (only fields)
#
#    W, L = size(J[2])
#    dig = zeros(Int, W)
#    sigma0 = zeros(W)
#    sigma1 = zeros(W)
#    h_vec = β * J[1][:,k]
#
#    if k < (L+1) # generic case
#        J_mat = β * Coupling_mat(J[2][:,k], J[3][:,k])
#        for i0 = 0:2^W-1
#            i2σ!(sigma0, dig, i0) 
#            h_sigma0 = fastdot(h_vec, sigma0)
#            J_mat_sigma0 = J_mat*sigma0
#            for i1 = 0:2^W-1
#                i2σ!(sigma1, dig, i1) 
#                Mat[i1+1,i0+1] = exp(fastdot(sigma1, J_mat_sigma0) + h_sigma0) # indices shifted by 1
#            end
#        end
#    else # special case (only fields). ie k = L+1
#        for i0 = 0:2^W-1
#            Mat[:,i0+1] .= 0
#            i2σ!(sigma0, dig, i0)
#            Mat[i0+1,i0+1] = exp(fastdot(h_vec, sigma0))
#        end
#    end
#    return nothing
#end
#
#
#
#function i2σ!(σ::Vector{Float64}, dig::Vector{Int}, i::Int)
#    digits!(dig, i, base=2)
#    σ .= 2 .* (.5 .- dig)
#    return nothing
#end

#function computeFreeEnergy(J::Vector{Matrix{Float64}},
#                           conf_left::Vector{Float64},
#                           conf_right::Vector{Float64}, β::Number)
#    # Calculation of free energy by transfer matrices
#    W, L = size(J[2])
#    Mat = zeros(2^W,2^W)
#    Z = conf_right
#    for layer in 1:L+1
#        MatLayer!(Mat, J, layer, β)
#        Z = Mat * Z
#    end
#    return -log(conf_left'*Z)[1,1]/β
#end

#function fastdot(a::AbstractVector, b::AbstractVector)
#    # A faster dot product.
#    s = zero(eltype(a))
#    @inbounds for i in eachindex(a)
#        s += a[i] * b[i]
#    end
#    return s
#end


#function Coupling_mat(Jleft, Jright)
#    # Matrix of coupling for a given layer
#    W = size(Jleft)[1]
#    M = Array(Bidiagonal(Jleft, Jright[1:(W-1)], :U)) # upper bidiagonal
#    M[W,1] = Jright[W]
#    return M
#end


#function evolvePop(seqs::Vector{Sequence},
#                   K::Table,
#                   Q::Settings,
#                   fieldsList::Vector{Fields},
#                   assayList::Vector{String},
#                   P::Int,
#                   a::Number, # selection threshold eg a=0.3 means remove bottom 30%
#                   N::Int, # number of generations
#                   μ::Number, # mutation rate
#                   τ::Int; # enviromental timescale
#                   sampleLastHalf::Bool=false,
#                   numSeqs2Keep::Int=length(seqs)) 
#
#    # Evolve a popluation of sequences under a fluctuating selection
#    @assert 0 <= μ <= 1
#    @assert 0 <= a <= 1
#    Q = deepcopy(Q) # detach from outside function
#    P = length(seqs)
#    τ ==0 && (τ=N)
#    k = 1 # environment index
#    numWinners = Int(floor((1-a)*P) + ((1-a)<1) * 1)
#    
#    # compute initial fitnesses
#    fits = map(x -> computeFitness(x, K, Q), seqs)
#    randIndicies = sortperm(fits, rev=true) # intialize vector randIndicies
#    seqs0 = deepcopy(seqs)
#    previousSelection = "single"
#    
#    # some book keeping for sampling sequence from the population
#    numGen2Sample = Int(floor(N / 2))
#    maxSeqsSampleSize = P * numGen2Sample
#    if numSeqs2Keep > maxSeqsSampleSize
#        println("You are trying to keep too many sequences\n"*
#                "Setting numSeqs2Keep = maxSeqsSampleSize")
#        numSeqs2Keep = maxSeqsSampleSize
#    end
#    numSeqs2KeepPerGen = Int(ceil(numSeqs2Keep / numGen2Sample))
#    gen2StartSampling = Int(ceil(N/2))
#    seqBucket = Sequence[]
#
#    for t in 0:N-1
#
#        env = updateEnvironment(t, τ, env)
#        # update fitness
#        applySelection!(seqs, fits, rankedIndicies, P, numWinners)
#        mutatedRecord = mutateAtRate!.(seqs, [μ], [Q])
#
#        # compute fitness
#        if t % τ == 0 # if time to switch
#            tempAssay = "DoubleBinding"
#            tempFields = [fieldsList[1][1], fieldsList[1][2], fieldsList[2][2]] # should be [hs,hr,hw]
#            tempSites = [Q.sites..., Q.sites[1]]
#            if previousSelection == "single"
#                fits = map(x -> computeFitness(x, K, Q, tempSites, tempFields,tempAssay), seqs)
#            elseif previousSelection == "both"
#                fits[mutatedRecord] = map(x -> computeFitness(x, K, Q, tempSites,
#                                          tempFields, tempAssay), seqs[mutatedRecord])
#            end
#            k +=1
#            Q.fields = fieldsList[(k-1)%length(fieldsList)+1]
#            previousSelection = "both"
#        else
#            if previousSelection == "both"
#                fits = map(x -> computeFitness(x, K, Q), seqs)
#            elseif previousSelection == "single"
#                fits[mutatedRecord] = map(x -> computeFitness(x, K, Q), seqs[mutatedRecord])
#            end
#            previousSelection = "single"
#        end
#        applySelection!(seqs, fits, rankedIndicies, P, numWinners)
#        mutatedRecord = mutateAtRate!.(seqs, [μ], [Q])
#
#        # sampling sequences
#        if sampleLastHalf && t >= gen2StartSampling
#            selectedSeqs = seqs[randperm(P)[1:numSeqs2KeepPerGen]]
#            append!(seqBucket, deepcopy(selectedSeqs))
#        end
#    end
#
#    if sampleLastHalf
#        seqs1 = seqBucket[randperm(length(seqBucket))[1:numSeqs2Keep]]
#    else
#        seqs1 = deepcopy(seqs)
#    end
#
#    return seqs0, seqs1
#end
#function evolvePop(seqs::Vector{Sequence},
#                   K::Table,
#                   Q::Settings,
#                   fieldsList::Vector{Fields},
#                   P::Int,
#                   a::Number, # selection threshold eg a=0.3 means remove bottom 30%
#                   N::Int, # number of generations
#                   μ::Number, # mutation rate
#                   τ::Int; # enviromental timescale
#                   sampleLastHalf::Bool=false,
#                   numSeqs2Keep::Int=length(seqs)) 
#    # Evolve a popluation of sequences under a fluctuating selection
#    Q = deepcopy(Q)
#    P = length(seqs)
#    if τ==0; (τ=N); end
#    k = 1 # environment index
#
#    numWinners = Int(floor((1-a)*P) + ((1-a)<1) * 1)
#    
#    # compute initial fitnesses
#    fits = map(x -> computeFitness(x, K, Q), seqs)
#    seqs0 = deepcopy(seqs)
#    
#    # some book keeping for sampling sequence from the population
#    numGen2Sample = Int(floor(N / 2))
#    maxSeqsSampleSize = P * numGen2Sample
#    if numSeqs2Keep > maxSeqsSampleSize
#        println("You are trying to keep too many sequences\n"*
#                "Setting numSeqs2Keep = maxSeqsSampleSize")
#        numSeqs2Keep = maxSeqsSampleSize
#    end
#
#    numSeqs2KeepPerGen = Int(ceil(numSeqs2Keep / numGen2Sample))
#    gen2StartSampling = Int(ceil(N/2))
#    seqBucket = Sequence[]
#
#    for t in 0:N-1
#        
#        # selection
#        rankedIndicies = sortperm(fits, rev=true)
#        winners = rankedIndicies[1:numWinners]
#        losers = rankedIndicies[numWinners+1:end]
#        replicates = rand(winners, P-numWinners)
#       
#        seqs[losers] = copy.(seqs[replicates])
#        #fits[losers] = copy(fits[replicates])
#
#        # mutation
#        mutatedRecord = mutateAtRate!.(seqs, [μ], [Q])
#        #mutateAtRate1!.(seqs, [μ], [Q])
#        
#        # compute fitness
#        if t % τ == 0
#            tempAssay = "DoubleBinding"
#            # should be [hs,hr,hw]
#            tempFields = [fieldsList[1][1], fieldsList[1][2], fieldsList[2][2]] 
#            tempSites = [Q.sites..., Q.sites[1]]
#            fits = map(x -> computeFitness(x, K, Q, tempSites, tempFields,
#                                           tempAssay), seqs)
#            k +=1
#            Q.fields = fieldsList[(k-1)%length(fieldsList)+1]
#        else
#            fits = map(x -> computeFitness(x, K, Q), seqs)
#        end
#
#        # sampling sequences
#        if sampleLastHalf && t >= gen2StartSampling
#            selectedSeqs = seqs[randperm(P)[1:numSeqs2KeepPerGen]]
#            append!(seqBucket, deepcopy(selectedSeqs))
#        end
#    end
#
#    if sampleLastHalf
#        seqs1 = seqBucket[randperm(length(seqBucket))[1:numSeqs2Keep]]
#    else
#        seqs1 = deepcopy(seqs)
#    end
#    
#    return seqs0, seqs1
#end



#function computeEvolvedEvolvabilityBindSpace(bind, dms, Q; limits = [)
#    
#    # return the evolvabilities of all N_collect sequences for each bin 
#    # in binding space.
#    @assert length(Q.fields)==3 
#    # gen rand seqs and compute functions
#    randseqs = map( x -> randSeq(Q), 1:N_randseq)
#    energies = pmap( x -> computeEnergies(x, K, Q), randseqs)
#    e1 = getindex.(energies, 1)
#    e2 = getindex.(energies, 2)
#    e3 = getindex.(energies, 3)
#    bindEnergies = [e2 - e1 e3 - e1]
#    
#    # bin binding space and get indices of seqs in each bin.
#    binSeqInds = binBindingSpace(bindEnergies, N_bins; N_collect);
#    temp(x) = computeEvolvability.(randseqs[binSeqInds[x]], [K], [Q])
#    evo =  pmap( temp, 1:N_bins);
#    return evo
#end


#################################################
### Old Code ####################################
#################################################

#function computeMagnetization(seq::Sequence, K::Table, sites::Sites, fields::Fields, Q::Settings)
#    # magnetization at each site
#    W, L = Q.W, Q.L 
#    J, β = seq2J(seq, K), Q.β 
#    conf_left, conf_right = ones(2^W), ones(2^W)
#    Mmat = [zeros(W,L+1) for i in eachindex(sites)]
#    H = copy(J[1])
#    for j in eachindex(sites)
#        site, field = sites[j], fields[j]
#        #h00 = J[1][site] # save old feild
#        J[1][site] += field # bind ligand by adding fields
#        for k = 1:(L+1) # compute magnetization
#            Z_right = conf_right
#            for layer in 1:(k-1)
#                Z_right = MatLayer(J, layer, β) * Z_right
#            end
#            Z_left = conf_left'
#            for layer in (L+1):-1:k
#                Z_left = Z_left * MatLayer(J, layer, β)
#            end 
#            Z_den = Z_left*Z_right
#            for i = 1:W 
#                seq = zeros(W)
#                seq[i] = 1 
#                #Z_num = Z_left * (2*diagm(seq2conf(seq)) - eye(2^W)) * Z_right
#                Z_num = Z_left * (2*diagm(seq2conf(seq)) - I) * Z_right
#                Mmat[j][i,k] = Z_num / Z_den
#            end 
#        end 
#        J[1] .= H # unbind ligand
#    end 
#    return Mmat
#end
#
#
#

#function saturated_mutagenesis(seq::Sequence,
#                               K::Table,
#                               Q::Settings; option="mean")
#    # fitness differences for all single point mutations
#    # option = mean or worst to have mean or worst effect over aa at a position
#    # option = all to have all mutational effects
#    f_wt = fitness(seq, K, Q) 
#    delta_f, delta_f_all = zeros(length(seq)), zeros(length(seq)*(Q.q-1))
#    k_all = 0
#    for i = 1:length(seq)
#        awt = seq[i]
#        df, k = zeros(Q.q-1), 0
#        for a = 1:Q.q
#            if a != awt
#                seq[i] = a
#                df_value = fitness(seq, Q) - f_wt
#                k += 1
#                df[k] = df_value
#                k_all += 1
#                delta_f_all[k_all] = df_value
#            end
#        end
#        seq[i] = awt
#        if option == "mean"
#            delta_f[i] = mean(df)
#        elseif option == "worst"
#            delta_f[i] = min(df)
#        elseif option == "least"
#            delta_f[i] = max(df)
#        end
#    end
#    if option == "all"
#        return delta_f_all
#    else
#        return delta_f
#    end
#end
#
#
#function computeMagnetization(seq::Sequence, K::Table, sites::Sites, fields::Fields, Q::Settings)
#    # magnetization at each site
#    W, L = Q.W, Q.L
#    J, β = seq2J(seq, K), Q.β
#    conf_left, conf_right = ones(2^W), ones(2^W)
#    Mmat = [zeros(W,L+1) for i in eachindex(sites)]
#    H = copy(J[1])
#    Mat = zeros(2^W,2^W)
#    states = getSpinStates(W)
#    
#    for j in eachindex(sites) # for each perturbation
#        site, field = sites[j], fields[j]
#        #h00 = J[1][site] # save old feild
#        J[1][site] += field # bind ligand by adding fields
#        for k = 1:(L+1) # for each layer
#            Z_right = conf_right
#            for layer in 1:(k-1)
#                MatLayer!(Mat, states, J, layer, β)
#                Z_right = Mat * Z_right
#            end
#            Z_left = conf_left'
#            for layer in (L+1):-1:k
#                MatLayer!(Mat, states, J, layer, β)
#                Z_left = Z_left * Mat
#            end
#            Z_den = Z_left*Z_right
#            for i = 1:W
#                seq = zeros(W)
#                seq[i] = 1
#                #Z_num = Z_left * (2*diagm(seq2conf(seq)) - eye(2^W)) * Z_right
#                Z_num = Z_left * (2*diagm(seq2conf(seq)) - I) * Z_right
#                Mmat[j][i,k] = Z_num / Z_den
#            end
#        end
#        J[1] .= H # unbind ligand
#    end
#    return Mmat
#end


