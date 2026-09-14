#function computeCorrelation(J::Array{Float64, 3}, Q::Settings)
#    # return the spin Correlation matrix for each ligand perturbation.
#
#    W, L, β = Q.W, Q.L, Q.β
#    Z_right = Vector{Float64}(undef, 2^W)  
#    Z_left  = Vector{Float64}(undef, 2^W)
#    Z_middle = Matrix{Float64}(undef, 2^W, 2^W)
#    buff_vec = Vector{Float64}(undef, 2^W)
#    buff_mat = Matrix{Float64}(undef, 2^W, 2^W)
#    T= [Matrix{Float64}(undef, 2^W, 2^W) for i in 1:(L+1)]
#    C = Array{Float64,4}(undef, W, L+1, W, L+1)
#    J_mat = Matrix{Float64}(undef, W, W) #
#    J_mat_states = Matrix{Float64}(undef, W, 2^W) #
#    h_states = Vector{Float64}(undef, 2^W) #
#    states = getSpinStates(W)
#
#    function mul_update!(B::AbstractArray, buff::AbstractArray, A::AbstractArray)
#        # does the following opperation inplace: B <-- A*B
#        mul!(buff, A, B)
#        B .= buff
#        return nothing
#    end
#    function set2Identity!(A::Matrix)
#        # write the identity matrix to A
#        @inbounds for j in axes(A,2), i in axes(A,1)
#            A[i,j] = ifelse(i==j, one(eltype(A)), zero(eltype(A)))
#        end
#        return nothing
#    end
#
#    
#    function computeCoor!(C)
#        # compute the coorelations <σ_{i1,l1} σ_{i1,l1}>
#        # where i1 and i2 are the positions within a layer
#        # and l1, l2 are the layers.
#
#        # compute all transfer matrices and partition function
#        Z_right .= 1.0
#        for k in 1:L+1 # each layer k
#            updateTransferMatrix!(T[k], J_mat, J_mat_states, h_states, states, J, k, β)
#            mul_update!(Z_right, buff_vec, T[k])
#        end
#        Z = sum(Z_right)
#        
#        # compute transfer matrix product components 1...i1-1, i1...i2-1, i2...L+1
#        for l1 in 1:L+1, l2 in l1:L+1 # coorelation between layer i and j
#            Z_right .=  1.0
#            for k in 1:(l1-1)
#                mul_update!(Z_right, buff_vec, T[k])
#            end
#            
#            set2Identity!(Z_middle)
#            for k in l1:(l2-1)
#                mul_update!(Z_middle, buff_mat, T[k])
#            end
#
#            Z_left .= 1.0
#            for k in (L+1):-1:l2
#                mul_update!(Z_left, buff_vec, T[k]')
#            end
#            
#            for i2 in 1:W
#                @views S_i2 = states[i2,:]
#                h_states .= S_i2 .* Z_left # using h_states as buffer
#                for i1 in 1:W
#                    @views S_i1 = states[i1, :]
#                    # solve Z_left * S_i2 * Z_middle * S_i1 * Z_right
#                    buff_vec .= S_i1 .* Z_right
#                    coor = dot(h_states, Z_middle, buff_vec) / Z
#                    C[i1,l1,i2,l2] = C[i2, l2, i1, l1] = coor
#                end
#            end
#        end
#        return nothing
#    end
#    computeCoor!(C)
#    return C
#end
#
#function computeCorrelation(seq::Sequence, K::Table, Q::Settings)
#    J = seq2J(seq, K)
#    return computeCorrelation(J, Q)
#end
#
#
#function computeCorrelations(seq::Sequence, K::Table, Q::Settings, ligs::Ligands)
#    # compute the coorlation matrix for each ligand perturbation condition
#    J = seq2J(seq, K)
#    function kernel(i)
#        perturb = ligs[i]
#        J[perturb.sites] += perturb.fields # apply perturbation
#        C = computeCorrelation(J, Q)
#        J[perturb.sites] -= perturb.fields # remove perturbation
#        return C
#    end
#    return map(kernel, 1:length(ligs))
#end
#
#function computeOverlaps(seq::Sequence, K::Table,  Q::Settings, ligs::Ligands)
#    # return the otherlap of the conformational change upon ligand binding
#    # and the modes of the coupling matrix J and the correlation matrix C.
#    # also return the eigvalues of both matrices.
#
#    # compute magnetization
#    M_list = computeMagnetization(seq, K, Q, ligs)
#    m1 = M_list[1][:]
#    m2 = M_list[2][:]
#    dm = normalize(m2 - m1)
#
#    # compute Coorelation
#    C_tensor = computeCorrelation(seq, K, Q, ligs)
#    C = reshape(C_tensor, length(M_list[1]), length(M_list[1])) - m1 .* m1'
#    EC = eigen(C)
#    UC, vC = EC.vectors, EC.values
#    OC = UC'dm # overlaps
#
#  #  # get coupling matrix 
#  #  J = getJMatrix(seq, K)
#  #  EJ= eigen(J)
#  #  UJ, vJ = EJ.vectors, EJ.values
#  #  OJ = UJ'dm # overlaps
#    return OC, vC
#end


#function seq2conf(config)
#    # Return a vector in configuration space associated with a particular sequence
#    # if the input is a sequence of +/-
#    # if a position in the sequence is 0 it is ignored, i.e. all configurations are 
#    # considered that are consistent with the +/- irrespectively of positions with 0
#    # a sequence has size L leads to a vector of size 2^L
#    # usage: seq2conf([1,-1,0])
#    W = size(config, 1)
#    U = zeros(2^W)
#    for i = 0:2^W-1
#        sigma = round.(Int64, 2 .* (.5 .- digits(i, base=2, pad=W)))
#        if sigma[config.!=0] == config[config.!=0]
#            U[i+1] = 1
#        end
#    end
#    return U
#end


#function computeEffectiveField(seq::Sequence, K::Table, site::CartesianIndex, Q::Settings)
#    # computes the effective feild on a spin with 
#    # a give magentization.
#    β = Q.β
#    sites = [[site]]
#    fields = [[0.0]]
#    m = computeMagnetization(seq, K, ligands, Q)[1][site]
#    return atanh(m) / β
#end

#function overlaps(seq::Sequence, K::Table, Q::Settings)
#    # overlap of magnetization (h=0,h=hs) and (h=0,h=hp)
#    # scaled to between 0 and 1
#    M_mat = magnetization(seq, K, Q)
#    w0s = sum(abs.(M_mat[1]-M_mat[2]))/(2*Q.W*(Q.L+1))
#    w0p = sum(abs.(M_mat[1]-M_mat[3]))/(2*Q.W*(Q.L+1))
#    return w0s, w0p
#end
#
#
#
#





#function computeCorrelation(J::Array{Float64, 3}, Q::Settings)
#    # return the spin Correlation matrix for each ligand perturbation.
#
#    W, L, β = Q.W, Q.L, Q.β
#    Z_right = Vector{Float64}(undef, 2^W)  
#    Z_left  = Vector{Float64}(undef, 2^W)
#    Z_middle = Matrix{Float64}(undef, 2^W, 2^W)
#    buff_vec = Vector{Float64}(undef, 2^W)
#    buff_mat = Matrix{Float64}(undef, 2^W, 2^W)
#    T= [Matrix{Float64}(undef, 2^W, 2^W) for i in 1:(L+1)]
#    C = Array{Float64,4}(undef, W, L+1, W, L+1)
#    J_mat = Matrix{Float64}(undef, W, W) #
#    J_mat_states = Matrix{Float64}(undef, W, 2^W) #
#    h_states = Vector{Float64}(undef, 2^W) #
#    states = getSpinStates(W)
#
#    function mul_update!(B::AbstractArray, buff::AbstractArray, A::AbstractArray)
#        # does the following opperation inplace: B <-- A*B
#        mul!(buff, A, B)
#        B .= buff
#        return nothing
#    end
#    function set2Identity!(A::Matrix)
#        # write the identity matrix to A
#        @inbounds for j in axes(A,2), i in axes(A,1)
#            A[i,j] = ifelse(i==j, one(eltype(A)), zero(eltype(A)))
#        end
#        return nothing
#    end
#
#    
#    function computeCoor!(C)
#        # compute the coorelations <σ_{i1,l1} σ_{i1,l1}>
#        # where i1 and i2 are the positions within a layer
#        # and l1, l2 are the layers.
#
#        # compute all transfer matrices and partition function
#        Z_right .= 1.0
#        for k in 1:L+1 # each layer k
#            updateTransferMatrix!(T[k], J_mat, J_mat_states, h_states, states, J, k, β)
#            mul_update!(Z_right, buff_vec, T[k])
#        end
#        Z = sum(Z_right)
#        
#        # compute transfer matrix product components 1...i1-1, i1...i2-1, i2...L+1
#        for l1 in 1:L+1, l2 in l1:L+1 # coorelation between layer i and j
#            Z_right .=  1.0
#            for k in 1:(l1-1)
#                mul_update!(Z_right, buff_vec, T[k])
#            enda
#
#
#function getBins(bindingEnergies::Matrix; # has dims (popSize, N_ligs)
#                 N_bins=50,
#                 θ = pi/4,
#                 limits=[-2cos(θ), 2cos(θ)])
#    if maximum(abs, bindingEnergies) > 1
#        error("Current version of `getBins` only handles binding space in [-1, 1] × [-1, 1]")
#    end
#
#    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
#    rotBE = R * bindingEnergies'
#    coors = rotBE[1,:]
#    m, M = extrema(coors)
#    bins = LinRange(limits[1], limits[2], N_bins+1)
#    binvec = zeros(Int, size(bindingEnergies,1))
#    for i in eachindex(coors)
#        binvec[i] = findfirst(y -> y>=coors[i], bins)-1
#    end
#    return binvec
#end
