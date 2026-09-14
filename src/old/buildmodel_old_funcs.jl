# Commented-out code moved from src/buildmodel.jl.

#function makeLigandList(assay::String,
#                        hs::Number, hr::Number, hw::Number,
#                        actSite::Vector{CartesianIndex{3}},
#                        alloSite::Vector{CartesianIndex{3}})
#    # Given the solvent, right and wrong ligand fields, and 
#    # the active and allosteric sites, return a list of sites
#    # and feilds that is nessicary for the assay.
#    if assay == "Stability"
#        sites = fill(actSite, 1)
#        fields = [[hs]]
#    elseif assay == "Binding"
#        sites = fill(actSite, 2)
#        fields = [[hs], [hr]]
#    elseif assay in ["Specificity", "GeneralSpecificity", "DoubleBinding"]
#        sites = fill(actSite, 3)
#        fields = [[hs], [hr], [hw]]
#    elseif assay == "Allostery"
#        sites = fill([actSite; alloSite], 4)
#        fields = [[hs, hs], [hr, hs], [hs, hr], [hr, hr]]
#    end
#    return Ligands( [Perturbation(sites[i], fields[i]) for i in eachindex(sites)]  )
#end

#function seq2J(seq::Sequence, K::Table)
#    # from seq to J using K
#    # the mapping is such that the ligand corresponds to the last 
#    # part of the sequence: (W*L+1):W*(L+1)
#    # consider as alphabet 1,2,...,q
#    # seq should of length W*(L+1)
#    # K should be of dimension 3 x W x (L+1) x q^2
#    # (some of its entries are not used)
#    _, W, Lp1, q2 = size(K)
#    L, q = Lp1-1, Int(sqrt(q2))
#    J = uni_J(W, L, 0)
#    for i = 1:W
#        # fields (r=-1):
#        for k = 1:(L+1)
#             si = seq[W*(k-1)+i]
#             J[1][i,k] = K[1, i, k, si]
#        end
#        # couplings (r=0,1):
#        for k = 1:L  
#            for r = 2:3
#                if r == 2
#                    si, sj = seq[W*(k-1)+i], seq[W*k+i]
#                else
#                    si, sj = seq[W*(k-1)+i%W+1], seq[W*k+i]
#                end
#                J[r][i,k] = K[r, i, k, q*(si-1) + sj]
#            end
#        end
#    end
#    return J
#end

#function getJMatrix(seq::Sequence, K::Table)
#    # return the coupling matrix J_ij where i,j are the
#    # ith and jth sites in the latice. 
#    # THIS IS A VERY CONFUSING FRANKENSTEIN FUNCTION
#    # THAT IS A QUICK FIX OF OLIVIER'S FUNCTIONS.
#
#    function findmy(x0::Number, x1::Number, y0::Number, y1::Number)
#        # find the linear indices of sites on lattice that have
#        # positions x0, x1, y0, y1.
#        x0_inds = findall( v -> v≈x0, x)
#        y0_inds = findall( v -> v≈y0, y)
#        i = intersect(x0_inds,y0_inds)[1] # this is the ith site in model
#        x1_inds = findall( v -> v≈x1, x)
#        y1_inds = findall( v -> v≈y1, y)
#        j = intersect(x1_inds,y1_inds)[1] # this is the jth site in model
#        return i,j
#    end
#
#    J_temp = seq2J(seq, K)
#    W, L = size(J_temp, 1), size(J, 2)-1
#    J = Dict(-1 => J_temp[:,:,1], 0 => J_temp[1:W,1:L,2], 1 => J_temp[1:W,1:L,3])
#
#    # Get the x and y positions of sites at (layer, index) on square lattice.
#    x = zeros(W*(L+1))
#    y = similar(x)
#    k = 0
#    for layer = 1:(L+1), index = 1:W
#        k += 1
#        Jval = J[-1][index, layer]
#        x[k] = L+1-layer+.5
#        y[k] = mod(index-1+.5*(layer-L-1),W)+1
#    end
#
#    # Fill J_Matrix with interactions between sites. 
#    J_matrix = zeros(length(x), length(x))
#    for ori in 1:-1:0
#        for layer = 1:L, index = 1:W
#            Jval = J[ori][index, layer]
#            y0, y1 = mod(index-1+.5*(layer-L),W)+1, mod(index-1+.5*(layer-L-1)+ori,W)+1
#            x0, x1 = L-layer+.5, L-layer+1.5
#            i,j = findmy(x0, x1, y0, y1)
#            J_matrix[i,j] = J_matrix[j,i] = Jval
#        end
#    end
#    return J_matrix
#end
