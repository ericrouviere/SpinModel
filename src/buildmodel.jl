function convertSettings(S::Dict)
    return Settings(S["W"], S["L"], S["q"])
end

function convertEvoParams(S::Dict)
    return EvoParams(S["P"], S["N"], Float64(S["μ"]), S["τ"], Float64(S["α"]))
end

function rand_table(W::Int, L::Int, q::Int, sh::Number, sJ::Number; J_pos=false)
    # random table for mapping sequences to fields and couplings
    # W,L are the size of the system
    # q is the size of the alphabet of sequences
    # the fields and couplings are drawn randomly from a normal distribution
    # with standard deviation sh for the fields and sJ for the couplings
    # first index of first dimention is the feilds, the 2nd and 3rd are
    # couplings
    K = randn(3, W, L+1, q^2)
    K[1,:,:,:] .*= sh
    if J_pos
        K[2:3,:,:,:] .= sJ*abs.(K[2:3,:,:,:])
    else
        K[2:3,:,:,:] .*= sJ
    end
    return K
end

"""
return random sequence of integers between 1 and q
of correct length.
"""
randSeq(q::Integer, N_sites::Integer) = Int8.(rand(1:q, N_sites))
randSeq(Q::Settings) = randSeq(Q.q, Q.W*(Q.L+1))


"""
In-place `seq2J`: fill `J`, of size `(W, L+1, 3)`, from `seq` and `K`.
Scanning the mutants of one background sequence calls this once per mutant, so the
allocation `seq2J` would make each time is worth avoiding.
"""
function seq2J!(J::AbstractArray{Float64,3}, seq::Sequence, K::Table)
    _, W, Lp1, q2 = size(K)
    L, q = Lp1 - 1, Int(sqrt(q2))

    for i = 1:W
        # fields (r=-1):
        for k = 1:(L + 1)
            si = seq[W * (k - 1) + i]
            J[i, k, 1] = K[1, i, k, si]
        end
        # couplings (r=0,1):
        for k = 1:L
            for r = 2:3
                if r == 2
                    si, sj = seq[W * (k - 1) + i], seq[W * k + i]
                else
                    si, sj = seq[W * (k - 1) + i % W + 1], seq[W * k + i]
                end
                J[i, k, r] = K[r, i, k, q * (si - 1) + sj]
            end
        end
    end
    J[:,L+1, 2:3] .= 0.    
    return J
end

"""
from seq to J using K
the mapping is such that the ligand corresponds to the last 
part of the sequence: (W*L+1):W*(L+1)
consider as alphabet 1,2,...,q
seq should of length W*(L+1)
K should be of dimension 3 x W x (L+1) x q^2
(some of its entries are not used)
"""
function seq2J(seq::Sequence, K::Table)
    _, W, Lp1, _ = size(K)
    return seq2J!(Array{Float64,3}(undef, W, Lp1, 3), seq, K)
end
