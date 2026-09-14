# Commented-out code moved from src/mutate.jl.

#function mutate!(seq::Sequence, Q::Settings)
#    # Mutate 1 postion of the sequence randomly by switching to new type.
#    imut, amut = rand(1:Q.W*(Q.L+1)), rand(1:(Q.q-1))
#    awt = seq[imut]
#    seq[imut] = ((awt-1 + amut) % Q.q) +1
#    return nothing
#end

#mutatePosition!(seq::Sequence, position::Integer, newType::Integer) = seq[position]=newType
