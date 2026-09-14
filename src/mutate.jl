


#function mutate!(seq::Sequence, Q::Settings)
#    # Mutate 1 postion of the sequence randomly by switching to new type.
#    imut, amut = rand(1:Q.W*(Q.L+1)), rand(1:(Q.q-1))
#    awt = seq[imut]
#    seq[imut] = ((awt-1 + amut) % Q.q) +1
#    return nothing
#end

"""
Mutate each position of sequence with probability μ.
"""
function mutateAtRate!(seq::Sequence, μ::Number, q::Integer)
    didMutate = false
    for i in eachindex(seq)
        if rand() < μ
            didMutate = true
            awt = seq[i]
            amut = rand(1:(q-1))
            seq[i] = ((awt-1 + amut) % q) +1
        end
    end
    return didMutate
end
mutateAtRate!(seq::Sequence, μ::Number, Q::Settings) = mutateAtRate!(seq, μ, Q.q)

#mutatePosition!(seq::Sequence, position::Integer, newType::Integer) = seq[position]=newType
