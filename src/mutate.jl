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
