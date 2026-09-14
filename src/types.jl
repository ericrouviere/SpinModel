# All of the structs are save in this file.


mutable struct Settings
    # Physical settings of the Spin Model that are invariant over the ensemble.
    W::Int # "width" of model
    L::Int # "Length" -1 of model
    q::Int # alphabet size
end

mutable struct EvoParams
    # Parameters of the evolutionary process.
    P::Int # population size
    N::Int # number of generations
    μ::Float64 # mutation rate
    τ::Int # environmental timescale
    α::Float64 # selection strength
end

mutable struct Perturbation
    # perturbations to the couplings and feilds (in J) to model
    # ligand binding.
    sites::Vector{CartesianIndex{3}}
    fields::Vector{Float64}
end


mutable struct Ligands
    perturbs::Vector{Perturbation}
end
Base.length(ligs::Ligands) = length(ligs.perturbs)
Base.getindex(ligs::Ligands, i::Integer) = ligs.perturbs[i]


# Define alias
const Sites = Vector{Vector{CartesianIndex{2}}}
const Fields = Vector{Vector{Float64}}
const Table = Array{Float64,4}
const Sequence = Vector{Int8}
const PopArray = Array{Int8, 3}
const TableArray = Array{Float64, 5}

