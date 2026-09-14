# Commented-out code moved from src/plotfuncs.jl.

# From inside plotCouplings, right after `J = seq2J(seq, K)`:
  #  J[:,:,2:3] .= -0.1
  #  J[CartesianIndex(3,L,3)] = 1
  #  J[CartesianIndex(3,1,2)] = 1

#"""
#Plot Free energy as a function of active site ligand field, h.
#"""
#function plotLigScape(seq::Sequence, K::Table, Q::Settings, site::CartesianIndex; h_range=[-2, 4], n=50)
#    h_list = LinRange(h_range[1], h_range[2], n)
#    energies = computeLigScape(seq, K, Q, site; h_range, n)
#    fig, ax = subplots(figsize=(4,3))
#    fig.subplots_adjust(bottom = 0.2, left = 0.2)
#    ax.plot(h_list, energies, c="b")
#    ax.set_xlabel(L"h", fontsize=14)
#    ax.set_ylabel(L"F", fontsize=14)
#    return fig
#end
#
