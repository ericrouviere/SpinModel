"""
Plot the fitnesses and binding energies as a function of evolutionary time.
"""
function plotFitnessEnergies(K::Table, Q::Settings, ligs::Ligands, assay::Assay,
                             times::Vector, fitnesses::Vector, sequences::Vector{Sequence})
    function vecOfTups2Matrix(vot)
        energies = zeros(length(vot), )
        [energies[i,:] .= vot[i] for i in eachindex(vot)]
        return energies
    end

    out= map(x -> computeFreeEnergies(x, K, Q, ligs), sequences )
    energies = zeros(length(out), length(ligs))
    [ energies[i,:] .= out[i] for i in eachindex(out)]
    
    fig, ax = subplots(1, 2, figsize=(10,3))
    fig.subplots_adjust(wspace=0.4)
    ax[1].plot(times, fitnesses)
    ax[1].set_xlabel("time", fontsize=16)
    ax[1].set_ylabel(L"\phi", fontsize=16)

    if assay isa Stability
        ax[2].plot(times, energies[:,1], c="b")
    elseif assay isa Binding
        ax[2].plot(times, energies[:,1], c="b")
        ax[2].plot(times, energies[:,2], c="g")
    elseif assay isa Union{Specificity, DoubleBinding}
        ax[2].plot(times, energies[:,1], c="b")
        ax[2].plot(times, energies[:,2], c="g")
        ax[2].plot(times, energies[:,3], c="r")
    elseif assay isa Union{Allostery, NegativeAllostery}
        ax[2].plot(times, energies[:,2] - energies[:,1], c="k")
        ax[2].plot(times, energies[:,4] - energies[:,3], c="r")
        ax[2].legend([L"\Delta F_0",L"\Delta F_1"])
    else
        error("That selective pressure is not accepted at this time.")
    end
    ax[2].set_xlabel("time", fontsize=16)
    ax[2].set_ylabel("Energy", fontsize=16)
    return fig, energies
end


"""
Visualization of fields and couplings
Note that y axis does not correspond to the indices in the matrix
(there is a translation that layer dependent to have a nice representation)
"""
function see_J!(ax,
                J::Dict;
                color_couplings=true,
                color_sites=true,
                maxCouplingSize=2,
                maxSiteSize=5,
                col_dict = Dict(true=>"r", false=>"b"),
                grey="0.8")

    W, L = size(J[0])
    for ori = 1:-1:-1
        for layer = 1:(L+(ori==-1))
            for index = 1:W
                Jval = J[ori][index, layer]
                col = col_dict[Jval>0]
                if ori < 0 # plot at lattice points
                    y0 = mod(index-1+.5*(layer-L-1),W)+1
                    ax.plot(L+1-layer+.5, y0, c="w", marker="o", mec="w", markersize=1.5*maxSiteSize)
                    !color_sites && (col=grey; Jval=1.5)
                    ax.plot(L+1-layer+.5, y0, c=col, marker="o", mec="w", markersize=maxSiteSize*sqrt(abs(Jval)))
                else # plot connections
                    !color_couplings && (col = grey; Jval=0.8)
                    y0, y1 = mod(index-1+.5*(layer-L),W)+1, mod(index-1+.5*(layer-L-1)+ori,W)+1
                    if abs(y0-y1) < 1
                        ax.plot([L-layer+.5, L-layer+1.5], [y0, y1], c=col, lw=.8*maxCouplingSize*abs(Jval))
                    else
                        if y0 == 1
                            ax.plot([L-layer+1, L-layer+1.5], [W+.75, y1], c=col, lw=.8*maxCouplingSize*abs(Jval))
                            ax.plot([L-layer+.5, L-layer+1], [y0, .75], c=col, lw=.8*maxCouplingSize*abs(Jval))
                        end
                        if y1 == 1
                            ax.plot([L-layer+.5, L-layer+1], [y0, W+.75], c=col, lw=.8*maxCouplingSize*abs(Jval))
                            ax.plot([L-layer+1, L-layer+1.5], [.75, y1], c=col, lw=.8*maxCouplingSize*abs(Jval))
                        end
                    end
                end
            end
        end
    end
    ax.axis("off");
    return nothing
end

"""
Plot the magnetization on the structure of the spin glass. 
"""
function plotMagnetization!(ax, M::Matrix; size_point=1)
    W, L = size(M,1), size(M,2)-1
    M_dict = Dict(-1 => size_point*M, 0 => zeros(W, L), 1 => zeros(W,L))
    see_J!(ax, M_dict, color_couplings=false)
    return nothing
end

"""
plot the magnetization and the difference in magnetization upon ligand binding for each ligand. 
"""
function plotMagnetizations(seq::Sequence, K::Table, Q::Settings, ligs::Ligands; size_point=1)
    M_list = computeMagnetization(seq, K, Q, ligs)
    fig, ax = subplots(2,3, figsize=(10,5))
    plotMagnetization!(ax[1,1], M_list[1]; size_point)
    ax[1,1].set_title(L"\langle \sigma_i \rangle_{s}")
    plotMagnetization!(ax[1,2], M_list[2]; size_point)
    ax[1,2].set_title(L"\langle \sigma_i \rangle_{r}")
    plotMagnetization!(ax[1,3], M_list[3]; size_point)
    ax[1,3].set_title(L"\langle \sigma_i \rangle_{w}")
    plotMagnetization!(ax[2,1], 0.5*(M_list[2] - M_list[1]); size_point)
    ax[2,1].set_title(L"\langle \sigma_i \rangle_{r} - \langle \sigma_i \rangle_{s}")
    plotMagnetization!(ax[2,2], 0.5*(M_list[3] - M_list[1]); size_point)
    ax[2,2].set_title(L"\langle \sigma_i \rangle_{w} -  \langle \sigma_i \rangle_{s} ")
    plotMagnetization!(ax[2,3], 0.5*(M_list[3] - M_list[2]); size_point)
    ax[2, 3].set_title(L"\langle \sigma_i \rangle_{w} - \langle \sigma_i \rangle_{r}")
    return fig
end


function plotMagnetizationsAllo(seq::Sequence, K::Table, Q::Settings, ligs::Ligands; size_point=1)
    # plot the magnetization and the difference in magnetization upon ligand binding for each ligand. 
    M_list = computeMagnetization(seq, K, Q, ligs)
    fig, ax = subplots(2,2, figsize=(7,5))
    plotMagnetization!(ax[1,1], M_list[1]; size_point)
    ax[1,1].set_title(L"\langle \sigma_i \rangle_{00}")
    plotMagnetization!(ax[2,1], M_list[2]; size_point)
    ax[2,1].set_title(L"\langle \sigma_i \rangle_{10}")
    plotMagnetization!(ax[1,2], M_list[3]; size_point)
    ax[1,2].set_title(L"\langle \sigma_i \rangle_{01}")
    plotMagnetization!(ax[2,2], M_list[4]; size_point)
    ax[2,2].set_title(L"\langle \sigma_i \rangle_{11}")
    return fig
end


"""
Plot magnitude of magnetization change between r and w ligands.
"""
function plotMagnetizationChange!(ax, seq::Sequence, K::Table, Q::Settings, ligs::Ligands;
                                  size_point=1,
                                  color_point="blue",
                                  maxCouplingSize=2)
    M_list = computeMagnetization(seq, K, Q, ligs)
    ΔM = abs.(M_list[3] .- M_list[2]) # |w - r|
    M_dict = Dict(-1 => size_point * ΔM, 0 => zeros(Q.W, Q.L ), 1 => zeros(Q.W, Q.L ))
    see_J!(ax, M_dict; color_couplings=false, col_dict=Dict(true=>color_point, false=>color_point),
          maxCouplingSize)
    ax.set_title("Magnetization Change", fontsize=16)
    return nothing
end

function plotMagnetizationChange(seq::Sequence, K::Table, Q::Settings, ligs::Ligands;
                                 size_point=1,
                                 maxCouplingSize=2)
    fig, ax = subplots(figsize=(3,2.3))
    plotMagnetizationChange!(ax, seq, K, Q, ligs; size_point, maxCouplingSize)
    return fig
end


"""
Plot magnitude of mean mutational sensitivity at each site.
"""
function plotMutationalSensitivity!(ax, seq::Sequence, K::Table, Q::Settings, ligs::Ligands; 
                                    size_point=1,
                                    maxCouplingSize=2)
    DMS = computeBindingDMS(seq, K, Q, ligs) # (q, L, ligs-1)
    data_1d = mean(abs, DMS, dims=(1,3))[1,:,1]
    data = reshape(data_1d, Q.W, Q.L+1) # Correct shape for see_J!
    M_dict = Dict(-1 => size_point * data, 0 => zeros(Q.W, Q.L), 1 => zeros(Q.W, Q.L))
    see_J!(ax, M_dict; color_couplings=false, col_dict=Dict(true=>"green", false=>"green"),
            maxCouplingSize)
    ax.set_title("Mutational Sensitivity", fontsize=16)
    return nothing
end
function plotMutationalSensitivity(seq::Sequence, K::Table, Q::Settings, ligs::Ligands;
                                   size_point=1,
                                   maxCouplingSize=2)
    fig, ax = subplots(figsize=(3,2.3))
    plotMutationalSensitivity!(ax, seq, K, Q, ligs; size_point, maxCouplingSize)
    return fig
end


"""
Plot the perturbation sensitivities for each site.
"""
function plotPerturbationScan(seq::Sequence, K::Table, Q::Settings, sites::Sites, fields::Fields; size_point=1)
    out = perturbationScan(seq, K, Q, sites, fields)
    data = mean(abs, out, dims=(1,4))[1,:,:,1] 
    M_dict = Dict(-1 => size_point*data, 0 => zeros(Q.W, Q.L), 1 => zeros(Q.W, Q.L))
    fig, ax = subplots(figsize=(4,3))
    see_J!(ax, M_dict; color_couplings=false, col_dict = Dict(true=>"g", false=>"r"))
    ax.set_title("Perturbation Sensitivity")
    return fig
end


"""
Plot the strength of the couplings on the lattice.
"""
function plotCouplings(seq::Sequence, K::Table, Q::Settings; maxCouplingSize=2, maxSiteSize=1)
    W,L = Q.W,Q.L
    J = seq2J(seq, K)
  #  J[:,:,2:3] .= -0.1
  #  J[CartesianIndex(3,L,3)] = 1
  #  J[CartesianIndex(3,1,2)] = 1
    J_dict = Dict(-1 => ones(W, L+1), 0 => J[:,1:L,2], 1 => J[:,1:L,3])
    fig, ax = subplots(figsize=(4,3))
    see_J!(ax, J_dict; color_sites=false, maxCouplingSize, maxSiteSize)

    return fig
end


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

"""
Plot change in magnetization due to ligand binding.
"""
function fig_states!(ax, Q::Settings, assay::Assay, seq::Sequence, K::Table; size_point = 1)
    # plot change in magnetization due to ligand binding.
    M_mat = computeMagnetization(seq, K, Q)

    if assay isa Specificity
        M_dict = Dict(-1 => size_point*.5*abs.(M_mat[2]-M_mat[1]),
                        0 => zeros(Q.W, Q.L), 1 => zeros(Q.W, Q.L))
        see_J!(ax, M_dict, col_J=false)
    elseif assay isa Allostery
        M_dict = Dict(-1 => size_point*.5*abs.(M_mat[4]-M_mat[1]), 
                        0 => zeros(Q.W, Q.L), 1 => zeros(Q.W, Q.L))
        see_J!(ax, M_dict, col_J=false)
    end
    return nothing
end

"""
Plot mutational sensitivity 
"""
function fig_mutagenesis!(ax, model, h0_list, seq; size_point = 1)
    delta_f = saturated_mutagenesis(model, h0_list, seq);
    J_dict = uni_J(model.W, model.L, 0)
    for i = 1:model.W
        for k = 1:(model.L+1)
            J_dict[-1][i,k] = size_point*delta_f[model.W*(k-1)+i]
        end
    end
    see_J!(ax, J_dict, col_J=false)
    return nothing
end


"""
Max effect on fitness of adding a field +/- h_add at each site
"""
function fig_perturb_field!(ax, model, h0_list, seq, h_add; size_point = 1)
    f_wt = fitness(model, h0_list, seq)
    J_dict = uni_J(model.W, model.L, 0)
    for i = 1:model.W
        for k = 1:(model.L+1)
            f_p = fitness(model, h0_list, seq, h_add=h_add, i=i, k=k)
            f_m = fitness(model, h0_list, seq, h_add=-h_add, i=i, k=k)
            if abs(f_m-f_wt) > abs(f_p-f_wt)
                f_p = f_m
            end
            J_dict[-1][i,k] = size_point*(f_p - f_wt)
        end
    end
    see_J!(ax, J_dict, col_J=false, col_dict = Dict(true=>"r", false=>"g"))
    return nothing
end


"""
Plot change in magnetization due to ligand binding.
"""
function computeWilsonCI(n::Int, p::Number, α::Number)
    z = cquantile(Normal(), α/2)
    d = (1 + z^2/n)
    center = (p + z^2/(2n)) / d
    spread = sqrt( p*(1-p)/n + z^2/(4n^2) )  / d
    lower = center - spread
    upper = center + spread
    return lower, upper
end


"""
this function takes the raw phenotypes data and 
returns the confidence intervals curves.
"""
function wilsonCIHelper(df::Vector, α::Number, ind::Int)
    CI = zeros(length(df), 2)
    for i in 1:length(df)
        data = df[i][:,ind]
        CI[i,:] .= computeWilsonCI(length(data), mean(float.(data)), α)
    end
    return CI
end


"""
n = size of sample
μ = mean of sample
σ = std of sample
α = confidence
"""
function computeWaldCI(n::Int, μ::Number, σ::Number, α::Number)
    z = cquantile(Normal(), α/2)
    lower = μ - z * σ / sqrt(n)
    upper = μ + z * σ / sqrt(n)
    return lower, upper
end


"""
 this function takes the raw phenotypes data and 
 returns the confidence intervals curves.
"""
function waldCIHelper(df::Vector, α::Number, ind::Int)
    CI = zeros(length(df), 2)
    for i in 1:length(df)
        data = df[i][:,ind]
        CI[i,:] .= computeWaldCI(length(data), mean(data), std(data), α)
    end
    return CI
end


"""
Plot the binding free energies points in negative binding space for many networks.

ΔF stores binding energies to ligand 1 in col 1 and ligand 2 in col 2.
"""

function plotDensityBindingSpaceScatter!(ax, ΔF::Matrix; buffer=0.04, alpha=0.1, markerSize=1) 
    mm,MM = extrema( [-ΔF[:]..., 0.0] )
    ax.scatter(-ΔF[:,1], -ΔF[:,2], color="k", s=markerSize,  alpha=alpha)
    ax.set_xlim([mm-buffer, MM+buffer])
    ax.set_ylim([mm-buffer, MM+buffer])
    ax.set_xlabel(L"-\Delta F_1", fontsize=16)
    ax.set_ylabel(L"-\Delta F_2", fontsize=16)
    ax.axvline(x=0, color="k", linestyle="--")
    ax.axhline(y=0, color="k", linestyle="--")
    ax.set_aspect("equal")
    return nothing
end
function plotDensityBindingSpaceScatter(ΔF::Matrix; buffer=0.04, alpha=0.1, markerSize=1) 
    fig, ax = subplots()
    plotDensityBindingSpaceScatter!(ax, ΔF; buffer, alpha, markerSize)
    return fig
end


"""
Plot population average negative binding energies over time.
ΔF1 and ΔF2 should be binding energies, ie negative for tight binding
"""
function plotBindingEvoTraces(ΔF1::Vector, ΔF2::Vector;
                              τ=0, μ=0, P=0, α=0,
                              color1="black", color2="lightgrey")
    fig, ax = subplots( figsize=(10,2))
    plotBindingEvoTraces!(ax, ΔF1, ΔF2; τ, μ, P, α, color1, color2)
    return fig
end

function plotBindingEvoTraces!(ax, ΔF1::Vector, ΔF2::Vector;
                              τ=0, μ=0, P=0, α=0, xlabel="generations",
                              plotTitle=true, plotLegend=true,
                              color1="black", color2="lightgrey")
    # Plotting the traces
    ax.plot(-ΔF1, lw=1.5, c=color1, label="Ligand 1")
    ax.plot(-ΔF2, lw=1.5, c=color2, label="Ligand 2")
    ax.set_ylim([-1.1, 1.1])
    
    # Adding vertical dashed lines every τ generations, behind the traces
    if τ > 0
        for i in 1:floor(length(ΔF1) / τ)
            ax.axvline(x=i*τ, color="black", linestyle="--", linewidth=0.5, zorder=-1)
        end
    end
    
    # Labels, limits, and other settings
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel("Binding Energy    \n"*L"-\langle \Delta F \rangle", fontsize=14, labelpad=-5)
    ax.set_xlim([0, length(ΔF1)])
    plotLegend && ax.legend()
    ax.axhline(y=0, color="grey", linestyle="-", linewidth=1, zorder=0)
    if plotTitle
        ax.set_title(L"\tau=" * @sprintf("%.3g", τ) * ", " *
                     L"\mu=" * @sprintf("%.3g", μ) * ", " *
                     L"P=" * @sprintf("%.3g", P) * ", " *
                     L"\alpha=" * @sprintf("%.3g", α))
    end
    return nothing
end



"""
    makepcoloraxes(x::AbstractVector)

Given a 1D vector `x` (e.g., bin centers), returns a new vector of edges suitable for use in pseudocolor plots (e.g., `pcolormesh` in matplotlib or `heatmap` in Julia).

Assumes `x` is approximately logarithmically spaced. The returned vector has length `length(x) + 1`, with inferred edges extrapolated at both ends to align with `x`.
"""
function makepcoloraxes(x)
    d = diff(x)
    p = x[1:end-1] .+ d ./ 2 
    p_first = exp(-diff(log.(p))[1])*p[1]
    p_last = exp(diff(log.(p))[1])*p[end]
    return [p_first; p ; p_last]
end


function plotHeatMap(x, y, z;
            xlabel="mutation rate "*L"\mu",
            ylabel="period "*L"\tau",
            title="Relative evolvability",
            fs=20, cmap="bwr_r", symCB=false)
    fig, ax = subplots()
    plotter!(fig, ax, x, y, z; xlabel, ylabel, title,fs, cmap, symCB)
    return fig
end

function plotHeatMap!(fig, ax, x, y, z;
            xlabel="mutation rate "*L"\mu",
            ylabel="period "*L"\tau",
            title="Relative evolvability",
            fs=20, cmap="bwr_r", symCB=false,
            vmin=nothing, vmax=nothing)
    px = makepcoloraxes(x)
    py = makepcoloraxes(y)
    if !isnothing(vmin) || !isnothing(vmax)
        im = ax.pcolor(py, px, z, cmap=cmap, vmin=vmin, vmax=vmax)
    elseif symCB # symetric colorbar range
        v = maximum(abs, z)
        im = ax.pcolor(py, px, z, cmap=cmap, vmin=-v, vmax=v)
    else
        im = ax.pcolor(py, px, z, cmap=cmap)
    end
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_yticks(x)
    ax.set_yticklabels(x)
    # Label every decade spanned by the x axis (e.g. 10^-4, 10^-3, 10^-2, 10^-1)
    # rather than letting matplotlib thin them out.
    ns = ceil(Int, log10(minimum(py))):floor(Int, log10(maximum(py)))
    ax.set_xticks([10.0^n for n in ns])
    ax.set_xticklabels([L"10^{%$n}" for n in ns])
    ax.set_xlabel(xlabel, fontsize=fs)
    ax.set_ylabel(ylabel, fontsize=fs)
    ax.set_title(title, fontsize=fs)
    clb = fig.colorbar(im, ax=ax)  # Ensure colorbar is associated with specific axes
    return nothing
end


