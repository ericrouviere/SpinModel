# Olivier Rivoire's code rewritten by Eric Rouviere in julia 1.7.
# April 2022

module SpinModel

using LinearAlgebra
using Distributions
using StatsBase
using PyPlot
using Random
using LaTeXStrings
using Distributed
using Printf
using JLD2
using Base: Threads

include("types.jl")
include("buildmodel.jl")
include("energies.jl")
include("fitness.jl")
include("evolve.jl")
include("mutate.jl")
include("analysis.jl")
include("ensemble.jl")
include("plotfuncs.jl")

export 

    # from types.jl
    Settings,
    EvoParams,
    Perturbation,
    Ligands,
    length,
    Sequence,

    # from buildmodel.jl
    convertSettings,
    convertEvoParams,
    rand_settings,
    rand_table,
    seq2J,
    seq2J!,
    TransferCache,
    mutantFreeEnergies!,
    randSeq,

    # from energies.jl
    computeFreeEnergies,
    computeFreeEnergiesFast,
    computeFreeEnergiesGeneral,
    computeFreeEnergy,
    getSpinStates,
    updateTransferMatrix!,
    computeBindingRandSeqs,

    # from fitness
    computeFitness,
    Assay,
    Stability,
    Binding,
    DoubleBinding,
    Specificity,
    Allostery,
    NegativeAllostery,
    Binding1,
    Binding2,

    # from evolve.jl
    evolve,
    evolvePop,
    evolvePopStatic,
    evolvePopSecondary,
    applySelection!,
    pickNumGens,
    evolvePop_saveBinding,
    evolvePop_saveSeqs,

    # from mutate.jl
    mutateAtRate!,

    # from analysis.jl
    computeMagnetization,
    computePartRatio,
    computeBindingDMS,
    computeBindingDoubleDMS,
    computeSingleMutantEvolvabilities,
    classifyMutantState,
    countMutantStates,
    computeAllMutantStates,
    computeSeqPhenotype,
    computeTrajPhenotypes,
    classifyTrajStates,
    computeStateOccupancy,
    foldPeriods,
    perturbationScan,
    scanAllostericSurface,
    computeAllostery,
    computeLigScape,
    analyseBindingDMS, 
    binBindingSpace,
    getBinnedSeqs,
    getBins,
    makeHistogram,
    binding2Histogram,

    # from ensemble.jl
    evolveEnsemble,
    evolve2Ensemble,      
    evolvePopEnsemble,
    evolvePopStaticEnsemble,
    evolvePop_saveBindingEnsemble,
    evolvePop_saveSeqsEnsemble,
    computeTrajPhenotypesEnsemble,
    computePartRatioEnsemble,
    computeAllosteryEnsemble,
    computeBindEnergiesEnsemble,
    perturbationScanEnsemble,
    convertPopulations2PopArray,
    convertPopArray2Populations,
    convertTables2TableArray,
    convertTableArray2Tables,
    computeBindingDMSEnsemble,
    analyseBindingDMSEnsemble,
    computeEffectiveFieldEnsemble,
    computeOverlapsEnsemble,
    buildAndAnalyze,

    # from plotfuncs.jl
    plotFitnessEnergies,
    see!,
    see_J!,
    fig_states!,
    fig_perturb_field!,
    fig_mutagenesis!,
    plotMagnetization!,
    plotMagnetizations,
    plotMagnetizationsAllo,
    plotMagnetizationChange,
    plotMagnetizationChange!,
    plotMutationalSensitivity,
    plotMutationalSensitivity!,
    plotSignSwitches,
    plotPerturbationScan,
    plotCouplings,
    plotLigScape,
    waldCIHelper,
    wilsonCIHelper,
    plotDensityBindingSpaceScatter!,
    plotDensityBindingSpaceScatter,
    plotBindingEvoTraces,
    plotBindingEvoTraces!,
    makepcoloraxes,
    plotHeatMap,
    plotHeatMap!

end
