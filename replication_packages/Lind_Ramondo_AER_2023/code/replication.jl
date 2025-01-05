#=

    This file replicates the results in Lind and Ramondo (2023)

    Nels Lind, January 27, 2023

=#

# set current directory to location of the folder containing this file
cd("C:/Users/nelsl/Dropbox/projects/Project Correlation/LR-replication-files")

# activate and instantiate the environment defined by the Manifest.toml and Project.toml files
using Pkg; Pkg.activate("."); Pkg.instantiate()

# load packages
using Queryverse, LinearAlgebra, StatsBase, Statistics, Random, Distributions
using GLFixedEffectModels, GLM, Flux
using Plots, Plots.PlotMeasures, StatsPlots, FreqTables, RegressionTables, Latexify

# directory structure
struct ProjDirectories
    proj::String
    code::String
    data::String
    results::String

    function ProjDirectories(proj)
        code = proj*"/code/"
        data = proj*"/data/"
        results = proj*"/results/"
        return new(proj,code,data,results)
    end
end
dirs = ProjDirectories(pwd())

# load tools
include(dirs.code*"tools/RegressionTable.jl")
include(dirs.code*"tools/rankCountries.jl")
include(dirs.code*"tools/loadSecData.jl")
include(dirs.code*"tools/loadEstimates.jl")
include(dirs.code*"tools/utilities.jl")
include(dirs.code*"tools/secNames.jl")

# A1: Prepare data for estimation
include(dirs.code*"steps/A1_prepareData.jl")
prepareDataFiles(dirs.data)
data = ProjData(dirs);

# A2: Prepare input-output shares
include(dirs.code*"steps/A2_prepareUseShares.jl")
prepareUseShares(dirs.data,data.countries,data.years)

# B1: gravity based CES/SGM estimation, IIA tests
include(dirs.code*"steps/B1_reducedFormEvidence.jl")
reducedFormEvidence(dirs,data)

# C1: LFM estimation, change "false" on the following line to "true" in order to do one repetition of algorithm
include(dirs.code*"steps/C1_estimateLFM.jl")
if false
    map(K->estimateLFM(dirs,data,K),vcat(1:8,14))
    estimateSGM(dirs,data)
end

# C2: Select LFM model order and create table witih specification tests
include(dirs.code*"steps/C2_selectOrderAndSave.jl")
selectOrderAndSave(dirs,data)

# C3: Compute standard errors on LFM estimates
include(dirs.code*"steps/C3_standardErrors.jl")
standardErrors(dirs,data,7)

# C4: Alternative estimation for θ
include(dirs.code*"steps/C4_secondStep.jl")
estimateSecondStep(dirs,data)

# D1: Calculate results from estimates
include(dirs.code*"tools/loadSitcDescriptions.jl")
include(dirs.code*"steps/D1_results.jl")
computeResults(dirs,data)

# D2: Evaluate goodness of fit
include(dirs.code*"steps/D2_goodnessOfFit.jl")
goodnessOfFit(dirs,data)

# D3: Calculate implied substitution elasticities
include(dirs.code*"steps/D3_elasticities.jl")
computeElasticities(dirs,data)

# E1: Counterfactual analysis
include(dirs.code*"steps/E1_counterfactuals.jl")
computeCounterfactuals(dirs,data)

# E2: Calculate gains from trade
include(dirs.code*"steps/E2_gainsFromTrade.jl")
computeGainsFromTrade(dirs,data)
