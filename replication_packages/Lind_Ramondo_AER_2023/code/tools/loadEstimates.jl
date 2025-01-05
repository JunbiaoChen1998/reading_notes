
function loadEstimates(dirs::ProjDirectories)
    # load estimates
    dfσ = DataFrame(load(dirs.results*"estimates/LFM/elasticitiesWithSEs.csv"))
    σ = dfσ.σ̂
    σse = dfσ.se
    dfΛ = load(dirs.results*"estimates/LFM/factorWeights.csv") |> DataFrame
    dfΦ = load(dirs.results*"estimates/LFM/factors.csv") |> DataFrame

    K = length(σ)
    K != length(unique(dfΛ.k)) && error("The number of factors is not consistent.")
    K != length(unique(dfΦ.k)) && error("The number of factors is not consistent.")

    sort(dfΛ,[:sitc,:k]) != dfΛ && error("The saved factor weights are not sorted.")

    # CES/SG< benchmark
    dfCES = load(dirs.results*"estimates/CES.csv") |> DataFrame
    dfSGM = load(dirs.results*"estimates/SGM.csv") |> DataFrame

    return σ,σse,dfΛ,dfΦ,K,dfCES,dfSGM
end
