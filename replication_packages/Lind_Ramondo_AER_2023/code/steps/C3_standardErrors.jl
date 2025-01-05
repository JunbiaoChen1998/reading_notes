
function calculateScoreMatrix(lf::LatentFactors,tdShares::TradeData)

    X̂k = [ predictXk(lf,k) for k=1:length(lf.η) ]
    X̂Filled = zeros(size(X̂k[1]))
    for x in X̂k
        X̂Filled += x
    end
    fillMissing!(tdShares,X̂Filled)
    R = predictRatios(tdShares,X̂Filled)
    miss = tdShares.sitcValuesMiss
    out = []
    for (k,x) in enumerate(X̂k)
        xScr = tdShares.logSitcTariffs .* x .*( R .- 1 )
        resids = []
        for m=1:size(miss,2)
            push!(resids,tdShares.secMap[miss[:,m],:]'*xScr[miss[:,m],m])
        end
        push!(out,vcat(xScr[.!tdShares.sitcValuesMiss],resids...))
    end
    return hcat(out...)
end


poissonNormalization(x) = mean(x)/var(x)
function vcovMatLFM(lf::LatentFactors,data::ProjData)
    tdShares,total = shares(data.trade)
    Ysitcshares = getObservations(tdShares)
    scaling = poissonNormalization(Ysitcshares)
    poissonNormalization(scaling.*Ysitcshares)
    scoreMatrix = scaling*calculateScoreMatrix(lf,tdShares)

    fisherInfo = (scoreMatrix'*scoreMatrix)
    vcovMat = inv(fisherInfo)
    return vcovMat
end

function standardErrors(dirs::ProjDirectories,data::ProjData,K::Int64)
    bestReps = load(dirs.results*"estimates/LFM/bestReps.csv") |> DataFrame
    name = bestReps.rep[findfirst(bestReps.k .== "$K")]
    lf = loadLatentFactors(dirs.results*"estimates/LFM/$K/",name,data.trade)
    vcovMat = vcovMatLFM(lf,data)

    ids = sortperm(lf.η,rev=true)
    save(dirs.results*"estimates/LFM/vcovMatLFM.csv",DataFrame(vcovMat))
    save(dirs.results*"estimates/LFM/elasticitiesWithSEs.csv",DataFrame(k=1:K,σ̂ = lf.η[ids],se=sqrt.(diag(vcovMat))[ids]))

end
