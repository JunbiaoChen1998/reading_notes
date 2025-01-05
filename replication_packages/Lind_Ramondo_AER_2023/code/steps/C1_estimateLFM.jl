
include(dirs.code*"tools/fitLatentFactorModelPPML.jl")

function shares(td::TradeData)
    tdShares = copy(data.trade)
    total = sum(data.trade.secValues,dims=1)
    tdShares.sitcValues ./= total
    tdShares.secValues ./= total
    tdShares.secResiduals ./= total
    tdShares.sitcValuesFilled ./= total
    return tdShares,total
end


function estimateLFM(dirs::ProjDirectories,data::ProjData,K::Int64)

    tdShares,total = shares(data.trade)

    doneReps = map(x->parse(Int64,x[4:end]),filter(x->x[1:3] .== "rep",readdir(dirs.results*"estimates/LFM/$K")))
    isempty(doneReps) ? rep = 1 : rep = minimum(setdiff(1:(1+maximum(doneReps)),doneReps))
    name = "rep$rep"
    in(name,readdir(dirs.results*"estimates/LFM/$K")) ? error("Attempting to use an existing repetition name.") : nothing

    dir = dirs.results*"estimates/LFM/"
    S,M = size(tdShares.sitcValues)
    data.J != size(tdShares.secValues,1) ? error("Dimension mismatch.") : nothing

    η0 = fill(10. * rand(),K)
    Λ0 = rand(S,K)
    Λ0 = Λ0 ./ sum(Λ0,dims=1)
    Φ0 = rand(K,M)

    lf = LatentFactors(η0,[tdShares.sitcTariffs.^-η0[k] for k=1:K],[Λ0[:,k] for k=1:K],[Φ0[k,:] for k=1:K])
    Y = getObservations(tdShares)
    Ŷ = getPredictions(tdShares,predictX(lf))

    history = DataFrame(iter=0,pd=poissonDeviance(Y,getPredictions(tdShares,predictX(lf))),R2=R2(Y,Ŷ))
    fitLatentFactorModel!(tdShares,lf,history,5,0)

    history = history[[nrow(history)],:]
    fitLatentFactorModel!(tdShares,lf,history,350,10)

    saveLatentFactors(dir*"$K/",name,tdShares,lf)
    save(dir*"$K/$name/history.csv",history)

    return nothing
end

function estimateSGM(dirs::ProjDirectories,data::ProjData)

    tdShares,total = shares(data.trade)

    doneReps = map(x->parse(Int64,x[4:end]),filter(x->x[1:3] .== "rep",readdir(dirs.results*"estimates/LFM/14SGM")))
    isempty(doneReps) ? rep = 1 : rep = minimum(setdiff(1:(1+maximum(doneReps)),doneReps))
    name = "rep$rep"
    in(name,readdir(dirs.results*"estimates/LFM/14SGM")) ? error("Attempting to use an existing repetition name.") : nothing

    K = 14
    dir = dirs.results*"estimates/LFM/"
    name = "rep$rep"
    S,M = size(tdShares.sitcValues)
    η0 = fill(10. * rand(),K)
    Λ0 = tdShares.secMap .* rand(S,K)
    Λ0 = Λ0 ./ sum(Λ0,dims=1)
    Φ0 = rand(K,M)

    lf = LatentFactors(η0,[tdShares.sitcTariffs.^-η0[k] for k=1:K],[Λ0[:,k] for k=1:K],[Φ0[k,:] for k=1:K])
    Y = getObservations(tdShares)
    Ŷ = getPredictions(tdShares,predictX(lf))

    history = DataFrame(iter=0,pd=poissonDeviance(Y,getPredictions(tdShares,predictX(lf))),R2=R2(Y,Ŷ))
    fitSectoralGravityModel!(tdShares,lf,history,5,0)
    
    history = history[[nrow(history)],:]
    fitSectoralGravityModel!(tdShares,lf,history,350,10)

    saveLatentFactors(dir*"14SGM/",name,tdShares,lf)
    save(dir*"14SGM/$name/history.csv",history)

    return nothing
end
