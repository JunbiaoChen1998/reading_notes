
using Base.Threads
import Base.copy
function copy(td::TradeData)
    values = []
    for name in fieldnames(TradeData)
        push!(values,copy(getfield(td,name)))
    end
    return TradeData(values...)
end

struct LatentFactors
    η::Vector{Float64}
    W::Vector{Matrix{Float64}}
    Λ::Vector{Vector{Float64}}
    Φ::Vector{Vector{Float64}}
end
function copy(lf::LatentFactors)
    η = []
    W = []
    Λ = []
    Φ = []
    for k=1:length(lf.η)
        push!(η,copy(lf.η[k]))
        push!(W,copy(lf.W[k]))
        push!(Λ,copy(lf.Λ[k]))
        push!(Φ,copy(lf.Φ[k]))
    end
    return LatentFactors(η,W,Λ,Φ)
end
function saveLatentFactors(dir::String,name::String,td::TradeData,lf::LatentFactors)
    K = length(lf.η)
    isdir(dir*name) || mkpath(dir*name)
    save(dir*"$(name)/eta.csv",DataFrame(η=lf.η))
    mValues = DataFrame(o=String[],d=String[],t=Int64[])
    for x in td.mValues
        o,d,t = x
        push!(mValues,[o d t])
    end
    save(dir*"$(name)/mValues.csv",mValues)
    save(dir*"$(name)/Lambda.csv",DataFrame(hcat(lf.Λ...),[Symbol("Λ"*string(k)) for k=1:K]))
    save(dir*"$(name)/Phi.csv",DataFrame(hcat(lf.Φ...),[Symbol("Φ"*string(k)) for k=1:K]))
    return nothing
end
function loadLatentFactors(dir::String,name::String,td::TradeData)
    η = DataFrame(load(dir*"$(name)/eta.csv")).η
    K = length(η)
    mValues = DataFrame(o=String[],d=String[],t=Int64[])
    for x in td.mValues
        o,d,t = x
        push!(mValues,[o d t])
    end
    mValues0 = load(dir*"$(name)/mValues.csv") |> DataFrame

    Λ = DataFrame(load(dir*"$(name)/Lambda.csv"))
    Φ = join(mValues,hcat(mValues,DataFrame(load(dir*"$(name)/Phi.csv"))),on=[:o,:d,:t],kind=:left)

    return LatentFactors(η,[td.sitcTariffs.^-η[k] for k=1:K],[Λ[:,Symbol("Λ$k")] for k=1:K],[Φ[:,Symbol("Φ$k")] for k=1:K])
end


predictXk(lf::LatentFactors,k::Int64) = lf.W[k] .* (lf.Λ[k] * lf.Φ[k]')
function predictX(lf::LatentFactors)
    out = zeros(size(lf.W[1]))
    for k=1:length(lf.η)
        out += predictXk(lf,k)
    end
    return out
end

function fillMissing!(td::TradeData,X::Matrix{Float64})
    @threads for m=1:size(X,2)
        miss = td.sitcValuesMiss[:,m]
        residuals = td.secMap[miss,:]'*X[miss,m]
        X[miss,m] = residuals[td.secInd[miss]]
    end
    return nothing
end
function fillMissing!(td::TradeData,X::Matrix{Float64})
    @threads for m=1:size(X,2)
        miss = td.sitcValuesMiss[:,m]
        residuals = td.secMap[miss,:]'*X[miss,m]
        X[miss,m] = residuals[td.secInd[miss]]
    end
    return nothing
end
function predictRatios(td::TradeData,X̂Filled::Matrix{Float64})
    out = td.sitcValuesFilled ./ X̂Filled
    z = X̂Filled .<= eps()
    out[td.sitcValuesZero .& (.!z) ] .= 0
    out[td.sitcValuesZero .& z] .= 1
    return out
end
function updateΦ!(td::TradeData,lf::LatentFactors)
    X̂Filled = predictX(lf)
    fillMissing!(td,X̂Filled)
    R = predictRatios(td,X̂Filled)
    for k=1:length(lf.Φ)
        lf.Φ[k] .*= ( (lf.W[k].*R)'*lf.Λ[k] ) ./ (lf.W[k]'*lf.Λ[k])
    end
    return nothing
end
function updateΦToMachinePrecision!(td::TradeData,lf::LatentFactors)
    X̂Filled = predictX(lf)
    fillMissing!(td,X̂Filled)
    R = predictRatios(td,X̂Filled)
    for k=1:length(lf.Φ)
        lf.Φ[k] .*= max.(eps(),( (lf.W[k].*R)'*lf.Λ[k] ) ./ (lf.W[k]'*lf.Λ[k]))
    end
    return nothing
end
function updateΛ!(td::TradeData,lf::LatentFactors)
    X̂Filled = predictX(lf)
    fillMissing!(td,X̂Filled)
    R = predictRatios(td,X̂Filled)
    for k=1:length(lf.Λ)
        Λnew = lf.Λ[k] .* ( ((lf.W[k].*R)*lf.Φ[k]) ./ (lf.W[k]*lf.Φ[k]) )
        scaling = sum(Λnew,dims=1)
        lf.Λ[k] .= Λnew ./ scaling
        lf.Φ[k] .*= scaling
    end
    return nothing
end
function updateη!(td::TradeData,lf::LatentFactors)
    X̂k = [ predictXk(lf,k) for k=1:length(lf.η) ]
    X̂Filled = zeros(size(X̂k[1]))
    for x in X̂k
        X̂Filled += x
    end
    fillMissing!(td,X̂Filled)
    R = predictRatios(td,X̂Filled)
    for (k,x) in enumerate(X̂k)
        z = td.logSitcTariffs .* x
        lf.η[k] *= ( sum( z ) / sum( z .* R ) )
        lf.W[k] .= td.sitcTariffs.^-lf.η[k]
    end
    return nothing
end

function getObservations(td::TradeData)
    vals = td.sitcValues
    resids = td.secResiduals
    return vcat(vals[.!td.sitcValuesMiss],vec(resids))
end
function getResiduals(td::TradeData,X::Matrix{Float64})
    miss = td.sitcValuesMiss
    resids = []
    for m=1:size(miss,2)
        push!(resids,td.secMap[miss[:,m],:]'*X[miss[:,m],m])
    end
    return hcat(resids...)
end
function getPredictions(td::TradeData,X::Matrix{Float64})
    return vcat(X[.!td.sitcValuesMiss],vec(getResiduals(td,X)))
end

poissonLogLik(Y,Ŷ) = sum( -Ŷ + X .* log.(Ŷ) - log.(factorial.(X)) )
pseudoPoissonLogLik(Y,Ŷ) = sum( -Ŷ + Y .* log.(Ŷ)  )
function poissonDeviance(Y,Ŷ)
    yly = Y .* log.(Y ./ Ŷ)
    yly[Y .== 0] .= 0
    return 2*sum( yly .- ( Y .- Ŷ ) )
end
R2(Y,Ŷ) = 1-sum((Y-Ŷ).^2)/sum(Y.^2)

function Δ(lf0::LatentFactors,lf1::LatentFactors)
    Δη = []
    ΔΛ = []
    ΔΦ = []
    for k=1:length(lf0.η)
        push!(Δη,lf0.η[k] - lf1.η[k])
        push!(ΔΛ,lf0.Λ[k] - lf1.Λ[k])
        push!(ΔΦ,lf0.Φ[k] - lf1.Φ[k])
    end
    return Δη,hcat(ΔΛ...),hcat(ΔΦ...)
end
function fitLatentFactorModel!(td::TradeData,lf::LatentFactors,history::DataFrame,maxit::Int64,Δt=1::Number)
    iter = history.iter[end]
    ηHistory = permutedims(lf.η)
    Y = getObservations(td)
    Ỹ = Y.-mean(Y)
    # lf_old = copy(lf)

    tStart = time()
    nextTime = Δt
    for i = 1:maxit
        iter += 1
        iterTime = time()-tStart
        @show iter, iterTime
        updateΦ!(td,lf)
        updateΛ!(td,lf)
        updateη!(td,lf)
        if iterTime > nextTime # mod(iter-1,checkRate) == 0 && iter > plotBurn
            nextTime += Δt
            Ŷ = getPredictions(td,predictX(lf))
            ηHistory = vcat(ηHistory,permutedims(lf.η))
            # Δη,ΔΛ,ΔΦ = Δ(lf,lf_old)
            push!(history,[ iter poissonDeviance(Y,Ŷ) R2(Ỹ,Ŷ) ])
            plots = []
            push!(plots,plot(history.iter,history.pd,legend=false,title="Poisson Deviance"))
            push!(plots,plot(history.iter,history.R2,legend=false,title="R Squared"))
            push!(plots,plot(ηHistory,legend=false,title="Eta"))
            # push!(plots,plot(history.iter,history.Δη,legend=false,title="Maximum Change in Eta",ylim=(0,maximum(history.Δη))))
            # push!(plots,plot(history.iter,history.ΔΛ,legend=false,title="Maximum Change in Lambda",ylim=(0,maximum(history.ΔΛ))))
            # push!(plots,plot(history.iter,history.ΔΦ,legend=false,title="Maximum Change in Phi",ylim=(0,maximum(history.ΔΦ))))
            display(plot(plots...,layout=(length(plots),1),size=(600,1200)))
        end

        # lf_old = copy(lf)
    end
    if iter > history.iter[end]
        Ŷ = getPredictions(td,predictX(lf))
        push!(history,[ iter poissonDeviance(Y,Ŷ) R2(Ỹ,Ŷ) ])
    end
    return nothing
end

function fitSectoralGravityModel!(td::TradeData,lf::LatentFactors,history::DataFrame,maxit::Int64,Δt=1::Number)
    iter = history.iter[end]
    ηHistory = permutedims(lf.η)
    Y = getObservations(td)
    Ỹ = Y.-mean(Y)
    # lf_old = copy(lf)

    tStart = time()
    nextTime = Δt
    for i = 1:maxit
        iter += 1
        iterTime = time()-tStart
        @show iter, iterTime
        updateΦToMachinePrecision!(td,lf)
        updateΛ!(td,lf)
        updateη!(td,lf)
        if iterTime > nextTime # mod(iter-1,checkRate) == 0 && iter > plotBurn
            nextTime += Δt
            Ŷ = getPredictions(td,predictX(lf))
            ηHistory = vcat(ηHistory,permutedims(lf.η))
            # Δη,ΔΛ,ΔΦ = Δ(lf,lf_old)
            push!(history,[ iter poissonDeviance(Y,Ŷ) R2(Ỹ,Ŷ) ])
            plots = []
            push!(plots,plot(history.iter,history.pd,legend=false,title="Poisson Deviance"))
            push!(plots,plot(history.iter,history.R2,legend=false,title="R Squared"))
            push!(plots,plot(ηHistory,legend=false,title="Eta"))
            # push!(plots,plot(history.iter,history.Δη,legend=false,title="Maximum Change in Eta",ylim=(0,maximum(history.Δη))))
            # push!(plots,plot(history.iter,history.ΔΛ,legend=false,title="Maximum Change in Lambda",ylim=(0,maximum(history.ΔΛ))))
            # push!(plots,plot(history.iter,history.ΔΦ,legend=false,title="Maximum Change in Phi",ylim=(0,maximum(history.ΔΦ))))
            display(plot(plots...,layout=(length(plots),1),size=(600,1200)))
        end

        # lf_old = copy(lf)
    end
    if iter > history.iter[end]
        Ŷ = getPredictions(td,predictX(lf))
        push!(history,[ iter poissonDeviance(Y,Ŷ) R2(Ỹ,Ŷ) ])
    end
    return nothing
end
