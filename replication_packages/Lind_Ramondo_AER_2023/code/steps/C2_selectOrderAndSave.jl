
### 1. Order Selection and Fit
function getBestRep(K,dir::String)
    # dir = dirs.results*"estimates/LFM/"
    pds = []
    for rep in readdir(dir*"$K")
        @show K, rep
        historyDf = load(dir*"$K/$rep/history.csv") |> DataFrame
        pd = historyDf.pd[end]
        ismissing(pd) || isnan(pd)  ? nothing : push!(pds,(rep,pd))
    end
    pdvals = map(x->x[2],pds)
    return pds[findfirst(pdvals .== minimum(pdvals))]
end

# table
# getPredictionsSec(td::TradeData,X::Matrix{Float64}) = vec(td.secMap'*X̂)
# getObservationsSec(td::TradeData) = vec(td.secValues)

function selectOrderAndSave(dirs::ProjDirectories,data::ProjData)
    Knames = vcat(string.(1:8),"14SGM","14")
    reps = map(K->getBestRep(K,dirs.results*"estimates/LFM/"),Knames)
    bestReps = DataFrame(k=Knames,rep=map(x->x[1],reps))
    save(dirs.results*"estimates/LFM/bestReps.csv",bestReps)

    selectKTable = DataFrame(K=Int64[],name=String[],pdsitc=Float64[],shareZeroLambda=Float64[],shareZeroPhi=Float64[],
                                    R2_sitc=Float64[],R2_sitc_w_odt=Float64[],
                                    R2_sec=Float64[],R2_sec_w_dt=Float64[],R2_sec_w_jdt=Float64[],R2_sec_w_odt=Float64[])

    tdShares,total = shares(data.trade)

    Ysitc = getObservations(data.trade)
    Nobs = length(Ysitc)
    Ysitcshares = getObservations(tdShares)
    poissonNormalization(x) = mean(x)/var(x)
    scaling = poissonNormalization(Ysitcshares)
    poissonNormalization(scaling.*Ysitcshares)

    dir = dirs.results*"estimates/LFM/"
    S,M = size(tdShares.sitcValues)
    data.J != size(tdShares.secValues,1) ? error("Dimension mismatch.") : nothing



    specs = vcat(map(k->(k,"LFM","$(bestReps.rep[k])"),1:8),(14,"SGM","$(bestReps.rep[9])"),(14,"LFM","$(bestReps.rep[10])"))
    for (K,name,file) in specs

        # (K,name,file) = specs[end]
        subdir = dir*"$K"
        name == "SGM" ? subdir *= "SGM/" : subdir *= "/"
        @show K,name,file

        lf = loadLatentFactors(subdir,file,data.trade)
        X̂sitcshares = predictX(lf)
        Ŷsitcshares = getPredictions(data.trade,X̂sitcshares)

        pdsitc = poissonDeviance(scaling.*Ysitcshares,scaling.*Ŷsitcshares)
        shareZeroΛ = mean(hcat(lf.Λ...) .<= eps())
        shareZeroΦ = mean(hcat(lf.Φ...) .<= eps())


        X̂sitc = X̂sitcshares.*total
        Ŷsitc = getPredictions(data.trade,X̂sitc)

        R2_sitc = R2(Ysitc,Ŷsitc)
        R2_sitc_w_odt = R2(Ysitcshares,Ŷsitcshares)


        sitcPredictions = DataFrame[]
        for (m,(o,d,t)) in enumerate(data.trade.mValues)
            sdf = DataFrame(X̂iodt = X̂sitc[:,m])
            sdf.sitc = 1:S
            sdf.o = o
            sdf.d = d
            sdf.t = t
            push!(sitcPredictions,sdf)
        end
        sitcPredictions = vcat(sitcPredictions...)

        sitcPredictions.sec = map(sitc->data.trade.secInd[sitc],sitcPredictions.sitc)
        secPredictions = by(sitcPredictions,[:sec,:o,:d,:t],Ŷjodt = :X̂iodt => sum)
        secPredictions = join(data.sec[:,Not(:residual)],secPredictions,on=[:sec,:o,:d,:t],kind=:left)
        secPredictions.value /= 1000
        rename!(secPredictions,:value => :Yjodt)

        secPredictions = join(secPredictions,by(secPredictions,[:o,:d,:t],Yodt = :Yjodt => sum,Ŷodt = :Ŷjodt => sum),on=[:o,:d,:t],kind=:left)
        secPredictions = join(secPredictions,by(secPredictions,[:sec,:d,:t],Yjdt = :Yjodt => sum,Ŷjdt = :Ŷjodt => sum),on=[:sec,:d,:t],kind=:left)
        secPredictions = join(secPredictions,by(secPredictions,[:d,:t],Ydt = :Yjodt => sum,Ŷdt = :Ŷjodt => sum),on=[:d,:t],kind=:left)

        secPredictions.πjodt = secPredictions.Yjodt ./ secPredictions.Ydt
        secPredictions.π̂jodt = secPredictions.Ŷjodt ./ secPredictions.Ŷdt

        Y = secPredictions.Yjodt
        Ŷ = secPredictions.Ŷjodt
        R2_sec = R2(Y,Ŷ)

        Y = secPredictions.πjodt
        Ŷ = secPredictions.π̂jodt
        R2_sec_w_dt = R2(Y,Ŷ)

        Y = secPredictions.Yjodt ./ secPredictions.Yodt
        Ŷ = secPredictions.Ŷjodt ./ secPredictions.Ŷodt
        R2_sec_w_odt = R2(Y,Ŷ)

        Y = secPredictions.Yjodt ./ secPredictions.Yjdt
        Ŷ = secPredictions.Ŷjodt ./ secPredictions.Ŷjdt
        R2_sec_w_jdt = R2(Y,Ŷ)

        push!(selectKTable,[K name pdsitc shareZeroΛ shareZeroΦ R2_sitc R2_sitc_w_odt R2_sec R2_sec_w_dt R2_sec_w_jdt R2_sec_w_odt])
    end
    selectKTableRaw = copy(selectKTable)




    selectKTable = copy(selectKTableRaw)
    selectKTable.dof = vcat([ (S+M)*K for K in 1:8 ],S - 14 + 14 + M*14,(S+M)*14)

    Δpd = selectKTable.pdsitc[1:7] - selectKTable.pdsitc[2:8]
    Δpd = vcat(Δpd,missing)
    Δpd = vcat(Δpd,selectKTable.pdsitc[7]-selectKTable.pdsitc[10])
    Δpd = vcat(Δpd,selectKTable.pdsitc[9]-selectKTable.pdsitc[10])
    Δdof = selectKTable.dof[2:8] - selectKTable.dof[1:7]
    Δdof = vcat(Δdof,missing)
    Δdof = vcat(Δdof,selectKTable.dof[10]-selectKTable.dof[7])
    Δdof = vcat(Δdof,selectKTable.dof[10]-selectKTable.dof[9])

    push!(selectKTable,selectKTable[end,:])
    selectKTable.name[end] = "LFM vs SGM"
    selectKTable.LLR = vcat(missing,Δpd)
    selectKTable.Δdof = vcat(missing,Δdof)
    selectKTable.pval = vcat(missing,map(n->ismissing(Δpd[n]) ? missing : ccdf(Chisq(Δdof[n]),Δpd[n]),1:length(Δpd)))

    for v in names(selectKTable)
        if v != :name
            selectKTable[:,v] = round.(selectKTable[:,v],digits=3)
        end
    end
    selectKTable.pdsitc = Int64.(round.(selectKTable.pdsitc,digits=0))
    selectKTable.dof = Int64.(selectKTable.dof)
    selectKTable.null = null = [missing,map(x->"LFM $x",1:7)...,missing,"LFM 7","SGM 14"]


    header = [Symbol("$(selectKTable.name[n]) $(selectKTable.K[n])") for n=1:nrow(selectKTable) ]


    selectKTable = selectKTable[:,[:R2_sitc,:R2_sitc_w_odt,:R2_sec,:R2_sec_w_dt,:R2_sec_w_jdt,:R2_sec_w_odt,:pdsitc,:dof,:null,:LLR,:Δdof,:pval]]

    @show names(selectKTable)
    selectKTableOut = DataFrame(permutedims(Array{Any}(selectKTable)),header)
    selectKTableOut.name = string.(names(selectKTable))
    selectKTableOut = selectKTableOut[:,[12,1:11...]]
    for v in names(selectKTableOut)
        selectKTableOut[ismissing.(selectKTableOut[:,v]),v] .= ""
    end

    @show selectKTableOut

    save(dirs.results*"onlineAppendix/TableO2.csv",selectKTableOut)
    io = open(dirs.results*"onlineAppendix/TableO2.tex","w")
    write(io,latexify(selectKTableOut,env=:table,latex=false))
    close(io)

    table1 = selectKTableOut[vcat(1,2,7:12...),vcat(1:9)]

    @show table1

    save(dirs.results*"paper/Table1.csv",table1)
    io = open(dirs.results*"paper/Table1.tex","w")
    write(io,latexify(table1,env=:table,latex=false))
    close(io)


    ### 2. 14-SGM Estimates
    dir = dirs.results*"estimates/LFM/14SGM/"
    name = bestReps.rep[findfirst(bestReps.k .== "14SGM")]
    lf = loadLatentFactors(dir,name,tdShares)
    save(dirs.results*"estimates/LFM/elasticities14SGM.csv",DataFrame(k=1:14,σ=lf.η))

    dfΛ = DataFrame(sitc=String[],k=Int64[],Λ=Float64[])
    for s=1:S,k = 1:14
        push!(dfΛ,[tdShares.sitcCodes[s] k lf.Λ[k][s]])
    end
    save(dirs.results*"estimates/LFM/factorWeights14SGM.csv",dfΛ)


    ### 3. Estimates
    K = 7
    dir = dirs.results*"estimates/LFM/7/"
    name = bestReps.rep[findfirst(bestReps.k .== "7")]
    lf = loadLatentFactors(dir,name,tdShares)

    krank = sortperm(lf.η,rev=true)
    lf.η .= lf.η[krank]
    lf.W .= lf.W[krank]
    lf.Λ .= lf.Λ[krank]
    lf.Φ .= lf.Φ[krank]

    save(dirs.results*"estimates/LFM/elasticities.csv",DataFrame(k=1:K,σ=lf.η))

    dfΛ = DataFrame(sitc=String[],k=Int64[],Λ=Float64[])
    for s=1:S,k = 1:K
        push!(dfΛ,[tdShares.sitcCodes[s] k lf.Λ[k][s]])
    end
    save(dirs.results*"estimates/LFM/factorWeights.csv",dfΛ)

    ###### Compute latent expenditure

    T̃k = [ (data.trade.sitcTariffs.^-lf.η[k])'*lf.Λ[k] for k=1:K ]
    tkodt = permutedims(hcat([ T̃k[k].^(-1/lf.η[k]) for k=1:K ]...))
    Xkodt =  permutedims(hcat([ T̃k[k].*lf.Φ[k] for k=1:K]...))
    Xkodt ./= sum(Xkodt,dims=1)
    Xkodt = Xkodt.*sum(data.trade.secValues,dims=1)

    df = DataFrame(k=Int64[],o=String[],d=String[],t=Int64[],Xkodt=Float64[],tkodt=Float64[])
    for k=1:K, m=1:M
        o,d,t = data.trade.mValues[m]
        push!(df,[k o d t Xkodt[k,m] tkodt[k,m]])
    end
    save(dirs.results*"estimates/LFM/factors.csv",df)

    return nothing
end
