
function computeResults(dirs::ProjDirectories,data::ProjData)

    # load estimates and labels
    σ,σse,dfΛ,dfΦ,K,dfCES,dfSGM = loadEstimates(dirs)
    #
    # df = copy(data.sec)
    # df = rankCountries(dirs,df)

    df = copy(dfΦ)

	J = data.J
	N = data.N
	T = data.T
	K = length(σ)
	S = length(data.trade.sitcCodes)

	sitcCodes = load(dirs.data*"sitcCodes.csv") |> DataFrame
	sitcToSec = load(dirs.data*"sitcToSec.csv") |> DataFrame
    dfΛ = join(dfΛ,sitcToSec,on=[:sitc])
    secNames = data.secNames
    secLabels = secNames.name
    dfΛ = join(dfΛ,secNames,on=[:sec])

	sitcDescriptions,firstCodes = loadSitcDescriptions(dirs,data)


	df = join(df,by(df,[:d,:t],Xdt = :Xkodt => sum),on=[:d,:t],kind=:left)
	df.πkodt = df.Xkodt ./ df.Xdt

    sort!(df,reverse([:k,:o,:d,:t]))
    Xkodt = reshape(df.Xkodt,(K,N,N,T))
    πkodt = reshape(df.πkodt,(K,N,N,T))
    tkodt = reshape(df.tkodt,(K,N,N,T))

    @show data.trade.sitcCodes
    sort!(dfΛ,reverse([:k,:sitc]))
    Λ = permutedims(reshape(dfΛ.Λ,(K,S)))
    sum(Λ,dims=1)


    for k in 1:K
        @show k,σ[k]
        srank = sortperm(Λ[1:773,k],rev=true)
        sitcNames = join(sitcCodes,rename(sitcDescriptions,:code=>:sitc),on=[:sitc],kind=:left)
        @show first(sitcNames[srank,:],10)
    end

    dfΛ.sitc2 = map(x->x[1:2],dfΛ.sitc)
    dfΛ = by(dfΛ,[:sitc2,:k,:sec,:name],Λ = :Λ => sum)
    dfΛ.sitc1 = map(x->"$(x[1])",dfΛ.sitc2)
    dfΛ = join(dfΛ,rename(sitcDescriptions,:code=>:sitc1,:description=>:d1),on=[:sitc1],kind=:left)
    dfΛ = join(dfΛ,rename(sitcDescriptions,:code=>:sitc2,:description=>:d2),on=[:sitc2],kind=:left)
    for k in 1:K
        @show k,σ[k]
        sdf = dfΛ[dfΛ.k .== k,[:Λ,:name,:d2]]
        sdf = sort(sdf,[:Λ],rev=true)
        @show first(sdf,10)
    end



    # Compute ρ
    θ = minimum(σ)
    ρk = min.(max.(0,1 .- θ./σ),1)


    dfSec = copy(data.sec)
    rename!(dfSec,:value=>:Xjodt,:sec=>:j)
    dfSec = join(dfSec,by(dfSec,[:d,:t],Xdt = :Xjodt => sum),on=[:d,:t],kind=:left)
    dfSec.πjodt = dfSec.Xjodt ./ dfSec.Xdt

    sort!(dfSec,reverse([:j,:o,:d,:t]))
    Xjodt = reshape(dfSec.Xjodt,(J,N,N,T))
    πjodt = reshape(dfSec.πjodt,(J,N,N,T))

    θSGM = copy(θ) #minimum(-dfSGM.σ[Not(2)])
    ρSGM = -dfSGM.σ
    # ρSGM[2] = θSGM
    ρSGM = 1 .- θSGM./ρSGM

    π̃jodt = (πjodt.^(1 .- reshape(ρSGM,(J,1,1,1)))) .* (sum(πjodt,dims=2).^reshape(ρSGM,(J,1,1,1)))
    ωjodt = π̃jodt ./ sum(π̃jodt,dims=1)


    # properties of infered expenditure
    df.trade = df.o .!= df.d

    dfAgg = rename(by(df[.!df.trade,:],[:k,:d,:t],STknt = :Xkodt => sum),:d=>:n)
    dfAgg = join(dfAgg,rename(by(df[df.trade,:],[:k,:o,:t],EXknt = :Xkodt => sum),:o=>:n),on=[:k,:n,:t])
    dfAgg = join(dfAgg,rename(by(df[df.trade,:],[:k,:d,:t],IMknt = :Xkodt => sum),:d=>:n),on=[:k,:n,:t])
    dfAgg = by(dfAgg,[:k],STk = :STknt => sum,EXk = :EXknt => sum,IMk = :IMknt => sum)
    @show DataFrame(k=1:K,  Xk_X = (dfAgg.STk + dfAgg.EXk)./sum((dfAgg.STk + dfAgg.EXk)),
                            STk_Xk = dfAgg.STk ./ (dfAgg.STk + dfAgg.EXk),
                            STk_ST = dfAgg.STk./sum(dfAgg.STk),
                            EXk_EX = dfAgg.EXk./sum(dfAgg.EXk))

    dfAgg.Xk_X = (dfAgg.STk + dfAgg.EXk)./sum((dfAgg.STk + dfAgg.EXk))
    dfAgg.STk_Xk = dfAgg.STk ./ (dfAgg.STk + dfAgg.EXk)
    dfAgg.STk_ST = dfAgg.STk./sum(dfAgg.STk)
    dfAgg.EXk_EX = dfAgg.EXk./sum(dfAgg.EXk)
    dfAgg0 = copy(dfAgg)


    dfAgg = rename(by(df[.!df.trade,:],[:k,:d,:t],STknt = :Xkodt => sum),:d=>:n)
    dfAgg = join(dfAgg,rename(by(df[df.trade,:],[:k,:o,:t],EXknt = :Xkodt => sum),:o=>:n),on=[:k,:n,:t])
    dfAgg = join(dfAgg,rename(by(df[df.trade,:],[:k,:d,:t],IMknt = :Xkodt => sum),:d=>:n),on=[:k,:n,:t])
    dfAgg = by(dfAgg,[:k,:t],STkt = :STknt => sum,EXkt = :EXknt => sum)
    dfAgg.Xkt = dfAgg.STkt .+ dfAgg.EXkt
    dfAgg = join(dfAgg,by(dfAgg,[:t],STt = :STkt => sum,EXt = :EXkt => sum,Xt = :Xkt => sum),on=[:t])

    dfAgg.STkt_Xkt = dfAgg.STkt ./ dfAgg.Xkt
    dfAgg.EXkt_Xt = dfAgg.EXkt ./ dfAgg.Xt
    dfAgg.STkt_Xkt = dfAgg.STkt ./ dfAgg.Xkt
    dfAgg.Xkt_Xt = dfAgg.Xkt ./ dfAgg.Xt
    dfAgg.STkt_STt = dfAgg.STkt ./ dfAgg.STt
    dfAgg.EXkt_EXt = dfAgg.EXkt ./ dfAgg.EXt

    vars = [:Xkt_Xt,:STkt_Xkt,:STkt_STt,:EXkt_EXt]
    plots = []
    for var in vars, k=1:K
        push!(plots,plot(data.years,dfAgg[dfAgg.k .== k,var],legend=false,title="Factor $k: $var",ylim=(0,1)))
    end
    plot(plots...,layout=(length(vars),K),size=(2000,300*length(vars)))

    dfAgg = rename(by(df[.!df.trade,:],[:k,:d,:t],STknt = :Xkodt => sum),:d=>:n)
    dfAgg = join(dfAgg,rename(by(df[df.trade,:],[:k,:o,:t],EXknt = :Xkodt => sum),:o=>:n),on=[:k,:n,:t])
    dfAgg = join(dfAgg,rename(by(df[df.trade,:],[:k,:d,:t],IMknt = :Xkodt => sum),:d=>:n),on=[:k,:n,:t])

    dfAgg.Yknt = dfAgg.STknt + dfAgg.EXknt
    dfAgg.Xknt = dfAgg.STknt + dfAgg.IMknt
    dfAgg.TBknt = dfAgg.EXknt - dfAgg.IMknt

    plots = []
    for k=1:K
        push!(plots,histogram(dfAgg.TBknt[dfAgg.k .== k],legend=false,title="Factor $k"))
    end
    plot(plots...,layout=(K,1),size=(400,1000))

    dfAgg = join(dfAgg,by(dfAgg,[:n,:t],Ynt = :Yknt => sum),on=[:n,:t])
    dfAgg = join(dfAgg,by(dfAgg,[:n,:t],Xnt = :Xknt => sum),on=[:n,:t])
    dfAgg = join(dfAgg,by(dfAgg,[:n,:t],STnt = :STknt => sum),on=[:n,:t])
    dfAgg = join(dfAgg,by(dfAgg,[:n,:t],EXnt = :EXknt => sum),on=[:n,:t])
    dfAgg = join(dfAgg,by(dfAgg,[:n,:t],IMnt = :IMknt => sum),on=[:n,:t])

    dfAgg.STknt_Xnt = dfAgg.STknt ./ dfAgg.Xnt
    dfAgg.IMknt_Xnt = dfAgg.IMknt ./ dfAgg.Xnt

    dfAgg.Xknt_Xnt = dfAgg.Xknt ./ dfAgg.Xnt
    dfAgg.STknt_Xknt = dfAgg.STknt ./ dfAgg.Xknt

    dfAgg.Yknt_Ynt = dfAgg.Yknt ./ dfAgg.Ynt
    dfAgg.STknt_Yknt = dfAgg.STknt ./ dfAgg.Yknt

    dfAgg.STknt_STnt = dfAgg.STknt ./ dfAgg.STnt
    dfAgg.EXknt_EXnt = dfAgg.EXknt ./ dfAgg.EXnt
    dfAgg.IMknt_IMnt = dfAgg.IMknt ./ dfAgg.IMnt

    vars = [:Xknt_Xnt,:STknt_Xknt,:Yknt_Ynt,:STknt_Yknt,:STknt_STnt,:EXknt_EXnt,:IMknt_IMnt]
    plots=[]
    for var in vars, k=1:K
            push!(plots,histogram(dfAgg[dfAgg.k .== k,var],title="Factor $k: $var"))
    end
    plot(plots...,layout=(length(vars),K),legend=false,size=(2000,300*length(vars)),xlim=(0,1))

    dfAgg = rename(by(df[.!df.trade,:],[:k,:d,:t],STknt = :Xkodt => sum),:d=>:n)
    dfAgg = join(dfAgg,rename(by(df[df.trade,:],[:k,:o,:t],EXknt = :Xkodt => sum),:o=>:n),on=[:k,:n,:t])
    dfAgg = join(dfAgg,rename(by(df[df.trade,:],[:k,:d,:t],IMknt = :Xkodt => sum),:d=>:n),on=[:k,:n,:t])
    dfAgg = join(dfAgg,by(dfAgg,[:k,:t],EXkt = :EXknt => sum),on=[:k,:t])
    dfAgg.EXknt_EXkt = dfAgg.EXknt ./ dfAgg.EXkt

    # create Table 3
    estimatesTableColumns = []
    for k=1:K
        column = [
            round(σ[k],digits=3),
            "("*"$(round(σse[k],digits=3))"*")",
            round(ρk[k],digits=3),
            # round(mean(Λ[:,k]),digits=3),
            # round(std(Λ[:,k]),digits=3),
            round(mean(Λ[:,k] .<= eps()),digits=3),
            # round(quantile(Λ[:,k],.5),digits=8),
            round(quantile(Λ[:,k],.9),digits=3),
            round(quantile(Λ[:,k],.99),digits=3),
            round(maximum(Λ[:,k]),digits=3),
            # round(sum(Λ[:,k].^2),digits=3),
            round.(sort(by(dfΛ[dfΛ.k .== k,:],[:sec,:name],Λ = :Λ => sum),:sec).Λ,digits=3)...,
            round.(dfAgg0.Xk_X[k],digits=3),
            round.(dfAgg0.STk_Xk[k],digits=3),
            round.(dfAgg0.STk_ST[k],digits=3),
            round.(dfAgg0.EXk_EX[k],digits=3),
            sort(dfAgg[(dfAgg.k .== k).&(dfAgg.t .== 1999),:],[:EXknt_EXkt],rev=true).n[1:3]...,
            sort(dfAgg[(dfAgg.k .== k).&(dfAgg.t .== 2007),:],[:EXknt_EXkt],rev=true).n[1:3]...]
        push!(estimatesTableColumns,column)
    end
    estimatesTable = DataFrame(hcat(estimatesTableColumns...),[Symbol("$k") for k=1:K])
    estimatesTableRows = [  "Sigma",
                            "SE",
                            "Rho",
                            # "Mean",
                            # "Standard Deviation",
                            "Zero Share",
                            # "50th Percentile",
                            "90th Percentile",
                            "99th Percentile",
                            "Maximum",
                            # "Herfindahl",
                            sort(unique(dfΛ[:,[:sec,:name]]),:sec).name...,
                            "Share of Expenditure",
                            "Self-Trade Share",
                            "Share of Total Self-Trade",
                            "Share of Total Exports",
                            "Rank 1 Exporter in 1999",
                            "Rank 2 Exporter in 1999",
                            "Rank 3 Exporter in 1999",
                            "Rank 1 Exporter in 2007",
                            "Rank 2 Exporter in 2007",
                            "Rank 3 Exporter in 2007"]
    estimatesTable.K = estimatesTableRows
    estimatesTable = estimatesTable[:,[8,1:7...]]

    @show estimatesTable

    save(dirs.results*"paper/Table2.csv",estimatesTable)
    io = open(dirs.results*"paper/Table2.tex","w")
    write(io,latexify(estimatesTable,env=:table,latex=false))
    close(io)

    Ntop=5
    top2Digit = DataFrame(k=String[],rank=Int64[],sitc2=String[],d2=String[],Weight=Float64[])
    for k = 1:K
        sdf = dfΛ[(dfΛ.k .== k).&(dfΛ.sitc2 .!= "R "),:]
        sdf = sort(sdf,[:Λ],rev=true)
        for rank = 1:Ntop
            push!(top2Digit,["F$k" rank sdf[rank,[:sitc2,:d2,:Λ]]... ])
        end
    end
    top2Digit.Weight = round.(top2Digit.Weight,digits=3)
    @show top2Digit

    save(dirs.results*"paper/Table3.csv",top2Digit)
    io = open(dirs.results*"paper/Table3.tex","w")
    write(io,latexify(top2Digit,env=:table,latex=false))
    close(io)

    dfΛbackup = copy(dfΛ)

    dfΛ = by(dfΛ,[:sitc2,:k,:d2],Λ = :Λ => sum)
    sort!(dfΛ,[:k,:sitc2,:d2])
    dfΛ.s = repeat(1:length(unique(dfΛ.sitc2)),K)


    @show dfΛ

    sitc2Descriptions = unique(dfΛ[:,[:s,:sitc2,:d2]])
    plots = []
    θ = minimum(σ)
    ρ = min.(max.(0,1 .- θ./σ),1)

    for k in 1:K
        sdf = dfΛ[dfΛ.k .== k,:]
        push!(plots,plot(sdf.s,sdf.Λ,ylabel="Factor $k",legend=false,xtickfont=font(7),xticks=(sitc2Descriptions.s,sitc2Descriptions.sitc2)))
    end
    plot!(plots[end],xticks=(sitc2Descriptions.s,sitc2Descriptions.d2),bottom_margin=125px)
    plot(plots...,layout=(K,1),size=(1000,1200),legend=false,linewidth=3,fontfamily=:Times,ylim=(0,.45),left_margin=20px,xrotation=75)
    savefig(dirs.results*"onlineAppendix/FigureO2.pdf")

    Λsitc2 = [dfΛ.Λ[findfirst((dfΛ.sitc2 .== sitc2).*(dfΛ.k .== k))] for sitc2 in sitc2Descriptions.sitc2, k=1:K]
    sum(Λsitc2,dims=1)

    dfΛ = copy(dfΛbackup)

    cosineSimilarity(x,y) = dot(x,y)/(norm(x)*norm(y))

    shareOfSectorsUsingFactors = [ mean(Λ[:,k] .> eps()) for k=1:K ]
    plt = bar(1:K,shareOfSectorsUsingFactors,legend=false,ylabel="Share of Sectors",xlabel="Factors",ylim=(0,1),fontfamily=:Times,size=(400,300),tick_direction=:out,grid=false)
    savefig(plt,dirs.results*"paper/Figure1a.pdf")

    shareOfFactorsUsedBySectors = [ sum(Λ[s,:] .> eps()) for s=1:S ]
    sRank = sortperm(shareOfFactorsUsedBySectors,rev=false)
    plt = plot((1:S)./S,shareOfFactorsUsedBySectors[sRank],legend=false,ylabel="Number of Factors",xlabel="Share of Sectors",xlim=(0,1),fontfamily=:Times,linewidth=3,size=(400,300),tick_direction=:out,grid=false)
    savefig(plt,dirs.results*"paper/Figure1b.pdf")


    simExt(x,y) = mean((x .> eps()).&(y .> eps()))      # fraction of jointly positive entries
    simInt(x,y) = dot(x,y)/(norm(x)*norm(y))            # cosine(angle between x and y)

    factorSimilarity = DataFrame(i=[],j=[],simExt=[],simInt=[])
    for i=1:K, j=1:K
        push!(factorSimilarity,[i j simExt(Λ[:,i],Λ[:,j]) simInt(Λ[:,i],Λ[:,j])])
    end
    sectorSimilarity = DataFrame(i=[],j=[],icode=[],jcode=[],simExt=[],simInt=[])
    for i=1:773, j=1:773
        push!(sectorSimilarity,[i j sitcCodes.sitc[i] sitcCodes.sitc[j] simExt(Λ[i,:],Λ[j,:]) simInt(Λ[i,:],Λ[j,:])])
    end
    factorOffDiagonal = factorSimilarity.i .> factorSimilarity.j
    sectorOffDiagonal = sectorSimilarity.i .> sectorSimilarity.j

    similarityMoments = []
    for x in [   factorSimilarity.simExt[factorOffDiagonal],
                    factorSimilarity.simInt[factorOffDiagonal],
                    sectorSimilarity.simExt[sectorOffDiagonal],
                    sectorSimilarity.simInt[sectorOffDiagonal]]
        push!(similarityMoments,round.([mean(x),std(x),minimum(x),quantile(x,.1),quantile(x,.5),quantile(x,.9),maximum(x)],digits=3))
    end
    similarityMomentNames = [
        "Mean",
        "Standard Deviation",
        "Minimum",
        "10th Percentile",
        "Median",
        "90th Percentile",
        "Maximum"
    ]
    similarityTable = DataFrame(hcat(similarityMomentNames,similarityMoments...),
                                [:Moment,Symbol("Factor Extensive"),Symbol("Factor Intensive"),
                                         Symbol("Sector Extensive"),Symbol("Sector Intensive")])

    @show similarityTable

    save(dirs.results*"onlineAppendix/TableO6.csv",similarityTable)
    io = open(dirs.results*"onlineAppendix/TableO6.tex","w")
    write(io,latexify(similarityTable,env=:table,latex=false))
    close(io)

    plt = histogram(factorSimilarity.simInt[factorOffDiagonal],legend=false,xlabel="Similarity Between Pairs of Factors",ylabel="Count",bins=10,xlim=(0,1),ylim=(0,6),fontfamily=:Times,size=(400,300),tick_direction=:out,grid=false)
    savefig(dirs.results*"paper/Figure1c.pdf")
    plt = histogram(sectorSimilarity.simInt[sectorOffDiagonal],legend=false,xlabel="Similarity Between Pairs of 4-Digit Sectors",ylabel="Count",xlim=(0,1),ylim=(0,30000),fontfamily=:Times,size=(400,300),tick_direction=:out,grid=false)
    savefig(dirs.results*"paper/Figure1d.pdf")

    # Figure for similarity of factor weights: what sectors load on the same factors?
    sectorSimilarity = DataFrame(i=[],j=[],icode=[],jcode=[],sim=[])
    for i=1:nrow(sitc2Descriptions), j=1:nrow(sitc2Descriptions)
        push!(sectorSimilarity,[i j sitc2Descriptions.s[i] sitc2Descriptions.s[j] cosineSimilarity(Λsitc2[i,:],Λsitc2[j,:])])
    end
    plt = scatter(sectorSimilarity.icode,sectorSimilarity.jcode,
                    msize=20,shape=:rect,markerstrokewidth=0,malpha=sectorSimilarity.sim,
                    xticks=(sitc2Descriptions.s,sitc2Descriptions.d2),xrotation=75,
                    yticks=(sitc2Descriptions.s,sitc2Descriptions.d2),
                    size=(5000,5000),legend=false,grid=false,
                    left_margin=0px,bottom_margin=550px,
                    xtickfontsize=30,ytickfontsize=30,fontfamily=:Times)
    savefig(dirs.results*"onlineAppendix/FigureO3.pdf")

	# Table O.5.
	secData = load(dirs.data*"secData.csv") |> @filter(_.t == 2007) |> DataFrame
	secData = join(by(secData |> DataFrame,[:sec],total = :value => sum),
		by(secData |> @filter(_.o == _.d) |> DataFrame,[:sec],self = :value => sum),on=[:sec],kind=:left)
	secData[:,Symbol("Self-Trade Share")] = round.(secData.self ./ secData.total,digits=2)
	secData[:,Symbol("Expenditure Share")] = round.(secData.total ./ sum(secData.total),digits = 4)

	estimatesSGM = load(dirs.results*"estimates/SGM.csv") |> DataFrame
	estimatesSGM.σ .*= -1
	estimatesSGM[:,Symbol("Sectoral elasticity")] = string.(round.(estimatesSGM.σ,digits=2)) .* " (" .* string.(round.(estimatesSGM.se,digits=3)) .* ")"

	tdf = load(dirs.data*"secNames.csv") |> DataFrame
	tdf = join(tdf,estimatesSGM[:,[:sec,Symbol("Sectoral elasticity")]],on=[:sec],kind=:left)
	tdf = join(tdf,secData[:,[:sec,Symbol("Self-Trade Share"),Symbol("Expenditure Share")]],on=[:sec],kind=:right)

	rename!(tdf,[:sec => :Code,:name => :Name])

	save(dirs.results*"onlineAppendix/TableO5.csv",tdf)
	io = open(dirs.results*"onlineAppendix/TableO5.tex","w")
	write(io,latexify(tdf,env=:table,latex=false))
	close(io)

	# Table O.7.
	out = vcat(["σ$K" for K=1:8],"θ","Average σ_k")
    for K = 1:8
		name = getBestRep(K,dirs.results*"estimates/LFM/")[1]
		ηvals = sort(DataFrame(load(dirs.results*"estimates/LFM/$K/$name/eta.csv")).η,rev=true)
		θ = minimum(ηvals)
		avgη = mean(ηvals)
		out = hcat(out,map(x->isnan(x) ? "" : string(round(x,digits=4)),vcat(ηvals,repeat([NaN],8-K),θ,avgη)))
    end

	out = DataFrame(out,vcat(Symbol(""),Symbol.(1:8)))

	save(dirs.results*"onlineAppendix/TableO7.csv",out)
	io = open(dirs.results*"onlineAppendix/TableO7.tex","w")
	write(io,latexify(out,env=:table,latex=false))
	close(io)

    return nothing
end
