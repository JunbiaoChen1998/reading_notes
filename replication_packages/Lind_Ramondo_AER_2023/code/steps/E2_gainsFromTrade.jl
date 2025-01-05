function computeGainsFromTrade(dirs::ProjDirectories,data::ProjData)


    σ,σse,dfΛ,dfΦ,K,dfCES,dfSGM = loadEstimates(dirs)

	J = data.J
	N = data.N
	T = data.T
	K = length(σ)
	S = length(data.trade.sitcCodes)
	countries = data.countries
	years = data.years

	# load data
	finalUseShares = load(dirs.data*"finalUseShares.csv") |> DataFrame
	intermediateUseShares = load(dirs.data*"intermediateUseShares.csv") |> DataFrame

	# create arrays containing expenditure data at the WIOT level
	sort!(finalUseShares,reverse([:j,:d,:t]))
	μjdt = reshape(finalUseShares.value,(J,N,T))
	sort!(intermediateUseShares,reverse([:i,:j,:d,:t]))
	γijdt = reshape(intermediateUseShares.value,(J,J,N,T))

    # create arrays containing expenditure data at the WIOT level
    df = copy(data.sec)
    rename!(df,:value=>:Xjodt,:sec=>:j)
    sort!(df,reverse([:j,:o,:d,:t]))
    Xjodt = reshape(df.Xjodt,(J,N,N,T))
    πjodt = Xjodt ./ sum(Xjodt,dims=(1,2))
    πWjodt = Xjodt ./ sum(Xjodt,dims=(2))
    πBjdt = dropdims(sum(Xjodt,dims=(2)) ./ sum(Xjodt,dims=(1,2)),dims=2)
    maximum(abs,πjodt .- πWjodt .* reshape(πBjdt,(J,1,N,T)))

    # create arrays containing expenditure estimates at the factor level
    sort!(dfΦ,reverse([:k,:o,:d,:t]))
    Xkodt = reshape(dfΦ.Xkodt,(K,N,N,T))
    πkodt = Xkodt ./ sum(Xkodt,dims=(1,2))
    tkodt = reshape(dfΦ.tkodt,(K,N,N,T))
    πWkodt = Xkodt ./ sum(Xkodt,dims=(2))
    πBkdt = dropdims(sum(Xkodt,dims=(2)) ./ sum(Xkodt,dims=(1,2)),dims=2)
    maximum(abs,πkodt .- πWkodt .* reshape(πBkdt,(K,1,N,T)))

    # aggregate shares
    πodt = dropdims(sum(πkodt,dims=1),dims=1)
    πddt = [ πodt[d,d,t] for d=1:N, t=1:T ]

    # estimates
    θ = minimum(σ)
    ρk = 1 .- θ./σ
    π̃kodt = (πkodt.^(1 .- reshape(ρk,(K,1,1,1)))) .* (sum(πkodt,dims=2).^reshape(ρk,(K,1,1,1)))
    maximum(abs,[sum(π̃kodt[:,d,d,t]) for d=1:N, t=1:T ] .- [ sum(πWkodt[:,d,d,t].^(1 .- ρk) .* πBkdt[:,d,t]) for d=1:N, t=1:T ])

    θCES = -dfCES.σ[1]
    θSGM = θ
    σSGM = -dfSGM.σ
    ρSGM = 1 .- θSGM./σSGM
    π̃jodt = (πjodt.^(1 .- reshape(ρSGM,(J,1,1,1)))) .* (sum(πjodt,dims=2).^reshape(ρSGM,(J,1,1,1)))

    #### self-trade and correlation-adjusted self-trade shares
    πodt = dropdims(sum(πkodt,dims=1),dims=1)
    π̃SGModt = dropdims(sum(π̃jodt,dims=1),dims=1)
    π̃LFModt = dropdims(sum(π̃kodt,dims=1),dims=1)

    kvals = 1:K
    tmp = [ sum(πWkodt[kvals,d,d,t] .^ (1 .- ρk[kvals]) .* πBkdt[kvals,d,t]) for d=1:N, t=1:T ]
    maximum(abs,[ π̃LFModt[d,d,t] for d=1:N, t=1:T ] .- tmp )

    # calculate gains from trade across specifications
    GTCESdt = [ πodt[d,d,t]^(-1/θ) for d=1:N, t=1:T ]
    GTSGMdt = [ π̃SGModt[d,d,t]^(-1/θ) for d=1:N, t=1:T ]
    GTLFMdt = [ π̃LFModt[d,d,t]^(-1/θ) for d=1:N, t=1:T ]

    Ξdt = [ inv(diagm(ones(J)) - γijdt[:,:,d,t]') for d=1:N, t=1:T]
    GTSGMIOdt_0 = [ prod( (πWjodt[:,d,d,t]') .^ ( - Ξdt[d,t] .* (μjdt[:,d,t] ./ (σSGM')) ) ) for d=1:N, t=1:T ]
    GTSGMIOdt = [ sum([ prod( ( πWjodt[:,d,d,t].^(1 ./ σSGM) ) .^ Ξdt[d,t][j,:] ) for j=1:J ] .^ θ .* πBjdt[:,d,t]).^(-1/θ) for d=1:N, t=1:T ]
    GTSGMdt_0 = [ prod(πWjodt[:,d,d,t].^(-πBjdt[:,d,t]./σSGM)) for d=1:N, t=1:T ]
    GTLFMdt_0 = [ prod(πWkodt[:,d,d,t].^(-πBkdt[:,d,t]./σ)) for d=1:N, t=1:T ]

    # factor and WIOT sector level shares
    πkddt = [ πkodt[k,d,d,t] for k=1:K, d=1:N, t=1:T]
    πjddt = [ πjodt[j,d,d,t] for j=1:J, d=1:N, t=1:T]


    fsz = 8
    psz = (600,500)



    σ̄LFMdt = [ sum(σ .* πBkdt[:,d,t]) for d=1:N, t=1:T ]
    ρ̄LFM = 1-θ/mean(σ)

    σ̄SGMdt = [ sum(σSGM .* πBjdt[:,d,t]) for d=1:N, t=1:T ]
    ρ̄SGM = 1-θ/mean(σSGM)

    scatter(σ̄LFMdt[:,end],σ̄SGMdt[:,end],label="",xlabel="Expenditure-Weighted Average of Factor-Level Elasticity",ylabel="Expenditure-Weighted Average of Sector-Level Elasticity",alpha=0,series_annotations=Plots.series_annotations(countries,Plots.font(fsz,:black)),size=psz,fontfamily=:Times)
    plot!(sort(σ̄LFMdt[:,end]),sort(σ̄LFMdt[:,end]),style=:solid,label="45 degree line",legend=:bottomright,color=:black)
    hline!([θCES],style=:dashdot,label="CES",color=:black)

    savefig(dirs.results*"paper/Figure2a.pdf")


    ρcutoff = .25
    ticks = [.0,.2,.4,.6,.8,1.]
    πddtLFM_select = [ sum(πBkdt[ρk .>= ρcutoff,d,t]) for d=1:N, t=1:T ]
    πddtSGM_select = [ sum(πBjdt[ρSGM .>= ρcutoff,d,t]) for d=1:N, t=1:T ]
    scatter(πddtLFM_select[:,end],πddtSGM_select[:,end],xticks=ticks,yticks=ticks,label="",xlabel="Expenditure Share in High Correlation Factors",ylabel="Expenditure Share in High Correlation Sectors",alpha=0,series_annotations=Plots.series_annotations(countries,Plots.font(fsz,:darkorange)),size=psz,fontfamily=:Times)
    plot!([0,1],[0,1],style=:solid,label="45 degree line",legend=:bottomright,color=:black)

    ρcutoff = .7
    πddtLFM_select = [ sum(πBkdt[ρk .>= ρcutoff,d,t]) for d=1:N, t=1:T ]
    πddtSGM_select = [ sum(πBjdt[ρSGM .>= ρcutoff,d,t]) for d=1:N, t=1:T ]
    scatter!(πddtLFM_select[:,end],πddtSGM_select[:,end],label="",alpha=0,series_annotations=Plots.series_annotations(countries,Plots.font(fsz,:blue)))

    ρcutoff = .86
    πddtLFM_select = [ sum(πBkdt[ρk .>= ρcutoff,d,t]) for d=1:N, t=1:T ]
    πddtSGM_select = [ sum(πBjdt[ρSGM .>= ρcutoff,d,t]) for d=1:N, t=1:T ]
    scatter!(πddtLFM_select[:,end],πddtSGM_select[:,end],label="",alpha=0,series_annotations=Plots.series_annotations(countries,Plots.font(fsz,:black)))

    savefig(dirs.results*"paper/Figure2b.pdf")


    # Figure 5
    scatter(πddt[:,end],log.(GTLFMdt[:,end]),ylim=(0,3.5),legend=:false,xlabel="Self-Trade Share",ylabel="Log Gains From Trade",alpha=0,series_annotations=Plots.series_annotations(countries,Plots.font(fsz,:black)),size=psz,fontfamily=:Times)
    p = countries .!= "SVN"
    scatter!(πddt[p,end],log.(GTSGMIOdt[p,end]),ylim=(0,3.5),alpha=0,series_annotations=Plots.series_annotations(countries,Plots.font(fsz,:blue)),size=psz,fontfamily=:Times)
    scatter!(πddt[:,end],log.(GTSGMdt[:,end]),ylim=(0,3.5),alpha=0,series_annotations=Plots.series_annotations(countries,Plots.font(fsz,:darkorange)),size=psz,fontfamily=:Times)
    p = sortperm(πddt[:,end])
    plot!(πddt[p,end],log.(πddt[p,end].^(-1/θCES)),label="",color=:black,style=:dash)

    savefig(dirs.results*"paper/Figure5a.pdf")


    p = πddt[:,end] .> .5
    any(countries[p] .== "SVN")
    scatter(πddt[p,end],log.(GTLFMdt[p,end]),legend=:false,xlabel="Self-Trade Share",ylabel="Log Gains From Trade",alpha=0,series_annotations=Plots.series_annotations(countries[p],Plots.font(fsz,:black)),size=psz,fontfamily=:Times)
    scatter!(πddt[p,end],log.(GTSGMIOdt[p,end]),alpha=0,series_annotations=Plots.series_annotations(countries[p],Plots.font(fsz,:blue)),size=psz,fontfamily=:Times)
    scatter!(πddt[p,end],log.(GTSGMdt[p,end]),alpha=0,series_annotations=Plots.series_annotations(countries[p],Plots.font(fsz,:darkorange)),size=psz,fontfamily=:Times)
    x = sort(πddt[p,end])
    plot!(x,log.(x.^(-1/θCES)),label="",color=:black,style=:dash)

    savefig(dirs.results*"paper/Figure5b.pdf")

    tdf = DataFrame(hcat(countries,round.(πddt[:,end],digits=3)),[Symbol("Country Code"),Symbol("Domestic Share")])
    tdf[:,Symbol("CES")] = round.(πddt[:,end].^(-1/θCES),digits=3)
    tdf[:,Symbol("SGM")] = round.(GTSGMdt[:,end],digits=3)
    tdf[:,Symbol("SGM + IO")] = round.(GTSGMIOdt[:,end],digits=3)
    tdf[:,Symbol("LFM")] = round.(GTLFMdt[:,end],digits=3)

	save(dirs.results*"onlineAppendix/TableO8.csv",tdf)
	io = open(dirs.results*"onlineAppendix/TableO8.tex","w")
	write(io,latexify(tdf,env=:table,latex=false))
	close(io)

    GTLFM_sameρ = [ sum(πWkodt[:,d,d,t] .^ (1 .- ρ̄LFM ) .* πBkdt[:,d,t])^(-1/θ) for d=1:N, t=1:T ]
    GTSGM_sameρ = [ sum(πWjodt[:,d,d,t] .^ (1 .- ρ̄SGM) .* πBjdt[:,d,t])^(-1/θ) for d=1:N, t=1:T ]

    GTLFM_sameW = [ sum(πddt[d,t] .^ (1 .- ρk ) .* πBkdt[:,d,t])^(-1/θ) for d=1:N, t=1:T ]
    GTSGM_sameW = [ sum(πddt[d,t] .^ (1 .- ρSGM) .* πBjdt[:,d,t])^(-1/θ) for d=1:N, t=1:T ]

    GTLFM_sameB = [ sum((πkodt[:,d,d,t].*K) .^ (1 .- ρk ) ./ K)^(-1/θ) for d=1:N, t=1:T ]
    GTSGM_sameB = [ sum((πjodt[:,d,d,t].*J) .^ (1 .- ρSGM) ./ J)^(-1/θ) for d=1:N, t=1:T ]

    GTLFM_sameρB = [ sum((πkodt[:,d,d,t].*K) .^ (1 .- ρ̄LFM ) ./ K)^(-1/θ) for d=1:N, t=1:T ]
    GTSGM_sameρB = [ sum((πjodt[:,d,d,t].*J) .^ (1 .- ρ̄SGM) ./ J)^(-1/θ) for d=1:N, t=1:T ]

    GTLFM_total = [ (πddt[d,t] .^ (1 .- ρ̄LFM) )^(-1/θ) for d=1:N, t=1:T ]
    GTSGM_total = [ (πddt[d,t] .^ (1 .- ρ̄SGM) )^(-1/θ) for d=1:N, t=1:T ]


    scatter(πddt[:,end],log.(GTLFMdt[:,end]),fg_legend=:transparent,label="",xlabel="Self-Trade Share",ylabel="Log Gains From Trade",alpha=0,series_annotations=Plots.series_annotations(countries,Plots.font(fsz,:black)),size=psz,fontfamily=:Times,ylim=(0,3.5))
    scatter!(πddt[:,end],log.(GTLFM_sameW[:,end]),label="Common Within-Factor Shares",color=:black,marker=:star)
    # scatter!(πddt[:,end],log.(GTLFM_sameB[:,end]),label="Common Between-Factor Shares",color=:darkgrey,marker=:circle)
    scatter!(πddt[:,end],log.(GTLFM_sameρ[:,end]),label="Common Correlation Parameters",color=:gray,marker=:circle)
    # scatter!(πddt[:,end],log.(GTLFM_sameρB[:,end]),label="Common Between-Factor Shares and Correlation Parameters",color=:blue,marker=:circle)
    scatter!(πddt[:,end],log.(GTLFM_total[:,end]),label="Both",color=:white,marker=:circle)
    p = sortperm(πddt[:,end])
    plot!(πddt[p,end],log.(πddt[p,end].^(-1/θCES)),label="",color=:black,style=:dash)

    savefig(dirs.results*"paper/Figure5c.pdf")


    scatter(πddt[:,end],log.(GTSGMdt[:,end]),fg_legend=:transparent,label="",xlabel="Self-Trade Share",ylabel="Log Gains From Trade",alpha=0,series_annotations=Plots.series_annotations(countries,Plots.font(fsz,:black)),size=psz,fontfamily=:Times,ylim=(0,3.5))
    scatter!(πddt[:,end],log.(GTSGM_sameW[:,end]),label="Common Within-Sector Shares",color=:black,marker=:star)
    # scatter!(πddt[:,end],log.(GTSGM_sameB[:,end]),label="Common Between-Factor Shares",color=:darkgrey,marker=:circle)
    scatter!(πddt[:,end],log.(GTSGM_sameρ[:,end]),label="Common Correlation Parameters",color=:gray,marker=:circle)
    scatter!(πddt[:,end],log.(GTSGM_total[:,end]),label="Both",color=:white,marker=:circle)
    p = sortperm(πddt[:,end])
    plot!(πddt[p,end],log.(πddt[p,end].^(-1/θCES)),label="",color=:black,style=:dash)

    savefig(dirs.results*"paper/Figure5d.pdf")

    return nothing
end
