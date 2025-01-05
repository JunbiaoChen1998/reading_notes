function computeCounterfactuals(dirs::ProjDirectories,data::ProjData)

    σ,σse,dfΛ,dfΦ,K,dfCES,dfSGM = loadEstimates(dirs)

	J = data.J
	N = data.N
	T = data.T
	K = length(σ)
	S = length(data.trade.sitcCodes)
	countries = data.countries
	years = data.years

    # LFM
    df = copy(dfΦ)
    df = join(df,by(df,[:d,:t],Xdt = :Xkodt => sum),on=[:d,:t],kind=:left)
    df.πkodt = df.Xkodt ./ df.Xdt

    θ = minimum(σ)
    ρ = 1 .- θ./σ

    sort!(df,reverse([:k,:o,:d,:t]))
    Xkodt = reshape(df.Xkodt,(K,N,N,T))
    πkodt = reshape(df.πkodt,(K,N,N,T))
    tkodt = reshape(df.tkodt,(K,N,N,T))
    π̃kodt = (πkodt.^(1 .- reshape(ρ,(K,1,1,1)))) .* (sum(πkodt,dims=2).^reshape(ρ,(K,1,1,1)))
    π̃LFModt = sum(π̃kodt,dims=1)
    ωkodt = π̃kodt ./ π̃LFModt

    # SGM
    dfSec = copy(data.sec)
    rename!(dfSec,:value=>:Xjodt,:sec=>:j)
    dfSec.Xjodt ./= 1000    # make sure expenditure is in the same units
    dfSec = join(dfSec,by(dfSec,[:d,:t],Xdt = :Xjodt => sum),on=[:d,:t],kind=:left)
    dfSec.πjodt = dfSec.Xjodt ./ dfSec.Xdt

    θSGM = copy(θ) #minimum(-dfSGM.σ)
    σSGM = -dfSGM.σ
    ρSGM = 1 .- θSGM./σSGM

    sort!(dfSec,reverse([:j,:o,:d,:t]))
    Xjodt = reshape(dfSec.Xjodt,(J,N,N,T))
    πjodt = reshape(dfSec.πjodt,(J,N,N,T))
    π̃jodt = (πjodt.^(1 .- reshape(ρSGM,(J,1,1,1)))) .* (sum(πjodt,dims=2).^reshape(ρSGM,(J,1,1,1)))
    π̃SGModt = sum(π̃jodt,dims=1)
    ωjodt = π̃jodt ./ π̃SGModt

    # check that aggregate expenditure matches
    maximum(abs,sum(Xkodt,dims=1) - sum(Xjodt,dims=1))

    # calculate aggregates
    Xdt = dropdims(sum(Xkodt,dims=(1,2)),dims=(1,2))
    maximum(abs,Xdt - dropdims(sum(Xjodt,dims=(1,2)),dims=(1,2)))
    Yot = dropdims(sum(Xkodt,dims=(1,3)),dims=(1,3))
    maximum(abs,Yot - dropdims(sum(Xjodt,dims=(1,3)),dims=(1,3)))
    TBdt = Yot - Xdt
    σ̃k = reshape(σ./θ,(K,1,1))
    σ̃j = reshape(σSGM./θSGM,(J,1,1))


    function expenditure(π̃od,Y,TB,σ̃k,θ,ω′kod,T̂od,Ŵ)
        Y′ = Ŵ.*Y
        xP′kod_Pd = (ω′kod .* π̃od .* T̂od .* reshape(Ŵ,(1,N,1)).^-θ ).^σ̃k
        xP′kd_Pd = sum( xP′kod_Pd , dims=2 ).^(1 ./ σ̃k)
        return ( xP′kod_Pd .* xP′kd_Pd.^(1 .- σ̃k) ./ sum(xP′kd_Pd,dims=1) ) .* reshape(Y′-TB,(1,1,N))
    end
    function changeInRealWage(π̃od,Y,TB,σ̃k,θ,ω′kod,T̂od,Ŵ)
        Y′ = Ŵ.*Y
        xP′kod_Pd = (ω′kod .* π̃od .* T̂od .* reshape(Ŵ,(1,N,1)).^-θ ).^σ̃k
        xP′kd_Pd = sum( xP′kod_Pd , dims=2 ).^(1 ./ σ̃k)
        P′d_Pd = sum( xP′kd_Pd, dims=1 ).^(-1 ./ θ)
        return Ŵ./dropdims(P′d_Pd,dims=(1,2))
    end
    function excessDemand(π̃od,Y,TB,σ̃k,θ,ω′kod,T̂od,Ŵ)
        Y′ = Ŵ.*Y
        xP′kod_Pd = (ω′kod .* π̃od .* T̂od .* reshape(Ŵ,(1,N,1)).^-θ ).^σ̃k
        xP′kd_Pd = sum( xP′kod_Pd , dims=2 ).^(1 ./ σ̃k)
        return dropdims(sum( (xP′kod_Pd .* xP′kd_Pd.^(1 .- σ̃k) ./ sum(xP′kd_Pd,dims=1)) .* reshape(Y′-TB,(1,1,N)),dims=(1,3)),dims=(1,3))./Y′ .- 1
    end
    maximum(abs,excessDemand(π̃LFModt[:,:,:,1],Yot[:,1],TBdt[:,1],σ̃k,θ,ωkodt[:,:,:,1],ones(1,N,N),ones(N)))


    function solve!(π̃od,Y,TB,σ̃k,θ,ω′kod,T̂od,Ŵ;tol=1e-14,maxit=5000,λ = .1,report=100)
        iter = 0
        err = 1.
        while err > tol && iter < maxit
            iter += 1
            resid = excessDemand(π̃od,Y,TB,σ̃k,θ,ω′kod,T̂od,Ŵ)
            Ŵ .*= max.(100*eps(),1 .+ λ .* resid)
            Ŵ .*= (sum(Y) ./ sum(Ŵ.*Y))
            err = maximum(abs,resid)
            if mod(iter,report) == 1
                @show iter,err
            end
        end
        if iter == maxit
            @warn "Maximum iterations reached! Final error = $err"
        end
        return nothing
    end

    function decomposedImpacts(Xkod,TBd,ωkod,σ,θ)
        K = length(σ)
        Yo = dropdims(sum(Xkod,dims=(1,3)),dims=(1,3))
        πkod = Xkod ./ sum(Xkod,dims=(1,2))
        ρ = 1 .- θ./σ
        π̃kod = (πkod.^(1 .- reshape(ρ,(K,1,1)))) .* (sum(πkod,dims=2).^reshape(ρ,(K,1,1)))
        π̃od = sum(π̃kod,dims=1)
        Πod = reshape(sum(Xkod,dims=1) ./ sum(Xkod,dims=(1,2)),(N,N))
        σ̃ = reshape(σ ./ θ,(K,1,1))

        M1 = - permutedims(hcat([gradient(lnŴ->excessDemand(π̃od,Yo,TBd,σ̃,θ,ωkod,ones(1,N,N),exp.(lnŴ))[o],zeros(N))[1] for o=1:N ]...))

        decomp = DataFrame(o=String[],d=String[],o′=String[],d′=String[],dlnW_lnP=Float64[])
        for d′=1:N
            M2 = permutedims(hcat([gradient(lnτod->excessDemand(π̃od,Yo,TBd,σ̃,θ,ωkod,reshape(exp.(-θ.*(lnτod*(1:N .==d′)')),(1,N,N)),ones(N))[o],zeros(N))[1] for o=1:N]...))
            for o′=1:N
                dlnW = zeros(N)
                dlnW[Not(o′)] = M1[Not(o′),Not(o′)]\M2[Not(o′),o′]
                dlnW_lnP_decomp = (I - Πod').*dlnW'  .- ((1:N .== d′)*(1:N .== o′)').*Πod[o′,d′]
                dlnW_lnP = (I - Πod')*dlnW - (1:N .== d′)*Πod[o′,d′]
                for d=1:N
                    for o=1:N
                        push!(decomp,[countries[o] countries[d] countries[o′] countries[d′] dlnW_lnP_decomp[d,o]])
                    end
                end
            end
        end

        return decomp
    end

    # Decompose local real wage effects
    t = findfirst(years .== 2007)
    impactsdf = rename(decomposedImpacts(Xkodt[:,:,:,t],TBdt[:,t],ωkodt[:,:,:,t],σ,θ),:dlnW_lnP => :dlnW_lnP_LFM)
    impactsdf = join(impactsdf,rename(decomposedImpacts(Xjodt[:,:,:,t],TBdt[:,t],ωjodt[:,:,:,t],σSGM,θSGM),:dlnW_lnP => :dlnW_lnP_SGM),on=[:o,:d,:o′,:d′],kind=:left)

    impactsdfAgg = impactsdf[(impactsdf.d .== impactsdf.d′).&(impactsdf.o′ .!= impactsdf.d′),:]
    impactsdfAgg.thirdParty = (impactsdfAgg.o .!= impactsdfAgg.d).&(impactsdfAgg.o .!= impactsdfAgg.o′)
    impactsdfAgg.self = (impactsdfAgg.o .== impactsdfAgg.d)
    impactsdfAgg = by(impactsdfAgg,[:d,:o′,:thirdParty,:self],dlnW_lnP_LFM=:dlnW_lnP_LFM=>sum,dlnW_lnP_SGM=:dlnW_lnP_SGM=>sum)

    # Create Figure O.5 and Table O.9
    moments = []

    name = "Domestic"
    ids = impactsdfAgg.self
    y = 100 .* (impactsdfAgg.dlnW_lnP_LFM[ids] ./ impactsdfAgg.dlnW_lnP_SGM[ids] .- 1)
    plt = density(y,label=name,xlim=(-100,100),xlabel="% Difference Real Wage Elasticity",ylabel="Density",color=:blue)
    push!(moments,[name,round.([mean(y) ,std(y), skewness(y), [ percentile(y,p) for p in [25,50,75,90] ]...],digits=4)...])
    name = "3rd Party"
    ids = impactsdfAgg.thirdParty
    y = 100 .* (impactsdfAgg.dlnW_lnP_LFM[ids] ./ impactsdfAgg.dlnW_lnP_SGM[ids] .- 1)
    density!(plt,y,label=name,color=:orange)
    push!(moments,[name,round.([mean(y) ,std(y), skewness(y), [ percentile(y,p) for p in [25,50,75,90] ]...],digits=4)...])

    name = "Total"
    subdf = by(impactsdfAgg,[:d,:o′],dlnW_lnP_LFM=:dlnW_lnP_LFM=>sum,dlnW_lnP_SGM=:dlnW_lnP_SGM=>sum)
    y =  100 .* (subdf.dlnW_lnP_LFM ./ subdf.dlnW_lnP_SGM .- 1)
    # y = y[(y .>= percentile(y,2.5)).&(y .<= percentile(y,97.5))]
    density!(plt,y,label=name,color=:purple)
    push!(moments,[name,round.([mean(y) ,std(y), skewness(y), [ percentile(y,p) for p in [25,50,75,90] ]...],digits=4)...])

    moments = hcat(["name","Mean", "Standard Deviation","Skewness",["$(p)th Percentile" for p in [25,50,75,90]]...],moments...)
    moments = DataFrame(moments[2:end,:],Symbol.(moments[1,:]))
    @show moments
    plot(plt,fontfamily=:Times,size=(350,250),linewidth=3,right_margin=10px,legend=:topleft,title="LFM vs SGM",titlefontsize=10)
    savefig(dirs.results*"onlineAppendix/FigureO5.pdf")

    save(dirs.results*"onlineAppendix/TableO9.csv",moments)
    io = open(dirs.results*"onlineAppendix/TableO9.tex","w")
    write(io,latexify(moments,env=:table,latex=false))
    close(io)

    # destination-specific verision of decomposedImpacts
    function decomposedImpacts(Xkod,TBd,ωkod,σ,θ,d::Int64)
        K = length(σ)
        Yo = dropdims(sum(Xkod,dims=(1,3)),dims=(1,3))
        πkod = Xkod ./ sum(Xkod,dims=(1,2))
        ρ = 1 .- θ./σ
        π̃kod = (πkod.^(1 .- reshape(ρ,(K,1,1)))) .* (sum(πkod,dims=2).^reshape(ρ,(K,1,1)))
        π̃od = sum(π̃kod,dims=1)
        Πod = reshape(sum(Xkod,dims=1) ./ sum(Xkod,dims=(1,2)),(N,N))
        σ̃ = reshape(σ ./ θ,(K,1,1))

        M1 = - permutedims(hcat([gradient(lnŴ->excessDemand(π̃od,Yo,TBd,σ̃,θ,ωkod,ones(1,N,N),exp.(lnŴ))[o],zeros(N))[1] for o=1:N ]...))

        decomp = DataFrame(o=String[],d=String[],o′=String[],d′=String[],dlnW_lnP=Float64[])
        d′= d
        M2 = permutedims(hcat([gradient(lnτod->excessDemand(π̃od,Yo,TBd,σ̃,θ,ωkod,reshape(exp.(-θ.*(lnτod*(1:N .==d′)')),(1,N,N)),ones(N))[o],zeros(N))[1] for o=1:N]...))
        for o′=1:N
            dlnW = zeros(N)
            dlnW[Not(o′)] = M1[Not(o′),Not(o′)]\M2[Not(o′),o′]
            dlnW_lnP_decomp = (I - Πod').*dlnW'  .- ((1:N .== d′)*(1:N .== o′)').*Πod[o′,d′]
            dlnW_lnP = (I - Πod')*dlnW - (1:N .== d′)*Πod[o′,d′]

            for o=1:N
                push!(decomp,[countries[o] countries[d] countries[o′] countries[d′] dlnW_lnP_decomp[d,o]])
            end
        end

        return decomp
    end

    #### USA raises tariffs on China
    o1 = "CHN"
    o2 = "USA"
    id1 = findfirst(countries .== o1)
    id2 = findfirst(countries .== o2)
    t = findfirst(years .== 2007)

    tariffResultsLFM = []
    tariffResultsSGM = []
    xResultsLFM = []
    xResultsSGM = []
    dResultsLFM = []
    dResultsSGM = []
    tvals = exp.(LinRange(log(1),log(101),100)).-1
    ŴLFM = ones(N)
    ŴSGM = ones(N)
    for Δt in tvals
        # global ŴLFM, ŴSGM
        ω′kod = ωkodt[:,:,:,t]
        ω′jod = ωjodt[:,:,:,t]

        T̂od = ones(1,N,N)
        T̂od[:,id1,id2] .= (1+Δt/100)^-θ

        @time solve!(π̃LFModt[:,:,:,t],Yot[:,t],TBdt[:,t],σ̃k,θ,ω′kod,T̂od,ŴLFM)
        Ŵ_P̂ = changeInRealWage(π̃LFModt[:,:,:,t],Yot[:,t],TBdt[:,t],σ̃k,θ,ω′kod,T̂od,ŴLFM)
        X = expenditure(π̃LFModt[:,:,:,t],Yot[:,t],TBdt[:,t],σ̃k,θ,ω′kod,T̂od,ŴLFM)
        impacts = decomposedImpacts(X,TBdt[:,t],ω′kod,σ,θ,id2)
        push!(tariffResultsLFM,Ŵ_P̂)
        push!(xResultsLFM,X)
        push!(dResultsLFM,impacts)

        T̂od = ones(1,N,N)
        T̂od[:,id1,id2] .= (1+Δt/100)^-θSGM

        @time solve!(π̃SGModt[:,:,:,t],Yot[:,t],TBdt[:,t],σ̃j,θSGM,ω′jod,T̂od,ŴSGM)
        Ŵ_P̂ = changeInRealWage(π̃SGModt[:,:,:,t],Yot[:,t],TBdt[:,t],σ̃j,θSGM,ω′jod,T̂od,ŴSGM)
        X = expenditure(π̃SGModt[:,:,:,t],Yot[:,t],TBdt[:,t],σ̃j,θSGM,ω′jod,T̂od,ŴSGM)
        impacts = decomposedImpacts(X,TBdt[:,t],ω′jod,σSGM,θSGM,id2)
        push!(tariffResultsSGM,Ŵ_P̂)
        push!(xResultsSGM,X)
        push!(dResultsSGM,impacts)
    end

    for (i,Δt) in enumerate(tvals)
        dResultsLFM[i].Δt = fill(Δt,nrow(dResultsLFM[i]))
        dResultsLFM[i].model = fill("LFM",nrow(dResultsLFM[i]))
    end
    for (i,Δt) in enumerate(tvals)
        dResultsSGM[i].Δt = fill(Δt,nrow(dResultsSGM[i]))
        dResultsSGM[i].model = fill("SGM",nrow(dResultsSGM[i]))
    end

    components = vcat(dResultsLFM...,dResultsSGM...)
    components.self = ( components.o .== components.d )
    components.thirdParty = ( components.o .!= components.d ).&( components.o .!= components.o′ )
    components.direct = ( components.o .== components.o′ )
    components = by(components,[:d,:o′,:self,:thirdParty,:direct,:Δt,:model],dlnW_lnP = :dlnW_lnP => sum)
    components = components[components.o′ .== o1,:]
    sort!(components,[:self,:thirdParty,:direct,:model,:Δt])

    function cumulativeComponent(subdf)
        dlnW_lnP = subdf.dlnW_lnP
        lnTariff = log.(1 .+ subdf.Δt ./ 100)
        return vcat(0,cumsum( ((dlnW_lnP[1:end-1] + dlnW_lnP[2:end])./2) .* (lnTariff[2:end] - lnTariff[1:end-1])))
    end

    sz=(400,350)

    subdf = components[components.model .== "LFM",:]
    pltLFM = plot(tvals,cumulativeComponent(subdf[subdf.self,:]),label="Domestic",xlabel="Percentage Point Change in Tariffs",ylabel="Change in Log Real Wage",title="LFM")
    plot!(pltLFM,tvals,cumulativeComponent(subdf[subdf.thirdParty,:]),label="Third Party")
    plot!(pltLFM,tvals,cumulativeComponent(subdf[subdf.direct,:]),label="Direct")
    # plot!(pltLFM,tvals,cumulativeComponent(by(subdf,[:Δt],dlnW_lnP = :dlnW_lnP => sum)),label="Total")
    plot!(pltLFM,tvals,map(x->log.(x[id2]),tariffResultsLFM),label="Total")
    plot(pltLFM,ylim=(-.0275,.0325),yticks=-.03:.01:.03,linewidth=3,size=sz,legend=:topleft,fontfamily=:Times)
    savefig(dirs.results*"paper/Figure6A_LFM.pdf")

    subdf = components[components.model .== "SGM",:]
    pltSGM = plot(tvals,cumulativeComponent(subdf[subdf.self,:]),label="Domestic",size=(425,375),legend=:bottomleft,xlabel="Percentage Point Change in Tariffs",ylabel="Change in Log Real Wage",title="SGM")
    plot!(pltSGM,tvals,cumulativeComponent(subdf[subdf.thirdParty,:]),label="Third Party")
    plot!(pltSGM,tvals,cumulativeComponent(subdf[subdf.direct,:]),label="Direct")
    plot!(pltSGM,tvals,map(x->log.(x[id2]),tariffResultsSGM),label="Total")
    plot(pltSGM,ylim=(-.0275,.0325),yticks=-.03:.01:.03,linewidth=3,size=sz,legend=:topleft,fontfamily=:Times)
    savefig(dirs.results*"paper/Figure6A_SGM.pdf")


    πResultsLFM = map(x->sum(x,dims=(1))./sum(x,dims=(1,2)),xResultsLFM)
    πResultsSGM = map(x->sum(x,dims=(1))./sum(x,dims=(1,2)),xResultsSGM)

    pltLFM = plot()
    y = map(x->x[1,id1,id2],πResultsLFM)
    plot!(pltLFM,tvals,y.-y[1],label="CHN")
    y = map(x->x[1,id2,id2],πResultsLFM)
    plot!(pltLFM,tvals,y.-y[1],label="USA")
    y = map(x->sum(x[1,Not([id1,id2]),id2]),πResultsLFM)
    plot!(pltLFM,tvals,y.-y[1],label="Third Parties")
    plot!(pltLFM,xlabel="Percentage Point Change in Tariffs",ylabel="Change in USA Expenditure Share",title="LFM")
    plot(pltLFM,ylim=(-.045,.03),yticks=-.04:.01:.03,linewidth=3,size=sz,legend=:bottomleft,fontfamily=:Times)
    savefig(dirs.results*"paper/Figure6B_LFM.pdf")

    pltSGM = plot()
    y = map(x->x[1,id1,id2],πResultsSGM)
    plot!(pltSGM,tvals,y.-y[1],label="CHN")
    y = map(x->x[1,id2,id2],πResultsSGM)
    plot!(pltSGM,tvals,y.-y[1],label="USA")
    y = map(x->sum(x[1,Not([id1,id2]),id2]),πResultsSGM)
    plot!(pltSGM,tvals,y.-y[1],label="Third Parties")
    plot!(pltSGM,xlabel="Percentage Point Change in Tariffs",ylabel="Change in USA Expenditure Share",title="SGM")
    plot(pltSGM,ylim=(-.045,.03),yticks=-.04:.01:.03,linewidth=3,size=sz,legend=:right,fontfamily=:Times)
    savefig(dirs.results*"paper/Figure6B_SGM.pdf")

    πResultsLFM_2 = map(x->x./sum(x,dims=(1,2)),xResultsLFM)
	πResultsLFM_2 = hcat(map(x->x[:,id1,id2],πResultsLFM_2)...)

    πResultsSGM_2 = map(x->x./sum(x,dims=(1,2)),xResultsSGM)
	πResultsSGM_2 = hcat(map(x->x[:,id1,id2],πResultsSGM_2)...)

	K = length(σ)
    pltLFM_2 = plot()
    for k=1:K
        plot!(pltLFM_2,tvals,πResultsLFM_2[k,:],label="F$k")
    end
    plot!(pltLFM_2,xlabel="Percentage Point Change in Tariffs",ylabel="Expenditure Share of $o2 on $o1",title="LFM")
    plot(pltLFM_2,ylim=(0,.03),yticks=-0:.01:.03,linewidth=3,size=sz,legend=:bottomleft,fontfamily=:Times)
    savefig(dirs.results*"paper/Figure6C_LFM.pdf")

    pltSGM_2 = plot()
    for j=1:J
        plot!(pltSGM_2,tvals,πResultsSGM_2[j,:],label="S$j")
    end
    plot!(pltSGM_2,xlabel="Percentage Point Change in Tariffs",ylabel="Expenditure Share of $o2 on $o1",title="SGM")
    plot(pltSGM_2,ylim=(0,.03),yticks=-0:.01:.03,linewidth=3,size=sz,legend=:topright,fontfamily=:Times)
    savefig(dirs.results*"paper/Figure6C_SGM.pdf")

    return nothing
end
