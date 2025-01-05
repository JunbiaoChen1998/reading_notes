function computeElasticities(dirs::ProjDirectories,data::ProjData)


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
	ωkodt = π̃kodt ./ sum(π̃kodt,dims=1)

	# SGM
	dfSec = copy(data.sec)
	rename!(dfSec,:value=>:Xjodt,:sec=>:j)
	dfSec = join(dfSec,by(dfSec,[:d,:t],Xdt = :Xjodt => sum),on=[:d,:t],kind=:left)
	dfSec.πjodt = dfSec.Xjodt ./ dfSec.Xdt

	θSGM = minimum(-dfSGM.σ)
	σSGM = -dfSGM.σ
	ρSGM = 1 .- θSGM./σSGM

	sort!(dfSec,reverse([:j,:o,:d,:t]))
	Xjodt = reshape(dfSec.Xjodt,(J,N,N,T))
	πjodt = reshape(dfSec.πjodt,(J,N,N,T))
	π̃jodt = (πjodt.^(1 .- reshape(ρSGM,(J,1,1,1)))) .* (sum(πjodt,dims=2).^reshape(ρSGM,(J,1,1,1)))
	ωjodt = π̃jodt ./ sum(π̃jodt,dims=1)

	# elasticities
	εLFM = [ sum((Xkodt[:,o,d,t]./sum(Xkodt[:,o,d,t])) .* (
	                    (σ .- θ) .* (Xkodt[:,o′,d,t] ./ sum(Xkodt[:,:,d,t],dims=2))
	                    .- (o == o′).* σ )) for o=1:N, o′=1:N, d=1:N, t=1:T ]

	εSGM = [ sum((Xjodt[:,o,d,t]./sum(Xjodt[:,o,d,t])) .* (
	                    (σSGM .- θSGM) .* (Xjodt[:,o′,d,t] ./ sum(Xjodt[:,:,d,t],dims=2))
	                    .- (o == o′).* σSGM )) for o=1:N, o′=1:N, d=1:N, t=1:T ]

	myColorGradient = cgrad([:darkorange, :orange, :white, :blue, :darkblue])

	subElast = DataFrame(o=String[],o′=String[],d=String[],t=Int64[],εLFM=Float64[],εSGM=Float64[])
	for o=1:N, o′=1:N, d=1:N, t=1:T
	    push!(subElast,[countries[o] countries[o′] countries[d] years[t] εLFM[o,o′,d,t] εSGM[o,o′,d,t]])
	end
	save(dirs.results*"estimates/substitutionElasticities.csv",subElast)

	# Figure 3a
	ownPriceId = subElast.o .== subElast.o′
	plot(sort(subElast.εSGM[ownPriceId]),sort(subElast.εLFM[ownPriceId]),xlabel="Sectoral gravity model",ylabel="Latent factor model",color=:black,legend=false,linewidth=3)
	bounds = [extrema(subElast.εSGM[ownPriceId])...]
	plot!(bounds,bounds,ls=:dash,color=:black)
	savefig(dirs.results*"paper/Figure3a.pdf")

	# Figure 3b
	plot(sort(subElast.εSGM[.!ownPriceId]),sort(subElast.εLFM[.!ownPriceId]),xlabel="Sectoral gravity model",ylabel="Latent factor model",color=:black,legend=false,linewidth=3)
	bounds = [extrema(subElast.εSGM[.!ownPriceId])...]
	plot!(bounds,bounds,ls=:dash,color=:black)
	savefig(dirs.results*"paper/Figure3b.pdf")

	# Figure 4a
	ids = (subElast.o .!= "CHN").&(subElast.o .!= "USA").&(subElast.o′ .== "CHN").&(subElast.d .== "USA").&(subElast.t .== 2007)
	y = subElast.o[ids]
	xSGM = subElast.εSGM[ids]
	xLFM = subElast.εLFM[ids]
	sp = sortperm(xSGM,rev=true)
	scatter(xSGM[sp],y[sp],label="Sectoral gravity model")
	scatter!(xLFM[sp],y[sp],label="Latent factor model")
	plot!(legend=:outerbottom)
	savefig(dirs.results*"paper/Figure4a.pdf")
	
	# Figure 4b
	ids = (subElast.o .== subElast.o′).&(subElast.d .== "USA").&(subElast.t .== 2007)
	y = subElast.o[ids]
	xSGM = subElast.εSGM[ids]
	xLFM = subElast.εLFM[ids]
	sp = sortperm(xSGM,rev=true)
	scatter(xSGM[sp],y[sp],label="Sectoral gravity model")
	scatter!(xLFM[sp],y[sp],label="Latent factor model")
	plot!(legend=:outerbottom)
	savefig(dirs.results*"paper/Figure4b.pdf")

	εmin = min(minimum(εSGM),minimum(εLFM))
	εmax = max(maximum(εSGM),maximum(εLFM))
	εrange = max(abs(εmin),abs(εmax))

	destination = "USA"
	year = 2007
	d = findfirst(countries .== destination)
	t = findfirst(years .== year)

	sim(x,y) = dot(x,y)/sqrt(dot(x,x)*dot(y,y))
	νkodt = ωkodt.^(1 ./(1 .- reshape(ρ,(K,1,1,1))))
	ws = sum(νkodt,dims=2).^(1 .- reshape(ρ,(K,1,1,1)))
	ws = (reshape(ρ,(K,1,1,1)) ./(1 .- reshape(ρ,(K,1,1,1)))) .* ws ./ sum(ws,dims=1)
	νkodt = sqrt.(ws) .* νkodt ./ sum(νkodt,dims=2)
	similarityLFM = [ sim(νkodt[:,o,d,t],νkodt[:,o′,d,t]) for o=1:N, o′=1:N, d=1:N, t=1:T ]
	simMat = mean(similarityLFM[:,:,:,:],dims=(3,4))[:,:,1,1]
	ranking = sortperm([ sum(1 ./ -log.(simMat[o,Not(o)])) for o=1:N ],rev=true)

	plt = heatmap(1:N,1:N,εLFM[ranking,ranking,d,t],yticks=(1:N,countries[ranking]),xticks=(1:N,countries[ranking]),ylabel="Origin",xlabel="Competitor",size=(500,400),color=myColorGradient,clim=(-εrange,εrange),fontfamily=:Times,xrotation=75)
	savefig(dirs.results*"onlineAppendix/FigureO4a.png")

	plt = heatmap(1:N,1:N,εSGM[ranking,ranking,d,t],yticks=(1:N,countries[ranking]),xticks=(1:N,countries[ranking]),ylabel="Origin",xlabel="Competitor",size=(500,400),color=myColorGradient,clim=(-εrange,εrange),fontfamily=:Times,xrotation=75)
	savefig(dirs.results*"onlineAppendix/FigureO4b.png")

	return nothing
end
