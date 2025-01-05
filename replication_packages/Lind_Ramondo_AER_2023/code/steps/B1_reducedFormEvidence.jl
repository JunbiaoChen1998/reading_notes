
function reducedFormEvidence(dirs::ProjDirectories,data::ProjData)

	df = copy(data.grav)
	J = data.J
	countries = data.countries
	years = data.years
	t0 = minimum(years)
	t1 = maximum(years)
	N = data.N
	T = data.T
	secNames = data.secNames
	# secLabels = secNames.name

	cgsigned = cgrad([:darkorange,:orange,:white, :blue, :darkblue])

	# 1. Load the gravity data
	# df,J,countries,years,N,T,t0,t1,secNames,secLabels = loadGravityData(dirs.data)
	df,countries = rankCountries(dirs,df)

	rename!(df,:sec=>:j)
	rename!(df,:secStr=>:jStr)
	rename!(df,:lnIncome_o => :lnYot,:lnIncome_d => :lnYdt)
	df = join(df,by(df,[:t],lnȲt = :lnYot => mean),on=[:t],kind=:left)
	df.lnŶot = df.lnYot - df.lnȲt
	df.lnŶdt = df.lnYdt - df.lnȲt
	df.absΔlnYodt = abs.(df.lnYot - df.lnYdt)
	rename!(df,:value => :Xjodt)
	rename!(df,:share => :πjodt)
	rename!(df,:lnTariff => :lntjodt)
	rename!(df,:tariff => :tjodt)

	# Create arrays for expenditure data
	sort!(df,reverse([:j,:o,:d,:t]))
	Xjodt = reshape(df.Xjodt,(J,N,N,T))
	πjodt = reshape(df.πjodt,(J,N,N,T))
	lntjodt = reshape(df.lntjodt,(J,N,N,T))

	πjod = dropdims(mean(πjodt,dims=4),dims=4)
	πod = reshape(sum(πjod,dims=1),(N,N))
	πjo = reshape(mean(πjod,dims=3),(J,N))

	#	Distance
	distanceDf = rename(unique(df[:,[:o,:d,:ΔlnDistw]]),:d=>:o′)
	sort!(distanceDf,reverse([:o,:o′]))
	reshape(distanceDf.o,(N,N))
	distance = reshape(distanceDf.ΔlnDistw,(N,N))
	maximum(abs,distance - distance')
	w = copy(distance)
	w -= Diagonal(w)
	indexDist = w*reshape(permutedims(lntjodt,(2,1,3,4)),(N,J*N*T))
	indexDist = permutedims(reshape(indexDist,(N,J,N,T)),(2,1,3,4))
	df.indexDist = vec(indexDist)

	#	GDP Per Capita Difference
	incomeDf = rename(unique(df[:,[:o,:d,:t,:lnYdt,:lnYot]]),:d=>:o′,:lnYdt=>:lnYot′)
	incomeDf.absΔlnYoo′t = abs.(incomeDf.lnYot - incomeDf.lnYot′)
	sort!(incomeDf,reverse([:o,:o′,:t]))
	gdfDifference = reshape(incomeDf.absΔlnYoo′t,(N,N,T))
	indexGDPDist = fill(NaN,J,N,N,T)
	for t = 1:T
		w = copy(gdfDifference[:,:,t])
		w -= Diagonal(w)
		tariffs = reshape(permutedims(lntjodt[:,:,:,t],(2,1,3)),(N,J*N))
		indexGDPDist[:,:,:,t] = permutedims(reshape(w*tariffs,(N,J,N)),(2,1,3))
	end
	any(isnan.(indexGDPDist))
	df.indexGDPDist = vec(indexGDPDist)

	#	Correlation of sector-origin shares (within destination)
	shares = reshape(πjod,(J*N,N))
	plt = heatmap(1:J*N,1:J*N,corrMat(shares),size=(500,400),xticks=(7.5:J:(7.5+J*(N-1)),countries),clim=(-1,1),yticks=(7.5:J:(7.5+J*(N-1)),countries),color=cgsigned,fontfamily=:Times,xrotation=75)
	savefig(dirs.results*"onlineAppendix/FigureO1a.png")

	#	Correlation of aggregate origin shares (within destination)
	sim = corrMat(πod)
	heatmap(1:N,1:N,sim,size=(500,400),xticks=(1:N,countries),yticks=(1:N,countries),color=cgsigned,clim=(-1,1),fontfamily=:Times,xrotation=75)
	savefig(dirs.results*"onlineAppendix/FigureO1b.png")
	w = sim
	w -= Diagonal(w)
	tariffs = reshape(permutedims(lntjodt,(2,1,3,4)),(N,J*N*T))
	df.indexOriginCorrViaGeo = vec(permutedims(reshape(w*tariffs,(N,J,N,T)),(2,1,3,4)))

	#	Correlation of origin shares (within sector)
	sim = corrMat( πjo' ./ sum(πjo',dims=1) )
	heatmap(1:N,1:N,sim,size=(500,400),xticks=(1:N,countries),yticks=(1:N,countries),color=cgsigned,clim=(-1,1),fontfamily=:Times,xrotation=75)
	savefig(dirs.results*"onlineAppendix/FigureO1c.png")
	w = sim
	w -= Diagonal(w)
	tariffs = reshape(permutedims(lntjodt,(2,1,3,4)),(N,J*N*T))
	df.indexOriginCorrViaSectors = vec(permutedims(reshape(w*tariffs,(N,J,N,T)),(2,1,3,4)))

	#	Correlation of sectoral shares (within origin)
	# add names and fix y axis
	sim = corrMat(πjo  ./ sum(πjo,dims=1))

	heatmap(1:J,1:J,sim,size=(500,400),xticks=(1:J),yticks=(1:J),color=cgsigned,clim=(-1,1),fontfamily=:Times)
	savefig(dirs.results*"onlineAppendix/FigureO1d.png")
	w = sim
	w -= Diagonal(w)
	tariffs = reshape(lntjodt,(J,N*N*T))
	df.indexSectorCorr = vec(reshape(w*tariffs,(J,N,N,T)))


	df.jod = string.(df.j) .* "_" .* string.(df.o) .* "_" .* string.(df.d)

	df = join(df,by(df,[:j,:d,:t],lnt̄jdt = :lntjodt => mean),on=[:j,:d,:t],kind=:left)
	df.lnt̂jodt = df.lntjodt - df.lnt̄jdt

	# Estimation
	σNames = ["Elasticity / Mean Elasticity"]
	controlNames = [ "lnDistw & lnt̂jodt","absΔlnYodt & lnt̂jodt",
	 						  "indexOriginCorrViaGeo", "indexOriginCorrViaSectors", "indexSectorCorr" ]
	varNames = vcat(σNames,controlNames)
	feNames = ["lntjodt","j & lntjodt","absΔlnYodt","j & o & d","j & o & t","j & d & t"]
	# vcovSpec = Vcov.cluster(:jod)

	vcovSpec = Vcov.cluster(:jod)
	rt = RegressionTable(varNames,feNames)

	# Specification: CES
	m = @formula πjodt ~ lntjodt + fe(j)&fe(o)&fe(d) + fe(j)&fe(o)&fe(t) + fe(j)&fe(d)&fe(t)
	@show result = nlreg(df, m, Poisson(), LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8)
	w = coefnames(result) .== "lntjodt"
	transform = ParameterTransformation(x->w'*x,x->w,"Elasticity / Mean Elasticity")
	adjustedResult=addTransformation(result,transform)
	@show rt = RegressionTable(rt,adjustedResult,[true,false,false,true,true,true],spec="(1)")

	@show dfCES = DataFrame(σ=coef(result),se=stderror(result))
	save(dirs.results*"estimates/CES.csv",dfCES)

	# Specification: SGM
	m = @formula πjodt ~ jStr&lntjodt + fe(j)&fe(o)&fe(d) + fe(j)&fe(o)&fe(t) + fe(j)&fe(d)&fe(t)
	@show result = nlreg(df, m, Poisson(), LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8)
	w = map(x->in(x,[ "jStr: $j & lntjodt" for j in unique(df.jStr) ]),coefnames(result))
	transform = ParameterTransformation(x->w'*x/J,x->w/J,"Elasticity / Mean Elasticity")
	adjustedResult=addTransformation(result,transform)

	β̂ = result.coef
	V̂ = result.vcov

	Nparam = length(β̂)
	R = hcat(zeros(Nparam-1),I) - hcat(I,zeros(Nparam-1))
	r = zeros(size(R)[1])
	χ2 = (R*β̂.-r)'*((R*V̂*R')\(R*β̂.-r))
	testVal = ccdf(Chisq(length(r)),χ2)

	@show rt = RegressionTable(rt,adjustedResult,[false,true,false,true,true,true],testExternal=["CES",χ2,testVal],spec="(2)")

	@show dfSGM = DataFrame(sec=map(x->parse(Int64,x[7:8]),coefnames(result)),σ=coef(result),se=stderror(result))
	save(dirs.results*"estimates/SGM.csv",dfSGM)

	# Specification: SGM + absΔlnYodt + (lnDistw + absΔlnYodt)&lnt̂jodt
	m = @formula πjodt ~ jStr&lntjodt  + absΔlnYodt + (lnDistw + absΔlnYodt)&lnt̂jodt + fe(j)&fe(o)&fe(d) + fe(j)&fe(o)&fe(t) + fe(j)&fe(d)&fe(t)
	@show result = nlreg(df, m, Poisson(), LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8)
	w = map(x->in(x,[ "jStr: $j & lntjodt" for j in unique(df.jStr) ]),coefnames(result))
	transform = ParameterTransformation(x->w'*x/J,x->w/J,"Elasticity / Mean Elasticity")
	adjustedResult=addTransformation(result,transform)
	@show rt = RegressionTable(rt,adjustedResult,[false,true,true,true,true,true],test=("SGM",["lnDistw & lnt̂jodt","absΔlnYodt & lnt̂jodt"]),spec="(3)")

	# Specification: SGM + indexDist + indexGDPDist + indexOriginCorrViaGeo + indexOriginCorrViaSectors + indexSectorCorr
	m = @formula(πjodt ~ jStr&lntjodt + indexOriginCorrViaGeo + indexOriginCorrViaSectors + indexSectorCorr + fe(j)&fe(o)&fe(t) + fe(j)&fe(d)&fe(t) + fe(j)&fe(o)&fe(d))
	@show result = nlreg(df, m, Poisson(), LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8)
	w = map(x->in(x,[ "jStr: $j & lntjodt" for j in unique(df.jStr) ]),coefnames(result))
	transform = ParameterTransformation(x->w'*x/J,x->w/J,"Elasticity / Mean Elasticity")
	adjustedResult=addTransformation(result,transform)
	@show rt = RegressionTable(rt,adjustedResult,[false,true,false,true,true,true],test=("SGM",["indexOriginCorrViaGeo","indexOriginCorrViaSectors","indexSectorCorr"]),spec="(4)")

	# Specification: SGM + absΔlnYodt + (lnDistw + absΔlnYodt)&lnt̂jodt + indexDist + indexGDPDist + indexOriginCorrViaGeo + indexOriginCorrViaSectors + indexSectorCorr
	m = @formula(πjodt ~ jStr&lntjodt + absΔlnYodt + (lnDistw + absΔlnYodt)&lnt̂jodt + indexOriginCorrViaGeo + indexOriginCorrViaSectors + indexSectorCorr + fe(j)&fe(o)&fe(t) + fe(j)&fe(d)&fe(t) + fe(j)&fe(o)&fe(d))
	@show result = nlreg(df, m, Poisson(), LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8)
	w = map(x->in(x,[ "jStr: $j & lntjodt" for j in unique(df.jStr) ]),coefnames(result))
	transform = ParameterTransformation(x->w'*x/J,x->w/J,"Elasticity / Mean Elasticity")
	adjustedResult=addTransformation(result,transform)
	@show rt = RegressionTable(rt,adjustedResult,[false,true,true,true,true,true],test=("SGM",["lnDistw & lnt̂jodt","absΔlnYodt & lnt̂jodt","indexOriginCorrViaGeo","indexOriginCorrViaSectors","indexSectorCorr"]),spec="(5)")

	saveRegressionTable(rt,dirs.results*"onlineAppendix/","TableO1")

	save(dirs.data*"indices.csv",df[:,[:j,:o,:d,:t,Symbol.(["indexDist","indexGDPDist","indexOriginCorrViaGeo","indexOriginCorrViaSectors","indexSectorCorr"])...]])

	return nothing
end
