
function goodnessOfFit(dirs::ProjDirectories,data::ProjData)

	# # load data
	# td = TradeData(dirs.data)

	# include(dirs.code*"estimationTools/utilities.jl")
	# include(dirs.code*"estimationTools/RegressionTable.jl")
	# include(dirs.code*"estimationTools/loadData.jl")
	# include(dirs.code*"estimationTools/loadEstimates.jl")


	σ,σse,dfΛ,dfΦ,K,dfCES,dfSGM = loadEstimates(dirs)

	sitcToSec = load(dirs.data*"sitcToSec.csv") |> DataFrame
	secNames = load(dirs.data*"secNames.csv") |> DataFrame
	secLabels = secNames.name
	S = length(data.trade.sitcCodes)
	N = data.N
	T = data.T
	J = data.J

	tsodt = reshape(data.trade.sitcTariffs,(S,N,N,T))
	sort!(dfΛ,reverse([:sitc,:k]))
	λsk = reshape(dfΛ.Λ,(S,K))
	dfΦ = join(dfΦ,by(dfΦ,[:d,:t],Xdt = :Xkodt => sum),on=[:d,:t],kind=:left)
	dfΦ.πkodt = dfΦ.Xkodt ./ dfΦ.Xdt
	sort!(dfΦ,reverse([:k,:o,:d,:t]))
	Xkodt = reshape(dfΦ.Xkodt,(K,N,N,T))
	πkodt = reshape(dfΦ.πkodt,(K,N,N,T))
	tkodt = reshape(dfΦ.tkodt,(K,N,N,T))

	π̂LFMjodt = reshape(data.trade.secMap'*reshape([ sum( (tsodt[s,o,d,t] ./ tkodt[:,o,d,t]).^(.-σ) .* λsk[s,:] .* πkodt[:,o,d,t] )
	                    												for s=1:S, o=1:N, d=1:N, t=1:T ],(S,N*N*T)),(J,N,N,T))

	# include(dirs.code*"estimationTools/loadGravityData.jl")
	df = copy(data.grav)


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
	df.tjodt = 1 .+ df.tariff ./ 100
	df = join(df,by(df,[:j,:d,:t],lnt̄jdt = :lntjodt => mean),on=[:j,:d,:t],kind=:left)
	df.lnt̂jodt = df.lntjodt - df.lnt̄jdt
	df.jod = string.(df.j) .* "_" .* df.o .* "_" .* df.d

	dfΛ = join(dfΛ,sitcToSec,on=[:sitc])
	dfΛ = by(dfΛ,[:sec,:k],Λ = :Λ => sum)
	sort!(dfΛ,reverse([:sec,:k]))
	λjk = reshape(dfΛ.Λ,(J,K))
	sum(λjk,dims=1)

	# LFM
	θ = minimum(σ)
	ρ = 1 .- θ./σ

	sort!(dfΦ,reverse([:k,:o,:d,:t]))
	Xkodt = reshape(dfΦ.Xkodt,(K,N,N,T))
	πkodt = reshape(dfΦ.πkodt,(K,N,N,T))
	tkodt = reshape(dfΦ.tkodt,(K,N,N,T))
	π̃kodt = (πkodt.^(1 .- reshape(ρ,(K,1,1,1)))) .* (sum(πkodt,dims=2).^reshape(ρ,(K,1,1,1)))
	ωkodt = π̃kodt ./ sum(π̃kodt,dims=1)

	# SGM
	θSGM = copy(θ)
	σSGM = -dfSGM.σ
	ρSGM = 1 .- θSGM./σSGM

	sort!(df,reverse([:j,:o,:d,:t]))
	Xjodt = reshape(df.Xjodt,(J,N,N,T))
	πjodt = reshape(df.πjodt,(J,N,N,T))
	π̃jodt = (πjodt.^(1 .- reshape(ρSGM,(J,1,1,1)))) .* (sum(πjodt,dims=2).^reshape(ρSGM,(J,1,1,1)))
	ωjodt = π̃jodt ./ sum(π̃jodt,dims=1)
	tjodt = reshape(df.tjodt,(J,N,N,T))
	lntjodt = reshape(df.lntjodt,(J,N,N,T))

	πjod = dropdims(mean(πjodt,dims=4),dims=4)
	πjo = mean(πjod,dims=3)
	πod = sum(πjod,dims=1)
	πjo = reshape(πjo,(J,N))

	πjo = reshape(mean(πjodt,dims=(3,4)),(J,N))
	πodt = reshape(sum(πjodt,dims=1),(N,N,T))
	πjdt = reshape(sum(πjodt,dims=2),(J,N,T))

	# Create indices for testing IIA

	#	Distance
	distanceDf = rename(unique(df[:,[:o,:d,:ΔlnDistw]]),:d=>:o′)
	sort!(distanceDf,reverse([:o,:o′]))
	reshape(distanceDf.o,(N,N))
	distance = reshape(distanceDf.ΔlnDistw,(N,N))
	maximum(abs,distance - distance')
	w = copy(distance)
	w -= Diagonal(w)
	tariffs = reshape(permutedims(lntjodt,(2,1,3,4)),(N,J*N*T))
	indexDist = w*tariffs
	indexDist = permutedims(reshape(indexDist,(N,J,N,T)),(2,1,3,4))
	df.indexDist = vec(indexDist)

	#	GDP Per Capita Difference
	incomeDf = rename(unique(df[:,[:o,:d,:t,:lnYdt,:lnYot]]),:d=>:o′,:lnYdt=>:lnYo′t)
	incomeDf.ΔlnYoo′t = abs.(incomeDf.lnYot - incomeDf.lnYo′t)
	sort!(incomeDf,reverse([:o,:o′,:t]))
	gdfDifference = reshape(incomeDf.ΔlnYoo′t,(N,N,T))
	indexGDPDist = fill(NaN,J,N,N,T)
	for t = 1:T
		w = copy(gdfDifference[:,:,t])
		w -= Diagonal(w)
		tariffs = reshape(permutedims(lntjodt[:,:,:,t],(2,1,3)),(N,J*N))
		indexGDPDist[:,:,:,t] = permutedims(reshape(w*tariffs,(N,J,N)),(2,1,3))
	end
	any(isnan.(indexGDPDist))
	df.indexGDPDist = vec(indexGDPDist)

	# Correlation based measures:
	πod = reshape(sum(πjod,dims=1),(N,N))
	πjo = reshape(mean(πjod,dims=3),(J,N))

	#	Correlation of aggregate origin shares (within destination)
	sim = corrMat(πod)
	w = sim
	w -= Diagonal(w)
	tariffs = reshape(permutedims(lntjodt,(2,1,3,4)),(N,J*N*T))
	df.indexOriginCorrViaGeo = vec(permutedims(reshape(w*tariffs,(N,J,N,T)),(2,1,3,4)))

	#	Correlation of origin shares (within sector)
	sim = corrMat( πjo' ./ sum(πjo',dims=1) )
	w = sim
	w -= Diagonal(w)
	tariffs = reshape(permutedims(lntjodt,(2,1,3,4)),(N,J*N*T))
	df.indexOriginCorrViaSectors = vec(permutedims(reshape(w*tariffs,(N,J,N,T)),(2,1,3,4)))

	#	Correlation of sectoral shares (within origin)
	sim = corrMat(πjo  ./ sum(πjo,dims=1))
	w = sim
	w -= Diagonal(w)
	tariffs = reshape(lntjodt,(J,N*N*T))
	df.indexSectorCorr = vec(reshape(w*tariffs,(J,N,N,T)))

	df.π̂LFMjodt = vec(π̂LFMjodt)
	df.lnπ̂LFMjodt = log.(df.π̂LFMjodt)
	df.π̂LFMjodt_Center = vec([ sum(  λjk[j,:] .* πkodt[:,o,d,t] ) for j=1:J, o=1:N, d=1:N, t=1:T ])
	df.lnπ̂LFMjodt_Center = log.(df.π̂LFMjodt_Center)
	df.Δlnπ̂LFMjodt_Δt = df.lnπ̂LFMjodt .- df.lnπ̂LFMjodt_Center
	df.constant = ones(nrow(df))
	df.jo = string.(df.j) .* "_" .* df.o
	df.od = df.o .* "_" .* df.d
	df.jd = string.(df.j) .* "_" .* "_" .* df.d
	df.jod = string.(df.j) .* "_" .* df.o .* "_" .* df.d
	df.odt = df.o .* "_" .* df.d .* "_" .* string.(df.t)
	df.jodt = string.(df.j) .* "_" .* df.o .* "_" .* df.d .* "_" .* string.(df.t)
	df.ΔabslnYodt = abs.(df.lnYot-df.lnYdt)


	# Estimation
	varNames = [ "lnπ̂LFMjodt", "Δlnπ̂LFMjodt_Δt",
				 "lnDistw & lnt̂jodt","absΔlnYodt & lnt̂jodt",
				 "indexOriginCorrViaGeo", "indexOriginCorrViaSectors", "indexSectorCorr"]
	feNames = ["j & lntjodt","absΔlnYodt","j & o & d","j & o & t","j & d & t"]


	vcovSpec = Vcov.cluster(:jod)
	# vcovSpec = Vcov.cluster(:j,:od) # looks pretty good too

	rt = RegressionTable(varNames,feNames)


	m = @formula(πjodt ~ lnπ̂LFMjodt + fe(constant)  )
	@show result = nlreg(df, m, Poisson(), LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8)

	β̂ = result.coef
	V̂ = result.vcov

	Nparam = length(β̂)
	R = ones(1,1)
	r = ones(1)
	χ2 = (R*β̂.-r)'*((R*V̂*R')\(R*β̂.-r))
	testVal = ccdf(Chisq(length(r)),χ2)

	@show rt = RegressionTable(rt,result,[false,false,false,false,false],testExternal=["LFM",χ2,testVal],spec="LFM")


	m = @formula πjodt ~ absΔlnYodt + (lnDistw + absΔlnYodt)&lnt̂jodt + indexOriginCorrViaGeo + indexOriginCorrViaSectors + indexSectorCorr + jStr&lntjodt + fe(j)&fe(o)&fe(t) + fe(j)&fe(d)&fe(t) + fe(j)&fe(o)&fe(d)
	@show result = nlreg(df, m, Poisson(), LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8)
	@show rt = RegressionTable(rt,result,[true,true,true,true,true],test=("SGM",["lnDistw & lnt̂jodt","absΔlnYodt & lnt̂jodt","indexOriginCorrViaGeo","indexOriginCorrViaSectors","indexSectorCorr"]),spec="SGM Test")


	m = @formula πjodt ~ lnπ̂LFMjodt  + absΔlnYodt + (lnDistw + absΔlnYodt)&lnt̂jodt + jStr&lntjodt + fe(j)&fe(o)&fe(d) + fe(j)&fe(o)&fe(t) + fe(j)&fe(d)&fe(t)
	@show result = nlreg(df, m, Poisson(), LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8)
	@show rt = RegressionTable(rt,result,[true,true,true,true,true],test=("SGM Restriction",["lnDistw & lnt̂jodt","absΔlnYodt & lnt̂jodt"]),spec="Spec Test 1")

	m = @formula πjodt ~ lnπ̂LFMjodt + indexOriginCorrViaGeo + indexOriginCorrViaSectors + indexSectorCorr + jStr&lntjodt + fe(j)&fe(o)&fe(t) + fe(j)&fe(d)&fe(t) + fe(j)&fe(o)&fe(d)
	@show result = nlreg(df, m, Poisson(), LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8)
	@show rt = RegressionTable(rt,result,[true,false,true,true,true],test=("SGM Restriction",["indexOriginCorrViaGeo","indexOriginCorrViaSectors","indexSectorCorr"]),spec="Spec Test 2")


	m = @formula πjodt ~ lnπ̂LFMjodt + absΔlnYodt + (lnDistw + absΔlnYodt)&lnt̂jodt + indexOriginCorrViaGeo + indexOriginCorrViaSectors + indexSectorCorr + jStr&lntjodt + fe(j)&fe(o)&fe(t) + fe(j)&fe(d)&fe(t) + fe(j)&fe(o)&fe(d)
	@show result = nlreg(df, m, Poisson(), LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8)
	@show rt = RegressionTable(rt,result,[true,true,true,true,true],test=("SGM Restriction",["lnDistw & lnt̂jodt","absΔlnYodt & lnt̂jodt","indexOriginCorrViaGeo","indexOriginCorrViaSectors","indexSectorCorr"]),spec="Spec Test 3")

	m = @formula πjodt ~ lnπ̂LFMjodt + Δlnπ̂LFMjodt_Δt + absΔlnYodt + (lnDistw + absΔlnYodt)&lnt̂jodt + indexOriginCorrViaGeo + indexOriginCorrViaSectors + indexSectorCorr + jStr&lntjodt + fe(j)&fe(o)&fe(t) + fe(j)&fe(d)&fe(t) + fe(j)&fe(o)&fe(d)
	@show result = nlreg(df, m, Poisson(), LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8)
	@show rt = RegressionTable(rt,result,[true,true,true,true,true],test=("SGM Restriction",["lnDistw & lnt̂jodt","absΔlnYodt & lnt̂jodt","indexOriginCorrViaGeo","indexOriginCorrViaSectors","indexSectorCorr"]),spec="Spec Test 4")

	@show table = createTable(rt)



	save(dirs.results*"onlineAppendix/TableO3.csv",table)
	io = open(dirs.results*"onlineAppendix/TableO3.tex","w")
	write(io,latexify(table,env=:table,latex=false))
	close(io)

	return nothing
end
