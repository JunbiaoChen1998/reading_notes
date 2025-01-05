
struct ParameterTransformation
	f::Function
	grad::Function
	name::String
end
function addTransformation(result::GLFixedEffectModel,transform::ParameterTransformation)
	β̂ = copy(result.coef)
	V̂ = copy(result.vcov)
	grad = hcat(transform.grad(β̂),I)

	β̂ = vcat(transform.f(β̂),β̂)
	V̂ = grad'*V̂*grad

	return GLFixedEffectModel(β̂, # Vector of coefficients
	        V̂,    # Covariance matrix
	        result.vcov_type,           # CovarianceEstimator
	        result.nclusters,
	        result.iterations,
	        result.converged,
	        result.esample,
	        result.augmentdf,
	        result.distribution,
	        result.link,
	        vcat(transform.name,result.coefnames),
	        result.yname,
	        result.formula, # Original formula
	        result.formula_schema,
	        result.nobs,   # Number of observations
	        result.dof_residual,  # nobs - degrees of freedoms
	        result.deviance, # Deviance of the fitted model
	        result.nulldeviance, # null deviance
	        result.gradient,   # concentrated gradient
	        result.hessian  # concentrated hessian
	    )
end

struct RegressionTable
	β::DataFrame
	se::DataFrame
	ttest::DataFrame
	fe::DataFrame
	stats::DataFrame
	specs::Vector{String}
end
function RegressionTable(varNames::Vector{String},feNames::Vector{String})
	β = DataFrame(names=varNames)
	se = DataFrame(names=varNames)
	ttest = DataFrame(names=varNames)
	fe = DataFrame(names=feNames)
	stats = DataFrame(names=["nobs","dof","dev","null","chi2","pval"])
	specs = String[]
	RegressionTable(β,se,ttest,fe,stats,specs)
end



function RegressionTable(rt::RegressionTable,result::GLFixedEffectModel,feIncluded::Vector{Bool};test=("",String[]),spec=nothing,testExternal=[missing,missing,missing])
	if isnothing(spec)
		spec = "$(length(rt.specs)+1)"
	end
	typeof(test) == Tuple{String,Vector{String}} || error("Keyword 'test' must be a Tuple{String,Vector{String}}")

	β̂ = result.coef
	V̂ = result.vcov

	Nparam = length(β̂)

	tTests = ccdf.(Ref(FDist(1, dof_residual(result))), abs2.(coef(result)./stderror(result)))
	decoratedβ̂ = decoration.(string.(round.(coef(result),digits=3)),tTests)

	β = join(rt.β,rename(DataFrame(names=coefnames(result),S=decoratedβ̂),:S=>Symbol(spec)),on=:names,kind=:left)
	se = join(rt.se,rename(DataFrame(names=coefnames(result),S=round.(stderror(result),digits=3)),:S=>Symbol(spec)),on=:names,kind=:left)
	ttest  = join(rt.ttest,rename(DataFrame(names=coefnames(result),S=tTests),:S=>Symbol(spec)),on=:names,kind=:left)
	fe = join(rt.fe,rename(DataFrame(names=rt.fe.names,S=feIncluded),:S=>Symbol(spec)),on=:names,kind=:left)

	if length(test[2]) > 0
		R = permutedims(hcat([ result.coefnames .== name  for name in test[2] ]...))
		r = zeros(size(R)[1])
		χ2 = (R*β̂.-r)'*((R*V̂*R')\(R*β̂.-r))
		pval = ccdf(Chisq(length(r)),χ2)
		stats = join(rt.stats,rename(DataFrame(names=["nobs","dof","dev","null","chi2","pval"],S=[Int64(nobs(result)),Int64(dof_residual(result)),round(deviance(result),digits=3),test[1],round(χ2,digits=3),round(pval,digits=3)]),:S=>Symbol(spec)),on=:names,kind=:left)
	else
		stats = join(rt.stats,rename(DataFrame(names=["nobs","dof","dev","null","chi2","pval"],S=[Int64(nobs(result)),Int64(dof_residual(result)),round(deviance(result),digits=3),testExternal...]),:S=>Symbol(spec)),on=:names,kind=:left)
	end

	specs = vcat(rt.specs,spec)

	return RegressionTable(β,se,ttest,fe,stats,specs)
end


function createTable(rt::RegressionTable)
	table = DataFrame(names=[])
	for spec in rt.specs
		table[!,Symbol(spec)] = []
	end
	for i=1:(nrow(rt.β))
		push!(table,rt.β[i,:])
	    seRow = [rt.se[i,spec] for spec in Symbol.(rt.specs)]
	    push!(table,hcat(missing,map(x->ismissing(x) ? "" : "($x)",seRow)...))
	end
	table = vcat(table,rt.fe,rt.stats)
	for var in names(table)
	    table[ismissing.(table[:,var]),var] .= ""
	end
	return table
end

Base.copy(rt::RegressionTable) = RegressionTable(copy(rt.β),copy(rt.se),copy(rt.fe),copy(rt.ttest),copy(rt.stats),copy(rt.specs))
Base.show(io::IO,rt::RegressionTable) = show(io,createTable(rt))

function subRegressionTable(rt::RegressionTable,subSpecs::Vector{String})
	syms = vcat(:names,Symbol.(subSpecs))
	β = rt.β[:,syms]
	se = rt.se[:,syms]
	ttest = rt.ttest[:,syms]
	fe = rt.fe[:,syms]
	stats = rt.stats[:,syms]
	return RegressionTable(β,se,ttest,fe,stats,subSpecs)
end


function subRegressionTable(rt::RegressionTable,drop::InvertedIndex{String})
	β = rt.β[:,Not(Symbol.(drop.skip))]
	se = rt.se[:,Not(Symbol.(drop.skip))]
	ttest = rt.ttest[:,Not(Symbol.(drop.skip))]
	fe = rt.fe[:,Not(Symbol.(drop.skip))]
	stats = rt.stats[:,Not(Symbol.(drop.skip))]
	return RegressionTable(β,se,ttest,fe,stats,setdiff(rt.specs,x.skip))
end

import DataFrames.rename!
function rename!(rt::RegressionTable,p::Pair{String,String})
	rename!(rt.β,p)
	rename!(rt.se,p)
	rename!(rt.ttest,p)
	rename!(rt.fe,p)
	rename!(rt.stats,p)
	rt.specs[rt.specs .== p[1]] = [p[2]]
	return rt
end

function saveRegressionTable(rt::RegressionTable,dir::String,name::String)
	table = createTable(rt)
	save(dir*"$name.csv",table)
	io = open(dir*"$name.tex","w")
	write(io,latexify(table,env=:table,latex=false))
	close(io)
	return nothing
end
