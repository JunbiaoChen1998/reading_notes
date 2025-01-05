function timeReport(tstart)
	tend = time()
	minutes = Int64(floor((tend-tstart)/60))
	seconds = round(tend-tstart-60*minutes,digits=3)
    return "Total time: $minutes minutes and $seconds seconds."
end



function offDiagonal(X::Matrix)
	(M,N) = size(X)
	out = []
	for m=1:M, n=1:N
		m != n && push!(out,X[m,n])
	end
	return out
end


function plotHistory(history::DataFrame,starts::Vector{Float64};cutoff=0)
    plotData = history[history.iterTime.>=cutoff,:]
    p1 = plot(plotData.iterTime,plotData.mpd,legend=false,title="Mean poisson deviance")
    p2 = plot(plotData.iterTime,plotData.η,legend=false,title="Eta")
    p3 = plot(plotData.iterTime,plotData.maxΔT̃k,legend=false,title="Maximum change in inferred normalized tariffs")
    p4 = plot(plotData.iterTime,plotData.maxΔXk,legend=false,title="Maximum change in inferred expenditure shares")
    p5 = plot(plotData.iterTime,plotData.ΛZeroShare,legend=false,title="Share of zero entries in factor weights")
    p6 = plot(plotData.iterTime,plotData.ΦZeroShare,legend=false,title="Share of zero entries in factors")
    plt = plot(p1,p2,p3,p4,p5,p6,layout=(6,1),size=1.5 .* (700,800))
    if length(starts) >= 0
        plot!(plt,starts[starts.>cutoff],seriestype="vline",color=:red,alpha=.5)
    end
    display(plt)
    return plt
end


function sectorLabels(sitcCodes::DataFrame,sitcDescriptions::DataFrame)

    codes = DataFrame(i=1:length(sitcCodes.sitc),sitc=sitcCodes.sitc)
    codes.sitc2 = map(x->x[1:2],codes.sitc)
    firstCodes = by(codes,[:sitc2],i = :i => first)
    firstCodes = join(firstCodes,rename(sitcDescriptions,:code=>:sitc2),on=[:sitc2],kind=:left)

    gaps = firstCodes.i[2:end]-firstCodes.i[1:end-1]
    gapsize = 1
    id = vcat(false,gaps .<= gapsize) .| vcat(gaps .<= gapsize,false)
    # @show firstCodes[id,:]

    firstCodes[firstCodes.i .== 691,:description] .= "Bags, clothing accessories, and footwear"
    firstCodes = firstCodes[map(x->!in(x,[692,720]),firstCodes.i),:]

    firstCodes[firstCodes.i .== 770,:description] .= "Other animals, weapons, coins, and gold"
    firstCodes = firstCodes[map(x->!in(x,771:773),firstCodes.i),:]

    gaps = firstCodes.i[2:end]-firstCodes.i[1:end-1]
    gapsize = 2
    id = vcat(false,gaps .<= gapsize) .| vcat(gaps .<= gapsize,false)
    id = [findall(id)...,findall(id)[end]+1]
    # @show firstCodes[id,:]

    firstCodes[firstCodes.i .== 225,:description] .= "Gas, oils, and waxes"
    firstCodes = firstCodes[map(x->!in(x,[227,229,241]),firstCodes.i),:]



    gaps = firstCodes.i[2:end]-firstCodes.i[1:end-1]
    gapsize = 3
    windowsize = 2
    ids = gaps .<= gapsize
    tmp = copy(firstCodes)
    tmp.tooClose = vcat(false,ids) .| vcat(ids,false)
    # @show tmp[union([ id-windowsize:id+windowsize for id in findall(tmp.tooClose) ]...),:]

    firstCodes[firstCodes.i .== 91,:description] .= "Beverages and misc. edible products"
    firstCodes = firstCodes[map(x->!in(x,[94]),firstCodes.i),:]

    firstCodes[firstCodes.i .== 124,:description] .= "Cork, wood, and crude rubber"
    firstCodes = firstCodes[map(x->!in(x,[127]),firstCodes.i),:]

    firstCodes[firstCodes.i .== 299,:description] .= "Fertilizers, explosives, and chemicals"
    firstCodes = firstCodes[map(x->!in(x,[303,306,329]),firstCodes.i),:]

    firstCodes[firstCodes.i .== 685,:description] .= "Household durables and parts thereof"
    firstCodes = firstCodes[map(x->!in(x,[688]),firstCodes.i),:]

    gaps = firstCodes.i[2:end]-firstCodes.i[1:end-1]
    gapsize = 4
    windowsize = 2
    ids = gaps .<= gapsize
    tmp = copy(firstCodes)
    tmp.tooClose = vcat(false,ids) .| vcat(ids,false)
    ids = union([ id-windowsize:id+windowsize for id in findall(tmp.tooClose) ]...)
    # @show tmp[intersect(ids,1:nrow(tmp)),:]


    firstCodes[firstCodes.i .== 770,:description] .= "Other animals, weapons, coins, gold, and residual"
    firstCodes = firstCodes[map(x->!in(x,[774]),firstCodes.i),:]

    return firstCodes
end



function decoration(s::String, pval::Float64)
  if pval<0.0
      error("p value needs to be nonnegative.")
  end
  if (pval > 0.1)
      return "$s"
  elseif (pval > 0.05)
      return "$s"
  elseif (pval > 0.01)
      return "$s*"
  elseif (pval > 0.001)
      return "$s**"
  else
      return "$s***"
  end
end


R2(Y,Ŷ) = 1-sum((Y .- Ŷ).^2)/sum((Y .- mean(Y)).^2)

function corrMat(shares::Matrix{Float64})
	shares = shares .- mean(shares,dims=2)
	sim = shares*shares'
	return sim ./ (sqrt.(diag(sim))*sqrt.(diag(sim))')
end

function cumR2(shares::Matrix{Float64})
	shares = shares .- mean(shares,dims=2)
	(eigvals,F) = eigen(shares*shares')
	p = sortperm(eigvals,rev=true)
	eigvals = eigvals[p]
	F = F[:,p]
	y = cumsum(eigvals)./sum(eigvals)
end


function waldTest(result,R,r)
    β̂ = result.coef
    V̂ = result.vcov
    Nparam = length(β̂)
    χ2 = (R*β̂.-r)'*((R*V̂*R')\(R*β̂.-r))
    pval = ccdf(Chisq(length(r)),χ2)

    return [χ2,pval]
end
