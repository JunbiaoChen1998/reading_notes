
function prepareUseShares(dataDir::String,countries::Vector{String},years::Vector{Int64})

	print("\nPreparing savings intermediate and final use datasets...\n\n")
	tstart = time()

    wiot = load(dataDir*"wiot/wiot_full.dta") |> DataFrame
    rename!(wiot,:row_country => :o,:col_country => :d,:year => :t)

	wiodToSec = load(dataDir*"wiodToSec.csv") |> DataFrame
    wiodToSec.WIOD = Int64.(wiodToSec.WIOD)
    wiodToSec.sec = Int64.(wiodToSec.sec)
    rename!(wiodToSec,:sec=>:j)

    # get final use shares
    finalUse = by(wiot[(wiot.row_item .<= 16) .& (wiot.col_item .> 35) .& map(x->in(x,years),wiot.t) .& map(x->in(x,countries),wiot.d),:],[:row_item,:d,:t],value = :value => sum)
    finalUse = by(join(rename(finalUse,:row_item => :WIOD),wiodToSec,on=[:WIOD],kind=:left),[:j,:d,:t],Xjdt = :value => sum)

    sort!(finalUse,[:j,:d,:t])

    finalUse = join(finalUse,by(finalUse,[:d,:t],Xdt = :Xjdt => sum),on=[:d,:t],kind=:left)
    finalUse.value = finalUse.Xjdt ./ finalUse.Xdt

    save(dataDir*"finalUseShares.csv",finalUse[:,[:j,:d,:t,:value]])

    # get intermediate use
    intermediateUse = by(wiot[(wiot.row_item .<= 16) .& (wiot.col_item .<= 16) .& map(x->in(x,years),wiot.t) .& map(x->in(x,countries),wiot.d),:],[:row_item,:col_item,:d,:t],value = :value => sum)
    intermediateUse = join(rename(intermediateUse,:row_item => :WIOD),wiodToSec,on=[:WIOD],kind=:left)
    rename!(intermediateUse,:j => :i,:WIOD => :WIOD′)
    intermediateUse = join(rename(intermediateUse,:col_item => :WIOD),wiodToSec,on=[:WIOD],kind=:left)
    intermediateUse = by(intermediateUse,[:i,:j,:d,:t],Xijdt = :value => sum)

    # get revenue
    revenue = by(wiot[(wiot.row_item .<= 16) .& map(x->in(x,years),wiot.t) .& map(x->in(x,countries),wiot.o),:],[:row_item,:o,:t],value = :value => sum)
    revenue = join(rename(revenue,:row_item => :WIOD),wiodToSec,on=[:WIOD],kind=:left)
    revenue = by(revenue,[:j,:o,:t],Rjot = :value => sum)

    rename!(revenue,:o => :d)
    sort!(intermediateUse,[:i,:j,:d,:t])

    intermediateUse = join(intermediateUse,revenue,on=[:j,:d,:t],kind=:left)
    intermediateUse.value = intermediateUse.Xijdt ./ intermediateUse.Rjot

    save(dataDir*"intermediateUseShares.csv",intermediateUse[:,[:i,:j,:d,:t,:value]])

	print("Done. $(timeReport(tstart))\n")

    return nothing
end