
function prepareDataFiles(dataDir::String)

	sitcData,secData,sitcToSec = prepareFinalData(loadAndCombineData(dataDir))

	print("\nPreparing and saving data files...\n\n")
	tstart = time()

    # some sitc codes are completely zero or missing in all origin-destination-years,
	# 	but only residuals are completely missing
    missingOrZero(x) = ismissing(x) || x == 0
    tmp = by(sitcData,[:sitc],missOrZero=:value=>x->mean(missingOrZero.(x)),allMiss=:value=>x->all(ismissing.(x)))
	sort(unique(tmp[tmp.missOrZero .== 1,[:sitc]]))
	sort(unique(tmp[tmp.allMiss,[:sitc]]))
	sort(unique(tmp[(tmp.missOrZero .== 1).&(.!tmp.allMiss),[:sitc]]))

	# Drop fully zero/miss sitc sectors.
	#	In estimation, they will get zero factorweights to match their zeros,
	#	and their missing values will be captured by their sec's "sitc" residual.
	dropSitc = sort(unique(tmp[(tmp.missOrZero .== 1).&(.!tmp.allMiss),:sitc]))
	sitcData = sitcData[map(x->!in(x,dropSitc),sitcData.sitc),:]

    sitcCodes = sort(unique(sitcData.sitc))
    # secCodes = sort(unique(secData.sec))

	countries = sort(unique(sitcData.o))
	countries != sort(unique(sitcData.d)) && error("Origins and destinations are not the same with sitc data.")
	countries != sort(unique(secData.o)) && error("Origins and destinations are not consistent between sitc and sec data.")
	countries != sort(unique(secData.d)) && error("Origins and destinations are not consistent between sitc and sec data.")

	years = sort(unique(sitcData.t))
    years != sort(unique(secData.t)) && error("The two data sets must contain the same years.")

	# print("Sorting and reshaping trade data...\n")

    # sort the data
    sort!(secData,reverse([:sec,:o,:d,:t]))
    sort!(sitcData,reverse([:sitc,:o,:d,:t]))

	save(dataDir*"secData.csv",secData)
	save(dataDir*"sitcData.csv",sitcData)
	save(dataDir*"sitcToSec.csv",sitcToSec)
	save(dataDir*"sitcCodes.csv",DataFrame(sitc=sitcCodes))

	print("Done. $(timeReport(tstart))\n")
	
	return nothing
end

missingOrFirst(x) = any(.!ismissing.(x)) ? first(skipmissing(x)) : missing
missingOrMinimum(x) = any(.!ismissing.(x)) ? minimum(skipmissing(x)) : missing

function createSitcToSec(dataDir::String)

	# sitc descriptions
	dfNameSITC = load(dataDir*"sectors/S2_descriptions.csv") |> DataFrame
	rename!(dfNameSITC,Symbol("Commodity Code") => :sitc)

	# get HS to sitc rev 2 concordance
	dfSITC = load(dataDir*"sectors/HS_to_S2.csv",colparsers=[String,String,String,String]) |> DataFrame

	rename!(dfSITC,Symbol("HS - Combined  Product Description") => :name_hs,
				Symbol("HS - Combined  Product Code") => :hs,
				Symbol("SITC Revision 2 Product Description") => :name_sitc,
				Symbol("SITC Revision 2 Product Code") => :sitc)
	dfSITC.sitc4 = map(x->length(x) > 4 ? x[1:4] : x,dfSITC.sitc)

	# get list of all sitc4 codes in the comtrade data
	dfSITC4 = DataFrame(sitc4=String[])
	for year in 1984:2007
		colparsers = Any[nothing for v=1:13]
		colparsers[8] = String
		tdf = load(dataDir*"comtrade/$year.csv",colparsers=colparsers) |> @unique() |> DataFrame
		dfSITC4 = join(dfSITC4,tdf,on=[:sitc4],kind=:outer)
	end

	# 4-digit sitc codes to classify (expands by "334", petroleum)
	dfSITC4 = join(unique(dfSITC[:,[:sitc4]]),dfSITC4,on=[:sitc4],kind=:outer)

	#  get HS to isic rev 3 concordance
	dfISIC = load(dataDir*"sectors/HS_to_I3.csv",colparsers=[String,String,String,String]) |> DataFrame
	rename!(dfISIC,Symbol("HS - Combined  Product Description") => :name_hs,
				Symbol("HS - Combined  Product Code") => :hs,
				Symbol("ISIC Revision 3 Product Description") => :name_isic,
				Symbol("ISIC Revision 3 Product Code") => :isic)
	dfISIC.isic2 = map(x->length(x) > 2 ? x[1:2] : x,dfISIC.isic)

	# match sitc rev 2 and isic rev 3
	df = join(dfISIC,dfSITC,on=[:hs,:name_hs],kind=:inner)
	for v in names(df)
		@show v,any(ismissing.(df[:,v]))
	end

	length(unique(df.sitc4))
	length(unique(df.isic2))

	# map isic codes to WIOD sector aggregates
	# map from isic to WIOD aggregates
	isicToSec(x::Missing) = missing
	isicToSec(x::Integer) = string(x)
	isicToSec(x::String) = isicToSec(_isicToSec(parse(Int64,x)))
	function _isicToSec(x::Integer)
		if x >= 0 && x <= 5
			return 1
		elseif	x >= 10  && x <= 14
			return 2
		elseif x >= 15  && x <= 16
			return 3
		elseif x >= 17  && x <= 19
			return 4
		elseif x == 20
			return 5
		elseif x >= 21  && x <= 22
			return 6
		elseif x == 23
			return 7
		elseif x >= 24  && x <= 25
			return 8
		elseif x == 26
			return 9
		elseif x >= 27  && x <= 28
			return 10
		elseif x == 29
			return 11
		elseif x >= 30  && x <= 33
			return 12
		elseif x >= 34  && x <= 35
			return 13
		elseif x >= 36  && x <= 37
			return 14
		elseif x >= 40  && x <= 41
			return 15
		elseif x == 45
			return 16
		elseif x >= 50  && x <= 55
			return 17
		elseif x == 60 || x == 63
			return 18
		elseif x == 61
			return 19
		elseif x == 62
			return 20
		elseif x == 64
			return 21
		elseif x >= 65  && x <= 66
			return 22
		elseif x >= 70  && x <= 74
			return 23
		elseif x == 75 || x == 80 || x == 85
			return 24
		elseif x >= 90 && x <= 95
			return 25
		else
			return missing
		end
	end
	df.sec = map(isicToSec,df.isic2)

	unique(df[ismissing.(df.sec),:isic2])
	unique(map(isicToSec,df.isic2))

	# df = df[:,[:sitc4,:sec]]
	df = df[:,[:hs,:sitc4,:sec]]
	# df = join(df,by(df,[:sitc4],counter = :sec => length),on=[:sitc4],kind=:left)

	# How many products were not classified?
	df.miss = map(ismissing,df.sec)
	sum(df.miss)

	### refine the map using the most common classification of products within the 4 digit level

	# count the number of hs products within each sitc4-sec pair
	tdf = by(df,[:sitc4,:sec],counter = :hs => length)

	# remove sitc4-bin pairs that have fewer products than the most common (including no classification)
	tdf = join(tdf,by(tdf,[:sitc4],total = :counter => sum,max = :counter => maximum),on=[:sitc4],kind=:left)
	tdf = tdf[tdf.counter .== tdf.max,:]

	# manually classify those that are missing
	manual = tdf[ismissing.(tdf.sec),:sitc4]
	dfNameSITC[map(x->x ∈ manual,dfNameSITC.sitc),:]
	replacements = Dict(
		"0811" => "1",
		"0812" => "1",
		"0819" => "14",
		"1213" => "14",
		"2332" => "14",
		"2511" => "14",
		"2614" => "4",
		"2633" => "14",
		"2672" => "14",
		"2686" => "14",
		"2690" => "4",
		"2820" => "14",
		"2881" => "14",
		"2882" => "14",
		"2890" => "14",
		"2911" => "1",
		"2919" => "1",
		"7933" => "14"
	)
	tdf.sec = map(x->x.sitc4 ∈ manual ? replacements[x.sitc4] : x.sec,eachrow(tdf))

	tdf = unique(tdf)
	select!(tdf,Not(:max))

	# are there ties?
	ties = join(tdf,by(tdf,[:sitc4],num_sec = :sec => length),on=[:sitc4],kind=:left)
	sort(ties[ties.num_sec .> 1,:],[:sitc4])

	# break ties manually
	manual = sort(unique(ties[ties.num_sec .> 1,:sitc4]))
	dfNameSITC[map(x->x ∈ manual,dfNameSITC.sitc),:]
	replacements = Dict(
		"0421" => "3",
		"0741" => "3",
		"2450" => "1",
		"2471" => "1",
		"2640" => "4",
		"2652" => "4",
		"2654" => "4",
		"2655" => "4",
		"2681" => "4",
		"2771" => "2",
		"2873" => "2",
		"3414" => "2",
		"4314" => "1",
		"5983" => "8",
		"6642" => "11",
		"7148" => "12",
		"8212" => "11"
	)
	tdf.sec = map(x->x.sitc4 ∈ manual ? replacements[x.sitc4] : x.sec,eachrow(tdf))

	tdf = unique(tdf)

	# no ties
	all(by(tdf,[:sitc4],num_sec = :sec => length).num_sec .== 1)

	# create sectoral classification for each sitc4 that was classified
	dfSITC4 = join(dfSITC4,tdf[:,[:sitc4,:sec]],on=[:sitc4],kind=:left)

	# the count of hs products in the new bins is the sum, replace
	select!(tdf,Not(:counter))
	rename!(tdf,:total => :counter)

	## refine the map using the most common classification within the 3 digit level
	tdf.sitc3 = map(x->length(x) > 3 ? x[1:3] : x,tdf.sitc4)

	# count the number of hs products within each sitc4-sec pair
	tdf = by(tdf,[:sitc3,:sec],counter = :counter => length)

	# remove sitc3-bin pairs that have fewer products than the most common (including no classification)
	tdf = join(tdf,by(tdf,[:sitc3],total = :counter => sum,max = :counter => maximum),on=[:sitc3],kind=:left)
	tdf = tdf[tdf.counter .== tdf.max,:]

	# no missings
	any(ismissing.(tdf.sec))
	select!(tdf,Not(:max))

	# are there ties?
	ties = join(tdf,by(tdf,[:sitc3],num_sec = :sec => length),on=[:sitc3],kind=:left)
	sort(ties[ties.num_sec .> 1,:],[:sitc3])

	# break ties manually
	manual = sort(unique(ties[ties.num_sec .> 1,:sitc3]))
	dfNameSITC[map(x->x ∈ manual,dfNameSITC.sitc),:]
	replacements = Dict(
		"025" => "3",
		"074" => "3",
		"081" => "3",
		"233" => "8",
		"261" => "4",
		"263" => "1",
		"267" => "8",
		"268" => "4",
		"271" => "8",
		"323" => "7",
		"334" => "7",
		"335" => "7",
		"341" => "7",
		"524" => "7",
		"592" => "8",
		"667" => "2",
		"714" => "11",
		"718" => "11",
		"773" => "11",
		"812" => "13"
	)
	tdf.sec = map(x->x.sitc3 ∈ manual ? replacements[x.sitc3] : x.sec,eachrow(tdf))

	tdf = unique(tdf)

	# no ties
	all(by(tdf,[:sitc3],num_sec = :sec => length).num_sec .== 1)

	# fill in using this mapping
	dfSITC4.sitc3 = map(x->length(x) > 3 ? x[1:3] : x,dfSITC4.sitc4)
	dfSITC4 = join(dfSITC4,rename(tdf[:,[:sitc3,:sec]],:sec => :sec_alt),on=[:sitc3],kind=:left)
	dfSITC4.sec = map(x->ismissing(x.sec) ? x.sec_alt : x.sec,eachrow(dfSITC4))

	select!(dfSITC4,Not([:sitc3,:sec_alt]))

	# twelve remain
	dfSITC4[ismissing.(dfSITC4.sec),:]

	# manually classify remaining

	# break ties manually
	manual = sort(unique(dfSITC4[ismissing.(dfSITC4.sec),:sitc4]))
	dfNameSITC[map(x->x ∈ manual,dfNameSITC.sitc),:]
	replacements = Dict(
		"2200" => "1",
		"2400" => "5",
		"3200" => "7",
		"4200" => "1",
		"6000" => "13",
		"6500" => "4",
		"6750" => "10",
		"7000" => "11",
		"7200" => "11",
		"8000" => "14",
		"9110" => "25",
		"9310" => "25"
	)
	dfSITC4.sec = map(x->x.sitc4 ∈ manual ? replacements[x.sitc4] : x.sec,eachrow(dfSITC4))

	# no missing
	dfSITC4[ismissing.(dfSITC4.sec),:]

	# drop them
	dropmissing!(dfSITC4)
	sort!(dfSITC4,[:sitc4,:sec])

	# convert sec code to integer
	dfSITC4.sec = parse.(Int64,dfSITC4.sec)

	save(dataDir*"sectors/sitc_rev2_4digit_to_sec.csv",dfSITC4)

	return nothing
end

function createWIODtoSec(dataDir::String)
		
	# map from WIOD to WIOD aggregates
	wiodToSecDictionary = Dict(
		1 => 1,
		2 => 2,
		3 => 3,
		4 => 4,
		5 => 4,
		6 => 5,
		7 => 6,
		8 => 7,
		9 => 8,
		10 => 8,
		11 => 9,
		12 => 10,
		13 => 11,
		14 => 12,
		15 => 13,
		16 => 14,
		17 => 15,
		18 => 16,
		19 => 17,
		20 => 17,
		21 => 17,
		22 => 17,
		23 => 18,
		24 => 19,
		25 => 20,
		26 => 18,
		27 => 21,
		28 => 22,
		29 => 23,
		30 => 23,
		31 => 24,
		32 => 24,
		33 => 24,
		34 => 25,
		35 => 25
	)

	df = DataFrame(WIOD=Int64[],sec=Int64[])
	for k in keys(wiodToSecDictionary)
		push!(df,[k wiodToSecDictionary[k]])
	end
	sort!(df,[:WIOD])

	save(dataDir*"wiodToSec.csv",df)

	return nothing
end

function loadAndCombineData(dataDir::String)

	### 1. prepare WIOD data
	
	print("\nLoading and preparing WIOD data...\n\n")
	tstart = time()

	# create file for WIOD to WIOD aggregates mapping
	createWIODtoSec(dataDir)

	# load WIOT data
	wiot = load(dataDir*"wiot/wiot_full.dta") |> @filter(_.row_item <= 35) |> DataFrame
	rename!(wiot,:row_country => :o,:col_country => :d)
	rename!(wiot,:row_item => :sec, :year => :t)

	# aggregate away destination use information
	wiot = by(wiot,[:sec,:o,:d,:t], value = :value => sum)

	# aggregate countries with years with completely missing sector data into RoW
	countryDrop = ["RoW","ROU","TWN","CYP","EST","LTU","LUX","LVA","MLT","IDN"]
	wiot = wiot[map(x->!in(x,countryDrop),wiot.o),:]
	wiot = wiot[map(x->!in(x,countryDrop),wiot.d),:]
	wiot = wiot[map(x->in(x,1999:2007),wiot.t),:]

	# replace sector classification with WIOD aggregates and collapse
	wiodToSec = load(dataDir*"wiodToSec.csv") |> DataFrame
	rename!(wiot,:sec => :WIOD)
	wiot = by(join(wiot,wiodToSec,on=[:WIOD],kind=:left),[:sec,:o,:d,:t],value = :value => sum)

	# no missing sec
	print("Missing in WIOT:\n")
	for var in names(wiot)
		@show var, mean(ismissing.(wiot[!,var]))
	end
	disallowmissing!(wiot)

	# check that the data is square
	wiotCountries = unique(wiot.d)
	wiotSectors = sort(unique(wiot.sec))
	years = sort(unique(wiot.t))
	origins = sort(unique(wiot.o))
	destinations = sort(unique(wiot.d))
	origins == destinations || error("Origins and destinations do not match.")

	# the data is balanced
	all(by(wiot,[:sec,:o,:t],N = :value => length).N .== by(wiot,[:sec,:d,:t],N = :value => length).N)
	all(by(wiot,[:sec,:o,:t],N = :value => length).N .== by(wiot,[:sec,:d,:t],N = :value => length).N)

	tmp = by(wiot,[:o,:d,:t],N = :value => length).N
	minimum(tmp),maximum(tmp)

	minimum(by(wiot,[:sec,:o,:t],N = :value => length).N)
	maximum(by(wiot,[:sec,:o,:t],N = :value => length).N)
	minimum(by(wiot,[:sec,:d,:t],N = :value => length).N)
	maximum(by(wiot,[:sec,:d,:t],N = :value => length).N)

	wiotNegatives = wiot[wiot.value .< 0,:]

	selfTrade = wiotNegatives.o .== wiotNegatives.d
	rowImport = wiotNegatives.d .== "RoW"
	(sum(selfTrade)+sum(rowImport))==nrow(wiotNegatives)

	# variable for cleaned flows
	wiot.secValueWiot = deepcopy(wiot.value)

	# when self-trade is negative, assume that output is incorrect, set to zero (increase output)
	wiot.secValueWiot[(wiot.o .== wiot.d) .& (wiot.value .< 0)] .= 0

	# when RoW imports are negative, assume that output in the origin is incorrect, set to zero (increase output in the origin)
	wiot.secValueWiot[(wiot.d .== "RoW") .& (wiot.value .< 0)] .= 0

	print("WIOD dataset finished. $(timeReport(tstart))\n")
	
	### 2. create and load sitcToSec concordance
	print("\nCreating 4-digit SITC rev. 2 and aggregated WIOD sector concordance...\n\n")
	tstart = time()
	createSitcToSec(dataDir)
	sitcRev2ToSec = load(dataDir*"sectors/sitc_rev2_4digit_to_sec.csv",colparsers=[String]) |> DataFrame
	sitcRev2ToSec[!,:sec] = convert(Vector{Integer},sitcRev2ToSec[!,:sec])
	disallowmissing!(sitcRev2ToSec)
	rename!(sitcRev2ToSec,:sitc4 => :sitc)

	print("Finished. $(timeReport(tstart))\n")

	### 3. Combine trade flow and tariff data
	data = []
	for year in years 

		print("\nCreating 4-digit SITC dataset for $year...\n\n")
		tstart = time()

		# load trade flows
		print("Loading comtrade data...\n")
		@time trade = load(dataDir*"comtrade/$year.csv") |> DataFrame
		select!(trade,Not([:quantity_dot1,:quantity_dot2,:unit]))

		print("Loading tariff data...\n")
		@time tariff = load(dataDir*"tariff/$year.csv") |> DataFrame

		# merge the two data sets
		print("Merging comtrade and tariff data...\n")
		@time df = join(trade,tariff,on=[:year,:iiso,:iwits,:importer,:eiso,:ewits,:exporter,:sitc4],kind=:left)

		# create value variable using importer report first, and exporter report if the former is missing
		df.value = map(x->ismissing(x.value_dot1) ? x.value_dot2 : x.value_dot1,eachrow(df))

		# rename variables and collapse
		print("Collapsing merged data...\n")
		rename!(df,:eiso => :o, :exporter => :origin, :iiso => :d, :importer => :destination, :year => :t,:sitc4 => :sitc)
		@time df = by(df,[:sitc,:o,:d,:t],value = :value => sum,t1 = :t1 => missingOrFirst,t1_pref = :t1_pref => missingOrMinimum)
		sort!(df,[:sitc,:o,:d,:t])

		print("Aggregating RoW...\n")
		for r in eachrow(df)
			oin = in(r.o,wiotCountries)
			din = in(r.d,wiotCountries)
			if !oin || !din
				r.t1 = missing
				r.t1_pref = missing
				if !oin
					r.o = "RoW"
				end
				if !din
					r.d = "RoW"
				end
			end
		end

		# aggregate RoW
		df = by(df,[:sitc,:o,:d,:t],value = :value => sum, t1 = :t1 => first, t1_pref = :t1_pref => first)

		# fix leading zeros
		df.sitc = map(x->"0"^(4-length(x))*x,string.(df.sitc))

		# drop self trade instances
		df = df[df.o .!= df.d,:]

		push!(data,df)

		print("Finished. $(timeReport(tstart))\n")
	end
	data = vcat(data...)

	### 4. Combine 4-digit SITC data with WIOD data

	print("\nCreating combined 4-digit SITC and WIOD dataset.\n\n")
	tstart = time()

	data = join(data,sitcRev2ToSec,on=[:sitc],kind=:left)
	rename!(data,:value=>:valueComtrade)
	print("Missing in trade data:\n")
	for var in names(data)
		@show var, mean(ismissing.(data[!,var]))
	end

	# expand wiot trade data to have all cases for merge with trade
	df = join(wiot,sitcRev2ToSec,on=[:sec],kind=:left)

	all(map(x->in(x,wiotSectors),unique(sitcRev2ToSec.sec)))

	# scale up units to be in thousands of dollars
	df.secValueWiot .*= 1000

	# merge in trade data
	id = map(x->in(x,wiotCountries),data.o) .& map(x->in(x,wiotCountries),data.d)
	df = join(data[id,:],df,on=[:sitc,:sec,:o,:d,:t],kind=:right)

	sumMissing(x) = all(ismissing.(x)) ? missing : sum(skipmissing(x))
	df = join(df,by(df,[:sec,:o,:d,:t],secValueComtrade = :valueComtrade => sumMissing),on=[:sec,:o,:d,:t],kind=:left)

	# sitc level values are missing only if all of sec is missing in comtrade
	df.missingSecValue = ismissing.(df.secValueComtrade)
	df.missingSitcValue = ismissing.(df.valueComtrade)
	@show freqtable(df,:missingSitcValue,:missingSecValue)

	print("Combined 4-digit SITC and WIOD dataset finished. $(timeReport(tstart))\n")

	return df[:,[:sitc,:sec,:o,:d,:t,:valueComtrade,:secValueComtrade,:secValueWiot,:t1,:t1_pref]]
end

function prepareFinalData(df::DataFrame)

	print("\nPreparing final dataset...\n\n")
	tstart = time()

	print("   1. Restricting sample...\n")

    # restrict to first 14 sectors
    df = df[df.sec .<= 14,:]

	print("   2. Inferring missing values versus zeros...\n")

    # variables for infered sitc-level values and sec-level values
    df.value = Union{Missing,Float64}[missing for n=1:nrow(df)]
    df.secValue = Union{Missing,Float64}[missing for n=1:nrow(df)]

    # # conditional on being zero in WIOT, 20.6% has a value in comtrade
    id = df.secValueWiot .== 0
	mean(.!ismissing.(df[id,:valueComtrade]))
    # mean(ismissing.(df[id,:secValueComtrade]))

    # conditional on being zero in WIOT, 21.1% has a value in comtrade
	id = df.secValueWiot .== 0
    mean(.!ismissing.(df[id,:valueComtrade]))

    # if secValueWiot is zero and secValueComtrade is missing, fill with zeros
    subsample = id
    id = subsample .& ismissing.(df.secValueComtrade)
    df.value[id] .= 0
    df.secValue[id] .= 0

    # otherwise use information from comtrade, treating missings as zeros
    id = subsample .& .!ismissing.(df.secValueComtrade)
    df.value[id] = map(x->ismissing(x) ? 0 : x,df.valueComtrade[id])
    df.secValue[id] = df.secValueComtrade[id]

    # Among remaining, if secValueComtrade is missing, use secValueWiot and treat value as missing
    subsample = .!subsample
    id = subsample .& ismissing.(df.secValueComtrade)
    df.value[id] .= missing       # note this is makes no change and can be commented out
    df.secValue[id] = df.secValueWiot[id]

    # Among remaining, if secValueWiot <= secValueComtrade,
    #   then use secValueComtrade and treat missings in valueComtrade as zeros
    subsample = subsample .& .!ismissing.(df.secValueComtrade)
    id = subsample .& (df.secValueWiot .<= df.secValueComtrade)
    df.value[id] = map(x->ismissing(x) ? 0 : x,df.valueComtrade[id])
    df.secValue[id] = df.secValueComtrade[id]

    # Otherwise (if secValueWiot > secValueComtrade),
    #   use secValueWiot and treat missing valueComtrade as missing (no change to missings)
    id = subsample .& (df.secValueWiot .> df.secValueComtrade)
    df.value[id] = df.valueComtrade[id]
    df.secValue[id] =  df.secValueWiot[id]

	mean(ismissing.(df.value))
	mean(map(x->ismissing(x) ? false : x==0,df.value))

	ids = ismissing.(df.value)
	sum(ids .& (df.o .== df.d))/sum(ids)
	sum(ids .& (df.o .== df.d))/sum((df.o .== df.d))
	mean(df.value[ids])

	print("   3. Interpolating tariffs...\n")

	rename!(df,:t1_pref => :t1Pref)

	df.t1Miss = ismissing.(df.t1)
	df.missT1Pref = ismissing.(df.t1Pref)
	# @show freqtable(df,:t1Miss,:missT1Pref)
	# @show unique(df[df.o .== df.d,[:t1Pref]])

    # interpolate tariffs
    df.intTariff = deepcopy(df.t1Pref)
    mean(ismissing.(df.intTariff))
	mean(df.o .== df.d)

	# df[(df.sitc .== sitcCodes.sitc[100]).&(df.d .== "CHN").&(df.t .== 1999),[:sitc,:o,:d,:t,:t1Pref]]

	tmp = by(df,[:sitc,:d,:t],var = :t1Pref => x -> all(ismissing.(x)) ? missing : var(skipmissing(x)))
	tmp.var = map(x->ismissing(x) || isnan(x) ? missing : x,tmp.var)
	ids = map(x->ismissing(x),tmp.var)
	maximum(tmp.var[.!ids])

	sort(tmp.var[.!ids],rev=true)[1:1000]

	percentile(tmp.var[.!ids],40)


	df.missingTariff = ismissing.(df.t1Pref)
	df.missingValue = ismissing.(df.value)

	@show freqtable(df,:missingTariff,:missingValue)

	mean(df.missingValue)
	mean(df.missingTariff)

	mean(df.missingTariff .& df.missingValue)/mean(df.missingTariff)
	mean(df.missingTariff .& .!df.missingValue)/mean(.!df.missingValue)

	# interpolate using minimum within sitc4
    missOrMin(x) = all(ismissing.(x)) ? missing : minimum(skipmissing(x))
    df = join(df,by(df,[:sitc,:d,:t],minTariff4=:t1=>missOrMin),on=[:sitc,:d,:t],kind=:left)
    df.intTariff[ismissing.(df.intTariff)] = df.minTariff4[ismissing.(df.intTariff)]
    mean(ismissing.(df.intTariff))

	# interpolate using minimum within sitc3
    df.sitc3 = map(x->x[1:3],df.sitc)
    df = join(df,by(df,[:sitc3,:d,:t],minTariff3=:t1=>missOrMin),on=[:sitc3,:d,:t],kind=:left)
    df.intTariff[ismissing.(df.intTariff)] = df.minTariff3[ismissing.(df.intTariff)]
	mean(ismissing.(df.intTariff))

	# interpolate using minimum within sitc2
    df.sitc2 = map(x->x[1:2],df.sitc)
    df = join(df,by(df,[:sitc2,:d,:t],minTariff2=:t1=>missOrMin),on=[:sitc2,:d,:t],kind=:left)
    df.intTariff[ismissing.(df.intTariff)] = df.minTariff2[ismissing.(df.intTariff)]
    mean(ismissing.(df.intTariff))

	# interpolate using minimum within sitc1
	df.sitc1 = map(x->x[1:1],df.sitc)
	df = join(df,by(df,[:sitc1,:d,:t],minTariff1=:t1=>missOrMin),on=[:sitc1,:d,:t],kind=:left)
	df.intTariff[ismissing.(df.intTariff)] = df.minTariff1[ismissing.(df.intTariff)]
	mean(ismissing.(df.intTariff))

	# this procedure filled in self-trade tariffs, set to zero
	df.intTariff[df.o .== df.d] .= 0.

    # any(ismissing.(df.intTariff))
    disallowmissing!(df,:intTariff)


	print("   4. Computing sec-level values...\n")

    # create residuals capturing unobserved data
    zeroOrSum(x) = all(ismissing.(x)) ? 0 : sum(skipmissing(x))
    dfSec = by(df,[:sec,:o,:d,:t],secSumSitcValue=:value=>zeroOrSum,secValue=:secValue=>first,
									range=:secValue=>x->maximum(x)-minimum(x),intTariff = :intTariff => minimum)
    any(ismissing.(dfSec.secValue))
	any(ismissing.(dfSec.intTariff))
    any(ismissing.(dfSec.range))

    mean(dfSec.range .== 0)
    dfSec.secResidual = dfSec.secValue - dfSec.secSumSitcValue
	dfSec = dfSec[:,Not([:range,:secSumSitcValue])]
	dfSec.sitc = map(x->"R"*" "^(3-length(x))*x,string.(dfSec.sec))
	dfSec.value = fill(missing,nrow(dfSec))

	# create output
    secData = rename(dfSec[:,[:sec,:o,:d,:t,:secValue,:secResidual]], :secValue => :value, :secResidual => :residual)
    disallowmissing!(secData)
	secData.value# ./ 1000000
    sitcData = 	rename(vcat( df[:,[:sitc,:sec,:o,:d,:t,:value,:intTariff]],dfSec[:,[:sitc,:sec,:o,:d,:t,:value,:intTariff]] ), :intTariff => :tariff)
    disallowmissing!(sitcData,Not(:value))
	sitcData.value# ./ 1000000
    sitcToSec = sort(unique(sitcData[:,[:sitc,:sec]]),[:sitc,:sec])
	disallowmissing!(sitcToSec)

	print("Final data done. $(timeReport(tstart))\n")


    return sitcData,secData,sitcToSec
end


struct TradeData
    sitcCodes::Vector{String}
    secCodes::Vector{Int64}
    countries::Vector{String}
    years::Vector{Int64}
	mValues::Vector{Tuple{String,String,Int64}}
    secMap::BitArray
    secInd::Vector{Int64}
	sitcValuesMiss::BitArray
	sitcValuesZero::BitArray
	sitcValues::Matrix{Float64}
	sitcTariffs::Matrix{Float64}
	logSitcTariffs::Matrix{Float64}
    secValues::Matrix{Float64}
    secResiduals::Matrix{Float64}
    sitcValuesFilled::Matrix{Float64}
end

function TradeData(dataDir::String)

	print("\nPreparing TradeData structure...\n\n")
	tstart = time()

	secData = DataFrame(load(dataDir*"secData.csv"))
	sitcData = DataFrame(load(dataDir*"sitcData.csv"))
	sitcToSec = DataFrame(load(dataDir*"sitcToSec.csv"))
	sitcCodes = DataFrame(load(dataDir*"sitcCodes.csv"))
	sitcCodes = sitcCodes.sitc
	secCodes = sort(unique(secData.sec))
	countries = sort(unique(sitcData.o))
	years = sort(unique(sitcData.t))

    # dimensions of data
    S = length(sitcCodes)
    J = length(secCodes)
    N = length(countries)
    T = length(years)
    M = N*N*T

    # create sector correspondence
    foo = x->in(x[1],sitcToSec.sitc[sitcToSec.sec .== x[2]])
    secMap = foo.(Base.Iterators.product(sitcCodes,secCodes)) .== true
    secInd = [ sitcToSec.sec[sitcToSec.sitc.==sitc][1] for sitc in sitcCodes ]
    any(sum(secMap,dims=2) .!= 1) && error("Sector groups are not disjoint and onto.")
    any( .![ secMap[s,secInd[s]] for s in 1:S ] ) && error("secMap and secInd are inconsistent.")

	secMValues = reshape(map(x->tuple(x...),eachrow(hcat(secData.o,secData.d,secData.t))),(J,M))
	sitcMValues = reshape(map(x->tuple(x...),eachrow(hcat(sitcData.o,sitcData.d,sitcData.t))),(S,M))

	for i=2:size(secMValues,1)
		secMValues[1,:] != secMValues[i,:] && error("Values for o-d-t are not consistent in reshaped secData.")
	end
	for i=2:size(sitcMValues,1)
		secMValues[1,:] != sitcMValues[i,:] && error("Values for o-d-t are not consistent in reshaped sitcData.")
	end
	mValues = secMValues[1,:]

	scaling = 1000 # units in millions of dollars
    secValues = reshape(Array(secData.value),(J,M)) ./ scaling
    secResiduals = reshape(Array(secData.residual),(J,M)) ./ scaling
    sitcValues = reshape(Array(sitcData.value),(S,M)) ./ scaling
    sitcTariffs = reshape(Array(sitcData.tariff),(S,M))

	sitcValuesMiss = ismissing.(sitcValues)
	sitcValuesZero = map(x->ismissing(x) ? false : x==0,sitcValues) .== true
	sitcValues[sitcValuesMiss] .= 0
    sitcValues = convert(Array{Float64,2},sitcValues)
    sitcTariffs = 1 .+ sitcTariffs./ 100

	# # check that residual calculation matches
	# secResidualsAlt = secValues - reshape(secMap'*reshape(sitcValues,(S,:)),(J,N,N,T))
	# secResiduals == secResidualsAlt
	#
	# # check that things add up to true secValues
	# err = (reshape(secMap'*reshape(sitcValues,(S,:)),(J,N,N,T)) + secResiduals) .- secValues
	# maximum(abs,reshape(secMap'*reshape(sitcValues,(S,:)),(J,N,N,T)) + secResiduals - secValues)
	#
	# # errors are small and unrelated to observables
	# id = err .!= 0
	# err[id]
	# @show DataFrame(sec=sort(unique(reshape(Array(secData.sec),(J,N,N,T))[id])))
	# @show DataFrame(o=sort(unique(reshape(Array(secData.o),(J,N,N,T))[id])))
	# @show DataFrame(d=sort(unique(reshape(Array(secData.d),(J,N,N,T))[id])))
	# @show DataFrame(t=sort(unique(reshape(Array(secData.t),(J,N,N,T))[id])))

	sitcValuesFilled = copy(sitcValues)
	fillValues = secMap*secResiduals
	sitcValuesFilled[sitcValuesMiss] .= fillValues[sitcValuesMiss]

	print("Done. $(timeReport(tstart))\n")

    return TradeData(sitcCodes,secCodes,countries,years,mValues,secMap,secInd,sitcValuesMiss,sitcValuesZero,sitcValues,sitcTariffs,log.(sitcTariffs),secValues,secResiduals,sitcValuesFilled)
end


function prepareGravityData(dir::String)

    print("\nPreparing gravity data...\n\n")
    tstart = time()

    print("Loading trade data... \n")

    # load data
    secData = load(dir*"secData.csv") |> DataFrame
    sitcData = load(dir*"sitcData.csv") |> DataFrame
    sitcToSec = load(dir*"sitcToSec.csv") |> DataFrame
    sitcCodes = load(dir*"sitcCodes.csv") |> DataFrame

    S = length(unique(sitcData.sitc))
    J = length(unique(secData.sec))

    # labels and sample range
    countries = sort(unique(sitcData.o))
    years = sort(unique(sitcData.t))
    N = length(countries)
    T = length(years)
    t0 = minimum(years)
    t1 = maximum(years)

    print("Loading cepii data... \n")

    cepii = load(dir*"cepii/cepii.csv") |> DataFrame
    cepii = cepii[:,[:o,:d,:t,:distw,:tdiff,:gdp_o,:gdp_d,:pop_o,:pop_d]]

    print("Preparing gravity dataset... \n")

    sitcData.miss = ismissing.(sitcData.value)
    sitcData.value = map(x->ismissing(x) ? 0. : x,sitcData.value)

    sitcData = join(sitcData,by(sitcData,[:sitc,:sec,:t],sitcTotal = :value => sum),on=[:sitc,:sec,:t],kind=:left)
    sitcData = join(sitcData,by(sitcData,[:sec,:t],secTotal = :value => sum),on=[:sec,:t],kind=:left)
    sitcData.weight = sitcData.sitcTotal ./ sitcData.secTotal
    gravData = join(secData,by(sitcData,[:sec,:o,:d,:t],tariff = [:tariff,:weight] => x-> sum(x.tariff .* x.weight)),on=[:sec,:o,:d,:t],kind=:left)

    gravData = join(gravData,by(gravData,[:d,:t],total = :value => sum),on=[:d,:t],kind=:left)
    gravData.share = gravData.value ./ gravData.total
    gravData.lnValue = map(x->x==0 ? missing : log(x),gravData.value)
    gravData.lnShare = map(x->x==0 ? missing : log(x),gravData.share)
    gravData.lnTariff = log.(1 .+ gravData.tariff./100)

    gravData = join(gravData,cepii,on=[:o,:d,:t],kind=:left)

    # Create covariates
    gravData = join(gravData,rename(gravData[(gravData.o .== "USA").&(gravData.d .== "USA"),[:sec,:t,:distw]],:distw => :distw_USAUSA),on=[:sec,:t],kind=:left)

    gravData.lnDistw = log.(gravData.distw)
    gravData.ΔlnDistw = log.(gravData.distw) - log.(gravData.distw_USAUSA)
    gravData.tdiff2 = gravData.tdiff.^2

    gravData.income_o = gravData.gdp_o ./ gravData.pop_o
    gravData.income_d = gravData.gdp_d ./ gravData.pop_d
    gravData = join(gravData,rename(gravData[gravData.d .== "USA",[:sec,:o,:t,:income_d]],:income_d => :income_USA),on=[:sec,:o,:t],kind=:left)
    gravData.lnIncome_o = log.(gravData.income_o)
    gravData.lnIncome_d = log.(gravData.income_d)
    gravData.ΔlnIncome_o = log.(gravData.income_o) - log.(gravData.income_USA)
    gravData.ΔlnIncome_d = log.(gravData.income_d) - log.(gravData.income_USA)

	# save ordering of countries by income in 1999
	countryOrder = sort(rename(unique(gravData[gravData.t .== 1999,[:d,:income_d]]),:d => :country),:income_d)
	countryOrder.n = 1:N
	sort!(countryOrder,:country)
	countryOrder.income_d = map(x->round(x,digits=3),countryOrder.income_d)
	save(dir*"countryOrder.csv",countryOrder)

    # merge in descriptions and create vector for labeling sectors
    secNames = load(dir*"secNames.csv") |> DataFrame
    gravData = join(gravData,secNames,on=[:sec],kind=:left)


    print("Saving gravity data... \n")

    save(dir*"gravityData.csv",gravData)

	print("Done. Total time: $(timeReport(tstart))\n")

    return nothing
end

function loadGravityData(dir::String)
    secData = load(dir*"gravityData.csv") |> DataFrame
    J = length(unique(secData.sec))
    countries = sort(unique(secData.o))
    years = sort(unique(secData.t))
    N = length(countries)
    T = length(years)
    t0 = minimum(years)
    t1 = maximum(years)

    secNames = load(dir*"secNames.csv") |> DataFrame
    secLabels = secNames.name

    secData.secStr = [ " "^(2-length(x))*x for x in string.(secData.sec) ]

    return secData,J,countries,years,N,T,t0,t1,secNames,secLabels
end

struct ProjData
    trade::TradeData
    sec::DataFrame
    secNames::DataFrame
    grav::DataFrame
    J::Int64
    N::Int64
    T::Int64
    countries::Vector{String}
    years::Vector{Int64}

    function ProjData(dirs::ProjDirectories)

		print("\nPreparing project data structre...\n\n")

        # load data
        td = TradeData(dirs.data)
        secData,J,countries,years,N,T,t0,t1 = loadSecData(dirs.data)
        prepareGravityData(dirs.data)  # creates data/gravityData.csv
        gravData,J,countries,years,N,T,t0,t1,secNames,secLabels = loadGravityData(dirs.data)

        return new(td,secData,secNames,gravData,J,N,T,countries,years)
    end
end