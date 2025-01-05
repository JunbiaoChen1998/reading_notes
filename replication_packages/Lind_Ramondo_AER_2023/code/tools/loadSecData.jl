

function loadSecData(dir::String)
    secData = load(dir*"secData.csv") |> DataFrame
    # sitcData = load(dataDir*"sitcData.csv") |> DataFrame
    # sitcToSec = load(dataDir*"sitcToSec.csv") |> DataFrame
    # sitcCodes = load(dataDir*"sitcCodes.csv") |> DataFrame

    # S = length(unique(sitcData.sitc))
    J = length(unique(secData.sec))

    # labels and sample range
    countries = sort(unique(secData.o))
    years = sort(unique(secData.t))
    N = length(countries)
    T = length(years)
    t0 = minimum(years)
    t1 = maximum(years)

    return secData,J,countries,years,N,T,t0,t1
end
