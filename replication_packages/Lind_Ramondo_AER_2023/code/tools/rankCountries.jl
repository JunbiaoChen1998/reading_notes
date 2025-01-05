function rankCountries(dirs::ProjDirectories,df::DataFrame)
    countryOrder = load(dirs.data*"countryOrder.csv") |> DataFrame
    countryOrder = countryOrder[:,[:n,:country]]
    sort!(countryOrder,:n)

    out = rename(copy(df),:o => :origin,:d => :destination)
    out = join(out,rename(countryOrder,:n=>:o,:country=>:origin),on=[:origin],kind=:left)
    out = join(out,rename(countryOrder,:n=>:d,:country=>:destination),on=[:destination],kind=:left)
    return out,countryOrder.country
end
