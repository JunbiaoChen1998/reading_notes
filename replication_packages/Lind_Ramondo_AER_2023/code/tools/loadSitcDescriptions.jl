
function loadSitcDescriptions(dirs::ProjDirectories,data::ProjData)
    # sitc names
    sitcDescriptions = load(dirs.data*"sectors/S2_descriptions.csv") |> DataFrame
    rename!(sitcDescriptions,Symbol("Commodity Code")=>:code)
    rename!(sitcDescriptions,Symbol("Commodity description")=>:description)
    push!(sitcDescriptions,["R" "Residual"])
    push!(sitcDescriptions,["R " "Residual"])
    for j=1:data.J
        jstr = string(j)
        push!(sitcDescriptions,["R"*" "^(3-length(jstr))*jstr "Residual for sec = $jstr"])
    end
    firstCodes = sectorLabels(DataFrame(sitc=data.trade.sitcCodes),sitcDescriptions)
    return sitcDescriptions,firstCodes
end
