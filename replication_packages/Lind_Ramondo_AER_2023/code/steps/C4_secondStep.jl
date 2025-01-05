
function estimateSecondStep(dirs::ProjDirectories,data::ProjData)

    # load estimates and labels
    σ,σse,dfΛ,dfΦ,K,dfCES,dfSGM = loadEstimates(dirs)

    df = copy(dfΦ)

    # second step
    df = join(df,by(df,[:d,:t],Xdt = :Xkodt => sum),on=[:d,:t],kind=:left)
    df = join(df,by(df,[:k,:d,:t],Xkdt = :Xkodt => sum),on=[:k,:d,:t],kind=:left)
    df.πkodt = df.Xkodt ./ df.Xdt
    df = join(df,by(df,[:o,:d,:t],πodt = :πkodt => sum),on=[:o,:d,:t],kind=:left)
    df = join(df,by(df,[:k,:d,:t],πkdt = :πkodt => sum),on=[:k,:d,:t],kind=:left)

    # Estimate θ
    df.lnπkodt = map(x-> isinf(x) ? missing : x,log.(df.πkodt))
    df.lnπkdt = map(x-> isinf(x) ? missing : x,log.(df.πkdt))


    df.kdt = string.(df.k).*"_".*df.d.*"_".*string.(df.t)
    df.kd = string.(df.k).*"_".*df.d
    df.σk = map(k->σ[k],df.k)
    df.Zkodt = map(x->isinf(x) ? missing : x,.1*log.(df.Xkodt ./ df.Xkdt)./df.σk)

    θ = minimum(σ)

    varNames = ["log(tkodt)"]
    feNames = ["Zkodt","k & Zkodt","k & o & t","d & t","o & d","k & o & d"]
    vcovSpec = Vcov.cluster(:kd)
    rt = RegressionTable(varNames,feNames)

    m = @formula πkdt ~ log(tkodt) + fe(k)&fe(o)&fe(t) + fe(d)&fe(t) + fe(o)&fe(d) + Zkodt
    @show result = nlreg(df, m, Poisson(), GLFixedEffectModels.LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8,save=true)
    R = permutedims(hcat([ result.coefnames .== name  for name in ["log(tkodt)"] ]...))
    r = [-θ]
    rt = RegressionTable(rt,result,[true,false,true,true,true,false],testExternal=["Baseline",round.(waldTest(result,R,r),digits=3)...])

    m = @formula πkdt ~ log(tkodt) + fe(k)&fe(o)&fe(t) + fe(d)&fe(t) + fe(k)&fe(o)&fe(d) + Zkodt
    @show result = nlreg(df, m, Poisson(), GLFixedEffectModels.LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8,save=true)
    R = permutedims(hcat([ result.coefnames .== name  for name in ["log(tkodt)"] ]...))
    r = [-θ]
    rt = RegressionTable(rt,result,[true,false,true,true,false,true],testExternal=["Baseline",round.(waldTest(result,R,r),digits=3)...])

    m = @formula πkdt ~ log(tkodt) + fe(k)&fe(o)&fe(t) + fe(d)&fe(t) + fe(k)&fe(o)&fe(d) + fe(k)&Zkodt
    @show result = nlreg(df, m, Poisson(), GLFixedEffectModels.LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8,save=true)
    R = permutedims(hcat([ result.coefnames .== name  for name in ["log(tkodt)"] ]...))
    r = [-θ]
    @show result = nlreg(df, m, Poisson(), GLFixedEffectModels.LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8,save=true)
    rt = RegressionTable(rt,result,[false,true,true,true,false,true],testExternal=["Baseline",round.(waldTest(result,R,r),digits=3)...])

    @show table = createTable(rt)

    table = createTable(rt)
    save(dirs.results*"onlineAppendix/TableO4.csv",table)
    io = open(dirs.results*"onlineAppendix/TableO4.tex","w")
    write(io,latexify(table,env=:table,latex=false))
    close(io)



    # SGM second step
    secData = copy(data.grav)
    rename!(secData,:value => :Xjodt,:share => :πjodt)
    secData = join(secData,by(secData,[:sec,:d,:t],Xjdt = :Xjodt => sum),on=[:sec,:d,:t],kind=:left)
    secData = join(secData,by(secData,[:sec,:d,:t],πjdt = :πjodt => sum),on=[:sec,:d,:t],kind=:left)

    # Estimate θ
    secData.jdt = string.(secData.sec).*"_".*secData.d.*"_".*string.(secData.t)
    secData.jd = string.(secData.sec).*"_".*secData.d
    secData.σj = map(j->-dfSGM.σ[j],secData.sec)
    secData.Zjodt = map(x->isinf(x) ? missing : x,.1*log.(secData.Xjodt ./ secData.Xjdt)./secData.σj)


    varNames = ["lnTariff"]
    feNames = ["Zjodt","j & Zjodt","j & o & t","d & t","o & d","j & o & d"]
    vcovSpec = Vcov.cluster(:jd)
    rt = RegressionTable(varNames,feNames)

    m = @formula πjdt ~ lnTariff + fe(sec)&fe(o)&fe(t) + fe(d)&fe(t) + fe(o)&fe(d) + Zjodt
    @show result = nlreg(secData, m, Poisson(), GLFixedEffectModels.LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8,save=true)
    R = permutedims(hcat([ result.coefnames .== name  for name in ["lnTariff"] ]...))
    r = [-θ]
    rt = RegressionTable(rt,result,[true,false,true,true,true,false],testExternal=["Baseline",round.(waldTest(result,R,r),digits=3)...])

    m = @formula πjdt ~ lnTariff + fe(sec)&fe(o)&fe(t) + fe(d)&fe(t) + fe(sec)&fe(o)&fe(d) + Zjodt
    @show result = nlreg(secData, m, Poisson(), GLFixedEffectModels.LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8,save=true)
    R = permutedims(hcat([ result.coefnames .== name  for name in ["lnTariff"] ]...))
    r = [-θ]
    rt = RegressionTable(rt,result,[true,false,true,true,false,true],testExternal=["Baseline",round.(waldTest(result,R,r),digits=3)...])

    m = @formula πjdt ~ lnTariff + fe(sec)&fe(o)&fe(t) + fe(d)&fe(t) + fe(sec)&fe(o)&fe(d) + fe(sec)&Zjodt
    @show result = nlreg(secData, m, Poisson(), GLFixedEffectModels.LogLink(),vcovSpec,maxiter=2000,rho_tol=1e-8,save=true)
    R = permutedims(hcat([ result.coefnames .== name  for name in ["lnTariff"] ]...))
    r = [-θ]
    rt = RegressionTable(rt,result,[false,true,true,true,false,true],testExternal=["Baseline",round.(waldTest(result,R,r),digits=3)...])


    @show table = createTable(rt)

    return nothing
end
