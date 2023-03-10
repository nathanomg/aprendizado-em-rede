using JuMP
using Plots
using GLPK

#=
1) Modelagem de déficit 
2) Modelagem de eólicas com cenário único 
3) Modelagem de eólicas com 10 cenários (pode ser com cenários aleatórios) 
4) Modelagem de hidrelétricas a fio d'água com 10 cenários
=#


T          = 24

wind_scen   = 5 * (1 .+ 1/2 .* cos.(-2*pi .* [1:1:T;]/24))
hidr_scen   = 10 * (1 .+ 1/2 .* sin.(-2*pi .* [1:1:T;]/24))
#wind_scen  = 15*rand(Float64, T)

therm1_gen = fill(10, T)
therm2_gen = fill(10, T)
therm3_gen = fill(10, T)

capacity = [therm1_gen, therm2_gen, therm3_gen]
nT       = size(capacity)[1]
nW       = 1
nH       = 1
produtib = [0.7]

cost     = [10, 20, 30]
demand   = 5 * (1 .+ 1/2 .* sin.(-2*pi .* [1:1:T;]/24))
disc_rate = 0.06
deficit_cost = 7000

model = JuMP.Model(GLPK.Optimizer)


@variable(model, thermo_gen[1:T, 1:nT]   >= 0)
@variable(model, wind_gen[1:T, 1:nW]     >= 0)
@variable(model, hidr_gen[1:T, 1:nH]     >= 0)
@variable(model, hidr_curtail[1:T, 1:nW] >= 0)
@variable(model, wind_curtail[1:T, 1:nW] >= 0)
@variable(model, deficit[1:T]            >= 0)
@variable(model, gen_tot[1:T]            >= 0)

@objective(model, Min, sum(( (1+disc_rate)^(t) * sum(cost .* thermo_gen[t, :] ) + sum(deficit[t])*deficit_cost ) for t in 1:T) )

@constraint(model, balanco_energia[t = 1:T]     , gen_tot[t] + sum(deficit[t]) .== demand[t])
@constraint(model, limite_cap[t = 1:T, i = 1:nT], thermo_gen[t,i] <= capacity[i][t])
@constraint(model, max_variacao[t = 2:T]        , thermo_gen[t,3] - thermo_gen[t-1,3] <= 1.)
@constraint(model, wind_const[t=1:T, i = 1:nW]  , wind_gen[t,i] + wind_curtail[t,i] == wind_scen[t])
@constraint(model, hidr_const[t=1:T, i = 1:nH]  , hidr_gen[t,i] + hidr_curtail[t,i] == hidr_scen[t]*produtib[i])
@constraint(model, generation[t=1:T]            , sum(thermo_gen[t,:]) + sum(wind_gen[t,:]) + sum(hidr_gen[t,:]) == gen_tot[t])

#@constraint(model, wind_const[t=1:T]            , wind_gen[t,nW] + wind_curtail[t,nW] == wind_scen[t])

optimize!(model)
solution_summary(model)

#plot(JuMP.value.(generation[:,3]))
#objective_function(model)
#objective_value(model)
#plot(JuMP.dual.(balanco_demanda))
plot(JuMP.value.(thermo_gen[:,1:nT]))
#plot(JuMP.dual.(balanco_cap[:,2]))
plot!(demand)
plot!(JuMP.value.(wind_gen[:,1:nW]), label="wind_gen")
plot!(JuMP.value.(hidr_gen[:,1:nW]), label="hidr_gen")

plot(demand)
plot!(JuMP.value.(gen_tot), label="wind_gen")
plot!(JuMP.value.(deficit))

plot!(JuMP.value.(wind_curtail[:,1:nW]))
plot!(JuMP.value.(hidr_curtail[:,1:nW]))

#plot(JuMP.value.(deficit))
