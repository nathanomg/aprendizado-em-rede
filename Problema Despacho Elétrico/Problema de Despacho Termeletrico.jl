using JuMP
using Plots
using GLPK

T          = 24
S          = 10

wind_scen   = ones(24,10) .* transpose(rand(Float64, S)) .* 10 .* (1 .+ 1/2 .* cos.(-2*pi .* [1:1:T;]/24))
hidr_scen   = 10 * (1 .+ 1/2 .* sin.(-2*pi .* [1:1:T;]/24)) 

therm1_gen = fill(10, T) 
therm2_gen = fill(10, T) 
therm3_gen = fill(10, T) 

capacity = [therm1_gen, therm2_gen, therm3_gen] 
nT       = size(capacity)[1] 
nW       = 1 
nH       = 1 
produtib = [0.7] 
F = 2

cost     = [10, 20, 30] 
demand   = [(5 * (1 .+ 1/2 .* ones(24))), (5 * (1 .+ 1/2 .* sin.(-2*pi .* [1:1:T;]/24)))] 

disc_rate = 0.06
deficit_cost = 7000 

model = JuMP.Model(GLPK.Optimizer) 


@variable(model, thermo_gen[1:T, 1:nT, 1:S]   >= 0) 
@variable(model, wind_gen[1:T, 1:nW, 1:S]     >= 0) 
@variable(model, hidr_gen[1:T, 1:nH, 1:S]     >= 0) 
@variable(model, hidr_curtail[1:T, 1:nW, 1:S] >= 0) 
@variable(model, wind_curtail[1:T, 1:nW, 1:S] >= 0) 
@variable(model, deficit[1:T, 1:2, 1:S]       >= 0)
@variable(model, gen_tot[1:T, 1:S]            >= 0) 
@variable(model, flow_se_s[1:T, 1:S]          >= 0)
@variable(model, flow_s_se[1:T, 1:S]          >= 0)

#ToDo: custo da reserva
@objective(model, Min, sum(((1+disc_rate)^(t) * sum(cost .* thermo_gen[t, :, :]) + sum(deficit[t,:, :])*deficit_cost ) for t in 1:T))

@constraint(model, balanco_energia_se[t = 1:T, s=1:S]           , deficit[t,1, s] + flow_s_se[t, s]                                      .== demand[1][t] + flow_se_s[t, s]) 
@constraint(model,                   [t = 1:T, s=1:S]           , sum(thermo_gen[t,:, s]) + sum(wind_gen[t,:, s]) + sum(hidr_gen[t,:, s]) + deficit[t,2, s] + flow_se_s[t, s] .== demand[2][t] + flow_s_se[t, s]) 
@constraint(model,        max_flow_se[t = 1:T, s=1:S]           , flow_se_s[t, s]                                                         <= F) 
@constraint(model,         max_flow_s[t = 1:T, s=1:S]           , flow_s_se[t, s]                                                         <= F)
@constraint(model,         limite_cap[t = 1:T, i = 1:nT, s=1:S] , thermo_gen[t,i, s]  + thermo_res[t,i, s]                                <= capacity[i][t]) 
@constraint(model,         wind_const[t = 1:T, i = 1:nW, s=1:S] , wind_gen[t,i, s] + wind_curtail[t,i, s]                                 == wind_scen[t, s])
@constraint(model,         hidr_const[t = 1:T, i = 1:nH, s=1:S] , hidr_gen[t,i, s] + hidr_curtail[t,i, s]                                 == hidr_scen[t]*produtib[i])
@constraint(model,        max_flow_se[t = 1:T, s=1:S]           , sum(thermo_res[t,:,s]) >= sum(demand[2][t])*0.05 )
@constraint(model,        max_flow_se[t = 1:T+1, s=1:S]         , v[]  =  v + ena - geracao - vertimento)

#ToDo: modelagem reservatorio
#ToDo: flow[t, s, de, pa] 
#media reservatorios no tempo. v ini = max; vol <> infinito
#periodo de 60 meses e resultados em 1, 2, 5, 10



optimize!(model)
solution_summary(model)

plot(JuMP.dual.(balanco_energia_se[:,:]))
plot(JuMP.value.(thermo_gen[:,1:nT, 1]))
plot(JuMP.dual.(max_flow_s[:,:]))
plot!(JuMP.value.(thermo_gen[:,1:nT, 2]))
plot!(demand) 

#plot(JuMP.value.(generation[:,3]))
#objective_function(model)
#objective_value(model)
#plot(JuMP.dual.(balanco_cap[:,2]))
#plot!(JuMP.value.(wind_gen[:,1:nW]), label="wind_gen")
#plot!(JuMP.value.(hidr_gen[:,1:nW]), label="hidr_gen")
#plot(demand)
#plot!(JuMP.value.(gen_tot), label="wind_gen")
#plot!(JuMP.value.(deficit))
#plot!(JuMP.value.(wind_curtail[:,1:nW]))
#plot!(JuMP.value.(hidr_curtail[:,1:nW]))
#plot(JuMP.value.(deficit))
#plot(JuMP.value.(max_flow_s) .- demand[1])
#plot(JuMP.value.(deficit))
