using JuMP
using Plots
using GLPK

T          = 60
S          = 5
nBus       = 4

flow_limits       = [ 0 10 20 15
                     10  0  0  0
                     20  0  0 40
                     15  0 40  0]

hidr_reser_limits = [10, 5, 10, 2]
hidr_reser_ini    = [ 0, 0,  0, 0]

#hidr_reser_ini    = hidr_reser_limits


wind_scen   = ones(T,S, nBus) .* transpose(rand(Float64, S)) .* 10 .* (1 .+ 1/2 .* cos.(-2*pi .* [1:1:T;]/T))
hidr_scen   = ones(T,S, nBus) .* transpose(rand(Float64, S)) .* 10 .* (1 .+ 1/2 .* sin.(-2*pi .* [1:1:T;]/T))


capacity = [[fill(S, T), fill(S, T), fill(S, T)],
            [fill( 0, T),  fill( 0, T),  fill( 0, T)],
            [fill( 0, T),  fill( 0, T),  fill( 0, T)],
            [fill(S, T), fill(S, T), fill(S, T)]]
nT       = 3 
nW       = 1 
nH       = 1 
produtib = [0.7] 
F = 2

cost     = [10, 20, 30] 
demand   = [(20 * (1 .+ 1/2 .* ones(T))), 
            (30 * (1 .+ 1/2 .* sin.(-2*pi .* [1:1:T;]/T))), 
            (15 * (1 .+ 1/2 .* sin.(-2*pi .* [1:1:T;]/T))), 
            (5 * (1 .+ 1/2 .* sin.(-2*pi .* [1:1:T;]/T)))] 

disc_rate = 0.06
deficit_cost = 7000 

model = JuMP.Model(GLPK.Optimizer) 


@variable(model, thermo_gen[1:T, 1:nT, 1:S, 1:nBus]   >= 0) 
@variable(model, wind_gen[1:T, 1:nW, 1:S, 1:nBus]     >= 0) 
@variable(model, hidr_gen[1:T, 1:nH, 1:S, 1:nBus]     >= 0) 
@variable(model, hidr_res[1:T, 1:nH, 1:S, 1:nBus]     >= 0) 
@variable(model, hidr_curtail[1:T, 1:nW, 1:S, 1:nBus] >= 0) 
@variable(model, wind_curtail[1:T, 1:nW, 1:S, 1:nBus] >= 0) 
@variable(model, deficit[1:T, 1:nBus, 1:S]          >= 0)
@variable(model, flow[1:T, 1:S, 1:nBus, 1:nBus]       >= 0)

@objective(model, Min, sum(((1+disc_rate)^(t) * sum(cost .* thermo_gen[t, :, :,:]) + sum(deficit[t,:, :])*deficit_cost ) for t in 1:T) +
                       sum(((1+disc_rate)^(t) * sum(cost .* thermo_gen[t+1, :, :,:]) + sum(deficit[t+1,:, :])*deficit_cost ) for t in 1:T-1))

@constraint(model,    balanco_energia[t = 1:T, s=1:S, b=1:nBus] , sum(thermo_gen[t,:, s, b]) + sum(wind_gen[t,:, s, b]) + sum(hidr_gen[t,:, s, b]) + deficit[t,b,s] + sum(flow[t, s, :, b]) .== demand[b][t] + sum(flow[t, s, b, :])) 
@constraint(model, max_flow[t = 1:T, s=1:S, o=1:nBus, d=1:nBus] , flow[t, s, o, d]                                        <= flow_limits[o, d]) 
@constraint(model,         limite_cap[t = 1:T, i = 1:nT, s=1:S, b=1:nBus] , thermo_gen[t,i, s, b]                        <= capacity[b][i][t]) 
@constraint(model,         wind_const[t = 1:T, i = 1:nW, s=1:S, b=1:nBus] , wind_gen[t,i, s, b] + wind_curtail[t,i, s, b] == wind_scen[t, s, b])
#@constraint(model,         hidr_const[t = 1:T, i = 1:nH, s=1:S, b=1:nBus] , hidr_gen[t,i, s, b] + hidr_curtail[t,i, s, b] <= hidr_scen[t, s, b]*produtib[i])
@constraint(model,         hidr_reserv_ini[i = 1:nH, s=1:S, b=1:nBus]     , hidr_res[1, i, s, b]  ==  hidr_reser_ini[b])
@constraint(model,         max_hidr_reserv[t = 1:T, i = 1:nH, s=1:S, b=1:nBus]     , hidr_res[t, i, s, b]  <=  hidr_reser_limits[b])
@constraint(model,         hidr_const[t = 1:T-1, i = 1:nH, s=1:S, b=1:nBus]    , hidr_res[t+1, i, s, b]  ==  hidr_res[t, i, s, b] + hidr_scen[t, s, b]*produtib[i] - hidr_gen[t,i, s, b] - hidr_curtail[t,i, s, b])

#@constraint(model,        max_flow_se[t = 1:T, s=1:S]           , sum(thermo_res[t,:,s]) >= sum(demand[2][t])*0.05 )
#@constraint(model, balanco_energia[t = 1:T, s=1:S, b=1:nBus]   , deficit[t,b, s] + sum(flow[t, s, :, o])   .== demand[b][t] + sum(flow[t, s, o, :])) 

#ToDo: custo da reserva
#ToDo: modelagem reservatorio
#ToDo: flow[t, s, de, pa] 
#media reservatorios no tempo. v ini = max; vol <> infinito
#periodo de 60 meses e resultados em 1, 2, 5, 10



optimize!(model)
#solution_summary(model)

#plot(JuMP.dual.(balanco_energia[:,:, 1]))
#plot(JuMP.value.(thermo_gen[:,1:nT, 5,2]))
#plot(demand[2,:])

#plot(JuMP.value.(flow[:,3, :, 1]))

#plot(JuMP.value.(flow[:,3, 1, :]))
#plot(JuMP.value.(flow[:,1, :, 2]))


plot(JuMP.value.(hidr_res[:,1, :, 3]))

#plot!(JuMP.value.(thermo_gen[:,1:nT, 2]))
#plot!(demand) 

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
