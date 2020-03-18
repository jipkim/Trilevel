using Pkg, DataFrames, Distributions, DelimitedFiles, CSV
using JuMP, KNITRO

s = 2
regset = [s]
timeset = 1:24
peaktime = 13:21
PdownMax = 2400
UserGap = 1e-2
η = 0.03
RegGap = 1
σ0_pv = 20e-2
σ0_wind = 20e-2/2
υ0 = 0.0
ScaleObj = 1e-5
RPStune = 1
α0inv = 1/0.25
π0 = 20


include("StructType.jl")
# Dπ_LB, Dπ_UB, τe_LB, τe_UB, τc_LB, τc_UB
bound = bounds(0, 100, 0, 0, 0, 0)

include("MPini.jl")
include("MP.jl")
include("SP.jl")

include("NetworkDataType.jl")
include("NetworkLoad.jl")

include("Demand24Load.jl")
filename_Load24 = string(pwd(),"/data/load24data-PJM.csv")
loadcurve = Demand24Load(filename_Load24)

include("GenCurveLoad.jl")
filename_windcurve = string(pwd(),"/data/hourly_wind_gen_2018.csv")
filename_pvcurve = string(pwd(),"/data/hourly_solar_gen_2018.csv")

WTcf = 0.262
PVcf = 0.145

temp_pv = GenCurveLoad(filename_pvcurve)
pvcurve_year = temp_pv / mean(temp_pv) * PVcf
temp_wind = GenCurveLoad(filename_windcurve)
windcurve_year = temp_wind / mean(temp_wind) * WTcf

pvcurve = sum(pvcurve_year, dims=1)/365
windcurve = sum(windcurve_year, dims=1)/365

testsystem = "ISONE8busTN"

stateset = 1:6
loadzone = Dict()
loadzone[1] = 1
loadzone[2] = 2
loadzone[3] = 3
loadzone[4] = [4, 5, 8]
loadzone[5] = 6
loadzone[6] = 7


buses, lines, generators, datamat = NetworkLoad(testsystem)
lineset = 1:length(lines)
busset = 1:length(buses)
genset = 1:length(generators)
Dp = datamat["bus"][:,2]


windset = Int64[]
pvset = Int64[]
hydroset = Int64[]
nucset = Int64[]
coalset = Int64[]
oilset = Int64[]
gasset = Int64[]

for g in genset
   if datamat["gen"][g,2] == "Wind"
      push!(windset, generators[g].gindex)
   elseif datamat["gen"][g,2] == "PV"
      push!(pvset, generators[g].gindex)
   elseif datamat["gen"][g,2] == "Hydro"
      push!(hydroset, generators[g].gindex)
   elseif datamat["gen"][g,2] == "Nuclear"
      push!(nucset, generators[g].gindex)
   elseif datamat["gen"][g,2] == "Coal"
      push!(coalset, generators[g].gindex)
   elseif datamat["gen"][g,2] == "Oil"
      push!(oilset, generators[g].gindex)
   elseif datamat["gen"][g,2] == "Gas"
      push!(gasset, generators[g].gindex)
   end
end

for g in genset
   generators[g].Pmin = 0
end

genset_hat = genset[end]+1 : genset[end]+3*length(busset)
for i in 1:length(busset)
   push!(windset, genset[end]+3*(i-1)+1)
   push!(pvset, genset[end]+3*(i-1)+2)
   push!(gasset, genset[end]+3*(i-1)+3)
end

B_g = []
for g in genset
   push!(B_g, generators[g].location)
end


for i in busset
   for j in 1:3
      push!(B_g, i)
   end
end

B_gn = [Int64[] for ii in busset]
for g in union(genset,genset_hat)
   push!(B_gn[B_g[g]], g)
end

S_g = Dict()
for s in stateset
   for z in loadzone[s]
      for g in B_gn[z]
         S_g[g] = s
      end
   end
end


S_gs = [Int64[] for ii in stateset]
for s in stateset, n in loadzone[s], g in B_gn[n]
   push!(S_gs[s], g)
end


# Investment cost source:  https://www.eia.gov/electricity/generatorcosts/
Cinv = Dict()

DCRF = 3.5481e-4 * length(timeset) / 24 # h=10years / 5%

Cinv["Wind"] = DCRF * 1630 * 1e3
Cinv["PV"] = DCRF * 2434 * 1e3
Cinv["Hydro"] = DCRF * 5312 * 1e3
Cinv["Nuclear"] = DCRF * 1000 * 1e3
Cinv["Coal"] = DCRF * 1000 * 1e3
Cinv["Oil"] = DCRF * 1672 * 1e3
Cinv["Gas"] = DCRF * 895 * 1e3

Coper = Dict()
Coper["Wind"] = [0, 2, 0]
Coper["PV"] = [0, 3, 0]
Coper["Hydro"] = [0, 4, 0]
Coper["Nuclear"] = [0, 7, 1250]
Coper["Coal"] = [0.0, 19.98, 1535.0]
Coper["Oil"] = [0.0, 192.06, 5000.0]
Coper["Gas"] = [0.0, 20.0, 612.84]



Cinvg = Vector{Float64}(undef, length(union(genset,genset_hat)))
Coperg = Vector{Vector{Float64}}(undef, length(union(genset,genset_hat)))
for g in genset_hat
   if g in windset
      Cinvg[g] = Cinv["Wind"]
      Coperg[g] = Coper["Wind"]
   elseif g in pvset
      Cinvg[g] = Cinv["PV"]
      Coperg[g] = Coper["PV"]
   elseif g in hydroset
      Cinvg[g] = Cinv["Hydro"]
      Coperg[g] = Coper["Hydro"]
   elseif g in nucset
      Cinvg[g] = Cinv["Nuclear"]
      Coperg[g] = Coper["Nuclear"]
   elseif g in coalset
      Cinvg[g] = Cinv["Coal"]
      Coperg[g] = Coper["Coal"]
   elseif g in oilset
      Cinvg[g] = Cinv["Oil"]
      Coperg[g] = Coper["Oil"]
   elseif g in gasset
      Cinvg[g] = Cinv["Gas"]
      Coperg[g] = Coper["Gas"]
   end
end


for g in genset
   if g in windset
      Cinvg[g] = Cinv["Wind"]
      Coperg[g] = Coper["Wind"]
   elseif g in pvset
      Cinvg[g] = Cinv["PV"]
      Coperg[g] = Coper["PV"]
   elseif g in hydroset
      Cinvg[g] = Cinv["Hydro"]
      Coperg[g] = Coper["Hydro"]
   elseif g in nucset
      Cinvg[g] = Cinv["Nuclear"]
      Coperg[g] = generators[g].cost
   elseif g in coalset
      Cinvg[g] = Cinv["Coal"]
      Coperg[g] = generators[g].cost
   elseif g in oilset
      Cinvg[g] = Cinv["Oil"]
      Coperg[g] = generators[g].cost
   elseif g in gasset
      Cinvg[g] = Cinv["Gas"]
      Coperg[g] = generators[g].cost
   end
end

ρ = 1.0*ones(length(union(genset, genset_hat)),timeset[end])
for g in pvset
   for t in timeset
      ρ[g,t] = pvcurve[1,t]
   end
end
for g in windset
   for t in timeset
      ρ[g,t] = windcurve[1,t]
   end
end

σ = zeros(length(union(genset, genset_hat)))
for g in pvset
   σ[g] = σ0_pv
end
for g in windset
   σ[g] = σ0_wind
end

υ = υ0*ones(length(union(genset, genset_hat)), timeset[end])
inv_ϕ = quantile.(Normal(), 1-η)
#-------------------------------------#
kappa = [1; 0.252; 0.75; 0.411; 0.44; 0.385]
kappa[s] = RPStune * kappa[s]

InvB = 1e13  * ones(length(stateset)) * DCRF
PolicyB = 1e8  * ones(length(stateset))
#-------------------------------------#
Γ = 0
I_R = intersect(genset, union(windset,pvset))
I_R_hat = intersect(genset_hat, union(windset,pvset))
I_C = setdiff(genset, union(windset,pvset))
I_C_hat = setdiff(genset_hat,union(windset,pvset))

α = zeros(length(union(genset, genset_hat)))
for nnn in loadzone[s], i in intersect(S_gs[s], I_C)
   α[i] = generators[i].Pmax/sum(generators[j].Pmax for j in intersect(S_gs[s], I_C))
end

time_M = @elapsed tempM = MPini(testsystem, s, kappa, PolicyB, timeset, υ0, UserGap, bound, peaktime)
π0 = ((value.(tempM[:Dπ_on])) * length(peaktime) + (value.(tempM[:Dπ_off])) * length(setdiff(timeset,peaktime)))/length(timeset)

time_S = @elapsed global tempS = SP(testsystem, s, kappa, timeset, PdownMax, η, σ0_pv, σ0_wind, υ0, value.(tempM[:Dπ]), value.(tempM[:tau_e]), value.(tempM[:tau_c]), RegGap, UserGap, α0inv, π0)

Niter = 100
conv = 777*ones(Niter,4)
time_iter = zeros(Niter,2)
MPbound_iter = zeros(Niter)
SPbound_iter = zeros(Niter)

# stack all ML results as a cut variable
Dπ_iter = 777*ones(Niter, length(busset), length(timeset))
tau_e_iter = 777*ones(Niter, length(stateset))
tau_c_iter = 777*ones(Niter, length(stateset))
ObjUL_iter = 777*ones(Niter)
ObjML_iter = 777*ones(Niter)
ObjLL_iter = 777*ones(Niter)

gpmax_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:gpmax])}}(undef,Niter)
gpo_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:gpo])}}(undef,Niter)
gp_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:gp])}}(undef,Niter)
lambda_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:lambda])}}(undef,Niter)
dp_iter = Vector{JuMP.Containers.DenseAxisArray{Float64,ndims(tempS[:Dpfun])}}(undef,Niter)

ii = 1
time_iter[ii,1] = time_M
time_iter[ii,2] = time_S


Dπ_iter[ii,:,:] = value.(tempM[:Dπ])
tau_e_iter[ii,:] = value.(tempM[:tau_e])
tau_c_iter[ii,:] = value.(tempM[:tau_c])
ObjUL_iter[ii] = objective_value(tempM)
ObjML_iter[ii] = objective_value(tempS)
ObjLL_iter[ii] = value.(tempS[:Obj_KKT])

gpmax_iter[ii] = value.(tempS[:gpmax])
gpo_iter[ii] = value.(tempS[:gpo])
gp_iter[ii] = value.(tempS[:gp])
lambda_iter[ii] = value.(tempS[:lambda])
dp_iter[ii] = value.(tempS[:Dpfun])
conv_temp = (value.(tempS[:Obj_SP]) - value.(tempM[:Obj_SP]))
global conv[ii,:] = [round.(value.(tempS[:Obj_SP]), digits = 3), round.(value.(tempS[:Obj_SP]), digits = 3), round.(value.(tempS[:Obj_SP]) - value.(tempM[:Obj_SP]), digits = 3), (value.(tempS[:Obj_SP]) - value.(tempM[:Obj_SP])) / value.(tempS[:Obj_SP])]
global SPbound_iter[ii] = value.(tempS[:SPbound])
global MPbound_iter[ii] = value.(tempM[:MPbound])

for ii in 2:Niter
   time_M = @elapsed global tempM = MP(testsystem, s, kappa, PolicyB, timeset, υ0, gpmax_iter, gpo_iter, dp_iter, ii-1, RegGap, UserGap, bound, peaktime, α0inv, π0)
   global Dπ_iter[ii,:,:] = value.(tempM[:Dπ])
   global tau_e_iter[ii,:] = value.(tempM[:tau_e])
   global tau_c_iter[ii,:] = value.(tempM[:tau_c])
   global ObjUL_iter[ii] = objective_value(tempM)
   global ObjML_iter[ii] = objective_value(tempS)
   global ObjLL_iter[ii] = value.(tempS[:Obj_KKT])

   time_S = @elapsed global tempS = SP(testsystem, s, kappa, timeset, PdownMax, η, σ0_pv, σ0_wind, υ0, value.(tempM[:Dπ]), value.(tempM[:tau_e]), value.(tempM[:tau_c]), RegGap, UserGap, α0inv, π0)
   global gpmax_iter[ii] = value.(tempS[:gpmax])
   global gpo_iter[ii] = value.(tempS[:gpo])
   global gp_iter[ii] = value.(tempS[:gp])
   global lambda_iter[ii] = value.(tempS[:lambda])
   global dp_iter[ii] = value.(tempS[:Dpfun])

   conv_temp = (value.(tempS[:Obj_SP]) - value.(tempM[:Obj_SP]))
   global conv[ii,:] = [round.(value.(tempS[:Obj_SP]), digits = 3), round.(value.(tempS[:Obj_SP]), digits = 3), round.(value.(tempS[:Obj_SP]) - value.(tempM[:Obj_SP]), digits = 3), (value.(tempS[:Obj_SP]) - value.(tempM[:Obj_SP])) / value.(tempS[:Obj_SP])]
   global SPbound_iter[ii] = value.(tempS[:SPbound])
   global MPbound_iter[ii] = value.(tempM[:MPbound])

   global time_iter[ii,1] = time_M
   global time_iter[ii,2] = time_S

end
