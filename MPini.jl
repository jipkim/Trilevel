using Pkg, DataFrames
using JuMP, KNITRO
using Distributions
function MPini(testsystem, s, kappa, PolicyB, timeset, υ0, UserGap, bounds::bounds, peaktime)
   ScaleParam = 1e-3
   ScaleObj = 1e-5

   include("NetworkDataType.jl")
   include("NetworkLoad.jl")
   include("Demand24Load.jl")
   filename_Load24 = string(pwd(),"/data/load24data-PJM.csv")
   loadcurve = Demand24Load(filename_Load24)

   include("GenCurveLoad.jl")
   filename_windcurve = string(pwd(),"/data/hourly_wind_gen_2018.csv")
   filename_pvcurve = string(pwd(),"/data/hourly_solar_gen_2018.csv")
   pvcurve = GenCurveLoad(filename_pvcurve)[2:31:end,:]/2883.8  # pick one day a month over the year / nameplate capacity = 2883MW
   windcurve = GenCurveLoad(filename_windcurve)[2:31:end,:]/1300 # pick one day a month over the year / nameplate capacity = 1300MW

   stateset = 1:6
   loadzone = Dict()
   loadzone[1] = 1
   loadzone[2] = 2
   loadzone[3] = 3
   loadzone[4] = [4, 5, 8] # WCMA, NEMA, SEMA
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

   # Avoid MIP (Not implementing UC problem)
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

   B_gn = Array{Array{Int64},1}(undef,length(buses))
   for ii in 1:length(buses)
      B_gn[ii] = Int64[]
   end
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

   S_gs = Array{Array{Int64},1}(undef,length(stateset))
   for ii in 1:length(stateset)
      S_gs[ii] = Int64[]
   end
   for s in stateset
      for n in loadzone[s]
         for g in B_gn[n]
            push!(S_gs[s], g)
         end
      end
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


   Γ = 0
   I_R = intersect(genset, union(windset,pvset))
   I_R_hat = intersect(genset_hat, union(windset,pvset))
   I_C = setdiff(genset, union(windset,pvset))
   I_C_hat = setdiff(genset_hat,union(windset,pvset))


   #---------------------------------------------------------------------------------#
   m = Model(with_optimizer(KNITRO.Optimizer, ms_enable = 1, ms_maxsolves = 5000, ms_terminate = 1, par_blasnumthreads = 10, maxit = 100000))
   #---------------------------------------------------------------------------------#
   @variables(m, begin
      0 <= gpo[g in union(genset, genset_hat), t in timeset]
      0 <= gpmax[g in genset_hat] <= 9000
      -PdownMax <= pdown[n in busset, t in timeset] <= PdownMax
      0 <= x[s in stateset, t in timeset, j in intersect(S_gs[s], union(I_R, I_R_hat))]
      0 <= yLB[g in union(I_C, I_C_hat), t in timeset]
      0 <= yUB[g in union(I_C, I_C_hat), t in timeset]
      0 <= Cinv_total[n in busset]
      bounds.Dπ_LB <= Dπ[n in busset, t in timeset] <= bounds.Dπ_UB
      bounds.Dπ_LB <= Dπ_off <= bounds.Dπ_UB
      bounds.Dπ_LB <= Dπ_on <= bounds.Dπ_UB
      bounds.tau_e_LB <= tau_e[s in stateset] <= bounds.tau_e_UB
      bounds.tau_c_LB <= tau_c[s in stateset] <= bounds.tau_c_UB
   end)

   @constraint(m, Dπ_all[n in busset, t in peaktime], Dπ[n,t] == Dπ_on)
   @constraint(m, Dπ_all2[n in busset, t in setdiff(timeset, peaktime)], Dπ[n,t] == Dπ_off)
   @constraint(m, Dπ_ineq, Dπ_on >= Dπ_off)

   #---------------------------------------------------------------------------------#
   @variables(m, begin
      fp[l in lineset, t in timeset]
      0 <= gp[g in union(genset, genset_hat), t in timeset]
      theta[n in busset, t in timeset]
      xi[l in lineset, t in timeset]
      0 <= lambda[n in busset, t in timeset] <= 500
      0 <= gammaLB[g in union(genset, genset_hat), t in timeset] <= 1e6
      0 <= gammaUB[g in union(genset, genset_hat), t in timeset] <= 1e6
      0 <= deltaLB[l in lineset, t in timeset] <= 1e6
      0 <= deltaUB[l in lineset, t in timeset] <= 1e6
      OCgen[g in union(genset, genset_hat), t in timeset]
      Obj_KKT
   end)
   #---------------------------------------------------------------------------------#
   @variable(m, 0 <= α[i in intersect(S_gs[s], union(I_C, I_C_hat)), t in timeset] <= 1)
   @constraint(m, Affine[t in timeset], sum(α[i,t] for i in intersect(S_gs[s], union(I_C, I_C_hat))) == 1)
   #---------------------------------------------------------------------------------#
   @variable(m, 1.0 * Dp[n] * loadcurve[n,t] <=  Dpfun[n in busset, t in timeset] <= 1.0 * Dp[n] * loadcurve[n,t])

   @objective(m, Min,
      ScaleObj*(sum(
         sum(Dπ[n,t] * Dpfun[n,t] for n in loadzone[s])
         + tau_e[s] * sum(gp[i,t] for i in intersect(union(I_R, I_R_hat), S_gs[s]))
      for t in timeset)
      + DCRF * tau_c[s] * sum(gpmax[i] for i in intersect(I_R_hat, S_gs[s]))
      )
      )

   @variable(m, MPbound)
   @constraint(m, MPbound ==
         ScaleObj*(sum(
            sum(Dπ[n,t] * Dpfun[n,t] for n in loadzone[s])
            + tau_e[s] * sum(gp[i,t] for i in intersect(union(I_R, I_R_hat), S_gs[s]))
         for t in timeset)
         + DCRF * tau_c[s] * sum(gpmax[i] for i in intersect(I_R_hat, S_gs[s]))
         )
      )

   @constraint(m, RPS, sum(sum(gp[g,t] for g in intersect(S_gs[s], union(I_R, I_R_hat))) for t in timeset) >= sum(sum(kappa[s] * Dpfun[n,t] for t in timeset) for n in loadzone[s]))
   @constraint(m, PolicyBudget,
      ScaleParam*(sum(tau_e[s]*sum(gp[g,t] for g in intersect(S_gs[s], union(I_R, I_R_hat))) for t in timeset)
      + DCRF * tau_c[s]*sum(gpmax[g] for g in intersect(S_gs[s], I_R_hat)))
      <= ScaleParam*PolicyB[s])

   @variable(m, Obj_SP >= 0)
   @constraint(m, Obj_SP_def, Obj_SP ==
      ScaleObj*(sum(
         sum(Dπ[n,t] * Dpfun[n,t] for n in loadzone[s])
         - sum(lambda[n,t] * Dpfun[n,t] for n in loadzone[s])
         + sum(lambda[S_g[i],t] * gp[i,t] for i in S_gs[s])
         + tau_e[s] * sum(gp[i,t] for i in intersect(union(I_R, I_R_hat), S_gs[s]))
         - sum(Coperg[i][2] * gp[i,t] for i in intersect(union(I_R, I_R_hat), S_gs[s]))
         - sum(Coperg[i][2] * (gp[i,t] - α[i,t]*υ[i,t]) for i in intersect(union(I_C, I_C_hat), S_gs[s]))
      for t in timeset)
      + DCRF * tau_c[s] * sum(gpmax[i] for i in intersect(I_R_hat, S_gs[s]))
      - sum(Cinv_total[n] for n in loadzone[s])
      )
   )

   @constraint(m, Cinv_def[n in busset], Cinv_total[n] == sum(Cinvg[g] * gpmax[g] for g in intersect(genset_hat, B_gn[n])))


   if testsystem == "ISONE8busTN"
      gpmax_goal_R = [3922.0457029487675, 601.92789666824, 658.245905433098, 0.0, 0.0, 543.6165982508121, 437.4257356863675, 3076.6108418250856]
      gpmax_goal_C = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
   elseif testsystem == "ISONE8busTN_rCoal"
      gpmax_goal_R = [3922.0457029487675, 601.9278966682401, 658.245905433098, 0.0, 3076.6108418250856, 543.6165982508121, 437.4257356863675, 0.0]
      gpmax_goal_C = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
   elseif testsystem == "ISONE8busTN_rCoalrNuc"
      gpmax_goal_R = [3922.0457029487666, 601.92789666824, 658.245905433098, 0.0, 3076.610841825086, 543.6165982508121, 437.42573568636743, 0.0]
      gpmax_goal_C = [0.0, 0.0, 12.846366719999999, 0.0, 0.0, 0.0, 0.0, 0.0]
   end

   #---------------------------------------------------------------------------------#
   #-----------------#
   @constraint(m, gpoDef1[g in intersect(I_R, S_gs[s]), t in timeset], gpo[g,t] == ρ[g,t] * generators[g].Pmax + υ[g,t])
   @constraint(m, gpoDef2[g in intersect(I_R_hat, S_gs[s]), t in timeset], gpo[g,t] == ρ[g,t] * gpmax[g] + υ[g,t])

   @constraint(m, xdef1[t in timeset, j in intersect(S_gs[s], I_R)], x[S_g[j],t,j] == σ[j] * generators[j].Pmax)
   @constraint(m, xdef2[t in timeset, j in intersect(S_gs[s], I_R_hat)], x[S_g[j],t,j] == σ[j] * gpmax[j])

   @constraint(m, yLBdef1[i in intersect(I_C, S_gs[s]), t in timeset], yLB[i,t] == gpo[i,t] - generators[i].Pmin)
   @constraint(m, yUBdef1[i in intersect(I_C, S_gs[s]), t in timeset], yUB[i,t] == generators[i].Pmax - gpo[i,t])

   @NLconstraint(m, CCreform1[i in intersect(I_C, S_gs[s]), t in timeset], sum(( inv_ϕ * α[i,t] * x[S_g[j],t,j])^2 for j in intersect(S_gs[S_g[i]], union(I_R, I_R_hat))) <= yLB[i,t]^2)
   @NLconstraint(m, CCreform2[i in intersect(I_C, S_gs[s]), t in timeset], sum(( inv_ϕ * α[i,t] * x[S_g[j],t,j])^2 for j in intersect(S_gs[S_g[i]], union(I_R, I_R_hat))) <= yUB[i,t]^2)

   #---------------------------------------------#
   #:::::|||| Middle-level constraints |||| :::::#
   @constraint(m, gpoLB2[g in genset, t in timeset], gpo[g,t] >= generators[g].Pmin)
   @constraint(m, gpoUB2_C[g in I_C, t in timeset], gpo[g,t] <= generators[g].Pmax)
   @constraint(m, gpoUB2_R[g in I_R, t in timeset], gpo[g,t] <= ρ[g,t] * generators[g].Pmax + υ[g,t])
   @constraint(m, gpoLB3[g in genset_hat, t in timeset], gpo[g,t] >= Γ*gpmax[g])
   @constraint(m, gpoUB3_C[g in I_C_hat, t in timeset], gpo[g,t] <= gpmax[g])
   @constraint(m, gpoUB3_R[g in I_R_hat, t in timeset], gpo[g,t] <= ρ[g,t] * gpmax[g] + υ[g,t])

   @constraint(m, SupplyDemandBalance[t in timeset, n in loadzone[s]], sum(gpo[g,t] for g in B_gn[n]) + pdown[n,t] >= Dpfun[n,t])
   #---------------------------------------------------------------------------------#
   @constraint(m, sumI_R_hatUB[nn in setdiff(busset, loadzone[s])], sum(gpmax[g] for g in intersect(B_gn[nn], I_R_hat)) == gpmax_goal_R[nn])
   @constraint(m, sumI_C_hatUB[nn in setdiff(busset, loadzone[s])], sum(gpmax[g] for g in intersect(B_gn[nn], I_C_hat)) == gpmax_goal_C[nn])
   #---------------------------------------------------------------------------------#
   @constraint(m, gpoDefF_C[g in setdiff(I_C, S_gs[s]), t in timeset], gpo[g,t] == generators[g].Pmax)
   @constraint(m, gpoDefF_R[g in setdiff(I_R, S_gs[s]), t in timeset], gpo[g,t] == ρ[g,t] * generators[g].Pmax + υ[g,t])
   @constraint(m, gpoDefF_C_hat[g in setdiff(I_C_hat, S_gs[s]), t in timeset], gpo[g,t] == gpmax[g])
   @constraint(m, gpoDefF_R_hat[g in setdiff(I_R_hat, S_gs[s]), t in timeset], gpo[g,t] == ρ[g,t] * gpmax[g] + υ[g,t])
   #---------------------------------------------------------------------------------#
   #---------------------------------------------------------------------------------#
   #---------------------------------------------#
   #:::::|||| Lower-level Equality constraints |||| :::::#
   @constraint(m, DCPF[l in lineset, t in timeset], fp[l,t] == (1/lines[l].x) * (theta[lines[l].fbus,t] - theta[lines[l].tbus,t]))
   @constraint(m, slack[t in timeset], theta[1,t] == 0)
   @constraint(m, NodeBalance[n in busset, t in timeset], sum(gp[g,t] for g in B_gn[n]) + sum(fp[l,t] for l in buses[n].inline)
      - sum(fp[l,t] for l in buses[n].outline) == Dpfun[n,t])


   #----------------------#
   @constraint(m, gpLB[g in union(genset, genset_hat), t in timeset], gp[g,t] >= 0)
   @constraint(m, gpUB[g in union(genset, genset_hat), t in timeset], gp[g,t] <= gpo[g,t])
   @constraint(m, LineCapLB[l in lineset, t in timeset], fp[l,t] >= - lines[l].u)
   @constraint(m, LineCapUB[l in lineset, t in timeset], fp[l,t] <= + lines[l].u)

   @variable(m, Obj_prim[t in timeset])
   @variable(m, Obj_dual[t in timeset])
   @constraint(m, Oprim[t in timeset], Obj_prim[t] == sum(- Coperg[g][2] * gp[g,t] for g in union(genset, genset_hat)))
   @constraint(m, Odual[t in timeset], Obj_dual[t] ==
      sum(gpo[g,t]*gammaUB[g,t] - 0 * gammaLB[g,t] for g in union(genset, genset_hat))
      + sum(lines[l].u * deltaUB[l,t] + lines[l].u * deltaLB[l,t] for l in lineset)
      - sum(Dpfun[n,t] * lambda[n,t] for n in busset)
      )
   @constraint(m, const_gp[g in union(genset, genset_hat), t in timeset], - gammaLB[g,t] + gammaUB[g,t] - lambda[B_g[g],t] == - Coperg[g][2])
   @constraint(m, const_fp[l in lineset, t in timeset], - deltaLB[l,t] + deltaUB[l,t] + xi[l,t] - lambda[lines[l].tbus,t] + lambda[lines[l].fbus,t] == 0)
   @constraint(m, const_theta[n in busset, t in timeset], - sum(xi[l,t]/lines[l].x for l in buses[n].outline) + sum(xi[l,t]/lines[l].x for l in buses[n].inline) == 0)
   @constraint(m, SD_split[t in timeset], Obj_dual[t] - Obj_prim[t] == 0)
   #-----------------#
   ##########################################################################################
   optimize!(m)
   return m
end
