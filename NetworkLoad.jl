using MAT, DelimitedFiles
function NetworkLoad(testsystem)
   PUcorrection = 1
   if testsystem == "AP15busDN"
      usingMatPower = 0
   elseif testsystem == "IEEE33busDN" || testsystem == "IEEE33busDN2" || testsystem == "ISONE8busTN" || testsystem == "ISONE8busTN_rCoal" || testsystem == "ISONE8busTN_rCoalrNuc" || testsystem == "ISONE2busTN"
      usingMatPower = 0   
   else
      usingMatPower = 1
      PUcorrection = 100
   end

   if usingMatPower == 0
      filename_Node = string(pwd(),"/data/$testsystem/Node.csv")
      filename_Generator = string(pwd(),"/data/$testsystem/Generator.csv")
      filename_Line = string(pwd(),"/data/$testsystem/Line.csv")

      busmat = readdlm(filename_Node, ',', header=true)[1]
      buses = Bus[]
      for i in 1:size(busmat,1)
         bindex = busmat[i,1]
         btype = 0
         Pd = busmat[i,2]
         Qd = busmat[i,3]
         Gs = busmat[i,6]
         Bs = busmat[i,7]
         area = 0
         Vm = 0
         Va = 0
         baseKV = busmat[i,8]
         bzone = 0
         Vmax = busmat[i,4]
         Vmin = busmat[i,5]
         b = Bus(bindex, btype, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, bzone, Vmax, Vmin)
         push!(buses, b)
      end

      genmat = readdlm(filename_Generator, ',', header=true)[1]
      generators = Generator[]
      for i in 1:size(genmat,1)
         gindex = genmat[i,1]
         gtype = genmat[i,2]
         location = genmat[i,3]
         Pg = genmat[i,4]
         Qg = genmat[i,5]
         Qmax = genmat[i,6]
         Qmin = genmat[i,7]
         Vg = genmat[i,8]
         mBase = genmat[i,9]
         status = genmat[i,10]
         Pmax = genmat[i,11]
         Pmin = genmat[i,12]
         cost = [genmat[i,13], genmat[i,14], genmat[i,15]]
         SUcost = genmat[i,16]
         SDcost = genmat[i,17]
         RU = genmat[i,18]
         RD = genmat[i,19]
         UPtime = genmat[i,20]
         DNtime = genmat[i,21]
         g = Generator(gindex, gtype, location, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, cost, SUcost, SDcost, RU, RD, UPtime, DNtime)
         push!(generators, g)
      end

      branchmat = readdlm(filename_Line, ',', header=true)[1]
      lines = Line[]
      for i in 1:size(branchmat,1)
         lindex = i
         fbus = Int(branchmat[i,2])
         tbus = Int(branchmat[i,1])
         r = branchmat[i,3]
         x = branchmat[i,4]
         b = branchmat[i,5]
         u = branchmat[i,6]
         push!(buses[tbus].children, fbus)#children
         push!(buses[fbus].ancestor, tbus)#ancestor
         push!(buses[tbus].inline, lindex)
         push!(buses[fbus].outline, lindex)
         l = Line(lindex, fbus, tbus, r, x, b, u)
         push!(lines,l)
      end
      for l in 1:length(lines)
         for k in 1:length(lines)
            if lines[l].fbus == lines[k].tbus
               push!(lines[l].ancestor, k)
            elseif lines[l].tbus == lines[k].fbus
               push!(lines[l].children, k)
            end
         end
      end
      datamat = Dict()
      datamat["bus"] = busmat
      datamat["branch"] = branchmat
      datamat["gen"] = genmat
   elseif usingMatPower == 1
      mpc = matread(string(pwd(),"/data/","$testsystem.mat"))["mpc"]
      # datamat = mpc
      #Sbase = mpc["baseMVA"] * 1e6 # MVA to VA
      #Vbase = mpc["bus"][:,10] * 1e3 # kV to V
      #Zbase = Vbase^2 / Sbase

      buses = Bus[]
      for i in 1:size(mpc["bus"],1)
         bindex = mpc["bus"][i,1]
         btype = mpc["bus"][i,2]
         Pd = mpc["bus"][i,3]
         Qd = mpc["bus"][i,4]
         Gs = mpc["bus"][i,5]
         Bs = mpc["bus"][i,6]
         area = mpc["bus"][i,7]
         Vm = mpc["bus"][i,8]
         Va = mpc["bus"][i,9]
         baseKV = mpc["bus"][i,10]
         bzone = mpc["bus"][i,11]
         Vmax = mpc["bus"][i,12]
         Vmin = mpc["bus"][i,13]
         b = Bus(bindex, btype, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, bzone, Vmax, Vmin)
         push!(buses, b)
      end

      generators = Generator[]
      for i in 1:size(mpc["gen"],1)
         gindex = i
         gtype = "NotDefined"
         location = mpc["gen"][i,1]
         Pg = mpc["gen"][i,2]
         Qg = mpc["gen"][i,3]
         Qmax = mpc["gen"][i,4]
         Qmin = mpc["gen"][i,5]
         Vg = mpc["gen"][i,6]
         mBase = mpc["gen"][i,7]
         status = mpc["gen"][i,8]
         Pmax = mpc["gen"][i,9]
         Pmin = mpc["gen"][i,10]
         cost = mpc["gencost"][i,5:end]
         SUcost = mpc["gencost"][i,2]
         SDcost = mpc["gencost"][i,3]
         RU = 0
         RD = 0
         UPtime = 0
         DNtime = 0
         g = Generator(gindex, gtype, location, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, cost, SUcost, SDcost, RU, RD, UPtime, DNtime)
         push!(generators, g)
      end

      lines = Line[]
      for i in 1:size(mpc["branch"],1)
         lindex = i
         fbus = Int(mpc["branch"][i,2]) #definition of fbus is based on  {Low, 2013} which is opposite of Matpower
         tbus = Int(mpc["branch"][i,1])
         r = mpc["branch"][i,3]/PUcorrection
         x = mpc["branch"][i,4]/PUcorrection
         b = mpc["branch"][i,5]
         u = 10000
         push!(buses[tbus].children, fbus)#children
         push!(buses[fbus].ancestor, tbus)#ancestor
         push!(buses[tbus].inline, lindex)
         push!(buses[fbus].outline, lindex)
         l = Line(lindex, fbus, tbus, r, x, b, u)
         push!(lines,l)
      end

      for l in 1:length(lines)
         for k in 1:length(lines)
            if lines[l].fbus == lines[k].tbus
               push!(lines[l].ancestor, k)
            elseif lines[l].tbus == lines[k].fbus
               push!(lines[l].children, k)
            end
         end
      end
   end
   return buses, lines, generators, datamat
end
