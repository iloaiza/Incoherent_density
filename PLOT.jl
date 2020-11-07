function sinc(x,tol=1e-15)
	if abs(x) <= tol
		return 1
	else
		return sin(x) / x
	end
end

# compare time evolution with theoretical value for dynamic averaging procedure
# tf: final time for integration. res: resolution for integral discretization
# ObsAdi is an array, each element should be the matrix operator for an observable in the molecular eigenstate representation
function timeEvo(tf,res,H,ϕ0;expOp=false,E=false,U=false,ObsNames=false,ObsAdi=false)
	Tarr = range(0,stop=tf,length=res)
	numObs = length(ObsAdi)
	dt = Tarr[2] - Tarr[1]
	Hsts = length(ϕ0)
	ϕarr = zeros(Complex,Hsts,res)
	muAdi = U' * ϕ0
	ϕarr[:,1] = muAdi

	# =
	if expOp == false
		if E == false || U == false
			println("Diagonalizing...")
			@time E,U = adiabatize(H)
		end
		expArr = exp.((-1im*dt) .* E)
		println("Building exponential operator")
		@time expOp = dia2en(Diagonal(expArr),U)
	end
	# =#
	

	println("Starting time evolution")
	@time for i in 2:res
		println("step $(i-1) of $(res-1)")
		#@time ϕarr[:,i] = expOp * ϕarr[:,i-1]
		ϕarr[:,i] = exp.((-1im * Tarr[i]) .* E) .* muAdi
		@show abs2(dot(ϕarr[:,i],ϕarr[:,1]))
	end

	if ObsAdi == false
		println("No observables, calculating purity")
		ξ = 0.0
		@time for i in 1:res
			for j in 1:res
				ξ += dot(ϕarr[:,i],ϕarr[:,j])
			end
		end
		ξ /= res^2
		@show ξ

		println("Returning array of times of evolution")
		return ϕarr
	end

	#Routine includes treatment of observables 
	println("Starting observables analysis: calculating observables for each time and observable")
	ObsVals = zeros(numObs,res)
	for opnum in 1:numObs
		println("Starting observable $opnum for all time evolution")
		@time for tnum in 1:res
			ObsVals[opnum,tnum] = real(dot(ϕarr[:,tnum],ObsAdi[opnum] * ϕarr[:,tnum]))
		end
	end
	
	println("Starting time averaging and analysis")
	ρVals = zeros(numObs+1,res) #ρVals[1,:] is purity
	ρVals[1,1] = 1

	for opnum in 1:numObs
		ρVals[opnum+1,1] = ObsVals[opnum,1]
	end
	
	@time for i in 1:res-1
		ovlps = ϕarr[:,i+1]' * ϕarr[:,1:i]
		ovlps = sum(abs2.(ovlps))
		@show ovlps
		ρVals[1,i+1] = i^2 * ρVals[1,i] +2*ovlps + 1
		ρVals[1,i+1] /= (i+1)^2

		for opnum in 1:numObs
			ρVals[opnum+1,i+1] = i*ρVals[opnum+1,i] + ObsVals[opnum,i+1]
			ρVals[opnum+1,i+1] /= i+1
		end
	end

	println("Building theoretical values")
	# =
	TheoVals = zeros(numObs+1,res)
	
	ρmu = muAdi * muAdi'
	for i in 1:res
		println("Building theoretical at time $i of $res")
		ρt = zeros(Hsts,Hsts)
		@time for st1 in 1:Hsts
			for st2 in 1:Hsts
				sincval = sinc((E[st1] - E[st2]) * Tarr[i] / 2)
				ρt[st1,st2] = sincval * ρmu[st1,st2]
				#TheoVals[1,i] +=  abs2(sincval * muAdi[st1] * muAdi[st2])
				for opnum in 1:numObs
					TheoVals[opnum+1,i] += ObsAdi[opnum][st2,st1] * ρt[st1,st2]
				end
			end
		end
		TheoVals[1,i] = sum(abs2.(ρt))
	end
	
	@show TheoVals[1,1], TheoVals[1,end]
	@show TheoVals[2,1], TheoVals[2,end]
	# =#

	@show ρVals[1,1], ρVals[1,end]
	@show ρVals[2,1], ρVals[2,end]
	
	println("Starting plotting procedure")
	if ObsNames == false
		ObsNames = ["O$i" for i in 1:numObs]
	end

	P = scatter(Tarr,ρVals[1,:],label="ξ",marker=(10,:blue,:circle),xlabel="Time interval (a.u.", ylabel="Normalized varied values")
	for opnum in 1:numObs
		scatter!(Tarr,ρVals[opnum+1,:],label=ObsNames[opnum],marker=(10,:red,:rect))
	end

	# =
	plot!(Tarr,TheoVals[1,:],label="ξtheo",line=(2.5,:blue))
	for opnum in 1:numObs
		plot!(Tarr,TheoVals[opnum+1,:],label=ObsNames[opnum]*"theo",line=(2.5,:red,:dash))
	end
	# =#

	return expOp,ρVals,P
end

# same as last time evolution, but plots and compares two different resolutions
function timeDoubleEvo(tf,H,ϕ0,res1,res2;expOp=false,E=false,U=false,ObsNames=false,ObsAdi=false)
	Tarr1 = range(0,stop=tf,length=res1)
	Tarr2 = range(0,stop=tf,length=res2)
	numObs = length(ObsAdi)
	dt1 = Tarr1[2] - Tarr1[1]
	dt2 = Tarr2[2] - Tarr2[1]
	Hsts = length(ϕ0)
	ϕarr1 = zeros(Complex,Hsts,res1)
	ϕarr2 = zeros(Complex,Hsts,res2)

	if E == false && U == false
		E,U = adiabatize(H)
	end

	muAdi = U' * ϕ0
	ϕarr1[:,1] = muAdi
	ϕarr2[:,1] = muAdi


	println("Starting time evolution 1")
	@time for i in 2:res1
		println("step $(i-1) of $(res1-1)")
		#@time ϕarr[:,i] = expOp * ϕarr[:,i-1]
		ϕarr1[:,i] = exp.((-1im * Tarr1[i]) .* E) .* muAdi
		@show abs2(dot(ϕarr1[:,i],ϕarr1[:,1]))
	end

	println("Starting time evolution 2")
	@time for i in 2:res2
		println("step $(i-1) of $(res2-1)")
		#@time ϕarr[:,i] = expOp * ϕarr[:,i-1]
		ϕarr2[:,i] = exp.((-1im * Tarr2[i]) .* E) .* muAdi
		@show abs2(dot(ϕarr2[:,i],ϕarr2[:,1]))
	end

	#Routine includes treatment of observables 
	println("Starting observables analysis 1: calculating observables for each time and observable")
	ObsVals1 = zeros(numObs,res1)
	for opnum in 1:numObs
		println("Starting observable $opnum for all time evolution")
		@time for tnum in 1:res1
			ObsVals1[opnum,tnum] = real(dot(ϕarr1[:,tnum],ObsAdi[opnum] * ϕarr1[:,tnum]))
		end
	end

	println("Starting observables analysis 2: calculating observables for each time and observable")
	ObsVals2 = zeros(numObs,res2)
	for opnum in 1:numObs
		println("Starting observable $opnum for all time evolution")
		@time for tnum in 1:res2
			ObsVals2[opnum,tnum] = real(dot(ϕarr2[:,tnum],ObsAdi[opnum] * ϕarr2[:,tnum]))
		end
	end
	
	println("Starting time averaging and analysis")
	ρVals1 = zeros(numObs+1,res1) #ρVals[1,:] is purity
	ρVals1[1,1] = 1

	for opnum in 1:numObs
		ρVals1[opnum+1,1] = ObsVals1[opnum,1]
	end
	
	@time for i in 1:res1-1
		ovlps = ϕarr1[:,i+1]' * ϕarr1[:,1:i]
		ovlps = sum(abs2.(ovlps))
		@show ovlps
		ρVals1[1,i+1] = i^2 * ρVals1[1,i] +2*ovlps + 1
		ρVals1[1,i+1] /= (i+1)^2

		for opnum in 1:numObs
			ρVals1[opnum+1,i+1] = i*ρVals1[opnum+1,i] + ObsVals1[opnum,i+1]
			ρVals1[opnum+1,i+1] /= i+1
		end
	end

	ρVals2 = zeros(numObs+1,res2) 
	ρVals2[1,1] = 1

	for opnum in 1:numObs
		ρVals2[opnum+1,1] = ObsVals2[opnum,1]
	end
	
	@time for i in 1:res2-1
		ovlps = ϕarr2[:,i+1]' * ϕarr2[:,1:i]
		ovlps = sum(abs2.(ovlps))
		@show ovlps
		ρVals2[1,i+1] = i^2 * ρVals2[1,i] +2*ovlps + 1
		ρVals2[1,i+1] /= (i+1)^2

		for opnum in 1:numObs
			ρVals2[opnum+1,i+1] = i*ρVals2[opnum+1,i] + ObsVals2[opnum,i+1]
			ρVals2[opnum+1,i+1] /= i+1
		end
	end



	println("Building theoretical values")
	# =
	res = maximum([res1,res2])
	if res == res1
		Tarr = Tarr1
	else
		Tarr = Tarr2
	end
	TheoVals = zeros(numObs+1,res)
	@show size(E), size(Tarr)
	
	ρmu = muAdi * muAdi'
	for i in 1:res
		println("Building theoretical at time $i of $res")
		ρt = zeros(Hsts,Hsts)
		@time for st1 in 1:Hsts
			for st2 in 1:Hsts
				sincval = sinc((E[st1] - E[st2]) * Tarr[i] / 2)
				ρt[st1,st2] = sincval * ρmu[st1,st2]
				#TheoVals[1,i] +=  abs2(sincval * muAdi[st1] * muAdi[st2])
				for opnum in 1:numObs
					TheoVals[opnum+1,i] += ObsAdi[opnum][st2,st1] * ρt[st1,st2]
				end
			end
		end
		TheoVals[1,i] = sum(abs2.(ρt))
	end
	
	@show TheoVals[1,1], TheoVals[1,end]
	@show TheoVals[2,1], TheoVals[2,end]
	# =#
	@show ρVals1[1,1], ρVals1[1,end]
	@show ρVals1[2,1], ρVals1[2,end]

	@show ρVals2[1,1], ρVals2[1,end]
	@show ρVals2[2,1], ρVals2[2,end]
	
	println("Starting plotting procedure")
	if ObsNames == false
		ObsNames = ["O$i" for i in 1:numObs]
	end


	P = scatter(Tarr[1:20:end],TheoVals[1,1:20:end],label="ξtheo",mark=(10,:blue,:circle))
	for opnum in 1:numObs
		scatter!(Tarr[1:20:end],TheoVals[opnum+1,1:20:end],label=ObsNames[opnum]*"theo",mark=(10,:red,:circle))
	end

	ytext = "Purity and " * L"S_0" * " population"
	plot!(Tarr1,ρVals1[1,:],label="ξ",line=(2.5,:blue),xlabel="Total time (a.u.)", ylabel=ytext)
	for opnum in 1:numObs
		plot!(Tarr1,ρVals1[opnum+1,:],label=ObsNames[opnum],line=(2.5,:red))
	end

	plot!(Tarr2,ρVals2[1,:],label="ξ",line=(2.5,:blue,:dash))
	for opnum in 1:numObs
		plot!(Tarr2,ρVals2[opnum+1,:],label=ObsNames[opnum],line=(2.5,:red,:dash))
	end

	plot!(legend=false,xguidefont=FONT,yguidefont=FONT,xtickfont=FONT,ytickfont=FONT,size=SIZE)
	
	@show dt1,dt2

	return TheoVals,ρVals1,ρVals2,P
end

# build observables for Rhodopsin system
# Hamiltonian and Franck-Condon wavefunctions can be built as shown in fdebug or loaded using running 'loading("rhodopsin")'
function makeRhodObservables(U=Urhod)
	Ps0 = zeros(8000,8000)
	for i in 1:4000
		Ps0[i,i] = 1
	end

	#Pcis0,Pcis1,Ptrans0,Ptrans1 = projectorsHeavisideRhodopsinRedDiaBasis()
	Pcis0,Pcis1,Ptrans0,Ptrans1 = projectorsHeavisideRhodopsinDiabatic()
	#[dia2en(Ps0,U),dia2en(Pcis0,U),dia2en(Pcis1,U),dia2en(Ptrans0,U),dia2en(Ptrans1,U)], [Ps0,Pcis0,Pcis1,Ptrans0,Ptrans1]
	return dia2en(Ps0,U), dia2en(Ptrans1,U) 
end

function makeLVCObservables(H)
	E,U = adiabatize(H)
	Hsts = length(E)
	Ps0 = zeros(Hsts,Hsts)
	for i in 1:Int(Hsts/2)
		Ps0[i,i] = 1
	end

	return dia2en(Ps0,U), Ps0
end

# perform evolution for Lindblad dynamics and compare with theoretical convergence values
function LindbladPlots(H,ϕ0,ϵ,steps,ObsDia;E=false,U=false,ObsAdi=false,ObsNames=false)
	@time T,VALS = lindbladIntegration(H,ϵ,ϕ0,steps,ObsDia;AbsTol=1e-11,RelTol=1e-11,emin=1e-7,emax=1e-5,verbose=true)
	#@time T,VALS = lindbladRK4(H,ϵ,ϕ0,steps,ObsDia)
	numObs = length(ObsDia)
	Hsts = length(ϕ0)
	THEO = zeros(numObs+1,steps+1)

	if E == false || U == false
		@time E,U = adiabatize(H)
	end

	if ObsAdi == false
		ObsAdi = zeros(Hsts,Hsts,numObs)
		for opnum in 1:numObs
			ObsAdi[:,:,opnum] = dia2en(ObsDia[opnum],U)
		end
	end

	muAdi = U' * ϕ0
	
	@time for i in 1:steps+1
		α = T[i]
		ρα = zeros(Hsts,Hsts)
		for st1 in 1:Hsts
			for st2 in 1:Hsts
				ω = E[st1] - E[st2]
				ρα[st1,st2] = exp(-ω^2 * α) * muAdi[st1] * muAdi[st2]
				for opnum in 1:numObs
					THEO[opnum+1,i] += ρα[st1,st2] * ObsAdi[st2,st1,opnum]
				end
			end
		end
		THEO[1,i] = sum(abs2.(ρα))
	end
	Tarr = T

	if ObsNames == false
		ObsNames = ["O$i" for i in 1:numObs]
	end

	# =
	P = scatter(1:120:steps+1,THEO[1,1:120:end],label="ξtheo",marker=(10,:blue,:circle))
	for opnum in 1:numObs
		scatter!(1:120:steps+1,THEO[opnum+1,1:120:end],label=ObsNames[opnum]*"theo",marker=(10,:red,:cirlce))
	end

	ytext = "Purity and " * L"S_0" * " population"
	plot!(VALS[1,:],label="ξ",line=(2.5,:blue),xlabel="Number of iterations", ylabel=ytext)
	for opnum in 1:numObs
		plot!(VALS[opnum+1,:],label=ObsNames[opnum],line=(2.5,:red))
	end

	plot!(legend=false,xguidefont=FONT,yguidefont=FONT,xtickfont=FONT,ytickfont=FONT,size=SIZE)

	return P,T,VALS,THEO
end


# just obtain theoretical estimates for Lindbladian evolution
function LindbladTheoPlots(H,ϕ0,tf,res,ObsDia;E=false,U=false,ObsAdi=false,ObsNames=false)
	Tarr = range(0,stop=tf,length=res)
	numObs = length(ObsDia)
	Hsts = length(ϕ0)
	THEO = zeros(numObs+1,res)

	T = Tarr .^2

	if E == false || U == false
		E,U = adiabatize(H)
	end

	if ObsAdi == false
		ObsAdi = zeros(Hsts,Hsts,numObs)
		for opnum in 1:numObs
			ObsAdi[:,:,opnum] = dia2en(ObsDia[opnum],U)
		end
	end

	muAdi = U' * ϕ0
	
	for i in 1:res
		println("Starting step $i of $res")
		α = T[i]
		ρα = zeros(Hsts,Hsts)
		@time for st1 in 1:Hsts
			for st2 in 1:Hsts
				ω = E[st1] - E[st2]
				ρα[st1,st2] = exp(-ω^2 * α) * muAdi[st1] * muAdi[st2]
				for opnum in 1:numObs
					THEO[opnum+1,i] += ρα[st1,st2] * ObsAdi[opnum][st2,st1]
				end
			end
		end
		THEO[1,i] = sum(abs2.(ρα))
	end

	if ObsNames == false
		ObsNames = ["O$i" for i in 1:numObs]
	end

	# =
	P = plot(Tarr,THEO[1,:],label="ξtheo",line=(3,:blue))
	for opnum in 1:numObs
		plot!(Tarr,THEO[opnum+1,:],label=ObsNames[opnum]*"theo",line=(3,:red))
	end

	FONT = font(24,"Helvetica")
	plot!(legend=false,xguidefont=FONT,yguidefont=FONT,xtickfont=FONT,ytickfont=FONT)

	return P,Tarr,THEO
end

# same as timeEvo, but already preassumes ObsAdi are going to be S0 and Ptrans1 projectors
function timeEvo3(tf,res,H,ϕ0;E=false,U=false,ObsAdi=false)
	Tarr = range(0,stop=tf,length=res)
	numObs = length(ObsAdi)
	dt = Tarr[2] - Tarr[1]
	@show dt
	Hsts = length(ϕ0)
	ϕarr = zeros(Complex,Hsts,res)
	muAdi = U' * ϕ0
	ϕarr[:,1] = muAdi
	

	println("Starting time evolution")
	@time for i in 2:res
		println("step $(i-1) of $(res-1)")
		#@time ϕarr[:,i] = expOp * ϕarr[:,i-1]
		ϕarr[:,i] = exp.((-1im * Tarr[i]) .* E) .* muAdi
		@show abs2(dot(ϕarr[:,i],ϕarr[:,1]))
	end

	#Routine includes treatment of observables 
	println("Starting observables analysis: calculating observables for each time and observable")
	ObsVals = zeros(numObs,res)
	for opnum in 1:numObs
		println("Starting observable $opnum for all time evolution")
		@time for tnum in 1:res
			ObsVals[opnum,tnum] = real(dot(ϕarr[:,tnum],ObsAdi[opnum] * ϕarr[:,tnum]))
		end
	end
	
	println("Starting time averaging and analysis")
	ρVals = zeros(numObs+1,res) #ρVals[1,:] is purity
	ρVals[1,1] = 1

	for opnum in 1:numObs
		ρVals[opnum+1,1] = ObsVals[opnum,1]
	end
	
	@time for i in 1:res-1
		ovlps = ϕarr[:,i+1]' * ϕarr[:,1:i]
		ovlps = sum(abs2.(ovlps))
		@show ovlps
		ρVals[1,i+1] = i^2 * ρVals[1,i] +2*ovlps + 1
		ρVals[1,i+1] /= (i+1)^2

		for opnum in 1:numObs
			ρVals[opnum+1,i+1] = i*ρVals[opnum+1,i] + ObsVals[opnum,i+1]
			ρVals[opnum+1,i+1] /= i+1
		end
	end

	println("Building theoretical values")
	# =
	TheoValsT = zeros(numObs+1,res)
	#TheoValsL = zeros(numObs+1,res)
	
	ρmu = muAdi * muAdi'
	for i in 1:res
		println("Building theoretical at time $i of $res")
		ρtT = zeros(Hsts,Hsts)
		#ρtL = zeros(Hsts,Hsts)
		@time for st1 in 1:Hsts
			for st2 in 1:Hsts
				ω = E[st1] - E[st2]
				sincval = sinc(ω * Tarr[i] / 2)
				#gaussval = exp(-ω^2 * Tarr[i]^2)
				ρtT[st1,st2] = sincval * ρmu[st1,st2]
				#ρtL[st1,st2] = gaussval * ρmu[st1,st2]
				#TheoVals[1,i] +=  abs2(sincval * muAdi[st1] * muAdi[st2])
				for opnum in 1:numObs
					TheoValsT[opnum+1,i] += ObsAdi[opnum][st2,st1] * ρtT[st1,st2]
					#TheoValsL[opnum+1,i] += ObsAdi[opnum][st2,st1] * ρtL[st1,st2]
				end
			end
		end
		TheoValsT[1,i] = sum(abs2.(ρtT))
		#TheoValsL[1,i] = sum(abs2.(ρtL))
	end
	
	@show TheoValsT[1,1], TheoValsT[1,end]
	@show TheoValsT[2,1], TheoValsT[2,end]
	#@show TheoValsL[1,1], TheoValsL[1,end]
	#@show TheoValsL[2,1], TheoValsL[2,end]
	# =#

	@show ρVals[1,1], ρVals[1,end]
	@show ρVals[2,1], ρVals[2,end]
	
	println("Starting plotting procedure")
	ytext = "Purity and " * L"S_0" * " and " * L"P_{trans}^{(1)}" * " population"
	P = plot(ρVals[1,:],label="ξ",line=(3,:blue),xlabel="Number of iterations", ylabel=ytext)
	for opnum in 1:numObs
		plot!(ρVals[opnum+1,:],line=(3,:red))
	end

	# =
	plot!(TheoValsT[1,:],label="ξtheo",line=(3,:blue,:dash))
	for opnum in 1:numObs
		plot!(TheoValsT[opnum+1,:],line=(3,:red,:dash))
	end

	#=
	scatter!(Tarr,TheoValsL[1,:],label="ξtheo",mark=(10,:blue,:circle))
	for opnum in 1:numObs
		scatter!(Tarr,TheoValsL[opnum+1,:],mark=(10,:red,:circle))
	end
	# =#
	plot!(legend=false,xguidefont=FONT,yguidefont=FONT,xtickfont=FONT,ytickfont=FONT,size=SIZE)

	return P,ρVals
end