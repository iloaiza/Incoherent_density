#shift and invert variant for Lanczos algorithm

# Shift and invert Lanczos procedure with exact inversion
function siLanczosDebug(ϕ0,H,σ,maxIts,Hinv=false)
	#println("Starting shift and invert Lanczos routine...")
	Hsts = length(ϕ0)
	if Hinv == false
		#println("Building shifted Hamiltonian")
		@time Hshift = H - σ*id(Hsts)
		#println("Inverting shift operator")
		@time Hinv = inv(Hshift)
	end
	
	βs = zeros(maxIts+1)	#βs[i] = β_{i-1}
	αs = zeros(maxIts+1)	#αs[i] = α_{i-1}, α0=0
	Ws = zeros(Hsts,maxIts+1)	#Ws[:,i] = w_{i-1} and w0≡0
	r = copy(ϕ0)
	
	βs[1] = norm(r)
	#println("Starting iterations...")
	for itnum in 1:maxIts
		#println("Starting iteration $itnum")
		Ws[:,itnum+1] = r / βs[itnum]
		r = Hinv * Ws[:,itnum+1]
		r -= βs[itnum] * Ws[:,itnum]
		αs[itnum+1] = dot(Ws[:,itnum+1],r)
		r -= Ws[:,itnum+1]*αs[itnum+1]
		#println("GS orthogonalization")
		r,rflag = GrammSchmidtVector(r,Ws[:,1:itnum])
		if rflag == true
			error("WARNING! Residual vector is linearly dependent beyond tolerance at step $itnum")
		end
		βs[itnum+1] = norm(r)
	end

	#println("Finished Lanczos process, calculating Ritz pairs...")
	T = SymTridiagonal(αs[2:maxIts+1], βs[2:maxIts])
	Θ,S = adiabatize(T)
	#@show sort!(σ .+ 1 ./Θ)
	X = Ws[:,2:end] * S
	#@show sort!(diagEls(dia2en(H,X)))
	return Ws[:,2:end],αs,βs
end

struct pol
	coeffs :: Array{Float64,1}
	poltype :: String
	deg :: Int64
end
#length(coeffs) = deg+1 since coeffs includes 0th degree (as last coefficient! C[1] is of largest coefficient)

function polEval(pol,H,x0)
	Hsts = length(x0)
	xk2 = zeros(Hsts)
	xk1 = zeros(Hsts)
	xk0 = zeros(Hsts)
	if pol.poltype == "Chebyshev"
		#Horner ruler for Chebyshev (i.e. Clenshaw rule): Tn(H) = |ϕn>, then |ϕn+2> = coeffs[deg+1-k]*|ϕ0>+2H|ϕn+1>-|ϕn>
		for k in 0:pol.deg-1
			xk2 = pol.coeffs[k+1]*x0 + 2*(H*xk1) - xk0
			xk0 = xk1
			xk1 = xk2
		end
		return pol.coeffs[end]*x0+H*xk1-xk0
	else
		error("polEval not defined for pol.type=$(pol.poltype)")
	end
end

function polynomialLanczos(ϕ0,H,maxIts,pol)
	println("Starting polynomial Lanczos routine...")
	Hsts = length(ϕ0)
	
	βs = zeros(maxIts+1)	#βs[i] = β_{i-1}
	αs = zeros(maxIts+1)	#αs[i] = α_{i-1}, α0=0
	Ws = zeros(Hsts,maxIts+1)	#Ws[:,i] = w_{i-1} and w0≡0
	r = copy(ϕ0)
	
	βs[1] = norm(r)
	println("Starting iterations...")
	@time for itnum in 1:maxIts
		println("Starting iteration $itnum")
		Ws[:,itnum+1] = r / βs[itnum]
		r = polEval(pol,H,Ws[:,itnum+1])
		r -= βs[itnum] * Ws[:,itnum]
		αs[itnum+1] = dot(Ws[:,itnum+1],r)
		r -= Ws[:,itnum+1]*αs[itnum+1]
		#println("GS orthogonalization")
		r,rflag = GrammSchmidtVector(r,Ws[:,1:itnum])
		if rflag == true
			error("WARNING! Residual vector is linearly dependent beyond tolerance at step $itnum, the algebra closes! Use maxIts = itnum-1")
		end
		βs[itnum+1] = norm(r)
	end

	println("Finished Lanczos process, calculating Ritz pairs...")
	T = SymTridiagonal(αs[2:maxIts+1], βs[2:maxIts])
	Θ,S = adiabatize(T)
	@show sort!(Θ)
	X = Ws[:,2:end] * S
	@show sort!(diagEls(dia2en(H,X)))
	return X
end

# Builds a random seed
function randBuilder(ϕ0,H)
	Hsts = length(ϕ0)
	ϕ1 = H*ϕ0
	avg = dot(ϕ0,ϕ1)
	dev = sqrt(dot(ϕ1,ϕ1) - avg^2)
	Emin = avg - 3dev
	Emax = avg + 3dev
	Hdia = diagEls(H)

	x0 = 2*rand(Hsts) .-1
	for st in 1:Hsts
		if !(Emin <= Hdia[st] <= Emax)
			x0[st] = 0
		end
	end
	x0 /= norm(x0)
	return x0
end

#Pcis0,Pcis1,Ptrans0,Ptrans1 = projectorsHeavisideRhodopsinDiabatic()

function shiftDecohere(ϕ0,σ,Hinv,X0,Ptrans)
	Hsts = length(ϕ0)

	ρSI = ϕ0 * ϕ0'
	ρSI = dia2en(ρSI,X0)
	ρSI = Diagonal(ρSI)
	ρSI = dia2en(ρSI,X0')
	#@show tr(ρSI)
	ρSI /= tr(ρSI)

	#@time densityQualityTrack(ρSI,H,ρex)
	#dist0 = quickQuality(ρSI0,ρex)
	cutSts = Int(Hsts/2)
	s0pop = sum(diagEls(ρSI[1:cutSts,1:cutSts]))
	coh = sum(abs2.(ρSI))
	trans = tr(Ptrans * ρSI)
	
	return s0pop,coh,trans
end
#=
function shiftDecohereAnalysis(ϕ0,H,itsItvl,ρex=false,σ=false)
	Hsts = length(ϕ0)
	@show Eavg = dot(ϕ0,H*ϕ0)
	if σ == false
		println("Setting σ")
		σ = Eavg
	end
	println("Inverting Hamiltonian matrix")
	@time Hinv = inv(H - σ*id(Hsts))
	if ρex == false
		println("Building exact solution")
		ρex = ϕ0 * ϕ0'
		E,U = adiabatize(H)
		ρex = dia2en(ρex,U)
		ρex = Diagonal(ρex)
		ρex = dia2en(ρex,U')
	end
	
	itsTot = length(itsItvl)
	COHS = zeros(itsTot,2)
	DISTS = zeros(itsTot,2)
	CIS = zeros(itsTot,2)
	TRANS = zeros(itsTot,2)
	cisex = 0.0; cisex = 0.0; cohex = 0.0; transex = 0.0
	#@time X = siLanczos(ϕ0,H,Eavg,maxIts)
	for (itNum,maxIts) in enumerate(itsItvl)
		println("Starting with maxIts=$maxIts")
		@time cispop0,cispop1,cisex,transpop0,transpop1,transex,coh0,coh1,cohex,dist0,dist1 = shiftDecohere(ϕ0,H,maxIts,ρex,σ,Hinv)
		CIS[itNum,1] = cispop0; CIS[itNum,2] = cispop1;
		TRANS[itNum,1] = transpop0; TRANS[itNum,2] = transpop1;
		COHS[itNum,1] = coh0; COHS[itNum,2] = coh1;
		DISTS[itNum,1] = dist0; DISTS[itNum,2] = dist1;
	end

	P = plot(itsItvl,DISTS[:,1],label="Dist R",line=(2,:green))
	plot!(itsItvl,DISTS[:,2],label="Dist UF",line=(2,:red))
	plot!(itsItvl,zeros(itsTot),label="Exact",line=(2,:black))
	dmax = maximum(DISTS)
	plot!(itsItvl,cohex*ones(itsTot) .+ (dmax*1.1),label="Exact",line=(2,:black))
	plot!(itsItvl,COHS[:,1] .+ dmax*1.1, label="Coh R",line=(2,:green))
	plot!(itsItvl,COHS[:,2] .+ dmax*1.1, label="Coh UF",line=(2,:red))
	cohsmax = 1.1*maximum(COHS) + dmax*1.1
	plot!(itsItvl,cisex*ones(itsTot) .+ cohsmax,label="Exact",line=(2,:black))
	plot!(itsItvl,CIS[:,1] .+ cohsmax, label="Cis R",line=(2,:green))
	plot!(itsItvl,CIS[:,2] .+ cohsmax, label="Cis UF",line=(2,:red))

	plot!(size=[1500,1000],xlabel="Iterations",ylabel="dist / cohs / cis")

	return P
end
# =#

function redTtoX(Ws,dims,αs,βs)
	#grabs tridiagonal output from Lanczos and Q array and builds X array at given dimension
	T = SymTridiagonal(αs[2:dims+1], βs[2:dims])
	Θ,S = adiabatize(T)
	#@show sort!(σ .+ 1 ./Θ)
	X = Ws[:,1:dims] * S
	#@show sort!(diagEls(dia2en(H,X)))
	return X
end

# Compares the Franck-Condon excitation, its modified version, and a random seed for Lanczos procedure (see Ref.[1])
function shiftDecohereSuperAnalysis(ϕ0,H,itsItvl,ρex=false,σ=false)
	Hsts = length(ϕ0)
	@show Eavg = dot(ϕ0,H*ϕ0)
	if σ == false
		println("Setting σ")
		σ = Eavg
	end
	Edev = sqrt(dot(H*ϕ0,H*ϕ0) - Eavg^2)
	println("Inverting Hamiltonian matrix")
	@time Hinv = inv(H - σ*id(Hsts))
	if ρex == false
		println("Building exact solution")
		ρex = ϕ0 * ϕ0'
		E,U = adiabatize(H)
		ρex = dia2en(ρex,U)
		ρex = Diagonal(ρex)
		ρex = dia2en(ρex,U')
	end
	
	itsTot = length(itsItvl)

	COHS = zeros(itsTot,3)
	CIS = zeros(itsTot,3)
	TRANS = zeros(itsTot,3)
	
	
	cutSts = Int(Hsts/2)
	println("Building PT stuff")
	PT1 = tiptO1(diagEls(H),H[1:cutSts,cutSts+1:end])
	x4 = PT1 * ϕ0
	@show norm(x4)
	x4 /= norm(x4)

	x3 = 2*rand(Hsts) .- 1; x3 /= norm(x3)
	
	#@time X = siLanczos(ϕ0,H,Eavg,maxIts)
	println("Doing Lanczos iterations...")
	@time Wuf,αsuf,βsuf = siLanczosDebug(ϕ0,H,σ,itsItvl[end],Hinv)
	@time Wrand,αsrand,βsrand = siLanczosDebug(x3,H,σ,itsItvl[end],Hinv)
	@time Wcorr,αscorr,βscorr = siLanczosDebug(x4,H,σ,itsItvl[end],Hinv)
	println("Starting loops")
	@time for (itNum,maxIts) in enumerate(itsItvl)
		println("Starting with maxIts=$maxIts")
		Xuf = redTtoX(Wuf,maxIts,αsuf,βsuf)
		@time cisuf,cohuf,transuf = shiftDecohere(ϕ0,σ,Hinv,Xuf)
		Xrand = redTtoX(Wrand,maxIts,αsrand,βsrand)
		@time cisrand,cohrand,transrand = shiftDecohere(ϕ0,σ,Hinv,Xrand)
		Xcorr = redTtoX(Wcorr,maxIts,αscorr,βscorr)
		@time ciscorr,cohcorr,transcorr = shiftDecohere(ϕ0,σ,Hinv,Xcorr)
		
		COHS[itNum,1] = cohuf; COHS[itNum,2] = cohrand; COHS[itNum,3] = cohcorr
		CIS[itNum,1] = cisuf; CIS[itNum,2] = cisrand; CIS[itNum,3] = ciscorr
		TRANS[itNum,1] = transuf; TRANS[itNum,2] = transrand; TRANS[itNum,3] = transcorr
	end

	P = plot(itsItvl,COHS[:,1],line=(3,:blue))
	plot!(itsItvl,COHS[:,3],line=(3,:blue,:dash))
	plot!(itsItvl,COHS[:,2],line=(3,:blue,:dot))

	plot!(itsItvl,CIS[:,1],line=(3,:red))
	plot!(itsItvl,CIS[:,3],line=(3,:red,:dash))
	plot!(itsItvl,CIS[:,2],line=(3,:red,:dot))

	plot!(itsItvl,TRANS[:,1],line=(2,:black))
	plot!(itsItvl,TRANS[:,3],line=(2,:black,:dash))
	plot!(itsItvl,TRANS[:,2],line=(2,:black,:dot))

	ytext = "Purity, " * L"S_0" * " and " * L"P^{(1)}_{trans}"
	plot!(xlabel="Number of iterations", ylabel=ytext)
	@show SIZE
	plot!(ylims=(0,1),legend=false,xguidefont=FONT,yguidefont=FONT,xtickfont=FONT,ytickfont=FONT,size=SIZE)


	TRex = 0.3358250803997233
	TRANS .-= TRex; TRANS ./= TRex; TRANS .= abs.(TRANS)
	display(TRANS)

	return P,COHS,CIS,TRANS
end


# Analizes Lanczos algorithm with Nrands different random seeds.
# Parallel shared memory implementation in parallel.jl and PARALLEL_RUN.jl
function randomSuperAnalysis(Nrands,ϕ0,H,itsItvl,ρex=false,σ=false)
	Hsts = length(ϕ0)
	@show Eavg = dot(ϕ0,H*ϕ0)
	if σ == false
		println("Setting σ")
		σ = Eavg
	end
	Edev = sqrt(dot(H*ϕ0,H*ϕ0) - Eavg^2)
	println("Inverting Hamiltonian matrix")
	@time Hinv = inv(H - σ*id(Hsts))
	if ρex == false
		println("Building exact solution")
		ρex = ϕ0 * ϕ0'
		E,U = adiabatize(H)
		ρex = dia2en(ρex,U)
		ρex = Diagonal(ρex)
		ρex = dia2en(ρex,U')
	end
	
	itsTot = length(itsItvl)

	S0 = zeros(itsTot,Nrands)
	COHS = zeros(itsTot,Nrands)
	TRANS = zeros(itsTot,Nrands)
	
	cutSts = Int(Hsts/2)
	Xvecs = 2*rand(Hsts,Nrands) .- 1
	for i in 1:Nrands
		Xvecs[:,i] /= norm(Xvecs[:,i])
	end

	totIts = itsItvl[end]
	Ws = zeros(Hsts,totIts,Nrands)
	αs = zeros(totIts+1,Nrands)
	βs = zeros(totIts+1,Nrands)
	println("Doing Lanczos iterations...")
	for i in 1:Nrands
		println("Starting $i-th random vector")
		@time Ws[:,:,i],αs[:,i],βs[:,i] = siLanczosDebug(Xvecs[:,i],H,σ,itsItvl[end],Hinv)
		println("Starting loops")
		@time for (itNum,maxIts) in enumerate(itsItvl)
			println("Starting with maxIts=$maxIts")
			X = redTtoX(Ws[:,:,i],maxIts,αs[:,i],βs[:,i])
			@time S0[itNum,i],COHS[itNum,i],TRANS[itNum,i] = shiftDecohere(ϕ0,σ,Hinv,X)
		end
	end

	COHex = 0.11252708949290838
	TRex = 0.3358250803997233
	S0ex = 0.26742651972223325

	TRANS .-= TRex; TRANS ./= TRex; TRANS .= abs.(TRANS)
	COHS .-= COHex; COHS ./= COHex; COHS .= abs.(COHS)
	S0 .-= S0ex; S0 ./= S0ex; S0 .= abs.(S0)
	display(TRANS)
	display(COHS)
	display(S0)
	
	TRANSmean = zeros(itsTot)
	for i in 1:itsTot
		TRANSmean[i] = sum(TRANS[i,:])
	end
	TRANSmean ./= Nrands
	S0mean = zeros(itsTot)
	for i in 1:itsTot
		S0mean[i] = sum(S0[i,:])
	end
	S0mean ./= Nrands
	COHSmean = zeros(itsTot)
	for i in 1:itsTot
		COHSmean[i] = sum(COHS[i,:])
	end
	COHSmean ./= Nrands
	@show TRANSmean
	@show S0mean
	@show COHSmean

	return TRANSmean,S0mean,COHSmean
end


function spectralDecohere(ϕ0,H,maxIts,Emin,Emax,degPol,Gcoeff,ρex=false,method="Gauss")
	Hsts = length(ϕ0)
	ϕ1 = H*ϕ0
	@show Eavg = dot(ϕ0,ϕ1)
	@show Edev = dot(ϕ1,ϕ1) - Eavg^2

	#Build Gaussian polynomial
	Emid = (Emax+Emin)/2
	ΔE = (Emax-Emin)/2
	Htilde = (H -Emid*id(Hsts))/ΔE
	e0 = (Eavg-Emid)/ΔE
	nodes = [cos((i+0.5)*π/(degPol+1)) for i in 1:degPol+1]
	if method == "Gauss"
		G = Gcoeff * exp.(-((nodes .- e0) .^2)/(2*Edev))
	elseif method == "Inv"
		G = Gcoeff * 1 ./ (nodes .- e0)
	end
	polcoeffs = chInterpolCoeffs(degPol,nodes,G)
	#=
	pdbg = chBuild(polcoeffs)
	plot(nodes,G,label="Gaussian")
	X = range(-1,stop=1,length=1000)
	Gvals = zeros(1000)
	for i in 1:1000
		Gvals[i] = polyval(pdbg,X[i])
	end
	return plot!(X,Gvals,label="Interpol")
	# =#
	
	p = pol(polcoeffs,"Chebyshev",degPol)


	x0 = rand(Hsts); x0 /= norm(x0)
	@time X = polynomialLanczos(x0,Htilde,maxIts,p)
	if ρex == false
		ρex = ϕ0 * ϕ0'
		E,U = adiabatize(H)
		ρex = dia2en(ρex,U)
		ρex = Diagonal(ρex)
		ρex = dia2en(ρex,U')
	end

	ρSI = ϕ0 * ϕ0'
	ρSI = dia2en(ρSI,X)
	ρSI = Diagonal(ρSI)
	ρSI = dia2en(ρSI,X')
	@show tr(ρSI)
	ρSI /= tr(ρSI)

	#@time densityQualityTrack(ρSI,H,ρex)
	@time quickQuality(ρSI,ρex)
	cutSts = Int(Hsts/2)
	cispop = sum(diagEls(ρSI[1:cutSts,1:cutSts]))
	transpop = sum(diagEls(ρSI[cutSts+1:end,cutSts+1:end]))
	cisex = sum(diagEls(ρex[1:cutSts,1:cutSts]))
	transex = sum(diagEls(ρex[cutSts+1:end,cutSts+1:end]))
	@time coh = sum(abs2.(ρSI))
	@time cohex = sum(abs2.(ρex))

	@show cispop,cisex
	@show transpop,transex
	@show coh,cohex

	return ρSI
end