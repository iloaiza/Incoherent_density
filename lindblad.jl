function lindbladODE(H,Pt)
	return commutator(commutator(H,Pt,true),H,true)
end

struct krylovRho
	k :: Int #iteration number, krylov space will go all the way to 2k
	aArr :: Array{Real} #coefficients of krylov structure, size is [2k+1,2k+1]
end

function δ(j,k)
	if j == k
		return 1
	else
		return 0
	end
end

function krylovNext(ρk::krylovRho,ϵ::Real)
	k = ρk.k
	aArrplus = zeros(Float64,2k+3,2k+3)
	aArrplus[1:2k+1,1:2k+1] = ρk.aArr


	for μ in 0:2k
		for λ in 0:2k
			for j in 0:2
				aArrplus[μ+j+1,λ+3-j] -= ϵ * binomial(2,j) * ρk.aArr[μ+1,λ+1] * (-1)^j
			end
		end
	end

	return krylovRho(k+1,aArrplus)
end

function algebraicIntegration(ϕ0,H,ϵ,steps)
	Hsts = length(ϕ0)
	ϕs = zeros(Hsts,2*steps+1)
	ϕs[:,1] = ϕ0
	ρk = krylovRho(0,[1])

	for k in 1:steps
		ρk = krylovNext(ρk,ϵ)
		ϕs[:,2k] = H*ϕs[:,2k-1]
		ϕs[:,2k+1] = H*ϕs[:,2k]
	end
	println("Finished algebraic integration, starting density calculation...")
	display(ρk.aArr)

	ρst = zeros(Hsts,Hsts)
	for α in 1:2*steps+1
		for β in 1:2*steps+1
			if ρk.aArr[α,β] != 0
				ρst += ρk.aArr[α,β] * ϕs[:,α] * ϕs[:,β]'
			end
		end
	end

	println("Finished rho building")
	@show tr(ρst),tr(ρst^2),tr(commutator(H,ρst)^2)
end

function runge4_step(y0,ϵ,H)
    k1=ϵ*lindbladODE(H,y0)
    y1=y0 + k1/2
    
    k2=ϵ*lindbladODE(H,y1)
    y2=y0 + k2/2
    
    k3=ϵ*lindbladODE(H,y2)
    y3=y0 + k3
   
    k4=ϵ*lindbladODE(H,y3)

    return y0 + k1/6 + k2/3 + k3/3 + k4/6
end

function runge_ks_54(y0,ϵ,H)
    k1=ϵ*lindbladODE(H,y0)

    y2=y0 + k1/5
    k2=ϵ*lindbladODE(H,y2)

    y3=y0 + 3k1/40 + 9k2/40
    k3=ϵ*lindbladODE(H,y3)

    y4=y0 + 44k1/45 - 56k2/15 + 32k3/9
    k4=ϵ*lindbladODE(H,y4)

    y5=y0 + 19372k1/6561 - 25360k2/2187 + 64448k3/6561 - 212k4/729
    k5=ϵ*lindbladODE(H,y5)

    y6=y0 + 9017k1/3168 - 355k2/33 + 46732k3/5247 + 49k4/176 - 5103k5/18656
    k6=ϵ*lindbladODE(H,y6)

    y7=y0 + 35k1/384 + 500k3/1113 + 125k4/192 - 2187k5/6784 + 11k6/84
    k7=ϵ*lindbladODE(H,y7)    
    
    return k1,k2,k3,k4,k5,k6,k7,y7
end

function rk54_step(y0,ϵ,H;force=false,ϵmax = 1e-2,ϵmin = 1e-5,AbsTol=1e-5,RelTol=1e-5,verbose=false)
    if ϵ>ϵmax
        ϵ=ϵmax
    end
    
    k1,_,k3,k4,k5,k6,k7,znew = runge_ks_54(y0,ϵ,H)
    
    ynew=y0 + 5179k1/57600 + 7571k3/16695 + 393k4/640 - 92097k5/339200 + 187k6/2100 + k7/40
    totnorm2=sum(abs2.(ynew-znew))
    R=sqrt(totnorm2)/(RelTol*sqrt(sum(abs2.(y0)))+AbsTol)
    ds=ϵ*0.84*sqrt(sqrt(1/R))
    if verbose
        @show R,ds
    end
    if R<=1 || force
        return ϵ,znew,ds
    else
        if ds<ϵmin #dt_min defined in code_config
            println("Runge-Kutta warning: convergence not achieved for minimum timestep")
            println("Continuing, be wary of results (specially if warning repeats!)...")
            @show ds
            @show ϵmin
            rk54_step(y0,ϵmin,H,force=true)
        else
            rk54_step(y0,ds,H)
        end
    end
end

function lindbladRK4(H,ϵ,ϕ0,steps,ObsDia)
	P0 = ϕ0 * ϕ0'
	Pold = copy(P0)
	COHS = zeros(steps+1)
	T = zeros(steps+1)
	numObs = length(ObsDia)
	Hsts = length(ϕ0)

	VALS = zeros(numObs+1,steps+1)
	VALS[1,1] = sum(abs2.(P0))
	for opnum in 1:numObs
		for st1 in 1:Hsts
			for st2 in 1:Hsts
				VALS[opnum+1,1] += ObsDia[opnum][st1,st2] * P0[st2,st1]
			end
		end
	end 

	for i in 1:steps
		#println("Starting step $i in Lindblad integration")
		#P0 += ϵ * lindbladODE(H,P0)
		P0 = runge4_step(P0,ϵ,H)
		#ϵeff,P0,ϵ = rk54_step(P0,ϵ,H,AbsTol=AbsTol,RelTol=RelTol,ϵmax=emax,ϵmin=emin,verbose=verbose)
		T[i+1] = T[i] + ϵ
		VALS[1,i+1] = sum(abs2.(P0))
		for opnum in 1:numObs
			for st1 in 1:Hsts
				for st2 in 1:Hsts
					VALS[opnum+1,i+1] += ObsDia[opnum][st1,st2] * P0[st2,st1]
				end
			end
		end
		Pold = P0
	end

	return T,VALS
end

function lindbladIntegration(H,ϵ,ϕ0,steps,ObsDia;AbsTol=1e-5,RelTol=1e-5,emin=1e-5,emax=1e-2,verbose=false)
	P0 = ϕ0 * ϕ0'
	Pold = copy(P0)
	COHS = zeros(steps+1)
	TDS = zeros(steps+1)
	T = zeros(steps+1)
	numObs = length(ObsDia)
	Hsts = length(ϕ0)

	VALS = zeros(numObs+1,steps+1)
	VALS[1,1] = sum(abs2.(P0))
	for opnum in 1:numObs
		for st1 in 1:Hsts
			for st2 in 1:Hsts
				VALS[opnum+1,1] += ObsDia[opnum][st1,st2] * P0[st2,st1]
			end
		end
	end 

	for i in 1:steps
		#println("Starting step $i in Lindblad integration")
		#P0 += ϵ * lindbladODE(H,P0)
		#P0 = runge4_step(P0,ϵ,H)
		ϵeff,P0,ϵ = rk54_step(P0,ϵ,H,AbsTol=AbsTol,RelTol=RelTol,ϵmax=emax,ϵmin=emin,verbose=verbose)
		T[i+1] = T[i] + ϵeff
		VALS[1,i+1] = sum(abs2.(P0))
		@show VALS[1,i+1]
		@show T[i+1]
		#=
		for opnum in 1:numObs
			for st1 in 1:Hsts
				for st2 in 1:Hsts
					VALS[opnum+1,i+1] += ObsDia[opnum][st1,st2] * P0[st2,st1]
				end
			end
		end
		# =#
		Pold = P0
	end

	return T,VALS
end

function lindbladCompare(H,ϵ,ϕ0,steps)
	P0 = ϕ0 * ϕ0'
	Pold = copy(P0)
	COHS = zeros(steps+1)
	TDS = zeros(steps+1)
	DIFF = zeros(steps+1)

	Pex = outerProduct(ϕ0,ϕ0)
	E,U = adiabatize(H)
	Pex = Diagonal(dia2en(Pex,U))
	Pex = dia2en(Pex,U')

	COHS[1] = tr(P0^2)
	TDS[1] = tr(commutator(H,P0)^2)
	DIFF[1] = sum(abs2.(Pex-outerProduct(ϕ0,ϕ0)))

	for i in 1:steps
		#println("Starting step $i in Lindblad integration")
		#P0 += ϵ * lindbladODE(H,P0)
		#P0 = runge4_step(P0,ϵ,H)
		ϵeff,P0,ϵ = rk54_step(P0,ϵ,H)
		COHS[i+1] = tr(P0^2)
		TDS[i+1] = tr(commutator(H,P0)^2)
		DIFF[i+1] = sum(abs2.(Pex-P0))

		if i > 1
			if abs(TDS[i]) > abs(TDS[i-1])
				println("Quality of iteration decreased at step $i in time dependence, returning last entry")
				@show TDS[i+1],TDS[i]
				return COHS[i],TDS[i]
			end
			if abs(COHS[i]) > abs(COHS[i-1])
				println("Quality of iteration decreased at step $i in coherence, returning last entry")
				@show COHS[i+1],COHS[i]
				return COHS[i],TDS[i]
			end
		end
		Pold = P0
	end

	return COHS[end],TDS[end],DIFF
end


function algebraicODE(γ,hk,hL = length(hk))
	hkPlus = zeros(hL+1)
	hkPlus[2:end] = 2*γ*hk
end

function algebraicLeapfrogStep(ak)
	expDeg = length(ak)
	aknew = zeros(expDeg+2)

	aknew[1:expDeg] = ak
	aknew[2:expDeg+1] += ak
	aknew[3:expDeg+2] += ak/2

	return aknew
end

function LanczosToKrylov(αvec,βvec,LKtype=Float64)
	Hkmat = SymTridiagonal(αvec, βvec[1:end-1])
	Ntot = length(αvec)
	HM = zeros(LKtype,Ntot,Ntot,Ntot)
	HM[:,:,1] = id(Ntot)
	for m in 1:Ntot-1
		HM[:,:,m+1] = Hkmat * HM[:,:,m]
	end

	Hk0m = HM[:,1,:]

	#= Debugging routine to check correct building of Hk0m matrix
	display(Hk0m)
	Mdebug = zeros(Ntot,Ntot)
	e0 = Qmat[:,1]
	for jnum in 1:Ntot
		ej = Qmat[:,jnum]
		Mdebug[jnum,1] = dot(ej,e0)
		for mnum in 1:Ntot-1
			ej = H * ej
			Mdebug[jnum,mnum+1] = dot(ej,e0)
		end
	end
	display(Mdebug)
	error("subspace error")	
	# =#
	return Float64.(Hk0m)
end

function leapfrogLindblad(H,ϵ,steps,ϕ0)
	α,β,Qvecs,flag = lanczosBasisRestarted(ϕ0,H,2*steps,BigFloat)
	Hk0m = LanczosToKrylov(α,β,BigFloat)
	display(Hk0m)
	Hsts = length(ϕ0)
	
	ak = zeros(2steps+1)
	ak[1] = 1
	for i in 1:steps
		ak[1:2i+1] = algebraicLeapfrogStep(ak[1:2i-1])
	end

	@show ak
	for i in 1:2steps+1
		ak[i] *= ϵ^(i-1)
	end
	ak = Hk0m*ak
	@show ak
	ρ = zeros(Hsts,Hsts)
	for i in 1:2*steps+1
		ρ += ak[i] * Qvecs[:,i] * Qvecs[:,i]'
	end
	@show tr(ρ),tr(ρ^2),tr(commutator(H,ρ)^2)

	krausExp = exp(-steps*ϵ/2*H*H)
	ρfin = krausExp*ρ*krausExp

	@show tr(ρfin),tr(ρfin^2),tr(commutator(H,ρfin)^2)
end

function TaylorKraus(H,Kmax,ϕ0,t)
	matExp = exp(-H^2*t)
	Hsts = length(ϕ0)	
	Jvecs = zeros(Hsts,Kmax+1)
	NORMS = zeros(Kmax+1)

	Jvecs[:,1] = matExp * ϕ0
	NORMS[1] = norm(Jvecs[:,1])
	@show norm(Jvecs[:,1])
	#ρ = Jvecs[:,1] * Jvecs[:,1]'
	for i in 2:Kmax
		Jvecs[:,i] = sqrt(2*t/(i-1))*H*Jvecs[:,i-1]
		#Jvecs[:,i] = H*Jvecs[:,i-1]
		NORMS[i] = norm(Jvecs[:,i])
		#Jvecs[:,i] /= NORMS[i]
		@show i,NORMS[i]
		#ρ += Jvecs[:,i] * Jvecs[:,i]'
	end

	#return tr(ρ),tr(ρ^2)/(tr(ρ)^2),tr(commutator(H,ρ)^2)/tr(ρ)

	return Jvecs,NORMS
end

function KrausAnalysis(H,Kmax,ϕ0,Tarr)
	TL = length(Tarr)
	TRS = zeros(TL)
	COHS = zeros(TL)
	TDS = zeros(TL)

	for (tnum,t) in enumerate(Tarr)
		TRS[tnum],COHS[tnum],TDS[tnum] = TaylorKraus(H,Kmax,ϕ0,t)
	end

	P = plot(Tarr,TRS,label="Traces",yscale=:log10)
	plot!(Tarr,COHS,label="Coherences")
	plot!(Tarr,abs.(TDS),label="Time dependence")



	return P
end

function TaylorKrausOptimizer(H,Kmax,ϕ0,t)
	Rmat = zeros(Kmax+1,Kmax+1)
	Hsts = length(ϕ0)

	Jvecs,NORMS = TaylorKraus(H,Kmax+2,ϕ0,t)
	for ki in 1:Kmax+1
		for kj in 1:ki
			Rmat[ki,kj] = 2*(NORMS[ki+1]*dot(Jvecs[:,ki+1],Jvecs[:,kj])*NORMS[kj+1]*dot(Jvecs[:,ki],Jvecs[:,kj+1]))
			Rmat[ki,kj] -= 2*(dot(Jvecs[:,ki],Jvecs[:,kj])*NORMS[kj+1]*NORMS[ki+1]*dot(Jvecs[:,ki+1],Jvecs[:,kj+1]))
			Rmat[kj,ki] = Rmat[ki,kj]
		end
	end

	display(Rmat)
	@show matConditioning(Rmat)
	βmat = -2*inv(Rmat)
	βmat /= sum(βmat)
	ωs = zeros(Kmax+1)
	ρ = zeros(Hsts,Hsts)
	for i in 1:Kmax+1
		ωs[i] = sum(βmat[i,:])
		ρ += ωs[i] * Jvecs[:,i] * Jvecs[:,i]'
	end

	@show ωs
	T = tr(ρ)
	ρ /= T
	er,_ = adiabatize(Symmetric(ρ))
	@show er

	@show T,tr(ρ^2),tr(commutator(H,ρ)^2)
end

function rescaledKraus(H,Kmax,ϕ0,t,Emax)
	Htilde = H / Emax
	Jvecs,Jnorms = TaylorKraus(Htilde,Kmax,ϕ0,t)

	ρ = Jvecs[:,1] * Jvecs[:,1]'
	for k in 1:Kmax
		ρ += Jvecs[:,k] * Jvecs[:,k]'
	end

	@show tr(ρ)
	ρ /= tr(ρ)
	@show tr(ρ^2),tr(commutator(H,ρ)^2),tr(ρ*H)
end

