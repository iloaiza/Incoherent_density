#parallel implementation for analysing many different initial seeds for Lanczos procedure

@everywhere using SharedArrays
@everywhere using LinearAlgebra
@everywhere using FastGaussQuadrature
@everywhere using SpecialFunctions
@everywhere using LaTeXStrings

@everywhere include("config.jl")
@everywhere include("general_functions.jl")
@everywhere include("rhodopsin.jl")
@everywhere include("heaviside.jl")
@everywhere include("shift.jl")



@everywhere function randomParallelAnalysis(Nrands,ϕ0,H,itsItvl,ρex=false)
	@show itsItvl
	_,_,_,Ptrans1 = projectorsHeavisideRhodopsinDiabatic()
	Ptrans1 = SharedArray(Ptrans1)
	itsItvlShared = SharedArray(itsItvl)
	ϕ0 = SharedArray(ϕ0)
	Hsts = length(ϕ0)
	σ = SharedArray([dot(ϕ0,H*ϕ0)])
	println("Inverting Hamiltonian matrix")
	Hinv = SharedArray{Float64}(Hsts,Hsts)
	@time Hinv .= inv(H - σ[1]*id(Hsts))
	if ρex == false
		println("Building exact solution")
		ρex = ϕ0 * ϕ0'
		E,U = adiabatize(H)
		ρex = dia2en(ρex,U)
		ρex = Diagonal(ρex)
		ρex = dia2en(ρex,U')
	end
	
	itsTot = length(itsItvl)

	S0 = SharedArray{Float64}(itsTot,Nrands)
	COHS = SharedArray{Float64}(itsTot,Nrands)
	TRANS = SharedArray{Float64}(itsTot,Nrands)
	
	cutSts = Int(Hsts/2)
	Xvecs = SharedArray{Float64}(Hsts,Nrands)
	Xvecs .= 2*rand(Hsts,Nrands) .- 1
	for i in 1:Nrands
		Xvecs[:,i] /= norm(Xvecs[:,i])
	end

	totIts = itsItvl[end]
	Ws = SharedArray{Float64}(Hsts,totIts,Nrands)
	αs = SharedArray{Float64}(totIts+1,Nrands)
	βs = SharedArray{Float64}(totIts+1,Nrands)
	println("Doing Lanczos iterations...")
	REMAINING = SharedArray([Nrands])
	@sync @distributed for i in 1:Nrands
		println("Starting random vector $i")
		time00 = time()
		Ws[:,:,i],αs[:,i],βs[:,i] = siLanczosDebug(Xvecs[:,i],H,σ,itsItvlShared[end],Hinv)
		#println("Starting loops")
		@time for (itNum,maxIts) in enumerate(itsItvlShared)
			#println("Starting with maxIts=$maxIts")
			X = redTtoX(Ws[:,:,i],maxIts,αs[:,i],βs[:,i])
			S0[itNum,i],COHS[itNum,i],TRANS[itNum,i] = shiftDecohere(ϕ0,σ,Hinv,X,Ptrans1)
		end
		REMAINING[1] -= 1
		println("Finished random seed $i for $(time()-time00) seconds, still remaining $(REMAINING[1]) seeds to do...")
	end

	# exact values obtained from diagonal, incoherent solution
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
	TRANSstd = zeros(itsTot)
	for i in 1:itsTot
		TRANSmean[i] = sum(TRANS[i,:]) / Nrands
		DEVS = TRANSmean[i] * ones(Nrands) - TRANS[i,:]
		TRANSstd[i] = sqrt(sum(DEVS .^2) / Nrands)
	end
	S0mean = zeros(itsTot)
	S0std = zeros(itsTot)
	for i in 1:itsTot
		S0mean[i] = sum(S0[i,:]) / Nrands
		DEVS = S0mean[i] * ones(Nrands) - S0[i,:]
		S0std[i] = sqrt(sum(DEVS .^2) / Nrands)
	end
	COHSmean = zeros(itsTot)
	COHSstd = zeros(itsTot)
	for i in 1:itsTot
		COHSmean[i] = sum(COHS[i,:]) / Nrands
		DEVS = COHSmean[i] * ones(Nrands) - COHS[i,:]
		COHSstd[i] = sqrt(sum(DEVS .^2) / Nrands)
	end
	@show TRANSmean
	@show S0mean
	@show COHSmean
	@show itsItvl
	@show TRANSstd + TRANSmean
	@show S0std + S0mean
	@show COHSstd + COHSmean

	return TRANSmean,S0mean,COHSmean
end