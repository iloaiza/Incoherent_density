#v0 is initial vector, H is matrix to generate Krylov space and N number of iterations
function krylovVecs(v0,H,N)
	nsts = length(v0)
	#Vs is matrix with all normalized Krylov vectors
	#NORM2 is vector with normalization coefficients to avoid numerical instabilities of high order vectors
	Vs = zeros(nsts,N+1)
	#NORM2 = zeros(N+1)
	#NORM2[1] = sum(abs2.(v0))
	Vs[:,1] = v0#/sqrt(NORM2[1])
	for i in 1:N
		Vs[:,i+1] = H*Vs[:,i]
		#NORM2[i+1] = sum(abs2.(Vs[:,i+1]))
		#Vs[:,i+1] /= sqrt(NORM2[i+1])
	end

	for i in 1:N
		#NORM2[i+1] *= NORM2[i]
		#Vs[:,i] *= sqrt(NORM2[i])
	end
	#Vs[:,end] *= sqrt(NORM2[end])

	return Vs
end

function krylovSubspaceBuilder(ϕ0,H,n,tol=1e-15)
	Hsts = length(ϕ0)
	Es = zeros(Hsts,n+1) #array of normalized subspace vectors
	NORMs = zeros(n+1) #array of norms of |ϕk>
	NORMs[1] = sqrt(sum(abs2.(ϕ0)))
	Es[:,1] = ϕ0 / NORMs[1]
	ϕkold = copy(ϕ0)
	ϕkplus = copy(ϕ0)
	PROJs = zeros(n,n)
	REMAINDERS = zeros(Hsts,n+1)
	tolflag = false
	for num in 1:n
		ϕkplus = H*ϕkplus
		ϕkold = ϕkplus
		for j in 1:num
			#println("Orthogonalizing ϕ$(num+1) and ϕ$j")
			PROJs[num,j] = sum(ϕkplus .* Es[:,j])
			#println("Projection value is $(PROJs[num,j])")
			ϕkplus -=  PROJs[num,j] * Es[:,j]
			#println("Dot product after projection is $(sum(ϕkplus .* Es[:,j]))")
			#println("Remaining norm of ϕ$(num+1) is $(sum(abs2.(ϕkplus)))")
		end
		REMAINDERS[:,num+1] = ϕkold - ϕkplus
		NORMs[num+1] = sqrt(sum(abs2.(ϕkplus)))
		if NORMs[num+1] < tol && !tolflag
			tolflag = true
			println("Warning, Krylov subspace operating trivially after $num repetitions")
		end
		Es[:,num+1] = ϕkplus / NORMs[num+1]
	end

	# = sanity checking subroutine
	println("Starting krylovSubspaceBuilder sanity check subroutine...")
	OVLPs = zeros(n+1,n+1)
	for n1 in 1:n+1
		for n2 in 1:n1
			OVLPs[n1,n2] = sum(Es[:,n1] .* Es[:,n2])
			OVLPs[n2,n1] = OVLPs[n1,n2]
		end
	end

	println("Overlap matrix (should be identity) is:")
	@show sum(abs2.(OVLPs - id(n+1)))
	if sum(abs2.(OVLPs - id(n+1))) > 1e-3
		#error("Krylov subspace not orthogonal anymore, not converging...")
	end
	#display(OVLPs)
	println("Norm of krylov subspace non-normalized vectors is:")
	@show maximum(NORMs)
	#display(NORMs)
	#println("Pojections matrix is:")
	#display(PROJs)
	println("End of sanity check subroutine")
	# =# # end of sanity check

	return NORMs,Es,PROJs,REMAINDERS
end

function lanczosBasis(ϕ0,H,n,FType = BigFloat,tol=1e-12)
	println("Finding Lanczos basis for $n iterations of H matrix on ϕ0")
	Hsts = length(ϕ0)
	@show Hsts
	#Qvecs are normalized vectors
	Qvecs = zeros(FType,Hsts,n+1)
	Norms = zeros(FType,n+1)
	Norms[1] = sqrt(dot(big.(ϕ0),big.(ϕ0)))
	Qvecs[:,1] = big.(ϕ0) / Norms[1]

	r = H*Qvecs[:,1]
	α = zeros(FType,n+1)	
	α[1] = dot(Qvecs[:,1],r)
	r -= α[1]*Qvecs[:,1]
	β = zeros(FType,n+1)
	β[1] = sqrt(dot(r,r))
	for k in 1:n
		qold = Qvecs[:,k]
		Qvecs[:,k+1] = r/β[k]
		r = H*Qvecs[:,k+1] - β[k]*qold
		α[k+1] = dot(Qvecs[:,k+1],r)
		r -= α[k+1]*Qvecs[:,k+1]
		β[k+1] = sqrt(dot(r,r))
		if β[k+1] < tol
			println("Lanczos invariant subspace found after $k iterations, returning reduced basis...")
			ηarr = zeros(n)
			for nnum in 1:k
				ηarr[nnum] = abs2(dot(Qvecs[:,1],Qvecs[:,nnum+1]))
			end
			println("Lanczos orthogonality test")
			@show sum(ηarr),maximum(ηarr),minimum(ηarr)
			return α[1:k+1],β[1:k+1],Qvecs[:,1:k+1],true
		end
	end

	# = SANITY CHECKING ROUTINE
	ηarr = zeros(n)
	for nnum in 1:n
		ηarr[nnum] = abs2(dot(Qvecs[:,1],Qvecs[:,nnum+1]))
	end
	println("Lanczos orthogonality test")
	@show sum(ηarr),maximum(ηarr),minimum(ηarr)
	#= Further debugging...
	S = zeros(n+1,n+1)
	for n1 in 1:n
		for n2 in 1:n
			S[n1,n2] = dot(Qvecs[:,n1],Qvecs[:,n2])
		end
	end
	@show sum(abs2.(S - id(n+1)))
	display(S)
	display(Qvecs)
	# =# # =#

	return α,β,Qvecs,false
end

function LanczosBasisRestarted(ϕ0,H,n;FType = BigFloat,βtol=1e-7,ηtol=1e-12,RestartTol=1e-7,GStol=1e-8)
	println("Finding Lanczos basis for $n iterations of H matrix on ϕ0")
	Hsts = length(ϕ0)
	@show Hsts
	#Qvecs are normalized vectors
	Qvecs = zeros(FType,Hsts,n+1)
	Norms = zeros(FType,n+1)
	Norms[1] = sqrt(dot(big.(ϕ0),big.(ϕ0)))
	Qvecs[:,1] = big.(ϕ0) / Norms[1]

	r = H*Qvecs[:,1]
	α = zeros(FType,n+1)	
	α[1] = dot(Qvecs[:,1],r)
	r -= α[1]*Qvecs[:,1]
	β = zeros(FType,n+1)
	β[1] = sqrt(dot(r,r))
	ηarr = zeros(n)
	Rcount = 0
	for k in 1:n
		qold = Qvecs[:,k]
		Qvecs[:,k+1] = r/β[k]
		ηarr[k] = abs2(dot(Qvecs[:,1],Qvecs[:,k+1]))
		if ηarr[k] > ηtol
			Rcount += 1
			#println("Reorthogonalizing at $k-th iteration...")
			println("GS routine")
			Qvecs[:,k+1],GSflag = GrammSchmidtVector(Qvecs[:,k+1],Qvecs[:,1:k],GStol)
			#for st in 1:k
			#	Qvecs[:,k+1] -= dot(Qvecs[:,k+1],Qvecs[:,st]) * Qvecs[:,st]
			#end
			if norm(Qvecs[:,k+1]) < RestartTol
				println("Lanczos invariant subspace found after $k iterations, returning reduced basis...")
				ηarr = zeros(n)
				for nnum in 1:k
					ηarr[nnum] = abs2(dot(Qvecs[:,1],Qvecs[:,nnum+1]))
				end
				println("Lanczos orthogonality test")
				@show sum(ηarr),maximum(ηarr),minimum(ηarr)
				return α[1:k+1],β[1:k+1],Qvecs[:,1:k+1],true

				ηarr[k] = abs2(dot(Qvecs[:,1],Qvecs[:,k+1]))
				Qvecs[:,k+1] /= sqrt(dot(Qvecs[:,k+1],Qvecs[:,k+1]))
			end
		end
		r = H*Qvecs[:,k+1] - β[k]*qold
		α[k+1] = dot(Qvecs[:,k+1],r)
		r -= α[k+1]*Qvecs[:,k+1]
		β[k+1] = sqrt(dot(r,r))
		if β[k+1] < βtol
			println("Lanczos invariant subspace found after $k iterations, returning reduced basis...")
			ηarr = zeros(n)
			for nnum in 1:k
				ηarr[nnum] = abs2(dot(Qvecs[:,1],Qvecs[:,nnum+1]))
			end
			println("Lanczos orthogonality test")
			@show sum(ηarr),maximum(ηarr),minimum(ηarr)
			return α[1:k+1],β[1:k+1],Qvecs[:,1:k+1],true
		end
	end
	println("Lanczos iteration finished, restarded $Rcount times")

	#= SANITY CHECKING ROUTINE
	println("Lanczos orthogonality test")
	@show sum(ηarr),maximum(ηarr),minimum(ηarr)
	# = Further debugging...
	S = zeros(n+1,n+1)
	for n1 in 1:n+1
		for n2 in 1:n+1
			S[n1,n2] = dot(Qvecs[:,n1],Qvecs[:,n2])
		end
	end
	@show sum(abs2.(S - id(n+1)))
	display(S)
	display(Qvecs)
	# =# # =#

	return α,β,Qvecs,false
end

function lanczosBlock(Q0,H,n,Ftype=Float64,dtol=1e-8)
	Hsts = length(Q0[:,1])
	p = length(Q0[1,:])
	Rj = zeros(Ftype,Hsts,p)
	Xmats = zeros(Ftype,Hsts,p,n+1)
	Mjs = zeros(Ftype,p,p,n+1)
	Bjs = zeros(Ftype,p,p,n+1)
	PCs = zeros(Int,n+1)
	PCs .= p
	ηs = zeros(n)
	Qkry = zeros(Ftype,Hsts,p*(n+1))

	Qkry[:,1:p] = Q0

	Xmats[:,:,1] = Q0

	pcount = p
	for j in 1:n
		Rj = H*Xmats[:,:,j]
		Mjs[:,:,j] = Xmats[:,:,j]' * Rj
		Rj -= Xmats[:,:,j] * Mjs[:,:,j]
		if j>1
			Rj -= Xmats[:,:,j-1] * Bjs[:,:,j]
		end
		F = qr(Rj)
		Xmats[:,:,j+1] = Matrix(F.Q)
		Bjs[:,:,j+1] = F.R
		#convergence stuff
		ηs[j] = sum(abs.(Xmats[:,:,j+1]' * Xmats[:,:,1]))
		if ηs[j] > dtol
			@show ηs[j]
			println("Orthogonality lost at j=$j, reorthogonalizing...")
			for pnum in 1:PCs[j]
				Xmats[:,pnum,j+1] = GrammSchmidtVector!(Xmats[:,pnum,j+1],Qkry[:,1:pcount])
			end
			Xmats[:,1:PCs[j],j+1] = GrammSchmidt(Xmats[:,1:PCs[j],j+1])
			ηs[j] = sum(abs.(Xmats[:,:,j+1]' * Xmats[:,:,1]))
			@show ηs[j]
		end
		Qkry[:,pcount+1:pcount+PCs[j]] = Xmats[:,1:PCs[j],j+1]
		pcount += PCs[j]
	end

	return Qkry,Mjs,Bjs
end 


function lanczosBand(Q0,H,n,Ftype=Float64,dtol=1e-8)
	bdims = length(Q0[1,:]) #initial band with
	Hsts = length(Q0[:,1])
	V = zeros(Ftype,Hsts,bdims*(n+1))
	#1
	V[:,1:bdims] = Q0
	#2
	pc = bdims
	Tmat = zeros(Ftype,bdims*n,bdims*n)
	tvals = zeros(Ftype,bdims*n,bdims*(n+1))
	svals = zeros(Ftype,bdims*n,bdims*(n+1))
	Ical = Int[]

	for j in 1:bdims*n
		println("Starting iteration $j")
		#3
		vnorm = sqrt(dot(V[:,j],V[:,j]))
		#4
		if vnorm < dtol
			println("Warning, deflation should happen, implement!")
			if j-pc > 0
				push!(Ical,j-pc)
				pc -= 1
				if pc == 0
					break
				end
				for k in j:j+pc-1
					V[:,k] = V[:,k+1]
				end
			end
		end #deflation routine
		#5
		tvals[j,j-pc+bdims] = sqrt(dot(V[:,j],V[:,j]))
		V[:,j] /= tvals[j,j-pc+bdims]
		#6
		for k in j+1:j+pc-1
			tvals[j,k-pc+bdims] = dot(V[:,j],V[:,k])
			V[:,k] -= tvals[j,k-pc+bdims] * V[:,j]
		end
		#7
		V[:,j+pc] = H*V[:,j]
		#8
		k0 = maximum([1,j-pc])
		for k in k0:j-1
			tvals[k,j] = conj(tvals[j,k])
			V[:,j+pc] -= tvals[k,j] * V[:,k]
		end
		#9
		usort = sort([Ical...,j])
		for k in usort
			tvals[k,j] = dot(V[:,k],V[:,j+pc])
			V[:,j+pc] -= tvals[k,j] * V[:,k]
		end
		#10
		for k in Ical
			svals[j,k] = conj(tvals[k,j])
		end
		#11
		Tmat = tvals[1:j,1+bdims:j+bdims] + svals[1:j,1+bdims:j+bdims]
	end

	return Tmat,V
end

function krylovTopN(Kmax,H,ϕ0,N,ltype=Float64)
	Dfacts = arrayTopN(abs2.(ϕ0),N)
	@show Dfacts,ϕ0[Dfacts]
	Hsts = length(ϕ0)
	ALPHAS = zeros(N,Kmax+1)
	BETAS = zeros(N,Kmax+1)
	QVECS = zeros(N,Hsts,Kmax+1)
	Qmat = zeros(Hsts,N*(Kmax+1))
	FLAGS = zeros(Bool,N)

	for nnum in 1:N
		vn = zeros(Hsts)
		vn[Dfacts[nnum]] = 1
		for n2 in 2:nnum
			vn = GrammSchmidtVector!(vn,QVECS[nnum-1,:,:])
		end
		@show vn
		ALPHAS[nnum,:],BETAS[nnum,:],QVECS[nnum,:,:],FLAGS[nnum] = lanczosBasisRestarted(vn,H,Kmax,ltype)
		Qmat[:,1+(nnum-1)*(Kmax+1):nnum*(Kmax+1)] = QVECS[nnum,:,:]
	end
	if sum(FLAGS) != 0
		println("Wargning, some Lanczos iteration is not converged...")
	end
	Qmat = GrammSchmidt(Qmat)

	Hred = Symmetric(dia2en(H,Qmat))
	Ered,Ured = adiabatize(Hred)
	@show Ered
	Qeig = Qmat * Ured
	Pμred = dia2en(outerProduct(ϕ0,ϕ0),Qeig)
	ρred = Diagonal(Pμred)
	@show coh = tr(ρred^2)
	ρfull = dia2en(ρred,Qeig')
	@show tr((ρfull-PμrhodDia)^2)

	return 0
end

function krylovAnalysis(Karr,H,ϕ0,ltype=Float64,E=false,U=false)
	if E == false || U == false
		println("Diagonalizing H...")
		@time E,U = adiabatize(H)
	end
	KL = length(Karr)
	println("Building Lanczos basis")
	@time αvec,βvec,Qvecs,Lflag = lanczosBasisRestarted(ϕ0,H,Karr[end],ltype)
	println("Building Pμ")
	@time Pμ = outerProduct(ϕ0,ϕ0)
	COHS = zeros(KL)
	TDS = zeros(KL)
	DISTS = zeros(KL)

	ρsol = Diagonal(dia2en(Pμ,U))
	cohreal = tr(ρsol^2)
	RHOS = Array{Float64}[]

	for (knum,k) in enumerate(Karr)
		print("Starting k=$k; ")
		time00 = time()
		Htri = SymTridiagonal(αvec[1:k+1],βvec[1:k])
		Etri,Utri = adiabatize(Htri)
		Qeig = Qvecs[:,1:k+1] * Utri
		Pμred = dia2en(Pμ,Qeig)
		ρtriLanc = Diagonal(Pμred)
		COHS[knum] = tr(ρtriLanc^2)
		ρtriFull = dia2en(ρtriLanc,Qeig'*U)
		TDS[knum] = tr(commutator(ρtriFull,Diagonal(E))^2)
		DISTS[knum] = tr((ρsol - ρtriFull)^2)
		push!(RHOS,ρtriLanc)
		println("finished after $(round(time() - time00,digits=3)) seconds")
	end
	#=
	P = plot(Karr,DISTS,yscale=:log10,ylabel="|exact-approx|^2",xlabel="Number of Lanczos iterations",title="Convergence to exact density vs NL",label="|ρApp-ρEx|^2")
	plot!(Karr,abs.(COHS .- cohreal),yscale=:log10,label="|cohApp-cohEx|")
	# =#

	return Karr,RHOS,αvec,βvec,Qvecs,DISTS,COHS,cohreal,ρsol
end



function rhoFromLanczos(Kmax,H,ϕ0,Ltype = BigFloat)
	@time αvec,βvec,Qvecs,Lflag = lanczosBasisRestarted(ϕ0,H,Kmax,Ltype)
	Htri = SymTridiagonal(αvec,βvec[1:end-1])

	if Lflag == true
		Kmax = length(αvec) - 1
		println("Changed Kmax to $Kmax since Lanczos basis converged")
	end

	println("Building Pμ")
	@time Pμ = outerProduct(ϕ0,ϕ0)
	println("Diagonalizing Lanzcos basis")
	@time Etri,Utri = adiabatize(Float64.(Htri))
	Qeig = Float64.(Qvecs) * Utri
	Pμ = dia2en(Pμ,Qeig)

	ρtriLanc = Diagonal(Pμ)
	#@show tr(ρtriLanc),tr(ρtriLanc^2),tr(commutator(Diagonal(Etri),ρtriLanc)^2),tr(ρtriLanc*Diagonal(Etri))

	#return 0

	ρtri = dia2en(Diagonal(Pμ),Qeig')
	#@show tr(ρeig),tr(ρeig^2),tr(commutator(H,ρeig)^2)

	return ρtri

	println("Finished ρtri calculation")

	#= Real solution routine
	E,U = adiabatize(H)
	ρtri = dia2en(ρtri,U)
	PμEn = dia2en(outerProduct(ϕ0,ϕ0),U)
	ρEn = Diagonal(PμEn)

	@show tr(ρtri),tr(ρEn)
	@show tr(ρtri^2),tr(ρEn^2)
	@show tr(commutator(Diagonal(E),ρtri)^2),tr(commutator(Diagonal(E),ρEn)^2),tr(commutator(Diagonal(E),PμEn)^2)
	@show tr(ρtri*Diagonal(E)),tr(ρEn*Diagonal(E))
	@show tr((ρtri - ρEn)^2),tr((PμEn - ρEn)^2)

	return 0
	# =#
end

function rhoFromHLanczosDebug(Kmax,H,E,U,ϕ0)
	@time αvec,βvec,Qvecs,Lflag = lanczosBasis(ϕ0,H,Kmax)
	Htri = SymTridiagonal(αvec,βvec[1:end-1])

	if Lflag == true
		Kmax = length(αvec) - 1
		println("Changed Kmax to $Kmax since Lanczos basis converged")
	end

	println("Building Pμ")
	@time Pμ = outerProduct(ϕ0,ϕ0)
	println("Diagonalizing Lanzcos basis")
	@time Etri,Utri = adiabatize(Float64.(Htri))
	Qeig = Float64.(Qvecs) * Utri
	Pμ = dia2en(Pμ,Qeig)

	ρtriLanc = Diagonal(Pμ)
	@show tr(ρtriLanc),tr(ρtriLanc^2),tr(commutator(Diagonal(Etri),ρtriLanc)^2),tr(ρtriLanc*Diagonal(Etri))


	ρtri = dia2en(Diagonal(Pμ),Qeig')
	#@show tr(ρeig),tr(ρeig^2),tr(commutator(H,ρeig)^2)

	#return ρtri

	println("Finished ρtri calculation")

	# = Real solution routine
	ρtri = dia2en(ρtri,U)
	PμEn = dia2en(outerProduct(ϕ0,ϕ0),U)
	ρEn = Diagonal(PμEn)

	@show tr(ρtri),tr(ρEn)
	@show tr(ρtri^2),tr(ρEn^2)
	@show tr(commutator(Diagonal(E),ρtri)^2),tr(commutator(Diagonal(E),ρEn)^2),tr(commutator(Diagonal(E),PμEn)^2)
	@show tr(ρtri*Diagonal(E)),tr(ρEn*Diagonal(E))
	@show tr((ρtri - ρEn)^2),tr((PμEn - ρEn)^2)

	return 0
	# =#
end



function kjDictBuilder(Nk)
	KJdict = Dict()

	kcount = 0
	for knum in 1:Nk
		for jnum in 1:knum
			kcount += 1
			KJdict[[knum,jnum]] = kcount
			KJdict[[jnum,knum]] = kcount
			KJdict[kcount] = [knum,jnum]
		end
	end

	return KJdict
end

function krylovSubspaceCoeffs(ϕ0,H,Nk,Nn)
	if iseven(Nn)
		println("Nn coefficients obtainable to one higher order, changing Nn to $(Nn+1)")
		Nn += 1
	end
	Nspace = Int((Nn-1)/2)
	Nkmod = Nk + Nspace
	println("Building Lanczos basis of order $Nkmod")
	@time αvec,βvec,Qvecs,_ = lanczosBasis(ϕ0,H,Nkmod)
	Kdict = kjDictBuilder(Nkmod+1)
	Ntot = Int(0.5*(Nkmod^2+Nkmod))
	NtotOrig = Int(0.5*(Nk^2+Nk))
	@show Ntot,NtotOrig
	Mkjn = zeros(Ntot,Nn+1)
	#we're only interested in filling up to the Ntot,Nn space. Each +2 iteration needs Mkj(n) = f({M(k+1)α}) for α⩽j+1.
	#Mkjn array will then be filled only to Ntot-Nk-(Nk-1)-...-(Nk-Nn iteration) in Nn coordinate

	#Mkj(0) and Mkj(1)
	println("Filling Mkj(0) and Mkj(1)")
	@time for kcount in 1:Ntot
		knum= Kdict[kcount][1]
		jnum= Kdict[kcount][2]
		if knum == jnum
			Mkjn[kcount,1] = 1
			Mkjn[kcount,2] = αvec[jnum]
		elseif knum == jnum+1
			Mkjn[kcount,2] = βvec[jnum]
		elseif knum == jnum-1
			Mkjn[kcount,2] = βvec[jnum-1]
		end
	end

	Ntoteff = Ntot
	Nkeff = Nk + Nspace
	println("Filling the rest of Mkj")
	@time for nnum in 1:Nspace
		Ntoteff -= Nkeff 
		Nkeff -= 1
		for nrep in 0:1
			for kcount in 1:Ntoteff
				knum = Kdict[kcount][1]
				jnum = Kdict[kcount][2]
				n = 2nnum-1+nrep

				Mkjn[kcount,n+2] += αvec[knum]*αvec[jnum]*Mkjn[kcount,n]
				
				kplusjplusnum = Kdict[[knum+1,jnum+1]]
				Mkjn[kcount,n+2] += βvec[knum]*βvec[jnum]*Mkjn[kplusjplusnum,n]

				kplusjnum = Kdict[[knum+1,jnum]]
				Mkjn[kcount,n+2] += βvec[knum]*αvec[jnum]*Mkjn[kplusjnum,n]

				kjplusnum = Kdict[[knum,jnum+1]]
				Mkjn[kcount,n+2] += αvec[knum]*βvec[jnum]*Mkjn[kjplusnum,n]

				if knum > 1
					kminusjplusnum = Kdict[[knum-1,jnum+1]]
					Mkjn[kcount,n+2] += βvec[knum-1]*βvec[jnum]*Mkjn[kminusjplusnum,n]

					kminusjnum = Kdict[[knum-1,jnum]]
					Mkjn[kcount,n+2] += βvec[knum-1]*αvec[jnum]*Mkjn[kminusjnum,n]
				end

				if jnum > 1
					kplusjminusnum = Kdict[[knum+1,jnum-1]]
					Mkjn[kcount,n+2] += βvec[knum]*βvec[jnum-1]*Mkjn[kplusjminusnum,n]

					kjminusnum = Kdict[[knum,jnum-1]]
					Mkjn[kcount,n+2] += αvec[knum]*βvec[jnum-1]*Mkjn[kjminusnum,n]
					
					if knum > 1
						kminusjminusnum = Kdict[[knum-1,jnum-1]]
						Mkjn[kcount,n+2] += βvec[knum-1]*βvec[jnum-1]*Mkjn[kminusjminusnum,n]
					end
				end
			end
		end
	end

	@show Ntoteff,NtotOrig
	
	#= Debug routine
	display(Mkjn)
	Mdebug = zeros(Ntoteff,Nn+1)
	kcount = 0
	for knum in 1:Nk
		for jnum in 1:knum
			kcount += 1
			ek = Es[:,knum]
			ej = Es[:,jnum]
			Mdebug[kcount,1] = dot(ek,ej)
			for nnum in 1:Nn
				ej = H * ej
				Mdebug[kcount,nnum+1] = dot(ek,ej)
			end
		end
	end
	display(Mdebug)
	display(Mkjn[1:Ntoteff,:])
	error("subspace error")

	# =#
	return Mkjn[1:Ntoteff,:],αvec,βvec,Qvecs
end

function krylovCoeffs(MU,H,N)
	mu0 = MU[:,1]
	Mvecs = krylovVecs(mu0,H,2N+1)
	M = zeros(2N+2)
	for i in 1:2N+2
		M[i] = mu0' * Mvecs[:,i]
	end

	return M
end


function vec2mat(Avec,Nk)
	if length(Avec) != Nk^2
		error("Wrong dimensions for vector to matrix conversion")
	end
	return reshape(Avec,(Nk,Nk))
end



function krylovMats(H,Nk,Mcoeffs)
	Ntot = Int((Nk^2+Nk)/2)
	
	Avec = zeros(Ntot)
	icount = 0
	for ik in 1:Nk
		for ij in 1:ik
			icount += 1
			Avec[icount] = Mcoeffs[ik]*Mcoeffs[ij]
		end
	end
	
	Rmat = zeros(Ntot,Ntot)
	kcount = 0
	for numk in 1:Nk
		for numj in 1:numk
			kcount += 1
			lcount = 0
			for numl in 1:Nk
				for numm in 1:numl
					lcount += 1
					Rmat[kcount,lcount] = Mcoeffs[numk+numm-1]*Mcoeffs[numj+numl-1] + Mcoeffs[numk+numl-1]*Mcoeffs[numj+numm-1]
				end
			end
		end
	end

	#=
	Imat = zeros(Ntot,Ntot)
	kcount = 0
	for numk in 1:Nk
		for numj in 1:numk
			kcount += 1
			lcount = 0
			for numl in 1:Nk
				for numm in 1:numl
					lcount += 1
					Imat[kcount,lcount] = Mcoeffs[numk+numl-1]*Mcoeffs[numj+numm-1] - Mcoeffs[numk+numm-1]*Mcoeffs[numj+numl-1]
				end
			end
		end
	end
	ImatRed = zeros(Ntot-Nk,Ntot-Nk)
	=#


	return Avec,Rmat
end

#krylovStructure has (k,j) coordinates and coefficient Ckj for each one. This corresponds to Ckj*Mk*Mj
import Base.+
import Base.zero

struct krylovStruct
    COORDS :: Array{Tuple{Int64,Int64},1}
    COEFFS :: Array{Int64,1}
    Ncords :: Int64
end

struct krylovOrthoStruct
	COORDS :: Array{Tuple{Int64,Int64,Int64,Int64,Int64,Int64},1}
	COEFFS :: Array{Float64,1}
	Ncords :: Int64
end

function zero(A::krylovStruct)
	return krylovStruct([(0,0)],[0],0)
end

function zero(A::krylovOrthoStruct)
	return krylovOrthoStruct([(0,0,0,0,0,0)],[0],0)
end

function zero(::Type{krylovStruct})
	return krylovStruct([(0,0)],[0],0)
end

function zero(::Type{krylovOrthoStruct})
	return krylovOrthoStruct([(0,0,0,0,0,0)],[0],0)
end

function krylovStructBuilder(COORDS,COEFFS)
	Ncoords = length(COORDS)
	if Ncoords != length(COEFFS)
		error("Different number of coefficients and coordinates when building Krylov structure")
	end

	coordCount = 0
	for coord1 in COORDS
		coordCount += 1
		if coord1[2] > coord1[1]
			#println("second coordinate is larger than first coordinate for Krylov structure in coordinate $coordCount, inverting...")
			COORDS[coordCount] = (coord1[2],coord1[1])
		end
	end

	for n1 in 1:Ncoords
		for n2 in n1+1:Ncoords
			if COORDS[n1] == COORDS[n2]
				#println("Repeating coordinate, merging coefficients n1=$n1 and n2=$n2...")
				COEFFS[n1] += COEFFS[n2]
				deleteat!(COORDS,n2)
				deleteat!(COEFFS,n2)
				return krylovStructBuilder(COORDS,COEFFS)
			end
		end
	end

	if Ncoords > 0
		for (n1,coeff) in enumerate(COEFFS)
			if coeff == 0
				deleteat!(COORDS,n1)
				deleteat!(COEFFS,n1)
				return krylovStructBuilder(COORDS,COEFFS)
			end
		end
	end

	return krylovStruct(COORDS,COEFFS,Ncoords)
end

function krylovOrthoStructBuilder(COORDS,COEFFS)
	Ncoords = length(COORDS)
	if Ncoords != length(COEFFS)
		@show COORDS
		@show COEFFS
		error("Different number of coefficients and coordinates when building Krylov structure")
	end

	coordCount = 0
	for coord1 in COORDS
		coordCount += 1
		if coord1[2] > coord1[1]
			#println("second coordinate is larger than first coordinate for Krylov structure in coordinate $coordCount, inverting...")
			COORDS[coordCount] = (coord1[2],coord1[1],coord1[3],coord1[4],coord1[5],coord1[6])
		end
		if coord1[5] > coord1[4]
			#println("second coordinate of second pair is larger than first coordinate for Krylov structure in coordinate $coordCount, inverting...")
			COORDS[coordCount] = (coord1[1],coord1[2],coord1[3],coord1[5],coord1[4],coord1[6])
		end
	end

	for n1 in 1:Ncoords
		for n2 in n1+1:Ncoords
			if COORDS[n1] == COORDS[n2]
				#println("Repeating coordinate, merging coefficients n1=$n1 and n2=$n2...")
				COEFFS[n1] += COEFFS[n2]
				deleteat!(COORDS,n2)
				deleteat!(COEFFS,n2)
				return krylovOrthoStructBuilder(COORDS,COEFFS)
			end
		end
	end

	if Ncoords > 0
		for (n1,coeff) in enumerate(COEFFS)
			if coeff == 0
				deleteat!(COORDS,n1)
				deleteat!(COEFFS,n1)
				return krylovOrthoStructBuilder(COORDS,COEFFS)
			end
		end
	end

	return krylovOrthoStruct(COORDS,COEFFS,Ncoords)
end


function +(A::krylovStruct,B::krylovStruct)
	nCoords = copy(A.COORDS)
	nCoeffs = copy(A.COEFFS)
	for i in 1:B.Ncords
		push!(nCoords,B.COORDS[i])
		push!(nCoeffs,B.COEFFS[i])
	end

	return krylovStructBuilder(nCoords,nCoeffs)
end

function +(A::krylovOrthoStruct,B::krylovOrthoStruct)
	nCoords = copy(A.COORDS)
	nCoeffs = copy(A.COEFFS)
	for i in 1:B.Ncords
		push!(nCoords,B.COORDS[i])
		push!(nCoeffs,B.COEFFS[i])
	end

	return krylovOrthoStructBuilder(nCoords,nCoeffs)
end

import Base.*

function *(A::krylovOrthoStruct,α::Number)
	return krylovOrthoStruct(A.COORDS,α .* A.COEFFS, A.Ncords)
end

function *(α::Number,A::krylovOrthoStruct)
	return krylovOrthoStruct(A.COORDS,α .* A.COEFFS, A.Ncords)
end

function krylovSuccessor(A::krylovStruct)
	newCOORDS = Tuple{Int64,Int64}[]
	newCOEFFS = Int64[]

	for i in 1:A.Ncords
		oldCoord = A.COORDS[i]
		oldCoeff = A.COEFFS[i]
		push!(newCOORDS,(oldCoord[1]+1,oldCoord[2]+1))
		push!(newCOORDS,(oldCoord[1]+2,oldCoord[2]))
		push!(newCOORDS,(oldCoord[1],oldCoord[2]+2))
		push!(newCOEFFS,oldCoeff*2)
		push!(newCOEFFS,-oldCoeff)
		push!(newCOEFFS,-oldCoeff)
	end

	return krylovStructBuilder(newCOORDS,newCOEFFS)
end

function krylovSuccessor(A::krylovOrthoStruct)
	newCOORDS = Tuple{Int64,Int64,Int64,Int64,Int64,Int64}[]
	newCOEFFS = Float64[]

	for i in 1:A.Ncords
		oldCoord = A.COORDS[i]
		oldCoeff = A.COEFFS[i]
		push!(newCOORDS,(oldCoord[1], oldCoord[2], oldCoord[3] + 1, oldCoord[4], oldCoord[5], oldCoord[6] + 1))
		push!(newCOEFFS,oldCoeff*2)

		push!(newCOORDS,(oldCoord[1], oldCoord[2], oldCoord[3] + 2, oldCoord[4], oldCoord[5], oldCoord[6]))
		push!(newCOEFFS,-oldCoeff)

		push!(newCOORDS,(oldCoord[1], oldCoord[2], oldCoord[3], oldCoord[4], oldCoord[5], oldCoord[6] + 2))
		push!(newCOEFFS,-oldCoeff)
	end

	return krylovOrthoStructBuilder(newCOORDS,newCOEFFS)
end

function krylovValueFromStructure(A::krylovStruct,Mcoeffs)
	val = 0.0
	for i in 1:A.Ncords
		ik = A.COORDS[i][1]
		ij = A.COORDS[i][2]
		val += A.COEFFS[i]*Mcoeffs[ik+1]*Mcoeffs[ij+1]
	end

	return val
end

function krylovValueFromStructure(A::krylovOrthoStruct,Mkjn,Kdict)
	val = 0.0
	for i in 1:A.Ncords
		iα = A.COORDS[i][1] + 1
		iβ = A.COORDS[i][2] + 1
		n1 = A.COORDS[i][3] + 1
		iγ = A.COORDS[i][4] + 1
		iδ = A.COORDS[i][5] + 1
		n2 = A.COORDS[i][6] + 1
		αβnum = Kdict[[iα,iβ]]
		γδnum = Kdict[[iγ,iδ]]
		val += A.COEFFS[i]*Mkjn[αβnum,n1]*Mkjn[γδnum,n2]
	end

	return val
end

function krylovAvecBuilder(Nk,Nad,Ntot = Int((Nk^2+Nk)/2))
	Astruct = zeros(krylovStruct,Ntot)
	Aold = zeros(krylovStruct,Ntot)

	#build A0 structure
	icount = 0
	for ik in 0:Nk-1
		for ij in 0:ik
			icount += 1
			Aold[icount] = krylovStructBuilder([(ik,ij)],[1])
		end
	end

	Astruct += Aold
	for nrec in 1:Nad
		for i in 1:Ntot
			Aold[i] = krylovSuccessor(Aold[i])
		end
		Astruct += Aold
	end

	return Astruct
end

function krylovAvecOrthoBuilder(Nk,Nad,Ntot = Int((Nk^2+Nk)/2))
	Astruct = zeros(krylovOrthoStruct,Ntot)
	Aold = zeros(krylovOrthoStruct,Ntot)

	#build A0 structure
	icount = 0
	Aold[1] = krylovOrthoStructBuilder([(0,0,0,0,0,0)],[1])

	Astruct += Aold
	for nrec in 1:Nad
		Aold[1] = krylovSuccessor(Aold[1])
		Astruct += Aold
	end

	return Astruct
end



function krylovRmatBuilder(Nk,Nad,Ntot = Int((Nk^2+Nk)/2))
	Rstruct = zeros(krylovStruct,Ntot,Ntot)
	Rold = zeros(krylovStruct,Ntot,Ntot)

	#build R0 structure
	kcount = 0
	for knum in 0:Nk-1
		for jnum in 0:knum
			kcount += 1
			lcount = 0
			for lnum in 0:Nk-1
				for mnum in 0:lnum
					lcount += 1
					Rold[kcount,lcount] = krylovStructBuilder([(knum+lnum,jnum+mnum),(knum+mnum,jnum+lnum)],[1,1])
				end
			end
		end
	end

	Rstruct += Rold
	for nrec in 1:Nad
		for n1 in 1:Ntot
			for n2 in 1:Ntot
				Rold[n1,n2] = krylovSuccessor(Rold[n1,n2])
			end
		end
		Rstruct += Rold
	end

	return Rstruct
end

function krylovRmatOrthoBuilder(Nk,Nad,Ntot = Int((Nk^2+Nk)/2))
	Rstruct = zeros(krylovOrthoStruct,Ntot,Ntot)
	Rold = zeros(krylovOrthoStruct,Ntot,Ntot)

	#build R0 structure
	kcount = 0
	for knum in 0:Nk-1
		for jnum in 0:knum
			kcount += 1
			lcount = 0
			for lnum in 0:knum
				for mnum in 0:lnum
					lcount += 1
					Rold[kcount,lcount] = krylovOrthoStructBuilder([(knum,lnum,0,jnum,mnum,0),(knum,mnum,0,jnum,lnum,0)],[1,1])
					Rold[lcount,kcount] = Rold[kcount,lcount]
				end
			end
		end
	end

	Rstruct += Rold
	for nrec in 1:Nad
		for n1 in 1:Ntot
			for n2 in 1:n1
				Rold[n1,n2] = krylovSuccessor(Rold[n1,n2])
				Rold[n2,n1] = Rold[n1,n2]
			end
		end
		Rstruct += Rold
	end

	return Rstruct
end

function krylovGeneralizedMats(H,Nk,Mcoeffs,Nad)
	Ntot = Int((Nk^2+Nk)/2)
	
	Astruct = krylovAvecBuilder(Nk,Nad,Ntot)
	Avec = zeros(Ntot)
	for i in 1:Ntot
		Avec[i] = krylovValueFromStructure(Astruct[i],Mcoeffs)
	end
	
	Rstruct = krylovRmatBuilder(Nk,Nad,Ntot)
	Rmat = zeros(Ntot,Ntot)
	for n1 in 1:Ntot
		for n2 in 1:Ntot
			Rmat[n1,n2] = krylovValueFromStructure(Rstruct[n1,n2],Mcoeffs)
		end
	end
	return Avec,Rmat
end

function krylovGeneralizedOrthoMats(H,Nk,Mkjn,Nad)
	Kdict = kjDictBuilder(Nk)
	Ntot = Int((Nk^2+Nk)/2)
	
	println("Building A structure...")
	@time Astruct = krylovAvecOrthoBuilder(Nk,Nad,Ntot)
	println("Evaluating A structure...")
	Avec = zeros(Ntot)
	@time for i in 1:Ntot
		Avec[i] = krylovValueFromStructure(Astruct[i],Mkjn,Kdict)
	end
	
	println("Building R structure...")
	@time Rstruct = krylovRmatOrthoBuilder(Nk,Nad,Ntot)
	@show Rstruct[1,1]
	@show Rstruct[end,end]
	println("Evaluating R structure...")
	Rmat = zeros(Ntot,Ntot)
	@time for n1 in 1:Ntot
		for n2 in 1:Ntot
			Rmat[n1,n2] = krylovValueFromStructure(Rstruct[n1,n2],Mkjn,Kdict)
		end
	end

	return Avec,Rmat
end


function krylovTraceConstraint(rmat,Avec,n,Mcoeffs,Nk)
	Ntot = Int((Nk^2+Nk)/2) 
	Bns = zeros(Ntot,n+1)
	for nnum in 0:n
		icount = 0
		for ik in 1:Nk
			for ij in 1:ik
				icount += 1
				Bns[icount,nnum+1] = Mcoeffs[ik+ij-1+nnum]
			end
		end
	end

	βmat = zeros(n+1,n+1)
	for n1 in 1:n+1
		for n2 in 1:n+1
			βmat[n2,n1] = dot(rmat*Bns[:,n1] , Bns[:,n2])
		end
	end

	cvec = zeros(n+1)
	rA = rmat * Avec
	for n1 in 1:n+1
		cvec[n1] = dot(rA,Bns[:,n1])
	end

	mu = -2 * inv(βmat) * cvec

	return mu,Bns
end

function krylovOrthoTraceConstraint(rmat,Avec,n,Mkjn,Nk)
	Ntot = Int((Nk^2+Nk)/2) 
	
	βmat = zeros(n+1,n+1)
	for n1 in 1:n+1
		for n2 in 1:n+1
			βmat[n2,n1] = dot(rmat*Mkjn[:,n1] , Mkjn[:,n2])
		end
	end

	cvec = zeros(n+1)
	rA = rmat * Avec
	for n1 in 1:n+1
		cvec[n1] = dot(rA,Mkjn[:,n1])
	end

	mu = -inv(βmat) * cvec

	return mu
end

function krylovSolution(MU,H,Nk,nmu,Nad)
	Mcoeffs = krylovCoeffs(MU,H,Nk+nmu+2*Nad)
	A,R = krylovGeneralizedMats(H,Nk,Mcoeffs,Nad)
	r = inv(R)
	mu,Bns = krylovTraceConstraint(r,A,nmu,Mcoeffs,Nk)
	a = -r*A
	for imu in 0:nmu
		a -= 0.5*mu[imu+1]*(r*Bns[:,imu+1])
	end

	#= #constraint checking that a⋅Bn=0 for all Bn
	for imu in 0:nmu
		@show dot(a,Bns[:,imu+1]) / sqrt(dot(a,a)*dot(Bns[:,imu+1],Bns[:,imu+1]))
	end
	# =#
	return a
end

function krylovOrthoSolution(MU,H,Nk,nmu,Nad)
	@show Nn = 2*maximum([nmu,Nad])
	ϕ0 = MU[:,1]
	Mkjn,αvecs,βvecs,Qvecs = krylovSubspaceCoeffs(ϕ0,H,Nk,Nn)
	println("Building generalized A vector and R matrix")
	@time A,R = krylovGeneralizedOrthoMats(H,Nk,Mkjn,Nad)
	#println("R=")
	#display(R)
	A *= 2
	R *= 2
	#@show A
	println("Inverting R matrix")
	@time r = inv(R)
	#display(r)
	println("Finding μ Lagrange multipliers")
	@time mu = krylovOrthoTraceConstraint(r,A,nmu,Mkjn,Nk)
	a = -r*A
	for imu in 0:nmu
		a -= mu[imu+1]*(r*Mkjn[:,imu+1])
	end
	#@show a

	#= #constraint checking that a⋅M(n)=0 for all M(n)
	for imu in 0:nmu
		@show dot(a,Mkjn[:,imu+1]) / sqrt(dot(a,a)*dot(Mkjn[:,imu+1],Mkjn[:,imu+1])),dot(a,Mkjn[:,imu+1])
	end
	# =#
	return a,βvecs,Qvecs
end


function rhoKrylovBuilder(MU,H,Nk,nmu,Nad,ortho=true,adj=true)
	if adj == false
		if ortho == false
			println("Non-orthogonal routine")
			Hsts = length(H[:,1])
			mu0 = MU[:,1]
			mu0 /= sqrt(sum(abs2.(mu0)))
			Pμ = outerProduct(mu0,mu0)
			Hkmu = krylovVecs(mu0,H,2Nk)

			a = krylovSolution(MU,H,Nk,nmu,Nad)
			#@show a
			ρδ = zeros(Hsts,Hsts)
			icount = 0
			for ik in 1:Nk
				for ij in 1:ik
					icount += 1
					Pkj = outerProduct(Hkmu[:,ik],Hkmu[:,ij])
					#@show ik,ij,Pkj == (H^(ik-1)) * Pμ * (H^(ij-1))
					#@show sum(abs.(Pkj - (H^(ik-1)) * Pμ * (H^(ij-1))))
					#@show sum(abs.(Pkj))
					ρδ += a[icount] * (Pkj + Pkj')
				end
			end

			ρst = Pμ + ρδ
			@show tr(ρδ)
			#@show ρδ == ρδ'
			@show tr(ρst)
			
			#@show tr(Pμ^2)
			#@show tr(ρst^2)
			#@show tr(commutator(H,ρst)^2)

			return ρst
		else
			println("Orthogonal routine")
			Hsts = length(H[:,1])
			a,βvecs,Qvecs = krylovOrthoSolution(MU,H,Nk,nmu,Nad)
			Pμ = outerProduct(Qvecs[:,1],Qvecs[:,1])
			#@show a
			ρδ = zeros(Hsts,Hsts)
			icount = 0
			println("Building ρδ from Lanczos basis")
			@time for ik in 1:Nk
				for ij in 1:ik
					icount += 1
					Pkj = outerProduct(Qvecs[:,ik],Qvecs[:,ij])
					ρδ += a[icount] * (Pkj + Pkj')
				end
			end

			@show sum(abs2.(a))
			ρst = Pμ + ρδ
			@show sum(abs2.(Pμ - ρδ)),sum(abs2.(Pμ)),sum(abs2.(ρδ)),sum(abs2.(ρst))
			@show tr(ρδ),tr(ρst)
			#@show ρδ == ρδ'
			
			#@show tr(Pμ^2)
			#@show tr(ρst^2)
			#@show tr(commutator(H,ρst)^2)

			return ρst
		end
	else #adj == true
		println("adjoint routine")
		Hsts = length(H[:,1])
		mu0 = MU[:,1]
		mu0 /= sqrt(sum(abs2.(mu0)))
		Hk0m,Qvecs = Hk0mBuilder(mu0,H,Nk,4Nk + 2nmu + 1)
		@show maximum(Hk0m)
		γvec = adjSolution(Nk,Hk0m,nmu)
		Pμ = outerProduct(Qvecs[:,1],Qvecs[:,1])
		ρst = zeros(Hsts,Hsts)
		adn = copy(Pμ)
		println("Building ρst from commutators")
		@time for ik in 1:Nk+1
			ρst += γvec[ik] * adn
			adn /= (2ik*(2ik-1))
			adn = commutator(H,commutator(H,adn)) 
		end
		@show γvec
		@show tr(ρst^2),tr(Pμ^2)
		@show tr(ρst)

		return ρst
	end
end
	
function Hk0mBuilder(ϕ0,H,Nj,Mmax)
	Ntot = Nj + Mmax + 1
	αvec,βvec,Qmat,Lflag = lanczosBasis(ϕ0,H,Ntot)
	Hkmat = SymTridiagonal(αvec, βvec[1:Ntot])
	if Lflag == true
		Ntot = length(αvec)
	end
	HM = zeros(Ntot+1,Ntot+1,Mmax+1)
	HM[:,:,1] = id(Ntot+1)
	for m in 1:Mmax
		HM[:,:,m+1] = Hkmat * HM[:,:,m]
	end

	Hk0m = HM[1:Nj+1,1,:]

	#= Debugging routine to check correct building of Hk0m matrix
	display(Hk0m)
	Mdebug = zeros(Nj+1,Mmax+1)
	e0 = Qmat[:,1]
	for jnum in 1:Nj+1
		ej = Qmat[:,jnum]
		Mdebug[jnum,1] = dot(ej,e0)
		for mnum in 1:Mmax
			ej = H * ej
			Mdebug[jnum,mnum+1] = dot(ej,e0)
		end
	end
	display(Mdebug)
	error("subspace error")	
	# =#
	return Hk0m,Qmat
end

function adjR0Factorial(Nα,Hk0m)
	R0 = zeros(Nα,Nα)

	for αnum in 1:Nα
		for βnum in 1:Nα
			for λ1 in 0:2αnum
				for λ2 in 0:2βnum
					FACT = factorial(big(λ1)) * factorial(big(λ2)) * factorial(big(2αnum-λ1))*factorial(big(2βnum-λ2))
					ffact = Float64(1/FACT)
					R0[βnum,αnum] += (-1)^(λ1+λ2) * Hk0m[1,2βnum+λ1-λ2+1] * Hk0m[1,2αnum+λ2-λ1+1] * ffact
				end
			end
		end
	end

	return R0
end

function adjRtotFactorial(Nα,Hk0m,n)
	R = zeros(BigFloat,Nα,Nα)

	for αnum in 1:Nα
		for βnum in 1:Nα
			for λ1 in 0:2αnum
				for λ2 in 0:2βnum
					for nnum in 0:n
						for σ1 in 0:nnum
							for σ2 in 0:nnum	
								FACT = factorial(big(λ1)) * factorial(big(λ2)) * factorial(big(2αnum-λ1))*factorial(big(2βnum-λ2))
								#ffact = Float64(1/FACT) * binomial(n,σ1) * binomial(n,σ2)
								ffact = binomial(n,σ1) * binomial(n,σ2) / FACT
								R[βnum,αnum] += (-1)^(λ1+λ2+σ1+σ2) * Hk0m[1,2βnum+λ1-λ2+1+nnum+σ1-σ2] * Hk0m[1,2αnum+λ2-λ1+1+nnum+σ2-σ1] * ffact
							end
						end
					end
				end
			end
		end
	end

	return R
end

function adjA0Factorial(Nα,Hk0m)
	A0 = zeros(Nα)

	for αnum in 1:Nα
		for λ in 0:2αnum
			ffact = Float64(1/(factorial(big(λ)) * factorial(big(2αnum-λ))))
			A0[αnum] += (-1)^λ * Hk0m[1,2αnum-λ+1] * Hk0m[1,λ+1] * ffact
		end
	end

	return A0
end

function adjAtotFactorial(Nα,Hk0m,n)
	A = zeros(Nα)

	for αnum in 1:Nα
		for λ in 0:2αnum
			for nnum in 0:n
				for σ1 in 0:nnum
					for σ2 in 0:nnum
						ffact = Float64(1/(factorial(big(λ)) * factorial(big(2αnum-λ)))) * binomial(n,σ1) * binomial(n,σ2)
						A[αnum] += (-1)^(λ+σ1+σ2) * Hk0m[1,2αnum-λ+1+nnum+σ1-σ2] * Hk0m[1,λ+1+nnum+σ2-σ1] * ffact
					end
				end
			end
		end
	end

	return A
end



function adjSolution(Nα,Hk0m,Mmax)
	NN = 1
	Rfact = adjRtotFactorial(Nα,Hk0m,NN)
	@show matConditioning(Float64.(Rfact))
	@show det(Rfact)
	Afact = adjAtotFactorial(Nα,Hk0m,NN)
	@show sum(abs2.(Afact))

	γfact = zeros(Nα+1)
	γfact[1] = 1
	γfact[2:end] = -inv(Rfact)*Afact

	@show γfact
	return γfact

end