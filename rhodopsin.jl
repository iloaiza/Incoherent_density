# Rhodopsin potential module

function diaRhod(x,ϕ,E0=E0rhod,E1=E1rhod,V0=V0rhod,V1=V1rhod,ω2=ω2rhod,κ=κrhod,λ=λrhod,m_inv=Iinvrhod)
	H = zeros(2,2)
	H[1,1] = E0 + 0.5*V0*(1-cos(ϕ))+0.5*ω2*x^2
	H[1,2] = λ*x
	H[2,1] = H[1,2]
	H[2,2] = E1 - 0.5*V1*(1-cos(ϕ))+0.5*ω2*x^2+κ*x

	return H
end

function rhodopsinBuild(Nx,Nϕ,E0=E0rhod,E1=E1rhod,V0=V0rhod,V1=V1rhod,ω2=ω2rhod,κ=κrhod,λ=λrhod,m_inv=Iinvrhod,cutConst=3)
	### Build harmonic oscillator states (over x coupling mode)
	# Basis: scaled Hermite polynomials (eigenfunctions of Harmonic oscillator centered in X=0 with frequency ω)
	#############################################################
	ω = sqrt(ω2)
	Idx = Diagonal(ones(Nx))
	# Ground state energy
	Exgs = zeros(Nx,Nx)
	for i in 1:Nx
		Exgs[i,i] = E0 + ω*(i-0.5)
	end

	# Build harmonic oscillator ladder operators
	aHat = zeros(Nx,Nx)
	adagHat = zeros(Nx,Nx)
	for i in 1:Nx-1
		aHat[i,i+1] = sqrt(i)
		adagHat[i+1,i] = sqrt(i)
	end
	nHat = zeros(Nx,Nx)
	for i in 1:Nx
		nHat[i,i] = i-1
	end

	# Get position and momentum operators. X2 and Px2 are squared operators to avoid inexact a*adag operation in last eigenvalue
	X = 1/sqrt(2*ω)*(adagHat + aHat)
	Px = 1im*sqrt(ω/2)*(adagHat - aHat)
	X2 = 1/(2*ω)*(adagHat^2 + aHat^2 + 2*nHat + Idx)
	Px2 = -ω/2*(adagHat^2 + aHat^2 - 2*nHat - Idx)
	##Exgs = 0.5*Px2 + 0.5*ω2*X2 #alternative definition for Exgs

	# Excited state energy
	Exes = E1*Idx + 0.5*Px2 + 0.5*ω2*X2 + κ*X
	##############################################################
	
	### Build planar wave basis states (over ϕ in [0:2π) angle)
	# Basis: planar wave functions exp(ikϕ) for k = 1,...,Nϕ
	##############################################################
	#Nϕeff will be effective number of planar wave functions, only lowest Nϕ eigenstates of each ground and excited states will be taken into account
	Nϕeff = Int(Nϕ * cutConst)
	kinvec = zeros(Nϕeff)

	# Kinetic energy contribution
	for i in 1:Nϕeff
		kinvec[i] = i^2*m_inv/2
	end

	Eϕeffgs = SymTridiagonal(kinvec .+ V0/2, -V0/4 .* ones(Nϕeff-1))
	Eϕeffes = SymTridiagonal(kinvec .- V1/2,  V1/4 .* ones(Nϕeff-1))
	
	##############################################################
	### Rotate planar wave basis to just use first ncutoff states
	eigenGS,Ugs = adiabatize(Eϕeffgs)
	eigenES,Ues = adiabatize(Eϕeffes)
	# only first Nϕ eigenvectors of GS and ES will be considered: basis consisting of 2Nϕ vectors
	# calculate energy of GS vectors in ES basis and vice-versa:
	# U*[1 0 0] -> PW decomp of [1 0 0] GS
	Smat = Ues'*Ugs #overlap matrix <ϕes|ϕgs>
	Sdag = Smat'

	Eϕgs = zeros(2Nϕ,2Nϕ)
	Eϕgs[1:Nϕ,1:Nϕ] = Diagonal(eigenGS[1:Nϕ])
	Eϕgs[1:Nϕ,Nϕ+1:end] = Diagonal(eigenGS[1:Nϕ]) * Sdag[1:Nϕ,1:Nϕ]
	Eϕgs[Nϕ+1:end,1:Nϕ] = Eϕgs[1:Nϕ,Nϕ+1:end]'
	Eϕgs[Nϕ+1:end,Nϕ+1:end] = (Ues'*Eϕeffgs*Ues)[1:Nϕ,1:Nϕ]
	
	Eϕes = zeros(2Nϕ,2Nϕ)
	Eϕes[1:Nϕ,1:Nϕ] = (Ugs'*Eϕeffes*Ugs)[1:Nϕ,1:Nϕ]
	Eϕes[1:Nϕ,Nϕ+1:end] = Sdag[1:Nϕ,1:Nϕ] * Diagonal(eigenES[1:Nϕ])
	Eϕes[Nϕ+1:end,1:Nϕ] = Eϕes[1:Nϕ,Nϕ+1:end]'
	Eϕes[Nϕ+1:end,Nϕ+1:end] = Diagonal(eigenES[1:Nϕ])

	### Electronic coupling
	##############################################################
	Exc = λ*X
	##############################################################

	### Overlap matrix
	Sϕ = zeros(2Nϕ,2Nϕ)
	for i in 1:2Nϕ
		Sϕ[i,i] = 1
	end
	#Sϕ[1:Nϕ,Nϕ+1:end] = Sdag[1:Nϕ,1:Nϕ]
	#Sϕ[Nϕ+1:end,1:Nϕ] = Smat[1:Nϕ,1:Nϕ]
	### Build complete Hamiltonian with Kronecker products
	Idϕ = id(2Nϕ)
	m = Nx*2Nϕ
	H = zeros(2m,2m)
	H[1:m,1:m] = kron(Idϕ, Exgs) + kron(Eϕgs, Idx)
	H[m+1:2m,m+1:2m] = kron(Idϕ, Exes) + kron(Eϕes, Idx)
	H[1:m,m+1:2m] = kron(Sϕ, Exc)
	H[m+1:2m,1:m] = H[1:m,m+1:2m]'

	return Symmetric(round.(H,digits=10))
end

function rhodopsinRedBuild(Nx,Nϕ,E0=E0rhod,E1=E1rhod,V0=V0rhod,V1=V1rhod,ω2=ω2rhod,κ=κrhod,λ=λrhod,m_inv=Iinvrhod,cutConst=3)
	### Build harmonic oscillator states (over x coupling mode)
	# Basis: scaled Hermite polynomials (eigenfunctions of Harmonic oscillator centered in X=0 with frequency ω)
	#############################################################
	ω = sqrt(ω2)
	Idx = Diagonal(ones(Nx))
	# Ground state energy
	Exgs = zeros(Nx,Nx)
	for i in 1:Nx
		Exgs[i,i] = E0 + ω*(i-0.5)
	end

	# Build harmonic oscillator ladder operators
	aHat = zeros(Nx,Nx)
	adagHat = zeros(Nx,Nx)
	for i in 1:Nx-1
		aHat[i,i+1] = sqrt(i)
		adagHat[i+1,i] = sqrt(i)
	end
	nHat = zeros(Nx,Nx)
	for i in 1:Nx
		nHat[i,i] = i-1
	end

	# Get position and momentum operators. X2 and Px2 are squared operators to avoid inexact a*adag operation in last eigenvalue
	X = 1/sqrt(2*ω)*(adagHat + aHat)
	Px = 1im*sqrt(ω/2)*(adagHat - aHat)
	X2 = 1/(2*ω)*(adagHat^2 + aHat^2 + 2*nHat + Idx)
	Px2 = -ω/2*(adagHat^2 + aHat^2 - 2*nHat - Idx)
	##Exgs = 0.5*Px2 + 0.5*ω2*X2 #alternative definition for Exgs

	# Excited state energy
	Exes = E1*Idx + 0.5*Px2 + 0.5*ω2*X2 + κ*X
	##############################################################
	
	### Build planar wave basis states (over ϕ in [0:2π) angle)
	# Basis: planar wave functions exp(ikϕ) for k = 1,...,Nϕ
	##############################################################
	#Nϕeff will be effective number of planar wave functions, only lowest Nϕ eigenstates of each ground and excited states will be taken into account
	Nϕeff = Int(Nϕ * cutConst)
	kinvec = zeros(Nϕeff)

	# Kinetic energy contribution
	for i in 1:Nϕeff
		kinvec[i] = i^2*m_inv/2
	end

	Eϕeffgs = SymTridiagonal(kinvec .+ V0/2, -V0/4 .* ones(Nϕeff-1))
	Eϕeffes = SymTridiagonal(kinvec .- V1/2,  V1/4 .* ones(Nϕeff-1))
	
	##############################################################
	### Rotate planar wave basis to just use first ncutoff states
	eigenGS,Ugs = adiabatize(Eϕeffgs)
	eigenES,Ues = adiabatize(Eϕeffes)
	# only first Nϕ eigenvectors of GS and ES will be considered: basis consisting of 2Nϕ vectors
	# calculate energy of GS vectors in ES basis and vice-versa:
	# U*[1 0 0] -> PW decomp of [1 0 0] GS

	Eϕgs = zeros(Nϕ,Nϕ)
	Eϕgs = Diagonal(eigenGS[1:Nϕ])
	
	Eϕes = zeros(Nϕ,Nϕ)
	Eϕes= Diagonal(eigenES[1:Nϕ])

	### Electronic coupling
	##############################################################
	Exc = λ*X
	##############################################################

	Smat = Ues'*Ugs
	Sdag = Smat'
	Sϕ = Smat[1:Nϕ,1:Nϕ]
	#Sϕ[Nϕ+1:end,1:Nϕ] = Smat[1:Nϕ,1:Nϕ]

	### Build complete Hamiltonian with Kronecker products
	Idϕ = id(Nϕ)
	m = Nx*Nϕ
	H = zeros(2m,2m)
	H[1:m,1:m] = kron(Idϕ, Exgs) + kron(Eϕgs, Idx)
	H[m+1:2m,m+1:2m] = kron(Idϕ, Exes) + kron(Eϕes, Idx)
	H[1:m,m+1:2m] = kron(Sϕ, Exc)
	H[m+1:2m,1:m] = H[1:m,m+1:2m]'

	return Symmetric(round.(H,digits=10))
end

function rhodopsinSimpleBuild(Nx,Nϕ,E0=E0rhod,E1=E1rhod,V0=V0rhod,V1=V1rhod,ω2=ω2rhod,κ=κrhod,λ=λrhod,m_inv=Iinvrhod)
	### Build harmonic oscillator states (over x coupling mode)
	# Basis: scaled Hermite polynomials (eigenfunctions of Harmonic oscillator centered in X=0 with frequency ω)
	#############################################################
	ω = sqrt(ω2)
	Idx = Diagonal(ones(Nx))
	# Ground state energy
	Exgs = zeros(Nx,Nx)
	for i in 1:Nx
		Exgs[i,i] = E0 + ω*(i-0.5)
	end

	# Build harmonic oscillator ladder operators
	aHat = zeros(Nx,Nx)
	adagHat = zeros(Nx,Nx)
	for i in 1:Nx-1
		aHat[i,i+1] = sqrt(i)
		adagHat[i+1,i] = sqrt(i)
	end
	nHat = zeros(Nx,Nx)
	for i in 1:Nx
		nHat[i,i] = i-1
	end

	# Get position and momentum operators. X2 and Px2 are squared operators to avoid inexact a*adag operation in last eigenvalue
	X = 1/sqrt(2*ω)*(adagHat + aHat)
	Px = 1im*sqrt(ω/2)*(adagHat - aHat)
	X2 = 1/(2*ω)*(adagHat^2 + aHat^2 + 2*nHat + Idx)
	Px2 = -ω/2*(adagHat^2 + aHat^2 - 2*nHat - Idx)
	##Exgs = 0.5*Px2 + 0.5*ω2*X2 #alternative definition for Exgs

	# Excited state energy
	Exes = E1*Idx + 0.5*Px2 + 0.5*ω2*X2 + κ*X
	##############################################################
	
	### Build planar wave basis states (over ϕ in [0:2π) angle)
	# Basis: planar wave functions exp(ikϕ) for k = 1,...,Nϕ
	##############################################################
	Idϕ = Diagonal(ones(2Nϕ))
	kinvec = zeros(Nϕ)

	# Kinetic energy contribution
	for i in 1:Nϕ
		kinvec[i] = i^2*m_inv/2
	end

	Eϕeffgs = SymTridiagonal(kinvec .+ V0/2, -V0/4 .* ones(Nϕ-1))
	Eϕeffes = SymTridiagonal(kinvec .- V1/2,  V1/4 .* ones(Nϕ-1))
	
	##############################################################
	eigenGS,Ugs = adiabatize(Eϕeffgs)
	eigenES,Ues = adiabatize(Eϕeffes)
	Smat = Ues'*Ugs #overlap matrix <ϕes|ϕgs>
	Sdag = Smat'

	Eϕgs = zeros(2Nϕ,2Nϕ)
	Eϕgs[1:Nϕ,1:Nϕ] = Diagonal(eigenGS)
	Eϕgs[1:Nϕ,Nϕ+1:end] = Diagonal(eigenGS) * Sdag
	Eϕgs[Nϕ+1:end,1:Nϕ] = Eϕgs[1:Nϕ,Nϕ+1:end]'
	Eϕgs[Nϕ+1:end,Nϕ+1:end] = Ues'*Eϕeffgs*Ues
	
	Eϕes = zeros(2Nϕ,2Nϕ)
	Eϕes[1:Nϕ,1:Nϕ] = Ugs'*Eϕeffes*Ugs
	Eϕes[1:Nϕ,Nϕ+1:end] = Sdag * Diagonal(eigenES)
	Eϕes[Nϕ+1:end,1:Nϕ] = Eϕes[1:Nϕ,Nϕ+1:end]'
	Eϕes[Nϕ+1:end,Nϕ+1:end] = Diagonal(eigenES)

	### Electronic coupling
	##############################################################
	Exc = λ*X
	##############################################################

	### Build complete Hamiltonian with Kronecker products
	m = Nx*2Nϕ
	H = zeros(2m,2m)
	H[1:m,1:m] = kron(Idϕ, Exgs) + kron(Eϕgs, Idx)
	H[m+1:2m,m+1:2m] = kron(Idϕ, Exes) + kron(Eϕes, Idx)
	H[1:m,m+1:2m] = kron(Idϕ, Exc)
	H[m+1:2m,1:m] = H[1:m,m+1:2m]'
	
	return Symmetric(H)
end

function rhodopsinEigenPlot(Nx,Nϕ,E0=E0rhod,E1=E1rhod,V0=V0rhod,V1=V1rhod,ω2=ω2rhod,κ=κrhod,λ=λrhod,m_inv=Iinvrhod)
	println("	Building Hamiltonian matrix")
	t00 = time()
	######### STARTING HAMILTONIAN BUILD
	ω = sqrt(ω2)
	Idx = Diagonal(ones(Nx))
	Exgs = zeros(Nx,Nx)
	for i in 1:Nx
		Exgs[i,i] = E0 + ω*(i-0.5)
	end
	aHat = zeros(Nx,Nx)
	adagHat = zeros(Nx,Nx)
	for i in 1:Nx-1
		aHat[i,i+1] = sqrt(i)
		adagHat[i+1,i] = sqrt(i)
	end
	nHat = zeros(Nx,Nx)
	for i in 1:Nx
		nHat[i,i] = i-1
	end
	X = 1/sqrt(2*ω)*(adagHat + aHat)
	Px = 1im*sqrt(ω/2)*(adagHat - aHat)
	X2 = 1/(2*ω)*(adagHat^2 + aHat^2 + 2*nHat + Idx)
	Px2 = -ω/2*(adagHat^2 + aHat^2 - 2*nHat - Idx)
	Exes = E1*Idx + 0.5*Px2 + 0.5*ω2*X2 + κ*X
	Eϕgs = zeros(Nϕ,Nϕ)
	Eϕes = zeros(Nϕ,Nϕ)
	Idϕ = Diagonal(ones(Nϕ))
	for i in 1:Nϕ
		Eϕgs[i,i] = i^2*m_inv/2
		Eϕes[i,i] = Eϕgs[i,i]
	end
	for i in 1:Nϕ
		Eϕgs[i,i] += V0/2
		Eϕes[i,i] += -V1/2
	end
	for i in 1:Nϕ-1
		Eϕgs[i,i+1] += -V0/4
		Eϕgs[i+1,i] += -V0/4
		Eϕes[i,i+1] += V1/4
		Eϕes[i+1,i] += V1/4
	end
	Exc = λ*X
	m = Nx*Nϕ
	H = zeros(2m,2m)
	H[1:m,1:m] = kron(Idϕ, Exgs) + kron(Eϕgs, Idx)
	H[m+1:2m,m+1:2m] = kron(Idϕ, Exes) + kron(Eϕes, Idx)
	H[1:m,m+1:2m] = kron(Idϕ, Exc)
	H[m+1:2m,1:m] = H[1:m,m+1:2m]'
	println("Time for Hamiltonian building is $(time()-t00) seconds")


	Xgrid = collect(-2:0.01:2)
	Φgrid = collect(-π/2:0.01:3π/2)
	XL = length(Xgrid)
	ΦL = length(Φgrid)
	PES = zeros(XL,ΦL,2,2)

	println("	Building PES")
	@time for (ix,x) in enumerate(Xgrid)
		for (iϕ,ϕ) in enumerate(Φgrid)
			PES[ix,iϕ,:,:] = diaRhod(x,ϕ,E0,E1,V0,V1,ω2,κ,λ,m_inv)
		end
	end

	println("	Diagonalizing ϕ basis and obtaining eigenfunctions")
	t00 = time()
	Eigenϕgs,Uϕgs = adiabatize(Eϕgs)
	Eigenϕes,Uϕes = adiabatize(Eϕes)
	t01 = time()
	println("Time for diagonalization is $(t01-t00) seconds")
	Ψϕ = zeros(Complex,Nϕ,ΦL)
	for k in 1:Nϕ
		Ψϕ[k,:] .= exp.(1im*k*Φgrid)/sqrt(2*π)
	end
	Ψϕgs = zeros(Complex,Nϕ,ΦL)
	Ψϕes = zeros(Complex,Nϕ,ΦL)
	for n1 in 1:Nϕ
		for n2 in 1:Nϕ
			Ψϕgs[n1,:] .+= Uϕgs[n2,n1] * Ψϕ[n2,:]
			Ψϕes[n1,:] .+= Uϕes[n2,n1] * Ψϕ[n2,:]
		end
	end
	println("Time spent in eigenfunction building and evaluation is $(time()-t01) seconds")
	ampGS = abs2.(Ψϕgs)
	ampES = abs2.(Ψϕes)

	return H,Xgrid,Φgrid,PES,Eigenϕgs,Eigenϕes,ampGS,ampES
end

function ϕPot(cosine=true,Nϕ=200,kmax=249,res=600,V0=V0rhod,V1=V1rhod,E0=E0rhod,E1=E1rhod,minv=Iinvrhod)
	# Build Hamiltonian matrices for CIS and TRANS ϕ states
	KIN = zeros(kmax+1,kmax+1) #H[k1+1,k2+1] corresponds to <k1|H|k2>
	POT = zeros(kmax+1,kmax+1)
	if cosine == true
		Karr = collect(0:kmax)
	else
		Karr = collect(1:kmax+1)
	end
	Φ = range(-π/2,stop=3π/2,length=res)

	for (ik1,k1) in enumerate(Karr)
		for (ik2,k2) in enumerate(Karr)
			ksum = k1+k2
			kdif = k1-k2

			prefactor = 0.5
			if k1 == 0; prefactor /= sqrt(2); end
			if k2 == 0; prefactor /= sqrt(2); end

			if k1 == k2
				KIN[ik1,ik2] = 0.5*minv*k1^2
				POT[ik1,ik2] = 0.5 
			end

			if abs(kdif) == 1
				POT[ik1,ik2] = -prefactor*0.5
			end

			if ksum == 1
				if cosine == true
					POT[ik1,ik2] += -prefactor*0.5
				else
					POT[ik1,ik2] += prefactor*0.5
				end
			end
		end
	end

	Hcis = KIN + V0*POT + E0*id(kmax+1)
	Htrans = KIN - V1*POT + E1*id(kmax+1)

	# Obtain energies and unitary transformations
	Ecis,Ucis = adiabatize(Hcis)
	Etrans,Utrans = adiabatize(Htrans)

	# Build wavefunctions
	Ψ = zeros(kmax+1,res)
	for (ik,k) in enumerate(Karr)
		if cosine == true
			Ψ[ik,:] .= cos.(k*Φ)
			if k == 0
				Ψ[ik,:] ./= sqrt(2)
			end
		else
			Ψ[ik,:] .= sin.(k*Φ)
		end
		Ψ[ik,:] ./= sqrt(π)
	end

	Ψcis = zeros(Nϕ,res)
	Ψtrans = zeros(Nϕ,res)
	for n1 in 1:Nϕ
		for n2 in 1:Nϕ
			Ψcis[n1,:] .+= Ucis[n2,n1] * Ψ[n2,:]
			Ψtrans[n1,:] .+= Utrans[n2,n1] * Ψ[n2,:]
		end
	end

	#get overlaps product matrices for calculations
	Scis = zeros(Nϕ,Nϕ,res)
	Strans = zeros(Nϕ,Nϕ,res)
	S = zeros(Nϕ,Nϕ)
	dϕ = (Φ[2]-Φ[1])

	for n1 in 1:Nϕ
		for n2 in 1:Nϕ
			Scis[n1,n2,:] = dϕ * Ψcis[n1,:] .* Ψcis[n2,:]
			Strans[n1,n2,:] = dϕ * Ψtrans[n1,:] .* Ψtrans[n2,:]
			S[n1,n2] = sum(dϕ * Ψcis[n1,:] .* Ψtrans[n2,:])
		end
	end

	return Ecis[1:Nϕ],Etrans[1:Nϕ],S,Hcis,Htrans,Scis,Strans,Ψcis,Ψtrans,Φ
end

function HOpot(ω=ω2rhod,NHO=20,Nquad=200,κ=κrhod,Nextra=10)
	# use larger N=NH0+Nextra for building, stay with desired NHO states in the end (to avoid inaccuracies from boundary values)
	N = NHO + Nextra
	# Build Hamiltonian for harmonic oscillator
	Hgs = zeros(N,N)
	for i in 1:N #harmonic oscillators energies
		Hgs[i,i] = ω*(i-0.5)
	end

	# =
	#find HO operators (generic forms, mode-independent)
	aHat = zeros(N,N)
	adagHat = zeros(N,N)
	for i in 1:N-1
		aHat[i,i+1] = sqrt(i)
		adagHat[i+1,i] = sqrt(i)
	end
	Xres = 1/sqrt(2*ω)*(adagHat + aHat)
	X = Xres .* sqrt(ω) #rescaled operator
	#= 
	#Px = 1im*sqrt(ω/2)*(adagHat - aHat)
	#Xei,Uxei = adiabatize(X)
	#Yei,Uyei = adiabatize(Y)
	=#
	nHat = zeros(N,N)
	for i in 1:N
		nHat[i,i] = i-1
	end

	nodes,weights = ghQuad(Nquad)
	Xgridres = nodes ./sqrt(ω)

	Hes = Hgs + κ*X
	
	herPols = hermitePol(N)
	#build N wavefunctions corresponding to each hermite polinomial (0 to N-1)
	Ψ = zeros(Nquad,N)
	for bst in 1:N
		for (i,x) in enumerate(nodes)
			Ψ[i,bst] = polyval(herPols[bst,1:bst],x)
		end
	end

	Ψgs = zeros(size(Ψ))
	for nbas in 1:N
		Ψgs[:,nbas] = Ψ[:,nbas] .* exp.(-nodes.^2/2)
	end

	Egs = [Hgs[st,st] for st in 1:N]
	Ees,Ues = adiabatize(Hes)
	Ψes = zeros(size(Ψ))
	Ψquad = zeros(size(Ψ))
	for n1 in 1:N
		for n2 in 1:N
			Ψes[:,n1] .+= Ues[n2,n1] * Ψgs[:,n2]
			Ψquad[:,n1] .+= Ues[n2,n1] * Ψ[:,n2]
		end
	end

	S = zeros(N,N,Nquad)
	for st1 in 1:N
		for st2 in 1:N
			G = Ψ[:,st1] .* Ψquad[:,st2]
			S[st1,st2,:] = G .* weights
		end
	end

	return Egs[1:NHO],Hes[1:NHO,1:NHO],Ues[1:NHO,1:NHO],Ees[1:NHO],Ψ[:,1:NHO],Ψgs[:,1:NHO],Ψes[:,1:NHO],weights,S[1:NHO,1:NHO,:],nodes
end

function rhodCosHamiltonian(NHO=20,Nϕ=200,kmax=249,Nquad=200,ϕres=600,V0=V0rhod,V1=V1rhod,E0=E0rhod,E1=E1rhod,minv=Iinvrhod,ω=ω2rhod,μ=1,κ=κrhod,λ=λrhod)
	EcisCOS,EtransCOS,SϕCOS,HcisCOS,HtransCOS,ScisCOS,StransCOS,ΨcisCOS,ΨtransCOS,ΦCOS = ϕPot(true,Nϕ,kmax,ϕres,V0,V1,E0,E1,minv)
	EgsHO,HesHO,UesHO,EesHO,ΨHO,ΨgsHO,ΨesHO,weightsHO,SHO,XnHO = HOpot(ω,NHO,Nquad,κ)

	# Eigenenergies printing
	#=
	println("Cis state eigenenergies (in cm-1)")
	@show EcisCOS .* au2cm
	println("Trans state eigenenergies (in cm-1)")
	@show EtransCOS .* au2cm
	println("HO ground state eigenenergies (in cm-1)")
	@show EgsHO .* au2cm
	println("HO excited state eigenenergies (in cm-1)")
	@show EesHO .* au2cm
	# =#


	m = NHO * Nϕ

	HOovlp = zeros(NHO,NHO)
	Ec = zeros(NHO,NHO)
	for n1 in 1:NHO
		for n2 in 1:NHO
			Ec[n1,n2] = λ*sum(SHO[n1,n2,:] .* XnHO)
			HOovlp[n1,n2] = sum(SHO[n1,n2,:])
		end
	end
	Stot = kron(HOovlp,SϕCOS)

	# Overlap printing
	#=
	println("Cis-trans states overlaps")
	display(SϕCOS)
	println("Harmonic oscillator overlap")
	HOovlp = zeros(NHO,NHO)
	for n1 in 1:NHO
		for n2 in 1:NHO
			HOovlp[n1,n2] = sum(SHO[n1,n2,:])
		end
	end
	display(HOovlp)
	println("Coupling matrix")
	display(Ec)
	# =#

	Hcos = zeros(2m,2m)
	Hcos[1:m,1:m] = kron(Diagonal(EgsHO),id(Nϕ)) + kron(id(NHO),Diagonal(EcisCOS))
	Hcos[m+1:2m,m+1:2m] = kron(Diagonal(EesHO),id(Nϕ)) + kron(id(NHO),Diagonal(EtransCOS))
	Hcos[1:m,m+1:2m] = kron(Ec,SϕCOS)
	Hcos[m+1:2m,1:m] = Hcos[1:m,m+1:2m]'
	Hsym = Symmetric(Hcos)

	Hwc = zeros(2m,2m)
	Hwc[1:m,1:m] = kron(Diagonal(EgsHO),id(Nϕ)) + kron(id(NHO),Diagonal(EcisCOS))
	Hwc[m+1:2m,m+1:2m] = kron(Diagonal(EesHO),id(Nϕ)) + kron(id(NHO),Diagonal(EtransCOS))
	HwcSYM = Symmetric(Hwc)
	#Ecos,Ucos = adiabatize(Hsym)
	#MUcosdia,MUcosadi = dipole4Plot(Ecos,Ucos,μ,m,Stot)

	# Full matrix eigenvalue printing
	#=
	println("Full system first 900 eigenenergies (in cm-1)")
	@show Ecos[1:900] .* au2cm
	=#

	#Wcm = (Ecos .- Ecos[1]) .* au2cm
	return Hsym,HwcSYM,Stot
end

function rhodHamiltonian(NHO=20,Nϕ=200,kmax=249,ϕres=600,V0=V0rhod,V1=V1rhod,E0=E0rhod,E1=E1rhod,minv=Iinvrhod,ω=ω2rhod,Nquad=40,xrange=[-11,11],μ=1,κ=κrhod,λ=λrhod)
	EcisCOS,EtransCOS,SϕCOS,HcisCOS,HtransCOS,ScisCOS,StransCOS,ΨcisCOS,ΨtransCOS,ΦCOS = ϕPot(true,Nϕ,kmax,ϕres,V0,V1,E0,E1,minv)
	EcisSIN,EtransSIN,SϕSIN,HcisSIN,HtransSIN,ScisSIN,StransSIN,ΨcisSIN,ΨtransSIN,ΦSIN = ϕPot(false,Nϕ,kmax,ϕres,V0,V1,E0,E1,minv)
	HOgs,HOes,ΨHO,ΨexpHO,weightsHO,SHO,XnHO = HOpot(sqrt(ω),NHO,Nquad,κ,E1,xrange)

	m = NHO * Nϕ

	Eyc = zeros(NHO,NHO)
	for n1 in 1:NHO
		for n2 in 1:NHO
			Eyc[n1,n2] = λ*sum(SHO[n1,n2,:] .* XnHO)
		end
	end

	Sϕ = zeros(2Nϕ,2Nϕ)
	Sϕ[1:Nϕ,1:Nϕ] = id(Nϕ)
	Sϕ[Nϕ+1:end,Nϕ+1:end] = id(Nϕ)

	H = zeros(4m,4m)
	H[1:2m,2m+1:4m] = kron(Eyc,id(2Nϕ))
	H[2m+1:4m,1:2m] = H[1:2m,2m+1:4m]'
	H[1:m,1:m] = kron(HOgs,Diagonal(EcisCOS))
	H[m+1:2m,m+1:2m] = kron(HOgs,Diagonal(EcisSIN))
	H[2m+1:3m,2m+1:3m] = kron(HOes,Diagonal(EtransCOS))
	H[3m+1:4m,3m+1:4m] = kron(HOes,Diagonal(EtransSIN))

	Ecos,Ucos = adiabatize(Hcos)
	MUcosdia,MUcosadi = dipole4Plot(Ecos,Ucos,μ)

	#E,U = adiabatize(H)
	#MUdia,MUadi = dipole4Plot(E,U,μ)

	return Ecos,MUcosdia,MUcosadi#,E,MUdia,MUadi
end