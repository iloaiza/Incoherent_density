# Linear Vibronic Coupling (LVC) potential module
function diaLVC(x,y,ω,a,Δ,c)
	V = zeros(2,2)
	Ω = ω.^2
	V[1,1] = 0.5*(Ω[1]*(x^2)+Ω[2]*(y^2))
	V[2,2] = (Ω[1]*(x-a)^2+Ω[2]*y^2)/2-Δ
	V[1,2] = c*y
	V[2,1] = V[1,2]

	return V
end
	
ωTest=[7e-3,5e-3]; aTest=25.0; ΔTest=0.021; cTest=5e-4;

LVCTest(x,y)=diaLVC(x,y,ωTest,aTest,ΔTest,cTest)

#return energies for harmonic oscillator basis for each minimum
function LVCenergies(N,Nquad,ω,a,Δ,c)
	#build ground state energies
	Exgs = zeros(N,N)
	Eygs = zeros(N,N)
	for i in 1:N #harmonic oscillators energies
		Exgs[i,i] = ω[1]*(i-0.5)
		Eygs[i,i] = ω[2]*(i-0.5)
	end

	nodes,weights = ghQuad(Nquad) #get nodes and weights for GH quadrature
	#build diabatic super matrix
	ngrid = length(nodes)
	Vx = zeros(ngrid)
	Vy = zeros(ngrid)
	Cy = zeros(ngrid)
	X = nodes ./sqrt(ω[1])
	Y = nodes ./sqrt(ω[2])
	for i in 1:ngrid
		Vx[i] = 0.5*(ω[1]*(X[i]-a))^2-Δ/2
		Vy[i] = 0.5*(ω[2]*Y[i])^2-Δ/2
		Cy[i] = c*Y[i]
	end
	
	#build hermite polinomials basis centered in 0
	herPols = hermitePol(N)
	#build N wavefunctions corresponding to each hermite polinomial (0 to N-1)
	Ψ = zeros(ngrid,N)
	for bst in 1:N
		for (i,x) in enumerate(nodes)
			Ψ[i,bst] = polyval(herPols[bst,1:bst],x)
		end
	end

	#build excited states and couplings
	Exes = zeros(N,N)
	Eyes = zeros(N,N)
	Eyc = zeros(N,N)
	for st1 in 1:N
		for st2 in 1:st1
			#first add potential contributions
			G = Ψ[:,st1] .* Ψ[:,st2]
			Exes[st1,st2] += sum(G .* Vx .* weights)
			Eyes[st1,st2] += sum(G .* Vy .* weights)
			Eyc[st1,st2] += sum(G .* Cy .* weights)
			#add kinetic energy contributions
			if st1 == st2
				Exes[st1,st2] += 0.5*ω[1]*(st1 - 0.5)
				Eyes[st1,st2] += 0.5*ω[2]*(st1 - 0.5)
			elseif st1-st2 == 2 #st1 is always bigger than st2
				Exes[st1,st2] += -0.25*ω[1]*sqrt(st2*(st2+1))
				Eyes[st1,st2] += -0.25*ω[2]*sqrt(st2*(st2+1))
			end
			#hermitize
			Exes[st2,st1] = Exes[st1,st2]
			Eyes[st2,st1] = Eyes[st1,st2]
			Eyc[st2,st1] = Eyc[st1,st2]
		end
	end

	#build full hamiltonian matrix
	m = N^2
	H = zeros(2m,2m)
	H[1:m,1:m] = kron(Diagonal(ones(N)), Exgs) + kron(Eygs, Diagonal(ones(N)))
	H[m+1:2m,m+1:2m] = kron(Diagonal(ones(N)), Exes) + kron(Eyes, Diagonal(ones(N)))
	H[1:m,m+1:2m] = kron(Eyc, Diagonal(ones(N)))
	H[m+1:2m,1:m] = H[1:m,m+1:2m]'


	for nbas in 1:N
		Ψ[:,nbas] .*= exp.(-nodes.^2/2)
	end

	return Symmetric(H),Ψ,nodes,weights,X,Y
end

#gets an approximate minimum number for the GH quadrature so that all expanded H-polynomials fit inside the nodes (for exact quadratures!)
function nmin(nhermite,a)
	hmax = sqrt(2nhermite+1)
	return round(Int,0.5*(hmax^2+a*hmax+a^2/4)) + 10
end

function LVCenergiesSeparate(Nx,Ny,ω,a,Δ,c,Nquad=nmin(Nx,a))
	#builds wells at -a/2 and a/2 for x (different from 0 and a!)
	@show Nquad

	nodes,weights = ghQuad(2Nquad) #get nodes and weights for GH quadrature
	#build diabatic super matrix
	ngrid = length(nodes)
	Vgs = zeros(ngrid)
	Ves = zeros(ngrid)
	Vy = zeros(ngrid)
	Cy = zeros(ngrid)
	X = nodes ./sqrt(ω[1])
	Y = nodes ./sqrt(ω[2])
	nodesxplus = nodes .+ sqrt(ω[1])*a/2
	nodesxminus = nodes .- sqrt(ω[1])*a/2
	for i in 1:ngrid
		Vgs[i] = 0.5*(ω[1]*(X[i]+a/2))^2
		Ves[i] = 0.5*(ω[1]*(X[i]-a/2))^2-Δ/2
		Vy[i] = 0.5*(ω[2]*Y[i])^2-Δ/2
		Cy[i] = c*Y[i]
	end
	
	#build hermite polinomials basis centered in 0
	herPols = hermitePol(maximum([Nx,Ny]))
	#build N wavefunctions corresponding to each hermite polinomial (0 to N-1)
	ΨxGS = zeros(ngrid,Nx)
	ΨxES = zeros(ngrid,Nx)
	Ψy = zeros(ngrid,Ny)
	for bst in 1:Nx
		for i in 1:ngrid
			ΨxGS[i,bst] = polyval(herPols[bst,1:bst],nodesxplus[i])
			ΨxES[i,bst] = polyval(herPols[bst,1:bst],nodesxminus[i])
		end
	end
	for nbas in 1:Nx
		ΨxGS[:,nbas] .*= exp.(-nodesxplus.^2/2) .* exp.(nodes.^2/2)
		ΨxES[:,nbas] .*= exp.(-nodesxminus.^2/2) .* exp.(nodes.^2/2)
	end
	for bst in 1:Ny
		for (i,y) in enumerate(nodes)
			Ψy[i,bst] = polyval(herPols[bst,1:bst],y)
		end
	end
	#build energies over x basis
	Exgs = zeros(Nx,Nx)
	Exes = zeros(Nx,Nx)
	S = zeros(Nx,Nx) #overlap matrix
	for st1 in 1:Nx
		for st2 in st1:-1:1
			# add contributions over corresponding centered basis
			if st1 == st2
				Exgs[st1,st1] = ω[1]*(st1-0.5)
				Exes[st1,st1] = ω[1]*(st1-0.5)
			end
		end

		for st2 in 1:Nx
			# get overlaps
			Gcross = ΨxES[:,st1] .* ΨxGS[:,st2] .* weights
			S[st1,st2] = sum(Gcross)
		end
	end	
	
	#build energies over y basis
	Eygs = zeros(Ny,Ny)
	for i in 1:Ny
		Eygs[i,i] = ω[2]*(i-0.5)
	end
	Eyes = zeros(Ny,Ny)
	Eyc = zeros(Ny,Ny)
	for st1 in 1:Ny
		for st2 in 1:st1
			#first add potential contributions
			G = Ψy[:,st1] .* Ψy[:,st2]
			Eyes[st1,st2] += sum(G .* Vy .* weights)
			Eyc[st1,st2] += sum(G .* Cy .* weights)
			#add kinetic energy contributions
			if st1 == st2
				Eyes[st1,st2] += 0.5*ω[2]*(st1 - 0.5)
			elseif st1-st2 == 2 #st1 is always bigger than st2
				Eyes[st1,st2] += -0.25*ω[2]*sqrt(st2*(st2+1))
			end
			#hermitize
			Eyes[st2,st1] = Eyes[st1,st2]
			Eyc[st2,st1] = Eyc[st1,st2]
		end
	end

	#build full hamiltonian matrix
	m = Ny*Nx
	H = zeros(2m,2m)
	H[1:m,1:m] = kron(id(Ny), Exgs) + kron(Eygs, id(Nx))
	H[m+1:2m,m+1:2m] = kron(id(Ny), Exes) + kron(Eyes, id(Nx))
	H[1:m,m+1:2m] = kron(Eyc, S)
	H[m+1:2m,1:m] = H[1:m,m+1:2m]'


	for nbas in 1:Ny
		Ψy[:,nbas] .*= exp.(-nodes.^2/2)
	end
	for nbas in 1:Nx
		ΨxES[:,nbas] .*= exp.(-nodes.^2/2)
		ΨxGS[:,nbas] .*= exp.(-nodes.^2/2)
	end

	return Symmetric(H),Ψy,nodes,weights,X,Y,Exgs,Exes,ΨxGS,ΨxES,nodesxplus,nodesxminus,S
end

function LVCenergiesDual(Nx,Ny,ω,a,Δ,c,Nquad=nmin(Nx,a),stol = 1e-2)
	#builds wells at -a/2 and a/2 for x (different from 0 and a!)
	@show Nquad

	nodes,weights = ghQuad(2Nquad) #get nodes and weights for GH quadrature
	#build diabatic super matrix
	ngrid = length(nodes)
	Vgs = zeros(ngrid)
	Ves = zeros(ngrid)
	Vy = zeros(ngrid)
	Cy = zeros(ngrid)
	X = nodes ./sqrt(ω[1])
	Y = nodes ./sqrt(ω[2])
	nodesxplus = nodes .+ sqrt(ω[1])*a/2
	nodesxminus = nodes .- sqrt(ω[1])*a/2
	for i in 1:ngrid
		Vgs[i] = 0.5*(ω[1]*(X[i]+a/2))^2
		Ves[i] = 0.5*(ω[1]*(X[i]-a/2))^2-Δ/2
		Vy[i] = 0.5*(ω[2]*Y[i])^2-Δ/2
		Cy[i] = c*Y[i]
	end
	
	#build hermite polinomials basis centered in 0
	herPols = hermitePol(maximum([Nx,Ny])+2)
	#build N wavefunctions corresponding to each hermite polinomial (0 to N-1)
	ΨxGS = zeros(ngrid,Nx+2)
	ΨxES = zeros(ngrid,Nx+2)
	Ψy = zeros(ngrid,Ny)
	for bst in 1:Nx+2
		for i in 1:ngrid
			ΨxGS[i,bst] = polyval(herPols[bst,1:bst],nodesxplus[i])
			ΨxES[i,bst] = polyval(herPols[bst,1:bst],nodesxminus[i])
		end
	end
	for nbas in 1:Nx+2
		ΨxGS[:,nbas] .*= exp.(-nodesxplus.^2/2) .* exp.(nodes.^2/2)
		ΨxES[:,nbas] .*= exp.(-nodesxminus.^2/2) .* exp.(nodes.^2/2)
	end
	for bst in 1:Ny
		for (i,y) in enumerate(nodes)
			Ψy[i,bst] = polyval(herPols[bst,1:bst],y)
		end
	end
	#build energies over x basis
	Exgs = zeros(2Nx,2Nx)
	Exes = zeros(2Nx,2Nx)
	S = zeros(2Nx,2Nx) #overlap matrix
	for st1 in 1:Nx
		for st2 in st1:-1:1
			# add contributions over corresponding centered basis
			if st1 == st2
				Exgs[st1,st1] = ω[1]*(st1-0.5)
				Exes[st1+Nx,st1+Nx] = ω[1]*(st1-0.5)
			end
				# add kinetic energy contributions
			# contributions for diagonal part of different basis
			if st1 == st2
				Exes[st1,st2] += 0.5*ω[1]*(st1 - 0.5)
				Exgs[st1+Nx,st2+Nx] += 0.5*ω[1]*(st1 - 0.5)
			elseif st1-st2 == 2 #st1 is always bigger than st2
				Exes[st1,st2] += -0.25*ω[1]*sqrt(st2*(st2+1))
				Exgs[st1+Nx,st2+Nx] += -0.25*ω[1]*sqrt(st2*(st2+1))
			end
		end

		for st2 in 1:Nx
			# get overlaps
			Gcross = ΨxES[:,st1] .* ΨxGS[:,st2] .* weights
			Ggs = ΨxGS[:,st1] .* ΨxGS[:,st2] .* weights #overlap of GS-centered wavefunctions
			Ges = ΨxES[:,st1] .* ΨxES[:,st2] .* weights #overlap of ES-centered wavefunctions
			S[st1,st2] = sum(Ggs)
			S[st1+Nx,st2+Nx] = sum(Ges)
			S[st1+Nx,st2] = sum(Gcross)
			S[st2,st1+Nx] = S[st1+Nx,st2]
			# calculate potential 
			Exes[st1,st2] += sum(Ggs .* Ves)
			Exgs[st1+Nx,st2+Nx] += sum(Ges .* Vgs)
			Exes[st1+Nx,st2] += sum(Gcross .* Ves)
			Exgs[st1+Nx,st2] += sum(Gcross .* Vgs)
			# kinetic contribution
			GPM = sum(Gcross)
			if st2 > 2
				GPMminus = sum(ΨxES[:,st1] .* ΨxGS[:,st2-2] .* weights)
			end
			GPMplus = sum(ΨxES[:,st1] .* ΨxGS[:,st2+2] .* weights)
			if st2 <= 2
				Exes[st1+Nx,st2] += GPM*0.5*ω[1]*(st2-0.5) - 0.25*ω[1]*sqrt(st2*(st2+1))*GPMplus
			else	
				Exes[st1+Nx,st2] += GPM*0.5*ω[1]*(st2-0.5) - 0.25*ω[1]*(sqrt(st2*(st2+1))*GPMplus + sqrt((st2-2)*(st2-1))*GPMminus)
			end
			Exgs[st1+Nx,st2] = Exes[st1+Nx,st2]
		end
	end
	#hermitize
	for st1 in 1:2Nx
		for st2 in st1:2Nx
			Exgs[st1,st2] = Exgs[st2,st1]
			Exes[st1,st2] = Exes[st2,st1]
		end
	end
	
	sbreak = 2Nx+1
	Svals,SU = adiabatize(S,false)
	for (scounter,s) in enumerate(Svals)
		if s < stol
			sbreak = scounter
			break
		end
	end

	##S = SU*Diagonal(Svals)*SU'

	#= # symmetricla orthogonalization
	S12 = SU*Diagonal(sqrt.(Svals))*SU'
	S12inv = SU*Diagonal(Svals .^(-0.5))*SU'
	Egstilde = S12inv*Exgs*S12inv
	egs,ugs = adiabatize(Egstilde)
	Egs = real.([(Ugs'*Exgs*Ugs)[st,st] for st in 1:2Nx])
	@show Egs
	# =#

	@show Svals
	Svals = Svals[1:sbreak-1]
	@show Svals
	S12inv = zeros(2Nx,sbreak-1)
	for (scount,s) in enumerate(Svals)
		S12inv[scount,scount] = 1/sqrt(s)
	end
	Xtilde = SU*S12inv
	#display(Xtilde)
	Egstilde = Xtilde'*Exgs*Xtilde
	egs,ugs = adiabatize(Egstilde)
	@show egs
	

	
	#build energies over y basis
	Eygs = zeros(Ny,Ny)
	for i in 1:Ny
		Eygs[i,i] = ω[2]*(i-0.5)
	end
	Eyes = zeros(Ny,Ny)
	Eyc = zeros(Ny,Ny)
	for st1 in 1:Ny
		for st2 in 1:st1
			#first add potential contributions
			G = Ψy[:,st1] .* Ψy[:,st2]
			Eyes[st1,st2] += sum(G .* Vy .* weights)
			Eyc[st1,st2] += sum(G .* Cy .* weights)
			#add kinetic energy contributions
			if st1 == st2
				Eyes[st1,st2] += 0.5*ω[2]*(st1 - 0.5)
			elseif st1-st2 == 2 #st1 is always bigger than st2
				Eyes[st1,st2] += -0.25*ω[2]*sqrt(st2*(st2+1))
			end
			#hermitize
			Eyes[st2,st1] = Eyes[st1,st2]
			Eyc[st2,st1] = Eyc[st1,st2]
		end
	end

	#build full hamiltonian matrix
	m = Ny*2Nx
	H = zeros(2m,2m)
	H[1:m,1:m] = kron(id(Ny), Exgs) + kron(Eygs, id(2Nx))
	H[m+1:2m,m+1:2m] = kron(id(Ny), Exes) + kron(Eyes, id(2Nx))
	H[1:m,m+1:2m] = kron(Eyc, id(2Nx))
	H[m+1:2m,1:m] = H[1:m,m+1:2m]'


	for nbas in 1:Ny
		Ψy[:,nbas] .*= exp.(-nodes.^2/2)
	end
	for nbas in 1:Nx+2
		ΨxES[:,nbas] .*= exp.(-nodes.^2/2)
		ΨxGS[:,nbas] .*= exp.(-nodes.^2/2)
	end

	return Symmetric(H),Ψy,nodes,weights,X,Y,Exgs,Exes,ΨxGS,ΨxES,nodesxplus,nodesxminus,S
end

function LVCbasisDebug(N,Nquad,ω,a,Δ,c)
	H,Ψ,nodes,weights,X,Y = LVCenergies(N,Nquad,ω,a,Δ,c)
	#remove exponential part (included in weights for GH quadrature)
	for i in 1:N
		Ψ[:,i] .*= exp.(nodes.^2/2)
	end
	for n1 in 1:N
		for n2 in 1:N
			if abs2(sum(Ψ[:,n1] .* Ψ[:,n2] .* weights)) > 1e-10
				@show n1,n2,sum(Ψ[:,n1] .* Ψ[:,n2] .* weights)
			end
		end
	end
end

function dia2adiFullBuild(X,Y,ω,a,Δ,c)
	xL = length(X)
	yL = length(Y)

	Earr = zeros(xL*yL,2)
	Uarr = zeros(xL*yL,2,2)
	counter = 0
	for iy in 1:yL
		for ix in 1:xL
			counter += 1
			V = diaLVC(X[ix],Y[iy],ω,a,Δ,c)
			E,U = adiabatize(V)
			Earr[counter,:] .= E
			Uarr[counter,:,:] = U
		end
	end

	return Earr,Uarr
end


#= Adiabatic building section
function dia2adiCleanBuild(X,Y,ω,a,Δ,c)
	xL = length(X)
	yL = length(Y)

	Earr = zeros(xL,yL,2)
	Uarr = zeros(xL,yL,2,2)
	counter = 0
	for iy in 1:yL
		for ix in 1:xL
			V = diaLVC(X[ix],Y[iy],ω,a,Δ,c)
			E,U = adiabatize(V)
			Earr[ix,iy,:] .= E
			Uarr[ix,iy,:,:] = U
		end
	end

	return Earr,Uarr
end

function LVCadi(N,Nquad,ω,a,Δ,c)
	nodes,weights = ghQuad(Nquad) #get nodes and weights for GH quadrature
	ngrid = length(nodes)
	
	#find points for GH quadrare and adiabatic energies of all points
	X = nodes ./sqrt(ω[1])
	Y = nodes ./sqrt(ω[2])
	Earr,Uarr = dia2adiCleanBuild(X,Y,ω,a,Δ,c)
	
	#build hermite polinomials basis centered in 0
	herPols = hermitePol(N)
	#build N wavefunctions corresponding to each hermite polinomial (0 to N-1)
	Ψ = zeros(ngrid,N)
	for bst in 1:N
		for (i,x) in enumerate(nodes)
			Ψ[i,bst] = polyval(herPols[bst,1:bst],x)
		end
	end

	Exgs = zeros(N,N)
	Eygs = zeros(N,N)
	Exes = zeros(N,N)
	Eyes = zeros(N,N)
	for st1 in 1:N
		for st2 in 1:st1
			#first add potential contributions
			G = Ψ[:,st1] .* Ψ[:,st2]
			Exes[st1,st2] += sum(G .* Vx .* weights)
			Eyes[st1,st2] += sum(G .* Vy .* weights)
			Eyc[st1,st2] += sum(G .* Cy .* weights)
			#add kinetic energy contributions
			if st1 == st2
				Exes[st1,st2] += 0.5*ω[1]*(st1 - 0.5)
				Eyes[st1,st2] += 0.5*ω[2]*(st1 - 0.5)
			elseif st1-st2 == 2 #st1 is always bigger than st2
				Exes[st1,st2] += -0.25*ω[1]*sqrt(st2*(st2+1))
				Eyes[st1,st2] += -0.25*ω[2]*sqrt(st2*(st2+1))
			end
			#hermitize
			Exes[st2,st1] = Exes[st1,st2]
			Eyes[st2,st1] = Eyes[st1,st2]
			Eyc[st2,st1] = Eyc[st1,st2]
		end
	end

	#build full hamiltonian matrix
	m = N^2
	H = zeros(2m,2m)
	H[1:m,1:m] = kron(Diagonal(ones(N)), Exgs) + kron(Eygs, Diagonal(ones(N)))
	H[m+1:2m,m+1:2m] = kron(Diagonal(ones(N)), Exes) + kron(Eyes, Diagonal(ones(N)))
	H[1:m,m+1:2m] = kron(Eyc, Diagonal(ones(N)))
	H[m+1:2m,1:m] = H[1:m,m+1:2m]'


	for nbas in 1:N
		Ψ[:,nbas] .*= exp.(-nodes.^2/2)
	end

	return H,Ψ,nodes,weights,X,Y
end
=#

function HObasisLVC(N,Nquad,ω,a,Δ,c)
	#build ground state energies
	Exgs = zeros(N,N)
	Eygs = zeros(N,N)
	for i in 1:N #harmonic oscillators energies
		Exgs[i,i] = ω[1]*(i-0.5)
		Eygs[i,i] = ω[2]*(i-0.5)
	end

	#find HO operators (generic forms, mode-independent)
	aHat = zeros(N,N)
	adagHat = zeros(N,N)
	for i in 1:N-1
		aHat[i,i+1] = sqrt(i)
		adagHat[i+1,i] = sqrt(i)
	end
	X = 1/sqrt(2*ω[1])*(adagHat + aHat)
	Px = 1im*sqrt(ω[1]/2)*(adagHat - aHat)
	Y = 1/sqrt(2*ω[2])*(adagHat + aHat)
	Py = 1im*sqrt(ω[2]/2)*(adagHat - aHat)
	#Xei,Uxei = adiabatize(X)
	#Yei,Uyei = adiabatize(Y)
	#nodes,weights = ghQuad(Nquad)
	#Xn = nodes ./sqrt(ω[1])
	#Yn = nodes ./sqrt(ω[2])
	### Xei and Xn coincide if Nquad = N
	
	#fill diagonals
	D = 0.5*(ω[1]^2)*(a^2)-Δ
	α = ω[1]^(3/2)*a/sqrt(2)
	Exes = Exgs + D*id(N)
	Eyes = Eygs
	Exes += -α*(aHat + adagHat)

	#coupling matrix
	Eyc = c/sqrt(2*ω[2])*Y

	### polaron transformation
	#γ = a/(ω[1]^2)
	#Ulog = kron(id,1im*a*Px)
	#Upol = exp(Ulog)


	m = N^2
	H = zeros(2m,2m)
	H[1:m,1:m] = kron(id(N), Exgs) + kron(Eygs, id(N))
	H[m+1:2m,m+1:2m] = kron(id(N), Exes) + kron(Eyes, id(N))
	H[1:m,m+1:2m] = kron(Eyc, id(N))
	H[m+1:2m,1:m] = H[1:m,m+1:2m]'

	#POL = real.(Upol*H[m+1:end,m+1:end]*Upol')
	#Epol,_ = adiabatize(POL)
	#@show Epol
	#Etest,_ = adiabatize(H[m+1:end,m+1:end])
	#@show Etest
	#@show [H[n,n] for n in 1:m]
	#display(Upol*H[m+1:end,m+1:end]*Upol')
	#display(Upol'*H[m+1:end,m+1:end]*Upol)
	return H
end

function HOdimensionalbasisLVC(Nx,Ny,Nquad,ω,a,Δ,c)
	#build ground state energies
	Exgs = zeros(Nx,Nx)
	Eygs = zeros(Ny,Ny)
	idx = Diagonal(ones(Nx))
	idy = Diagonal(ones(Ny))
	for i in 1:Nx #harmonic oscillators energies
		Exgs[i,i] = ω[1]*(i-0.5)
	end
	for i in 1:Ny
		Eygs[i,i] = ω[2]*(i-0.5)
	end

	#find HO operators (generic forms, mode-independent)
	aHatx = zeros(Nx,Nx)
	adagHatx = zeros(Nx,Nx)
	aHaty = zeros(Ny,Ny)
	adagHaty = zeros(Ny,Ny)
	for i in 1:Nx-1
		aHatx[i,i+1] = sqrt(i)
		adagHatx[i+1,i] = sqrt(i)
	end
	for i in 1:Ny-1
		aHaty[i,i+1] = sqrt(i)
		adagHaty[i+1,i] = sqrt(i)
	end
	X = 1/sqrt(2*ω[1])*(adagHatx + aHatx)
	Px = 1im*sqrt(ω[1]/2)*(adagHatx - aHatx)
	Xshift = X + idx*5
	Ulog = 1im*5*Px
	Upol = exp(Ulog)
	EX2,U2 = adiabatize(Xshift)
	EX1,U1 = adiabatize(X)
	@show EX1
	@show EX2
	@show maximum(abs.(U1-U2))
	display(U1*U2')
	@show sum(abs2.((Upol*Xshift*Upol')[1:Int(Nx/2),1:Int(Nx/2)] - X[1:Int(Nx/2),1:Int(Nx/2)]))

	return 0

	Y = 1/sqrt(2*ω[2])*(adagHat + aHat)
	Py = 1im*sqrt(ω[2]/2)*(adagHat - aHat)
	#Xei,Uxei = adiabatize(X)
	#Yei,Uyei = adiabatize(Y)
	#nodes,weights = ghQuad(Nquad)
	#Xn = nodes ./sqrt(ω[1])
	#Yn = nodes ./sqrt(ω[2])
	### Xei and Xn coincide if Nquad = N
	
	#fill diagonals
	D = 0.5*(ω[1]^2)*(a^2)-Δ
	α = ω[1]^(3/2)*a/sqrt(2)
	Exes = Exgs + D*id
	Eyes = Eygs
	Exes += -α*(aHat + adagHat)

	#coupling matrix
	Eyc = c/sqrt(2*ω[2])*Y

	### polaron transformation
	#γ = a/(ω[1]^2)
	Ulog = kron(id,1im*a*Px)
	Upol = exp(Ulog)


	m = N^2
	H = zeros(2m,2m)
	H[1:m,1:m] = kron(id, Exgs) + kron(Eygs, id)
	H[m+1:2m,m+1:2m] = kron(id, Exes) + kron(Eyes, id)
	H[1:m,m+1:2m] = kron(Eyc, id)
	H[m+1:2m,1:m] = H[1:m,m+1:2m]'

	POL = real.(Upol*H[m+1:end,m+1:end]*Upol')
	Epol,_ = adiabatize(POL)
	@show Epol
	Etest,_ = adiabatize(H[m+1:end,m+1:end])
	@show Etest
	@show [H[n,n] for n in 1:m]
	#display(Upol*H[m+1:end,m+1:end]*Upol')
	#display(Upol'*H[m+1:end,m+1:end]*Upol)
	return H
end

function singleHO(N,Nquad,ω)
	nodes,weights = ghQuad(Nquad) #get nodes and weights for GH quadrature
	Hho = zeros(N,N)
	Ψgh = zeros(Nquad,N) #wavefunctions for quadrature evaluation, don't include gaussian decaying term
	Q = nodes ./sqrt(ω) #coordinates for Q grid, given by GH quadrature and harmonic frequency	
	herPols = hermitePol(N)

	for bst in 1:N
		for i in 1:Nquad
			Ψgh[i,bst] = polyval(herPols[bst,1:bst],nodes[i])
		end
	end

	for st1 in 1:N
		for st2 in st1:-1:1
			# add contributions over corresponding centered basis
			if st1 == st2
				Hho[st1,st1] = ω*(st1-0.5)
			end
		end
	end	

	return Hho,Ψgh,Q,weights

	#usually don't need to return this part since just Ψgh is needed for overlap calculations
	Ψho = zeros(Nquad,N) #real wavefunctions, include decaying exponential term
	for nbas in 1:Nx
		Ψho[:,nbas] .= Ψho[:,nbas] .* exp.(-nodesxplus.^2/2) .* exp.(nodes.^2/2)
	end

end

#κ: linear potential term, accounts for displacement of acceptor HO states
#Δ: constant energy term added to all acceptor states
function donorAcceptorHO(N,Nquad,ω,κ,Δ,Nextra=10)
	# use larger Ntot=N+Nextra for building, stay with desired N states in the end (to avoid inaccuracies from boundary values)
	Ntot = N + Nextra

	#find HO operators (generic forms, mode-independent)
	aHat = zeros(Ntot,Ntot)
	adagHat = zeros(Ntot,Ntot)
	for i in 1:Ntot-1
		aHat[i,i+1] = sqrt(i)
		adagHat[i+1,i] = sqrt(i)
	end
	Xop = 1/sqrt(2*ω)*(adagHat + aHat)
	nHat = zeros(Ntot,Ntot)
	for i in 1:Ntot
		nHat[i,i] = i-1
	end

	nodes,weights = ghQuad(Nquad)
	Q = nodes ./sqrt(ω)

	Hgs = ω * (nHat + 0.5*id(Ntot))
	#Hes = Hgs + κ*Xop + Δ*id(Ntot)

	D = 0.5*(ω^2)*(κ^2)
	α = ω^(3/2)*κ/sqrt(2)
	Hes = Hgs + (D+Δ)*id(Ntot) - α*(aHat + adagHat)
	
	herPols = hermitePol(Ntot)
	#build Ntot wavefunctions corresponding to each hermite polinomial (0 to Ntot-1)
	ΨquadGS = zeros(Nquad,Ntot)
	for bst in 1:Ntot
		for (i,x) in enumerate(nodes)
			ΨquadGS[i,bst] = polyval(herPols[bst,1:bst],x)
		end
	end

	Egs = diagEls(Hgs)
	Ees,Ues = adiabatize(Hes)
	ΨquadES = zeros(size(ΨquadGS))
	for n1 in 1:Ntot
		for n2 in 1:Ntot
			ΨquadES[:,n1] .+= Ues[n2,n1] * ΨquadGS[:,n2]
		end
	end

	return Diagonal(Egs[1:N]),Diagonal(Ees[1:N]),ΨquadGS[:,1:N],ΨquadES[:,1:N],weights,Q

	#=
	#real wavefunction building, include exponential decay terms
	Ψgs = zeros(size(Ψ))
	Ψes = zeros(size(Ψ))
	for nbas in 1:N
		Ψgs[:,nbas] = ΨquadGS[:,nbas] .* exp.(-nodes.^2/2)
	end
	for n1 in 1:N
		for n2 in 1:N
			Ψes[:,n1] .+= Ues[n2,n1] * Ψgs[:,n2]
		end
	end

	S = zeros(N,N,Nquad)
	for st1 in 1:N
		for st2 in 1:N
			S[st1,st2,:] = ΨquadGS[:,st1] .* ΨquadES[:,st2] .* weights
		end
	end
	=#
end

function multiDimensionalLVC(Narr,Nquad,Ω,Ai,Bi,Ci,Δ,Aij,Bij,Cij,Nextra=10)
	M = length(Ω) #total number of dimensions
	if length(Narr) == 1 && M > 1
		Narr = Narr * ones(M) #use same number of basis functions for all dimensions if only one value is specified
	end

	HsGS = [zeros(n,n) for n in Narr]
	HsES = [zeros(n,n) for n in Narr]
	ΨarrGS = [zeros(Nquad,n) for n in Narr]
	ΨarrES = [zeros(Nquad,n) for n in Narr]
	Qarr = [zeros(Nquad) for n in Narr]
	WTarr = [zeros(Nquad) for n in Narr]
	Sarr = [zeros(n,n) for n in Narr]
	linCoupsArr = [zeros(n,n) for n in Narr]

	if sum(abs.(Aij)) != 0 || sum(abs.(Bij)) != 0 || sum(abs.(Cij)) != 0
		error("Not implemented for aij,bij or cij not 0!")
	end

	#Hamiltonian building routine over all dimensions
	for dim in 1:M
		if Ai[dim] != 0
			error("Ai is not zero, not implemented!")
		end

		if Bi[dim] == 0
			HsGS[dim],ΨarrGS[dim],Qarr[dim],WTarr[dim] = singleHO(Narr[dim],Nquad,Ω[dim])
			HsES[dim] = HsGS[dim]
			ΨarrES[dim] = ΨarrGS[dim]
			Sarr[dim] = id(Narr[dim])
		else
			HsGS[dim],HsES[dim],ΨarrGS[dim],ΨarrES[dim],WTarr[dim],Qarr[dim] = donorAcceptorHO(Narr[dim],Nquad,Ω[dim],Bi[dim],0,Nextra)
			for st1 in 1:Narr[dim]
				for st2 in 1:Narr[dim]
					for quad in 1:Nquad
						Sarr[dim][st1,st2] += ΨarrGS[dim][quad,st1] * ΨarrES[dim][quad,st2] * WTarr[dim][quad]
					end
				end
			end
		end

		#Linear coupling routine (over Ci's)
		if Ci[dim] != 0
			Cy = Ci[dim] * Qarr[dim]
			for st1 in 1:Narr[dim]
				for st2 in 1:Narr[dim]
					#first add potential contributions
					G = ΨarrGS[dim][:,st1] .* ΨarrES[dim][:,st2]
					linCoupsArr[dim][st1,st2] = sum(G .* Cy .* WTarr[dim])
				end
			end
		end
	end

	if M > 1
		HtotGS = kronsum(HsGS...)
		HtotES = kronsum(HsES...)
		COUPS = kronsum(linCoupsArr...)
		Stot = kron(Sarr...)
	else
		HtotGS = HsGS[1]
		HtotES = HsES[1]
		COUPS = linCoupsArr[1]
		Stot = Sarr[1]
	end



	krondim = length(HtotGS[:,1])
	Htot = zeros(2*krondim,2*krondim)
	Htot[1:krondim,1:krondim] = HtotGS - Δ*id(krondim)
	Htot[krondim+1:2*krondim,krondim+1:2*krondim] = HtotES + Δ*id(krondim)
	Htot[1:krondim,krondim+1:2*krondim] = COUPS
	Htot[krondim+1:2*krondim,1:krondim] = COUPS'

	MU = muBuild(id(2*krondim),1,krondim,Stot)

	return Htot,MU,Stot,ΨarrGS,ΨarrES,Qarr,WTarr
end

function MUtreat(Htot,MU)
	E = diagEls(Htot)
	eSts = length(E)
	ind = sortperm(E)
	if ind[1] != 1
		println("Warning! First index is not 1, ground diabatic state is not total ground state")
	end

	E .= E[ind]
	U = id(eSts)
	U = U[:,ind]
	MU = dia2en(MU,U)

	ρ0 = zeros(size(Htot))
	ρ0[1,1] = 1
	MUeffdia = zeros(eSts) #if ρ0 is in diabatic basis
	pop = diagEls(ρ0)
	
	for i in 1:eSts
		MUeffdia[i] = sum(abs2.(MU[i,:]) .* pop)
	end

	return E,MUeffdia
end

function quickcheck(Narr,Ω,Bi,numdims=length(Narr),Δ=0.1,Nquad=100,Nextra=10,Ai=zeros(numdims),Ci=zeros(numdims),Aij=zeros(numdims,numdims),Bij=zeros(numdims,numdims),Cij=zeros(numdims,numdims))
	Htot,MU,_ = multiDimensionalLVC(Narr,Nquad,Ω,Ai,Bi,Ci,Δ,Aij,Bij,Cij,Nextra)
	E,MUeff = MUtreat(Htot,MU)

	plot(E,MUeff)
end

