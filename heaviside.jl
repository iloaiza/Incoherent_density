# X: coordinate grid. Ψi and Ψj are numerical values of functions in that grid, x0 value. Returns matrix element of <Ψi|Θ(x0-x)|Ψj>
function Hθ(x0,x)
	if x0 < x
		return 1
	else
		return 0
	end 
end

function heavisideElement(Ψi,Ψj,x0,X)
	θ = 0.0
	dx = X[2] - X[1]
	for (k,x) in enumerate(X)
		θ += conj(Ψi[k]) * Ψj[k] * Hθ(x0,x)
	end

	return θ*dx
end

function heavisideInvElement(Ψi,Ψj,x0,X)
	θ = 0.0
	dx = X[2] - X[1]
	for (k,x) in enumerate(X)
		θ += conj(Ψi[k]) * Ψj[k] * Hθ(x,x0)
	end

	return θ*dx
end

function heavisideRhodopsinMats(Nx,Nϕ,res,cutConst=3,V0=V0rhod,V1=V1rhod,m_inv=Iinvrhod)
	Φgrid = range(-π/2,stop=3π/2,length=res)
	## Building ϕ eigenfunctions. First build Hamiltonian from plane wave basis, use cutConst as in energy building
	Nϕeff = Int(Nϕ * cutConst)
	kinvec = zeros(Nϕeff)
	for i in 1:Nϕeff
		kinvec[i] = i^2*m_inv/2
	end
	Eϕeffgs = SymTridiagonal(kinvec .+ V0/2, -V0/4 .* ones(Nϕeff-1))
	Eϕeffes = SymTridiagonal(kinvec .- V1/2,  V1/4 .* ones(Nϕeff-1))
	_,Ugs = adiabatize(Eϕeffgs)
	_,Ues = adiabatize(Eϕeffes)
	Ψϕ = zeros(Complex,Nϕeff,res)
	for k in 1:Nϕeff
		Ψϕ[k,:] .= exp.(1im*k*Φgrid)/sqrt(2*π)
	end
	Ψϕgs = zeros(Complex,Nϕ,res)
	Ψϕes = zeros(Complex,Nϕ,res)
	for n1 in 1:Nϕ
		for n2 in 1:Nϕeff
			Ψϕgs[n1,:] .+= Ugs[n2,n1] * Ψϕ[n2,:]
			Ψϕes[n1,:] .+= Ues[n2,n1] * Ψϕ[n2,:]
		end
	end

	Θtrans = zeros(2Nϕ,2Nϕ) #Capital Theta, not θ
	Θcis = zeros(2Nϕ,2Nϕ)
	for n1 in 1:Nϕ
		for n2 in 1:Nϕ
			Θtrans[n1,n2] = real(heavisideElement(Ψϕgs[n1,:],Ψϕgs[n2,:],π/2,Φgrid))
			Θtrans[n1+Nϕ,n2] = real(heavisideElement(Ψϕes[n1,:],Ψϕgs[n2,:],π/2,Φgrid))
			Θtrans[n1,n2+Nϕ] = real(heavisideElement(Ψϕgs[n1,:],Ψϕes[n2,:],π/2,Φgrid))
			Θtrans[n1+Nϕ,n2+Nϕ] = real(heavisideElement(Ψϕes[n1,:],Ψϕes[n2,:],π/2,Φgrid))

			Θcis[n1,n2] = real(heavisideInvElement(Ψϕgs[n1,:],Ψϕgs[n2,:],π/2,Φgrid))
			Θcis[n1+Nϕ,n2] = real(heavisideInvElement(Ψϕes[n1,:],Ψϕgs[n2,:],π/2,Φgrid))
			Θcis[n1,n2+Nϕ] = real(heavisideInvElement(Ψϕgs[n1,:],Ψϕes[n2,:],π/2,Φgrid))
			Θcis[n1+Nϕ,n2+Nϕ] = real(heavisideInvElement(Ψϕes[n1,:],Ψϕes[n2,:],π/2,Φgrid))
		end
	end

	return Θcis, Θtrans
end

function heavisideRhodopsinDiabaticMats(Nϕ=200,kmax=249,res=600,V0=V0rhod,V1=V1rhod,m_inv=Iinvrhod,E0=E0rhod,E1=E1rhod,cosine=true)
	_,_,_,_,_,_,_,Ψcis,Ψtrans,Φgrid = ϕPot(cosine,Nϕ,kmax,res,V0,V1,E0,E1,m_inv)

	ΘtransGS = zeros(Nϕ,Nϕ) #Capital Theta, not θ
	ΘcisGS = zeros(Nϕ,Nϕ)
	ΘtransES = zeros(Nϕ,Nϕ) #Capital Theta, not θ
	ΘcisES = zeros(Nϕ,Nϕ)
	for n1 in 1:Nϕ
		for n2 in 1:Nϕ
			ΘtransGS[n1,n2] = real(heavisideElement(Ψtrans[n1,:],Ψtrans[n2,:],π/2,Φgrid))
			ΘcisGS[n1,n2] = real(heavisideInvElement(Ψcis[n1,:],Ψcis[n2,:],π/2,Φgrid))

			ΘtransES[n1,n2] = real(heavisideInvElement(Ψtrans[n1,:],Ψtrans[n2,:],π/2,Φgrid))
			ΘcisES[n1,n2] = real(heavisideElement(Ψcis[n1,:],Ψcis[n2,:],π/2,Φgrid))			
		end
	end

	return ΘcisGS, ΘtransGS, ΘcisES, ΘtransES
end

function projectorsHeavisideRhodopsinDiabatic(Nx=20,Nϕ=200,kmax=249,res=600,V0=V0rhod,V1=V1rhod,m_inv=Iinvrhod,E0=E0rhod,E1=E1rhod,cosine=true)
	ΘcisGS,ΘtransGS, ΘcisES, ΘtransES = heavisideRhodopsinDiabaticMats(Nϕ,kmax,res,V0,V1,m_inv,E0,E1,cosine)
	diaGS = [1 0;0 0]
	diaES = [0 0;0 1]
	CISgs = kron(id(Nx),ΘcisGS)
	TRANSgs = kron(id(Nx),ΘtransGS)
	CISes = kron(id(Nx),ΘcisES)
	TRANSes = kron(id(Nx),ΘtransES)
	Pcis0 = kron(diaGS,CISgs)
	Pcis1 = kron(diaGS,CISes)
	Ptrans0 = kron(diaES,TRANSes)
	Ptrans1 = kron(diaES,TRANSgs)

	return Pcis0,Pcis1,Ptrans0,Ptrans1
end

function projectorsHeavisideRhodopsinAdiabatic(U,Nx=20,Nϕ=200,kmax=249,res=600,V0=V0rhod,V1=V1rhod,m_inv=Iinvrhod,E0=E0rhod,E1=E1rhod,cosine=true)
	ΘcisGS,ΘtransGS, ΘcisES, ΘtransES = heavisideRhodopsinDiabaticMats(Nϕ,kmax,res,V0,V1,m_inv,E0,E1,cosine)
	diaGS = [1 0;0 0]
	diaES = [0 0;0 1]
	CISgs = kron(id(Nx),ΘcisGS)
	TRANSgs = kron(id(Nx),ΘtransGS)
	CISes = kron(id(Nx),ΘcisES)
	TRANSes = kron(id(Nx),ΘtransES)
	Pcis0 = kron(diaGS,CISgs)
	Pcis1 = kron(diaGS,CISes)
	Ptrans0 = kron(diaES,TRANSes)
	Ptrans1 = kron(diaES,TRANSgs)

	return dia2en(Pcis0,U),dia2en(Pcis1,U),dia2en(Ptrans0,U),dia2en(Ptrans1,U)
end

function projectorsHeavisideRhodopsin(Nx,Nϕ,U,res=Int(1e4),cutConst=3,V0=V0rhod,V1=V1rhod,m_inv=Iinvrhod)
	Θcis,Θtrans = heavisideRhodopsinMats(Nx,Nϕ,res,cutConst,V0,V1,m_inv)
	diaGS = [1 0;0 0]
	diaES = [0 0;0 1]
	CISred = kron(Θcis,id(Nx))
	TRANSred = kron(Θtrans,id(Nx))
	Pcis0 = kron(diaGS,CISred)
	Pcis1 = kron(diaES,CISred)
	Ptrans0 = kron(diaGS,TRANSred)
	Ptrans1 = kron(diaES,TRANSred)

	return dia2en(Pcis0,U),dia2en(Pcis1,U),dia2en(Ptrans0,U),dia2en(Ptrans1,U)
end

function heavisideRhodopsinRedMats(Nx,Nϕ,res,cutConst=3,V0=V0rhod,V1=V1rhod,m_inv=Iinvrhod)
	Φgrid = range(-π/2,stop=3π/2,length=res)
	## Building ϕ eigenfunctions. First build Hamiltonian from plane wave basis, use cutConst as in energy building
	Nϕeff = Int(Nϕ * cutConst)
	kinvec = zeros(Nϕeff)
	for i in 1:Nϕeff
		kinvec[i] = i^2*m_inv/2
	end
	Eϕeffgs = SymTridiagonal(kinvec .+ V0/2, -V0/4 .* ones(Nϕeff-1))
	Eϕeffes = SymTridiagonal(kinvec .- V1/2,  V1/4 .* ones(Nϕeff-1))
	_,Ugs = adiabatize(Eϕeffgs)
	_,Ues = adiabatize(Eϕeffes)
	Ψϕ = zeros(Complex,Nϕeff,res)
	for k in 1:Nϕeff
		Ψϕ[k,:] .= exp.(1im*k*Φgrid)/sqrt(2*π)
	end
	Ψϕgs = zeros(Complex,Nϕ,res)
	Ψϕes = zeros(Complex,Nϕ,res)
	for n1 in 1:Nϕ
		for n2 in 1:Nϕeff
			Ψϕgs[n1,:] .+= Ugs[n2,n1] * Ψϕ[n2,:]
			Ψϕes[n1,:] .+= Ues[n2,n1] * Ψϕ[n2,:]
		end
	end

	Θtrans = zeros(Nϕ,Nϕ) #Capital Theta, not θ
	Θcis = zeros(Nϕ,Nϕ)
	for n1 in 1:Nϕ
		for n2 in 1:Nϕ
			Θtrans[n1,n2] = real(heavisideElement(Ψϕgs[n1,:],Ψϕgs[n2,:],π/2,Φgrid))
			
			Θcis[n1,n2] = real(heavisideInvElement(Ψϕes[n1,:],Ψϕes[n2,:],π/2,Φgrid))
		end
	end

	return Θcis, Θtrans
end

function projectorsHeavisideRhodopsinRed(Nx,Nϕ,U,res=Int(1e4),cutConst=3,V0=V0rhod,V1=V1rhod,m_inv=Iinvrhod)
	Θcis,Θtrans = heavisideRhodopsinRedMats(Nx,Nϕ,res,cutConst,V0,V1,m_inv)
	diaGS = [1 0;0 0]
	diaES = [0 0;0 1]
	CISred = kron(Θcis,id(Nx))
	TRANSred = kron(Θtrans,id(Nx))
	Pcis0 = kron(diaGS,CISred)
	Pcis1 = kron(diaES,CISred)
	Ptrans0 = kron(diaGS,TRANSred)
	Ptrans1 = kron(diaES,TRANSred)

	return dia2en(Pcis0,U),dia2en(Pcis1,U),dia2en(Ptrans0,U),dia2en(Ptrans1,U)
end

function projectorsHeavisideRhodopsinRedDiaBasis(Nx=20,Nϕ=200,res=800,cutConst=1.25,V0=V0rhod,V1=V1rhod,m_inv=Iinvrhod)
	println("Building heaviside ϕ projectors")
	@time Θcis,Θtrans = heavisideRhodopsinRedMats(Nx,Nϕ,res,cutConst,V0,V1,m_inv)
	diaGS = [1 0;0 0]
	diaES = [0 0;0 1]
	println("Starting Kronecker products for electronic states")
	CISred = kron(Θcis,id(Nx))
	TRANSred = kron(Θtrans,id(Nx))
	Pcis0 = kron(diaGS,CISred)
	Pcis1 = kron(diaES,CISred)
	Ptrans0 = kron(diaGS,TRANSred)
	Ptrans1 = kron(diaES,TRANSred)

	return Pcis0,Pcis1,Ptrans0,Ptrans1
end

function heavisideLVCMats(H,Ψ,nodes,weights,X,N)
	#H,Ψ,nodes,weights,X,Y = LVCenergies(N,Nquad,ω,a,Δ,c)
	#remove exponential part (included in weights for GH quadrature)
	for i in 1:N
		Ψ[:,i] .*= exp.(nodes.^2/2)
	end
	E,U = adiabatize(H)

	xsplit = a/2
	xsplitHeaviside = Hθ.(xsplit,X)
	xsplitHeavisideInv = Hθ.(X,xsplit)
	Θtrans = zeros(N,N) #Capital Theta, not θ
	Θcis = zeros(N,N)
	for n1 in 1:N
		for n2 in 1:N
			Θtrans[n1,n2] = sum(Ψ[:,n1] .* Ψ[:,n2] .* xsplitHeaviside .* weights)
			Θcis[n1,n2] = sum(Ψ[:,n1] .* Ψ[:,n2] .* xsplitHeavisideInv .* weights)
		end
	end

	return Θcis, Θtrans
end

function projectorsHeavisideLVC(N,Nquad,ω,a,Δ,c)
	H,Ψ,nodes,weights,X,Y = LVCenergies(N,Nquad,ω,a,Δ,c)
	E,U = adiabatize(H)
	Θcis,Θtrans = heavisideLVCMats(H,Ψ,nodes,weights,X,N)
	diaGS = [1 0;0 0]
	diaES = [0 0;0 1]
	CISred = kron(id(N),Θcis)
	TRANSred = kron(id(N),Θtrans)
	Pcis0 = kron(diaGS,CISred)
	Pcis1 = kron(diaES,CISred)
	Ptrans0 = kron(diaGS,TRANSred)
	Ptrans1 = kron(diaES,TRANSred)
	
	return dia2en(Pcis0,U),dia2en(Pcis1,U),dia2en(Ptrans0,U),dia2en(Ptrans1,U),U
end

function heavisideLVCMatsSeparate(H,ΨGS,ΨES,nodes,weights,X,a)
	#H,Ψ,nodes,weights,X,Y = LVCenergies(N,Nquad,ω,a,Δ,c)
	#remove exponential part (included in weights for GH quadrature)
	Nx = length(ΨGS[1,:])
	for i in 1:Nx
		ΨGS[:,i] .*= exp.(nodes.^2/2)
		ΨES[:,i] .*= exp.(nodes.^2/2)
	end
	E,U = adiabatize(H)

	xsplit = a/2
	xsplitHeaviside = Hθ.(xsplit,X)
	xsplitHeavisideInv = Hθ.(X,xsplit)
	Θtrans0 = zeros(Nx,Nx) #Capital Theta, not θ
	Θcis0 = zeros(Nx,Nx)
	Θtrans1 = zeros(Nx,Nx) #Capital Theta, not θ
	Θcis1 = zeros(Nx,Nx)
	for n1 in 1:Nx
		for n2 in 1:Nx
			Θtrans0[n1,n2] = sum(ΨES[:,n1] .* ΨES[:,n2] .* xsplitHeaviside .* weights)
			Θcis0[n1,n2] = sum(ΨGS[:,n1] .* ΨGS[:,n2] .* xsplitHeavisideInv .* weights)
			Θtrans1[n1,n2] = sum(ΨGS[:,n1] .* ΨGS[:,n2] .* xsplitHeaviside .* weights)
			Θcis1[n1,n2] = sum(ΨES[:,n1] .* ΨES[:,n2] .* xsplitHeavisideInv .* weights)
		end
	end

	return Θcis0, Θtrans0, Θcis1, Θtrans1
end

function projectorsHeavisideLVCSeparate(Nx,Ny,ω,a,Δ,c,Nquad)
	H,_,nodes,weights,X,_,_,_,ΨGS,ΨES,_,_,_ = LVCenergiesSeparate(Nx,Ny,ω,a,Δ,c,Nquad) #H,Ψy,nodes,weights,X,Y,Exgs,Exes,ΨxGS,ΨxES,nodesxplus,nodesxminus,S 
	E,U = adiabatize(H)
	#a=0 is where the wells are split
	Θcis0,Θtrans0,Θcis1,Θtrans1 = heavisideLVCMatsSeparate(H,ΨGS,ΨES,nodes,weights,X,0)
	diaGS = [1 0;0 0]
	diaES = [0 0;0 1]
	Pcis0 = kron(diaGS,id(Ny),Θcis0)
	Pcis1 = kron(diaGS,id(Ny),Θcis1)
	Ptrans0 = kron(diaGS,id(Ny),Θtrans0)
	Ptrans1 = kron(diaGS,id(Ny),Θtrans1)
	
	return dia2en(Pcis0,U),dia2en(Pcis1,U),dia2en(Ptrans0,U),dia2en(Ptrans1,U),U
end

function planeChooser(a,x)
	proj = sum(a .*(x - a))
	if proj < 0
		return 0
	else
		return 1
	end
end

function heavisideMultiDimLVCMats(Narr,ΨGS,ΨES,Q,weights,Bi)
	a = Bi ./ 2
	Quads = [length(warr) for warr in weigths]
	dims = length(Narr)

	SprodMatsCis0 = [zeros(Narr[i],Narr[i],Quads[i]) for i in 1:dims]
	SprodMatsCis1 = [zeros(Narr[i],Narr[i],Quads[i]) for i in 1:dims]
	SprodMatsTrans0 = [zeros(Narr[i],Narr[i],Quads[i]) for i in 1:dims]
	SprodMatsTrans1 = [zeros(Narr[i],Narr[i],Quads[i]) for i in 1:dims]
	#= UNDER CONSTRUCTION!
	for idim in 1:dims
		for n1 in 1:Narr[idim]
			for n2 in 1:Narr[idim]
				for (iquad,quadval) in enumerate(Q[dim])
					SprodMatsCis0[n1,n2,iquad] = ΨGS[dim][iquad,n1] * ΨGS[dim][iquad,n2] * weights[dim][]
					=#


	xsplitHeaviside = Hθ.(xsplit,X)
	xsplitHeavisideInv = Hθ.(X,xsplit)
	Θtrans0 = zeros(Nx,Nx) #Capital Theta, not θ
	Θcis0 = zeros(Nx,Nx)
	Θtrans1 = zeros(Nx,Nx) #Capital Theta, not θ
	Θcis1 = zeros(Nx,Nx)
	for n1 in 1:Nx
		for n2 in 1:Nx
			Θtrans0[n1,n2] = sum(ΨES[:,n1] .* ΨES[:,n2] .* xsplitHeaviside .* weights)
			Θcis0[n1,n2] = sum(ΨGS[:,n1] .* ΨGS[:,n2] .* xsplitHeavisideInv .* weights)
			Θtrans1[n1,n2] = sum(ΨGS[:,n1] .* ΨGS[:,n2] .* xsplitHeaviside .* weights)
			Θcis1[n1,n2] = sum(ΨES[:,n1] .* ΨES[:,n2] .* xsplitHeavisideInv .* weights)
		end
	end

	return Θcis0, Θtrans0, Θcis1, Θtrans1
end