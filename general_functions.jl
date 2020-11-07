function polyval(p,x)
	Px = [0.0]
	polDeg=length(p) #degree of polynome +1 (since p includes degree 0)
	for i in 1:polDeg
		Px[1] += p[i]*x^Int(polDeg-i)
	end
	return Px[1]
end


function hermitePol(n)
	p = zeros(n+1,n+1)
	p[1,1] = 1
	if n >= 1
		p[2,1:2] .= [2,0]
		if n >= 2
			for k in 3:n+1
				p[k,1:n] = 2*p[k-1,1:n]
				p[k,3:n+1] += -2*(k-2)*p[k-2,1:n-1]
			end
		end
	end
	#normalization
	for i in 1:n
		p[i,:] ./= Float64(sqrt(sqrt(π)*2^(i-1)*factorial(big(i-1))))
	end 

	return p
end


function ghQuad(n)
  #Gauss-Hermite quadrature obtained from GaussQuadrature package :)
  nodes,weights = gausshermite(n)
  return nodes,weights
end


function polyMult(p1,p2) #multiply polynomials p1 and p2
	n1 = length(p1)
	n2 = length(p2)
	p = zeros(n1+n2-1)

	for i1 in 1:n1
		for i2 in 1:n2
			p[i1+i2-1] += p1[i1]*p2[i2]
		end
	end

	return p
end

#transforms A matrix in diabatic representation to energy eigenbasis representation by transformation U
function dia2en(A,U)
	return U'*A*U
end

function id(n;diag=true,idtype=Float64)
	if diag == true
		return Diagonal(ones(idtype,n))
	elseif diag == false
		Id = zeros(idtype,n,n)
		Id .= id(n,idtype=idtype)
		return Id
	else
		error("Not specified if identity is Diagonal type array or not")
	end
end

#matrix polynomial multiplication by Clenshaw's rule, C is Chebyshev coefficients
# =
function matrixPoly(H,C,Hsts = length(H[:,1]))
	#=N = length(C) - 1
	ρ = zeros(size(H))
	ID = id(Hsts)
	for k in 0:N
		ρ = ID*C[k+1] + H*ρ
	end
	=#
	ρ = clenshawEval(C,H)

	return ρ/tr(ρ)
end
# =#

#= naive polynomial matrix multiplication
function matrixPoly(H,chCoeffs,Id=id(length(H[:,1])))
	Heval = Id
	ρs = zeros(size(H))
	for i in length(chCoeffs):-1:1
		ρs += chCoeffs[i] * Heval
		Heval *= H
	end

	return ρs/tr(ρs)
end
# =#

function adiabatize(H,ascending=true)
	E,U = eigen(H)
	ind = sortperm(E)
	if ascending == false #order in descending order
		ind = ind[end:-1:1]
	end
	E .= E[ind]
	U = U[:,ind]

	return E,U
end

function diagEls(A,L = length(A[:,1]))
	diag = zeros(typeof(A[1,1]),L)
	for i in 1:L
		diag[i] = A[i,i]
	end

	return diag
end

function kronsum(A...)
	Narrs = length(A)
	ArrDims = zeros(Int,Narrs)
	for n in 1:Narrs
		ArrDims[n] = length(A[n][:,1])
	end

	DimTot = prod(ArrDims)
	SUM = zeros(DimTot,DimTot)

	kronArr = [zeros(ArrDims[n],ArrDims[n]) for n in 1:Narrs]
	for n in 1:Narrs
		for n2 in 1:Narrs
			if n2 == n
				kronArr[n2] = A[n2]
			else
				kronArr[n2] = id(ArrDims[n2])
			end
		end
		SUM += kron(kronArr...)
	end

	return SUM
end

function commutator(A,B,Herm=false)
	if Herm == false
		return A*B-B*A
	else
		M = A*B
		return M - M'
	end
end

function matNorm(A)
	E,_ = eigen(A)

	return maximum(E)
end

function matConditioning(A)
	E,_ = eigen(A)

	return minimum(A)/maximum(A)
end

function GrammSchmidt(A,tol=1e-8)
	Ags = zeros(size(A))
	L = length(A[1,:])

	Dflag = Int[]
	for i in 1:L
		dotcum = 0.0
		Ags[:,i] = A[:,i]
		for j in 2:i
			dotcum += dot(Ags[:,i],Ags[:,j-1])
			Ags[:,i] -= dot(Ags[:,i],Ags[:,j-1]) * Ags[:,j-1]
		end
		agsnorm = sqrt(dot(Ags[:,i],Ags[:,i]))
		if agsnorm <= tol
			println("Warning, basis is not linearly dependent in Gramm-Schmidt orthogonalization at index $i")
			push!(Dflag,i)
			Ags[:,i] .= 0
		end
		Ags[:,i] /= agsnorm
		@show i,dotcum
	end

	INDS = collect(1:L)
	dcount = 0
	for dnum in Dflag
		deleteat!(INDS,dnum-dcount)
		dcount += 1
	end	

	#@show INDS

	return Ags[:,INDS]
end

function GrammSchmidtVector(v,A,tol=1e-8)
	L = length(A[1,:])

	for j in 1:L
		v -= dot(v,A[:,j]) * A[:,j]
	end
	vnorm = sqrt(dot(v,v))
	if vnorm < tol
		println("Warning, Gramm-Schmidt orthogonalization returns close to tolerance null vector!")
		v .= 0
		flag = true
	else
		flag = false
		#v /= vnorm
	end

	return v,flag
end



function GrammSchmidtBlock!(Qblock,Barr,tol=1e-8)
	#Qblock = [vecLength,blockSize]
	#Barr = [vecLength,blockSize,numBlocks]
	vecLength = length(Qblock[:,1])
	blockSize = length(Qblock[1,:])
	numBlocks = length(Barr[1,1,:])

	for bnum in 1:numBlocks
		for bpart1 in 1:blockSize
			for bpart2 in 1:blockSize
				Qblock[:,bpart1] -= dot(Qblock[:,bpart1],Barr[:,bnum,bpart2]) * Barr[:,bnum,bpart2]
			end
		end
	end

	for bpart1 in 1:blockSize
		for bpart2 in 1:bpart1 - 1
			dotpart = dot(Qblock[:,bpart1],Qblock[:,bpart2])
			Qblock[:,bpart1] -= dotpart * Qblock[:,bpart2]
			Qblock[:,bpart2] -= dotpart * Qblock[:,bpart1]
		end
	end

	for bpart in 1:blockSize
		qnorm = sqrt(dot(Qblock[:,bpart],Qblock[:,bpart]))
		if qnorm < tol
			println("Gramm-Schmidt block warning, norm is smaller than tolerance!")
		end
		Qblock[:,bpart] /= qnorm
	end

	return Qblock
end
			


#finds the N maximum elements and coordinates of array A
function arrayTopN(A,N)
	ind = sortperm(A)

	return ind[end:-1:end-N+1]
end

function densityQualityTrack(ρ,H,ρex;printout=true,tol=1e-10)
	Hsts = length(H[:,1])
	ρtrace = tr(ρ)
	if abs(ρtrace-1) > tol
		println("Density not normalized, returning relative normalized values")
		ρ /= ρtrace
	end
	dist = sum(abs2.(ρex - ρ))
	coh = tr(ρ^2)
	td = tr((commutator(H,ρ))^2)
	lind = tr((commutator(H,commutator(H,ρ)))^2)

	if printout
		@show ρtrace,dist,coh,td,lind,tr(H*ρ),tr(H^2*ρ)
	end

	return ρtrace,dist,coh,td,lind
end

function quickQuality(ρ,ρex)
	dist = sum(abs2.(ρex - ρ))
	@show dist

	return dist
end

function outerProduct(u,v)
	nsts = length(v)
	if nsts != length(u)
		error("Outer product vectors have different dimensions")
	end

	UV = zeros(nsts,nsts)
	for i in 1:nsts
		for j in 1:nsts
			UV[i,j] = u[i] * v[j]
		end
	end

	return UV
end

function outerProductComplex(u,v)
	nsts = length(v)
	if nsts != length(u)
		error("Outer product vectors have different dimensions")
	end

	UV = zeros(Complex,nsts,nsts)
	for i in 1:nsts
		for j in 1:nsts
			UV[i,j] = conj(u[i]) * v[j]
		end
	end

	return UV
end

function muBuild(U,μ=1,m = Int(length(U[:,1])/2),S=id(m))
	#this function must start with the dipole moment in the diabatic representation and return the dipole matrix in the energy representation
	#U is the diagonalizing transformation for the full energy matrix
	
	MU = zeros(2m,2m)

	MU[1:m,m+1:2m] = μ * S
	MU[m+1:2m,1:m] = MU[1:m,m+1:2m]'

	return dia2en(MU,U)
end
