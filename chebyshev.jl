function chPol(n::Int)
	p = zeros(n+1,n+1)
	p[1,1] = 1
	if n >= 1
		p[2,1:2] .= [1,0]
		for i in 3:n+1
			p[i,1:n] = 2*p[i-1,1:n]
			p[i,3:n+1] += -p[i-2,1:n-1]
		end
	end
	return p
end

function chRiemannCoeffs(N::Int,f,eps=1e-3,steps=Int(1e4)) #yields N coefficients of decomposition of f function in first N Chebyshev polynomials (T0,..,T(N-1)) in interval [-1+eps,1-eps]
	a = eps-1
	b = 1-eps
	ITVL = range(a,stop=b,length=steps)
	dx=ITVL[2]-ITVL[1]
	POL = chPol(N)

	C=zeros(N)
	for coeffNum in 1:N
		grid=zeros(steps)
		for (i,x) in enumerate(ITVL)
			grid[i]=polyval(POL[coeffNum,1:coeffNum],x)*f(x)/sqrt(1-x^2)
		end
		C[coeffNum] = 2*sum(grid)*dx/pi
		if coeffNum == 1
			C[1] /= 2
		end
	end

	return C
end

function chInterpolCoeffs(N::Int,nodes,F,nodecheck=false) #nodes of Chebyshev polynomial of degree N+1, F is array of function to be decomposed evaluated at said nodes
	if nodecheck
		for (i,x) in enumerate(nodes)
			println("checking node nodes[$i]-x$i=$(nodes[i]-cos((i+0.5)*π/(N+1)))")
		end
	end

	POLS = chPol(N)
	coeffs = zeros(N+1)
	for coeffnum in 1:N+1
		for i in 1:N+1
			coeffs[coeffnum] += F[i]*polyval(POLS[coeffnum,1:coeffnum],nodes[i])
		end
	end
	coeffs *= 2/(N+1)

	return coeffs
end

function chSetCoeffs(N,X,F,eps=1e-5) #yields N coefficients of decomposition of F set over X values in first N Chebyshev polynomials (T0,..,T(N-1))
	if X[1] <= -1
		X[1] = -1+eps
	end
	if X[end] >= 1
		X[end] = 1-eps
	end
	XL = length(X)
	DX=[X[i+1] - X[i] for i in 1:XL-1]
	POL = chPol(N)

	C=zeros(N)
	for coeffNum in 1:N
		grid=zeros(XL)
		for (i,x) in enumerate(X)
			grid[i]=polyval(POL[coeffNum,1:coeffNum],x)*F[i]/sqrt(1-x^2)
		end
		for i in 1:XL-1
			C[coeffNum] += (grid[i+1]+grid[i])*DX[i]
		end
		C[coeffNum] /= pi
		if coeffNum == 1
			C[1] /= 2
		end
	end

	return C
end

function chBuild(C) #builds polynomial coefficients from Chebyshev coefficients
	N = length(C)
	p = zeros(N)
	POL = chPol(N)

	for i in 1:N
		coeff = C[i].*POL[i,1:i]
		p[N+1-i:N] .+= coeff
	end

	return p
end

function chRoots(N) #returns roots of Chebyshev polynoial of degree N
	X=zeros(N)
	for i in 1:N
		X[i]=cos(π*(2i-1)/2/N)
	end

	return X
end

#C coefficients for Chebyshev polynomials. C[1] is coefficient for highest degree polynomial (of degree length(C) - 1)
function clenshawEval(C,H,Hsts = length(H[:,1]))
	N = length(C) - 1
	ID = id(Hsts)
	bkPlus = zeros(Hsts,Hsts)
	bk = zeros(Hsts,Hsts)
	ρ = zeros(Hsts,Hsts)

	for k in N+1:-1:2
		ρ = ID*C[k] + 2*H*bk - bkPlus
		bkPlus = bk
		bk = ρ
	end

	ρ = ID*C[1] + H*bk - bkPlus

	return ρ
end

function chebyshevKernel(coeffs,kernel="Jackson")
	N = length(coeffs)
	G = zeros(N)
	if kernel == "Jackson"
		for i in 0:N-1
			G[i+1] = ((N-i)*cos(π*i/N) + sin(π*i/N)*cot(π/N+1))/N
		end
	elseif kernel == "Fejer"
		for i in 0:N-1
			G[i+1] = 1 - i/(N-1)
		end
	else
		println("chebyshev Kernel $kernel not defined, returning original array")
		G = ones(N)
	end

	return coeffs .* G
end