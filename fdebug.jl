#fast debug module. Quickly creates and diagonalizes some Hamiltonians for several dimensions
#=
println("Debugging random Hamiltonians...")
for n in [2,5,10,20,50,100]
	eval(Meta.parse("H$n = Symmetric(rand($n,$n))"))
	eval(Meta.parse("ϕ$n = rand($n)"))
	eval(Meta.parse("ϕ$n /= sum(abs2.(ϕ$n))"))
	println("Created H$n and ϕ$n")
end
# =#

println("1D-LVC debugging")
Narr1Ddbg = [30]; Nquad1Ddbg = 50; Ω1Ddbg = [2]; Bi1Ddbg = [-3.0]; Ci1Ddbg = [1.7]; Δ1Ddbg = 2; Nextra1Ddbg=10; Ai1Ddbg=0; Aij=0;Bij=0;Cij=0;
H1D,MU1D,_,_,_,_,_ = multiDimensionalLVC(Narr1Ddbg,Nquad1Ddbg,Ω1Ddbg,Ai1Ddbg,Bi1Ddbg,Ci1Ddbg,Δ1Ddbg,Aij,Bij,Cij,Nextra1Ddbg)
Hsts1D = length(H1D[:,1])
ϕ01D = zeros(Hsts1D)
for st in 1:Hsts1D
	ϕ01D[st] = MU1D[st,1]
end
ϕ01D /= sqrt(dot(ϕ01D,ϕ01D))
E1D,U1D = adiabatize(H1D)
τ1D = E1D[3]
u1D = U1D[:,3] + 0.1*rand(60); u1D /= norm(u1D);


println("Created H1D, MU1D and ϕ01D of with Hsts = $Hsts1D")

#=
println("2D-BIG-LVC debugging")
Narr2Ddbg = [40,40]; Nquad2Ddbg = 60; Ω2Ddbg = [2,0.7]; Bi2Ddbg = [1.0,2.1]; Ci2Ddbg = [1.0,0.5]; Δ2Ddbg = 0.5; Nextra2Ddbg=10; Ai2Ddbg=0; Aij=0;Bij=0;Cij=0;
H2D,MU2D,_,_,_,_,_ = multiDimensionalLVC(Narr2Ddbg,Nquad2Ddbg,Ω2Ddbg,Ai2Ddbg,Bi2Ddbg,Ci2Ddbg,Δ2Ddbg,Aij,Bij,Cij,Nextra2Ddbg)
Hsts2D = length(H2D[:,1])
ϕ02D = zeros(Hsts2D)
for st in 1:Hsts2D
	ϕ02D[st] = MU2D[st,1]
end
ϕ02D /= sqrt(dot(ϕ02D,ϕ02D))

println("Created H2D, MU2D and ϕ02D of with Hsts = $Hsts2D")

println("MD-LVC debugging")
H3D,MU3D,_,_,_,_,_ = multiDimensionalLVC(N3,Nquad,Ω3,Ai3,Bi3,Ci3,Δ3,0,0,0,10)
Hsts3D = length(H3D[:,1])
ϕ03D = zeros(Hsts3D)
for st in 1:Hsts3D
	ϕ03D[st] = MU3D[st,1]
end
ϕ03D /= sqrt(dot(ϕ03D,ϕ03D))

H5D,MU5D,_,_,_,_,_ = multiDimensionalLVC(N5,Nquad,Ω5,Ai5,Bi5,Ci5,Δ5,0,0,0,10)
Hsts5D = length(H5D[:,1])
ϕ05D = zeros(Hsts5D)
for st in 1:Hsts5D
	ϕ05D[st] = MU5D[st,1]
end
ϕ05D /= sqrt(dot(ϕ05D,ϕ05D))

println("Created HnD, MUnD and ϕ0nD for n=2,3,5")

# =#
#=
println("Rhodopsin debug")
loading("rhodopsin")
# =#
# =
Hrhod,_,Srhod = rhodCosHamiltonian()
Hrhodsts = length(Hrhod[:,1])
Erhod,Urhod = adiabatize(Hrhod)
MUrhod = muBuild(Urhod,1,Int(Hrhodsts/2),Srhod) # returns MU in adiabatic basis
MUrhod = dia2en(MUrhod,Urhod')
ϕ0rhod = MUrhod[:,1]
ρrhod = outerProduct(ϕ0rhod,ϕ0rhod);
ρrhod = dia2en(ρrhod,Urhod);
ρrhod = Diagonal(ρrhod)
ρrhodDia = dia2en(ρrhod,Urhod')
@saving "rhodopsin" Hrhod Erhod Urhod ϕ0rhod MUrhod ρrhod ρrhodDia
rhodlist = loadnames("rhodopsin")
h5Arraywrite("rhodopsin",rhodlist,[Hrhod,Erhod,Urhod,ϕ0rhod,MUrhod,ρrhod,ρrhodDia])
# =#



