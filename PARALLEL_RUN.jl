include("parallel.jl")
using HDF5
include("saving.jl")
loading("rhodopsin")

# =
Nrands = 1000
itsItvl = collect(20:50)
#for i in 150:5:300
#	push!(itsItvl,i)
#end
# =#
@time TRANS,S0,COHS = randomParallelAnalysis(Nrands,ϕ0rhod,Hrhod,itsItvl,ρrhod)