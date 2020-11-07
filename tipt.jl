#time independent perturbation theory to obtain 1st-order correction to wavefunctions
#(more orders to come)

function tiptO1(E0,V)
	Hsts = length(E0)
	Vsts = Int(Hsts/2)
	Ψvecs = id(Hsts,diag=false) #perturbed basis in 0th-order representation
	Cvecs = zeros(Hsts,Hsts,2) #will hold coefficients for perturbation expansion
	#println("extended V system...")
	#@time Vext = kron([0 1;1 0],V)

	#first order energy correction is zero since Vkk0 is zero for energy transfer models
	println("Starting first order calculation...")
	#E1 = zeros(Hsts)
	println("First order energy correction is zero for energy transfer model")
	
	@time for st1 in 1:Vsts
		for st2 in 1:Vsts
			#|Ψm(i)> = sum(k) ckm(i)|Ψk(0)>
			#ckm(1) = Vkm(0)/ωmk
			Cvecs[st1,st2+Vsts,1] = V[st1,st2] / (E0[st2+Vsts] - E0[st1])
			Cvecs[st1+Vsts,st2,1] = V[st1,st2] / (E0[st2] - E0[st1+Vsts])
		end
	end
	Ψ1 = Ψvecs + Cvecs[:,:,1]
	# =
	println("Normalizing first order for answer...")
	@time for st in 1:Hsts
		Ψ1[:,st] /= norm(Ψ1[:,st])
	end
	return Ψ1
end

#E0 diagonal matrix with diagonal of energies, V matrix with diabatic coupling (size is Hsts/2). Full V matrix would be (0 1;1 0)⊗V
function tiptO2(E0,V)
	Hsts = length(E0)
	Vsts = Int(Hsts/2)
	Ψvecs = id(Hsts,diag=false) #perturbed basis in 0th-order representation
	Cvecs = zeros(Hsts,Hsts,2) #will hold coefficients for perturbation expansion
	#println("extended V system...")
	#@time Vext = kron([0 1;1 0],V)

	#first order energy correction is zero since Vkk0 is zero for energy transfer models
	println("Starting first order calculation...")
	#E1 = zeros(Hsts)
	println("First order energy correction is zero for energy transfer model")
	
	@time for st1 in 1:Vsts
		for st2 in 1:Vsts
			#|Ψm(i)> = sum(k) ckm(i)|Ψk(0)>
			#ckm(1) = Vkm(0)/ωmk
			Cvecs[st1,st2+Vsts,1] = V[st1,st2] / (E0[st2+Vsts] - E0[st1])
			Cvecs[st1+Vsts,st2,1] = V[st1,st2] / (E0[st2] - E0[st1+Vsts])
		end
	end
	Ψ1 = Ψvecs + Cvecs[:,:,1]
	# =
	println("Normalizing first order for answer...")
	@time for st in 1:Hsts
		Ψ1[:,st] /= norm(Ψ1[:,st])
	end
	display(Ψ1)
	# =#

	println("Starting second order calculation...")
	#E2 = zeros(Hsts)
	println("Second order energy correction is zero for energy transfer model")
	#= inneficient way, does not use symmetry of V
	@time for st1 in 1:Hsts
		for st2 in 1:Hsts
			if st1 == st2
				#cmm(2) = -sum(j) abs2(cjm(1))
				Cvecs[st1,st1,2] = -sum(abs2.(Cvecs[:,st1,1]))
			else
				#ckm(2) = Vkm(1)/ωmk
				#Vkm(i) = sum(j) cjm(i) Vkj(0)
				vkm1 = dot(Cvecs[:,st2,1],Vext[st1,:])
				Cvecs[st1,st2,2] = vkm1 / (E0[st2] - E0[st1])
			end
		end
	end
	# =#
	@time for st in 1:Hsts
		Cvecs[st,st,2] = -sum(abs2.(Cvecs[:,st,1]))
	end
	println("Finished diagonal correction, passing to off-diagonal corrections...")
	@time for vst1 in 1:Vsts
		println("vst1 = $vst1")
		@time for vst2 in 1:Vsts
			if vst1 != vst2
				#vkm1 = dot(Cvecs[:,st2,1],Vext[st1,:])
				#case 1: st2 ∈ Donor ≡ 1 (acceptor ≡ 2)
					#st2 ∈ 1 -> Cvecs[i,st2,1] = 0 ∀ i ∈ 1
					#not zero only if Vext[st1,i] ≠ 0 ∀ i ∈ 2
					#not zero only if st1 ∈ 1
				Cvecs[vst1,vst2,2] = dot(Cvecs[Vsts+1:end,vst2,1],V[vst1,:]) / (E0[vst2] - E0[vst1])
				#case 2: st2 ∈ 2
					#st2 ∈ 2 -> Cvecs[i,st2,1] = 0 ∀ i ∈ 2
					#not zero only if Vext[st1,i] ≠ 0 ∀ i ∈ 1
					#not zero only if st1 ∈ 2
				Cvecs[vst1+Vsts,vst2+Vsts,2] = dot(Cvecs[1:Vsts,vst2+Vsts,1],V[vst1,:]) / (E0[vst2+Vsts] - E0[vst1+Vsts])
			end
		end
	end

	Ψvecs += Cvecs[:,:,1] + Cvecs[:,:,2]
	# =
	NORMS = zeros(Hsts)
	@time for st in 1:Hsts
		NORMS[st] = norm(Ψvecs[:,st])
		Ψvecs[:,st] /= NORMS[st]
	end
	display(NORMS)
	# =#

	return Ψ1,Ψvecs
end

function tiptCheck(H,ϕ0,ρex)
	Hsts = length(ϕ0)
	Vsts = Int(Hsts/2)
	V = H[Vsts+1:end,1:Vsts]

	Ψ1,Ψ2 = tiptO2(diagEls(H),V)

	ϕ01 = Ψ1' * ϕ0
	ϕ02 = Ψ2' * ϕ0
	Pμ0 = Diagonal(diagEls(ϕ0 * ϕ0'))
	Pμ1 = zeros(Hsts,Hsts)
	Pμ2 = zeros(Hsts,Hsts)
	println("Building Pμ's")
	@time Pμ1 = dia2en(Diagonal(diagEls(ϕ01 * ϕ01')),Ψ1')
	@time Pμ2 = dia2en(Diagonal(diagEls(ϕ02 * ϕ02')),Ψ2')

	println("Franck-Condon")
	@time densityQualityTrack(ϕ0 * ϕ0',H,ρex)
	println("Zeroth order quality")
	@time densityQualityTrack(Pμ0,H,ρex)
	println("First order quality")
	@time densityQualityTrack(Pμ1,H,ρex)
	println("Second order quality")
	@time densityQualityTrack(Pμ2,H,ρex)
end