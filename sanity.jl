#sanity check functions

function unitaritySanityCheck(U,nsts=length(U[:,1]))
	return sum(abs2.(U'*U - id(nsts)))
end

#acceptComplex flag indicates if small imaginary components are accepted (come from numerical noise). If true, will
function densitySanity(ρ,nsts=length(ρ[:,1]),normTol=1e-10;acceptComplex=true,acceptNegative=true)
	flagLifted = false
	ρdiag = diagEls(ρ)
	if typeof(ρ[1,1]) <: Complex
		if acceptComplex == true
			imagTot = sum(abs.(imag.(ρdiag)))
			realTot = sum(abs.(real.(ρdiag)))
			if imagTot > 0
				println("Total imaginary part of density is $imagTot, and real part is $realTot")
				println("Changing diagonal density to real density")
				for st in 1:nsts
					ρ[st,st] = real(ρ[st,st])
				end
				flagLifted = true
			end
			ρdiag = zeros(nsts)
			for st in 1:nsts
				ρdiag[st] = real(ρ[st,st])
			end
		else
			imagTot = sum(abs2.(imag.(ρdiag)))
			if real(imagTot) > 0
				error("Imaginary component in diagonal element of density matrix")
			end 
		end
	end

	#positivity check
	for st in 1:nsts
		if ρdiag[st] < 0
			if acceptNegative == false
				error("Negative element in diagonal density matrix, terminating")
			end
			println("Warning, density matrix has negative elements in diagonal")
			posArr = ρdiag .> 0
			negArr = ρdiag .< 0
			posSum = sum(posArr .* ρdiag)
			negSum = sum(negArr .* ρdiag)
			println("Sum of positive elements is $posSum, while sum of negative is $negSum")
			println("Adding positive element to diagonal")
			negVal = -minimum(ρdiag)
			for st in 1:nsts
				ρ[st,st] += negVal
			end
			flagLifted = true
			break
		end
	end

	ρtrace = sum(ρdiag)
	if abs(ρtrace - 1) > normTol
		println("Trace of density is not one, normalizing...")
		if typeof(ρ) <: Symmetric
			println("Type of density matrix is Symmetric which is not mutable, remove trace elsewhere...")
			return true
		end
		ρ[:] ./= ρtrace
		flagLifted = true
	end

	if flagLifted == false
		println("Density passed all sanity checks")
	end

	return flagLifted
end	
