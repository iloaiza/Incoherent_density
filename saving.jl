#saving module, very practical :)
#use @saving "name" var1 var2 to save in DATAFOLDER ("./" unless defined elsewhere) a "name.h5" file (with movement to .old if preexisting) "var1", "var2", ... values
#use @loading "name" for automatic loading of all saved variables

### Define DATAFOLDER if it's not already defined
if !(@isdefined(DATAFOLDER))
	DATAFOLDER = "./"
end


### Section for implementing saving of "complicated" data types
import HDF5.h5write

function h5write(filename::String,varname::String,S::Symmetric)
	Scopy = zeros(size(S))
	Scopy .= S
	h5write(filename,varname,Scopy)
end

function h5write(filename::String,varname::String,D::Diagonal)
	Dcopy = zeros(size(D))
	Dcopy .= D
	h5write(filename,varname,Dcopy)
end
# =
function h5Arraywrite(filename::String,varnames,vars)
	h5name = DATAFOLDER*filename*".h5"
	for (ivar,var) in enumerate(vars)
		h5write(h5name,varnames[ivar],var)
		println("Saved var $(varnames[ivar])")
	end
end
# =#
### Auxiliary functions for saving macro

function stringSymbolSeparator(params,numparams=length(params))
	paramsNames = String[]
	for i in 1:numparams
		current = "$(params[i])"
		push!(paramsNames,current)
	end

	return paramsNames
end

function oldfile(funcname::String)
	filename = DATAFOLDER*funcname*".h5"
	if isfile(filename)
		oldname = filename*".old"
		if isfile(oldname)
			oldcount = 1
			while isfile(oldname*"$oldcount")
				oldcount += 1
			end
			oldcount -= 1
			for (i,oldnum) in enumerate(oldcount:-1:1)
				run(`mv -f $oldname$oldnum $oldname$(oldnum+1)`)
				println("Moved file $oldname$oldnum to $oldname$(oldnum+1)")
			end
			run(`mv -f $oldname $(oldname)1`)
			println("Moved file $oldname to $(oldname)1")
		end
		run(`mv -f $filename $filename.old`)
		println("Moved file $filename to $oldname")
	end
end

### Saving and loading macros


macro saving(funcname::String,params...)
	quote
		oldfile($funcname)
		h5write(DATAFOLDER*$funcname*".h5","paramsNames",stringSymbolSeparator($params))
		#=
		for (n,p) in enumerate($(esc(params)))
			h5write(DATAFOLDER*$funcname*".h5",stringSymbolSeparator($params)[n],:p)
			println("Saved :$p as $(stringSymbolSeparator($params)[n])")
		end
		=#
	end
end


function loadnames(funcname::String)
	h5name = DATAFOLDER*funcname*".h5"
	return h5read(h5name,"paramsNames")
end
# =
function loading(funcname::String,addname=false)
	h5name = DATAFOLDER*funcname*".h5"
	paramsNames = h5read(h5name,"paramsNames")
	if addname
		paramsNames .= funcname .* paramsNames
	end

	retArr = Any[]

	for (n,name) in enumerate(paramsNames)
		nstring = """$name = h5read("$h5name","$name")"""
		eval(Meta.parse(nstring))
		push!(retArr,eval(Meta.parse(name)))
	end

	println("Loaded $(length(paramsNames)) variables from function named $funcname with names ($paramsNames)")

	return paramsNames,retArr
end
# =#