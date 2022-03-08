function HDF5.h5write(filename::String,Group::String,Data::Dict)
    for (key,val) in Data
        h5write(filename,Group*"/"*key,val)
    end
end
HDF5.h5write(filename::String,Data::Dict) = HDF5.h5write(filename::String,"",Data::Dict)

function H5Merge(target::String,origin::String,Groups=h5keys(origin))
    h5open(origin,"r") do f
        for key in Groups
            data = read(f[key])
            try
                h5write(target,key,data)
            catch e
                @warn "Merging of field $key errored with exception :\n $e "
            end
        end
    end
end

function getSourceFilesWith(key::String,Dir::String)
    Dir .* filter!(x -> occursin(key,x), readdir(Dir))
end
Files = getSourceFilesWith("N=32.h5","/storage/niggeni/PMFRG_Results/Triangular/")
for f in Files
    newf = first(split(f,".h5"))*"_l_1.h5"
    H5Merge(newf,f)
end
##
rm.(Files)