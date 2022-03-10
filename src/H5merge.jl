function HDF5.h5write(filename::String,Group::String,Data::Dict)
    for (key,val) in Data
        h5write(filename,Group*"/"*key,val)
    end
end
HDF5.h5write(filename::String,Data::Dict) = HDF5.h5write(filename::String,"",Data::Dict)

function h5Merge(target::String,origin::String,Groups=h5keys(origin))
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