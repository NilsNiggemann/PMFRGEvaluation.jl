function H5Merge(target::String,origin::String)
    h5open(origin,"r") do f
        for i in keys(f)
            for j in keys(f[i])
                keystr = "$i/$j"
                h5write(target,keystr,Array(f[keystr]))
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