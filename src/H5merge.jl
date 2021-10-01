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
Files = getSourceFilesWith("(v_","/storage/niggeni/PMFRG_Results/Octochlore_NLen=10/")
for f in Files
    H5Merge("/storage/niggeni/PMFRG_Results/Octochlore_NLen=10/Octochlore_NLen=10_N=32_a_OneLoop.h5",f)
end
##
rm.(Files)