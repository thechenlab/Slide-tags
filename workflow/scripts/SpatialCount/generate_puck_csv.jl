using Pkg
project_path = ENV["JULIA_PROJECT_PATH"]
Pkg.activate(project_path)
Pkg.precompile()

using DataFrames
using CSV

inpath = ENV["LIB_PUCK_IN"]
outpath = ENV["LIB_PUCK_PATH"]

indirs = readdir(inpath)
outdirs = [replace(file, ".csv" => "") for file in readdir(outpath)]

indirs = filter(indir -> isfile(joinpath(inpath,indir,"BeadBarcodes.txt")) && isfile(joinpath(inpath,indir,"BeadLocations.txt")), indirs)
println("$(length(indirs)) total puck libraries")
indirs = filter(x -> !(x in outdirs), indirs)
println("$(length(indirs)) new libraries")

if length(indirs) == 0
    println("nothing to do, exiting...")
    exit(0)
end

for dir in indirs
    GC.gc(true)

    sbs = CSV.File(joinpath(inpath,dir,"BeadBarcodes.txt"), header=false, delim=',')
    sbs = [join(row) for row in sbs]
    file = open(joinpath(inpath,dir,"BeadLocations.txt"), "r")
    x = split(readline(file), ",")
    y = split(readline(file), ",")
    close(file)

    @assert length(x)==length(y)==length(sbs)

    df = DataFrame(sb=sbs, x=x, y=y)
    CSV.write(joinpath(outpath,"$(dir).csv"), df, header=false)
end

println("Done")
exit(0)
