using Pkg
project_path = ENV["JULIA_PROJECT_PATH"]
Pkg.activate(project_path)
Pkg.precompile()

using CSV
using HDF5
using FASTX
using CodecZlib
using IterTools: product
using StatsBase: countmap, sample
using DataFrames
using StringViews
using LinearAlgebra: dot
using Combinatorics: combinations


# Read the command-line arguments
if length(ARGS) != 2
    println("Usage: julia readfastq.jl fastqpath puckpath")
    @assert false
end
fastqpath = ARGS[1] ; println("FASTQS path: "*fastqpath)
puckpath = ARGS[2] ; println("Puck path: "*puckpath)
CBdictpath = ENV["CB_DIC_PATH"]

@assert isdir(fastqpath)
@assert isdir(puckpath)
@assert isfile(CBdictpath)


################################################################################
##### Load the data ############################################################
################################################################################

# Load the FASTQ paths
fastqs = readdir(fastqpath,join=true)
R1s = filter(s -> occursin("_R1_", s), fastqs) ; println("R1s: ", basename.(R1s))
R2s = filter(s -> occursin("_R2_", s), fastqs) ; println("R2s: ", basename.(R2s))
@assert length(R1s) == length(R2s) > 0
@assert [replace(R1, "_R1_"=>"", count=1) for R1 in R1s] == [replace(R2, "_R2_"=>"", count=1) for R2 in R2s]

# load the pucks
pucks = readdir(puckpath,join=true) ; println("Pucks: ",basename.(pucks))
puckdfs = [rename(CSV.read(puck,DataFrame,header=false), [:sb,:x,:y]) for puck in pucks]
SBwhitelist = Set{String15}()
for (puck,puckdf) in zip(pucks,puckdfs)
    println("Loaded $(basename(puck)): $(nrow(puckdf)) spatial barcodes found")
    # @assert all(length(s) == 14 for s in puckdf.sb) "Some SB in $puck do not have 14bp"
    @assert length(puckdf.sb) == length(Set(puckdf.sb)) "Some SB in $puck are repeated"
    for sb in puckdf.sb
        push!(SBwhitelist,sb)
    end
end
SBwhitelist = collect(SBwhitelist)
println("Total SB: $(length(SBwhitelist))")
@assert length(SBwhitelist) < 2^32-1 "must change sb_i from UInt32 to UInt64"

################################################################################
##### Learn R1 and R2 ##########################################################
################################################################################
charlist = ['A','C','G','T','N']
function listHDneighbors(str, hd)
    res = Set()
    for inds in combinations(1:length(str), hd)
        chars = [str[i] for i in inds]
        pools = [setdiff(charlist, [char]) for char in chars]
        prods = product(pools...)
        for prod in prods
            s = str
            for (i, c) in zip(inds, prod)
                s = s[1:i-1]*string(c)*s[i+1:end]
            end
            push!(res,s)
        end
    end
    return(res)
end
function list1delneighbors(str)
    return([string(str[1:i-1], str[i+1:end]) for i in 1:length(str)])
end
function expand_N(s::String15)::Vector{String15}
    if !occursin('N', s)
        return [s]
    end
    combins = String15[]
    for nucleotide in ['A','C','G','T']
        new_str = String15(replace(s, 'N' => nucleotide, count=1))
        append!(combins, expand_N(new_str))
    end
    return combins
end

# UP matching sets
UPseq = String31("TCTTCAGCGTTCCCGAGA")
UPseqHD1 = Set{String31}(listHDneighbors(UPseq, 1))
UPseqHD2 = Set{String31}(listHDneighbors(UPseq, 2))
UPseqLD1 = Set{String31}(list1delneighbors(UPseq))

# GG matching set
GGseq = "GGGGGGGGGGGGGGGGGG"
GGset = Set{String31}(reduce(union, [listHDneighbors(GGseq,i) for i in 0:3]))

# Switch R1 and R2 if needed
switch = false
R1len = R1s[1] |> open |> GzipDecompressorStream |> FASTQ.Reader |> first |> FASTQ.sequence |> length
R2len = R2s[1] |> open |> GzipDecompressorStream |> FASTQ.Reader |> first |> FASTQ.sequence |> length
if (R1len < 32) & (R2len < 32)
    error("Unexpected read structure: one fastq should have at least 32 characters")
elseif (R1len < 32) & (R2len >= 32)
    switch = false
elseif (R1len >= 32) & (R2len < 32)
    switch = true
else
    println("Learning the correct R1 and R2 assignment")
    s1 = 0 ; s2 = 0
    i1 = R1s[1] |> open |> GzipDecompressorStream |> FASTQ.Reader
    i2 = R2s[1] |> open |> GzipDecompressorStream |> FASTQ.Reader
    for (i,record) in enumerate(zip(i1, i2))
        i > 100000 ? break : nothing
        global s1 += FASTQ.sequence(record[1])[9:26] == UPseq
        global s2 += FASTQ.sequence(record[2])[9:26] == UPseq
    end
    println("R1: ", s1, " R2: ", s2)
    switch = s2 > s1 ? false : true
end

if switch == true
    println("Switching R1 and R2")
    temp = R1s
    R1s = R2s
    R2s = temp
end

################################################################################
##### Create Whitelists ########################################################
################################################################################

print("Creating matching dictionaries... ")

# store exact SB matches
SBtoindex = Dict{Tuple{String15,String15},UInt32}()
for (i,sb) in enumerate(SBwhitelist)
    if !occursin('N',sb)
        SBtoindex[(sb[1:8],sb[9:14])] = i
    end
end

# store fuzzy SB matches (==0 is ambiguous, >0 is unique, -1 is DNE)
SBtoindexHD1 = Dict{Tuple{String15,String15},UInt32}()
for (i,sb) in enumerate(SBwhitelist)
    if !occursin('N',sb)
        for sb_f in listHDneighbors(sb,1)
            sbtuple = (sb_f[1:8],sb_f[9:14]) # why 9:14
            get(SBtoindexHD1, sbtuple, i) == i ? SBtoindexHD1[sbtuple] = i : SBtoindexHD1[sbtuple] = 0
        end
        sb2 = sb[9:14]
        for sb1 in list1delneighbors(sb[1:8])
            sbtuple = (sb1,sb2)
            get(SBtoindexHD1, sbtuple, i) == i ? SBtoindexHD1[sbtuple] = i : SBtoindexHD1[sbtuple] = 0
        end
    elseif count(c -> c == 'N', sb) <= 2
        for sb_f in expand_N(sb)
            sbtuple = (sb_f[1:8],sb_f[9:14])
            get(SBtoindexHD1, sbtuple, i) == i ? SBtoindexHD1[sbtuple] = i : SBtoindexHD1[sbtuple] = 0
        end
    else
        println("Too many N: ",sb)
    end
end
println("Done")

# UMI compressing (between 0x00000000 and 0x00ffffff)
p12 = [convert(UInt32,4^i) for i in 0:11]
function UMItoindex(UMI::StringView{SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int64}}, true}})::UInt32
    return(dot(p12, (codeunits(UMI).>>1).&3))
end
# bases = ['A','C','T','G']
# function indextoUMI(i::UInt32)::String15
#     return(String15(String([bases[(i>>n)&3+1] for n in 0:2:22])))
# end

# Pass in R2, get back an index into SBwhitelist
function R2process(r2::StringView{SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int64}}, true}})::Int64
    if r2[9:26]==UPseq # exact match
        sb1=r2[1:8]; sb2=r2[27:32]; m["-"]+=1
    elseif in(r2[9:26],UPseqHD1) # one base flip (might be del at last base - could check cap seq)
        sb1=r2[1:8]; sb2=r2[27:32]; m["-1X"]+=1
    elseif in(r2[9:26], GGset) # discard
        m["GG"]+=1; return(-1)
    elseif r2[8:25]==UPseq # deletion in sb1, exact match
        sb1=r2[1:7]; sb2=r2[26:31]; m["1D-"]+=1
    elseif in(r2[8:25],UPseqHD1) # deletion in sb1, one base flip
        sb1=r2[1:7]; sb2=r2[26:31]; m["1D-1X"]+=1
    elseif in(r2[9:25],UPseqLD1) # deletion in UP
        sb1=r2[1:8]; sb2=r2[26:31]; m["-1D"]+=1
    elseif in(r2[9:26],UPseqHD2) # two base flips
        sb1=r2[1:8]; sb2=r2[27:32]; m["-2X"]+=1
    else # No detectable UP sequence
        m["none"]+=1; return(-1)
    end

    key = (sb1, sb2)
    res = get(SBtoindex, key, -1)
    if res > 0
        p["exact"]+=1
        return(res)
    end
    res = get(SBtoindexHD1, key, -1)
    if res > 0
        p["HD1"]+=1
    elseif res == 0
        p["HD1ambig"]+=1
    else
        p["none"]+=1
    end
    return(res)
end

################################################################################
##### Read the FASTQS ##########################################################
################################################################################

print("Reading FASTQS... ")

reads = 0
m = Dict("-"=>0,"-1X"=>0,"GG"=>0,"1D-"=>0,"1D-1X"=>0,"-1D"=>0,"-2X"=>0,"none"=>0,"R1lowQ"=>0)
p = Dict("exact"=>0,"HD1"=>0,"HD1ambig"=>0,"none"=>0)
CBwhitelist = Dict{String31,UInt32}() ; CBnum = 0
mat = Dict{Tuple{UInt32, UInt32, UInt32},UInt32}()
for fastqpair in zip(R1s,R2s)
    it1 = fastqpair[1] |> open |> GzipDecompressorStream |> FASTQ.Reader;
    it2 = fastqpair[2] |> open |> GzipDecompressorStream |> FASTQ.Reader;
    for record in zip(it1, it2)
        r1 = FASTQ.sequence(record[1])
        r2 = FASTQ.sequence(record[2], 1:32)
        global reads += 1

        occursin('N', r1[17:28]) ? (m["R1lowQ"] += 1 ; continue) : nothing

        sb_i = R2process(r2)
        sb_i>0 ? nothing : continue
        
        cb_i = get(CBwhitelist, r1[1:16], 0)
        if cb_i == 0
            global CBnum += 1
            CBwhitelist[r1[1:16]] = CBnum
            cb_i = CBnum
        end
        
        umi_i = UMItoindex(r1[17:28])
        key = (cb_i,umi_i,sb_i)
        mat[key] = get(mat,key,0) + 1
    end
end
print("D4")
@assert reads == sum(values(m))
@assert m["1D-"]+m["1D-1X"]+m["-"]+m["-1X"]+m["-1D"]+m["-2X"] == sum(values(p))
@assert p["exact"]+p["HD1"] == sum(values(mat))
empty!(SBtoindex)
empty!(SBtoindexHD1)
GC.gc()
println("Done")

################################################################################
##### Process the results ######################################################
################################################################################

print("Processing the results... ")

# Turn matrix into dataframe
df = DataFrame(cb = UInt32[], umi = UInt32[], sb = UInt32[], reads = UInt32[])
for (key,value) in zip(keys(mat),values(mat))
     push!(df, (key[1],key[2],key[3],value))
end
empty!(mat)
GC.gc()

# Turn CB whitelist into a dataframe and sort
CBwhitelist_df = DataFrame(cb = collect(String31,keys(CBwhitelist)), cb_i = collect(UInt32, values(CBwhitelist)))
empty!(CBwhitelist)
sort!(CBwhitelist_df, :cb_i)
GC.gc()

# Remap the cell barcodes
CBdf = CSV.File(CBdictpath, delim='\t', header=false) |> DataFrame
@assert size(CBdf) == (6794880, 2)
CBremap = Dict(CBdf[!, 2] .=> CBdf[!, 1])
empty!(CBdf)
CBwhitelist_df.cb_remap = [get(CBremap, cb, String31("")) for cb in CBwhitelist_df.cb]
empty!(CBremap)
GC.gc()

# Concatenate the puck dataframes
for (i,puckdf) in enumerate(puckdfs)
    puckdf[!, :puck_index] = fill(UInt8(i), nrow(puckdf))
end
puckdf = vcat(puckdfs...)
empty!(puckdfs)
GC.gc()

# Create a downsampling curve
downsampling = UInt32[]
table = countmap(df.reads)
for prob in 0:0.05:1
    s = [length(unique(floor.(sample(0:k*v-1, round(Int,k*v*prob), replace=false)/k))) for (k,v) in zip(keys(table),values(table))]
    append!(downsampling,sum(s))
end

println("Done")

################################################################################
##### Save results #############################################################
################################################################################

print("Saving Results... ")

h5open("SBcounts.h5", "w") do file
    create_group(file, "lists")
    file["lists/cb_list", compress=9] = CBwhitelist_df.cb
    file["lists/cb_list_remap", compress=9] = CBwhitelist_df.cb_remap
    file["lists/sb_list", compress=9] = SBwhitelist

    create_group(file, "matrix")
    file["matrix/cb_index", compress=9] = df.cb
    file["matrix/umi", compress=9] = df.umi
    file["matrix/sb_index", compress=9] = df.sb
    file["matrix/reads", compress=9] = df.reads

    create_group(file, "puck")
    file["puck/sb", compress=9] = puckdf.sb
    file["puck/x", compress=9] = puckdf.x
    file["puck/y", compress=9] = puckdf.y
    file["puck/puck_index", compress=9] = puckdf.puck_index
    file["puck/puck_list"] = basename.(pucks)

    create_group(file, "metadata")
    
    file["metadata/R1s"] = R1s
    file["metadata/R2s"] = R2s
    file["metadata/switch"] = convert(Int8, switch)
    file["metadata/num_reads"] = reads

    create_group(file, "metadata/UP_matching")
    file["metadata/UP_matching/type"] = keys(m) |> collect
    file["metadata/UP_matching/count"] = values(m) |> collect

    create_group(file, "metadata/SB_matching")
    file["metadata/SB_matching/type"] = keys(p) |> collect
    file["metadata/SB_matching/count"] = values(p) |> collect
    
    file["metadata/downsampling"] = downsampling
end;

println("Done")

################################################################################
##### Documentation ############################################################
################################################################################

# The workflow for processing the FASTQs is as follows:
#   1) throw out reads that have an N in the umi
#   2) throw out reads that don't have a detectable UP site
#      - UP site does fuzzy matching - allows 1 deletion or up to 2 mismatches
#      - this is implemented by a few Set() objects made in advance that contain all possible fuzzy matches
#   3) throw out reads whose spatial barcode (sb) doesn't match the puck whitelist
#      - allow fuzzy matching with either 1 mismatch or 1 deletion
#   4) dictionary update
#      - for reads that made it this far, the (cb,umi,sb) key of the dictionary will be incremented by 1
#      - umi is 2-bit encoded
#      - cb,sb is stored as an index into the CBwhitelist,SBwhitelist vector
# The final result of the FASTQ reading step is a dictionary, where the keys are encoded (cb, umi, sb) and the values are the number of reads
