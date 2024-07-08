using Pkg
project_path = ENV["JULIA_PROJECT_PATH"]
Pkg.activate(project_path)
Pkg.precompile()

# Pkg.add("CSV")
# Pkg.add("HDF5")
# Pkg.add("FASTX")
# Pkg.add("CodecZlib")
# Pkg.add("IterTools")
# Pkg.add("StatsBase")
# Pkg.add("DataFrames")
# Pkg.add("StringViews")
# Pkg.add("Combinatorics")

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


# fastqpath is the path to a directory containing all fastq files to process
# puckpath is the path to a directory containing all puck files to reference
# Running this script will process all files in both directories and produce one output
# - so be sure that these directories only contain the files for a specific run

# Read the command-line arguments
if length(ARGS) != 2
    error("Usage: julia spatial_count.jl fastqpath puckpath")
end
fastqpath = ARGS[1] 
puckpath = ARGS[2]

println("FASTQS path: "*fastqpath)
@assert isdir(fastqpath) "FASTQ path not found"
@assert !isempty(readdir(fastqpath)) "FASTQ path is empty"

println("Puck path: "*puckpath)
@assert isdir(puckpath) "Puck path not found"
@assert !isempty(readdir(puckpath)) "Puck path is empty"


##### Load the FASTQ data ######################################################

# Load the FASTQ paths
fastqs = readdir(fastqpath, join=true)
R1s = filter(s -> occursin("_R1_", s), fastqs) ; println("R1s: ", basename.(R1s))
R2s = filter(s -> occursin("_R2_", s), fastqs) ; println("R2s: ", basename.(R2s))
@assert length(R1s) == length(R2s) > 0
@assert [replace(R1, "_R1_"=>"", count=1) for R1 in R1s] == [replace(R2, "_R2_"=>"", count=1) for R2 in R2s]

# Validate the FASTQ sequence lengths
function fastq_seq_len(path)
    return(path |> open |> GzipDecompressorStream |> FASTQ.Reader |> first |> FASTQ.sequence |> length)
end
@assert length(unique([fastq_seq_len(R1) for R1 in R1s])) == 1 "WARNING: R1s have different FASTQ sequence lengths - proceed only if you are sure they have the same read structure"
@assert length(unique([fastq_seq_len(R2) for R2 in R2s])) == 1 "WARNING: R2s have different FASTQ sequence lengths - proceed only if you are sure they have the same read structure"

# Scan the first 100k records of both fastqs and identify R2 as the FASTQ with more UP
function learn_switch(R1, R2)
    @assert (fastq_seq_len(R1) >= 32 || fastq_seq_len(R2) >= 32) "ERROR: at least one FASTQ should have a read length of 32 bases"
    i1 = R1 |> open |> GzipDecompressorStream |> FASTQ.Reader
    i2 = R2 |> open |> GzipDecompressorStream |> FASTQ.Reader
    UPseq = String31("TCTTCAGCGTTCCCGAGA")
    s1 = 0 ; s2 = 0
    for (i, record) in enumerate(zip(i1, i2))
        i > 100000 ? break : nothing
        s1 += FASTQ.sequence(record[1])[9:26] == UPseq
        s2 += FASTQ.sequence(record[2])[9:26] == UPseq
    end
    println("R1: ", s1, " R2: ", s2)
    return(s2 >= s1 ? false : true)
end

# Switch R1 and R2 if needed
println("Learning the correct R1 and R2 assignment")
switch = false
res_list = [learn_switch(R1, R2) for (R1, R2) in zip(R1s, R2s)]
@assert all(res_list) || !any(res_list) "ERROR: the R1/R2 read assignment is not consistent"
switch = all(res_list)

if switch == true
    println("Switching R1 and R2")
    temp = R1s
    R1s = R2s
    R2s = temp
end

##### Load the puck data #######################################################

# Load the pucks
pucks = filter(x -> endswith(x, ".csv"), readdir(puckpath, join=true)) ; println("Pucks: ", basename.(pucks))
puckdfs = [rename(CSV.read(puck, DataFrame, header=false, types=[String15, Float64, Float64]), [:sb,:x,:y]) for puck in pucks]

# Validate the pucks
function assert_not_missing(df, col_name)
    @assert all(!ismissing, df[!, col_name]) "Puck column $(col_name) contains missing values, is the data type correct?"
end
for (puck, puckdf) in zip(pucks, puckdfs)
    println("Loaded $(basename(puck)): $(nrow(puckdf)) spatial barcodes found")
    assert_not_missing(puckdf, :sb)
    assert_not_missing(puckdf, :x)
    assert_not_missing(puckdf, :y)
    @assert all(length(s) == 14 for s in puckdf.sb) "Some spatial barcodes in $puck do not have 14bp"
    @assert length(puckdf.sb) == length(Set(puckdf.sb)) "Some spatial barcodes in $puck are repeated"
end

# Concatenate the puck dataframes
for (i, puckdf) in enumerate(puckdfs)
    puckdf[!, :puck_index] = fill(UInt8(i), nrow(puckdf))
end
puckdf = vcat(puckdfs...)
empty!(puckdfs)

# Create the sb_whitelist
sb_whitelist = sort(collect(Set{String15}(puckdf.sb)))
num_lowQbeads = count(str -> count(c -> c == 'N', str) >= 2, sb_whitelist)
if num_lowQbeads > 0
    println("Removed $(num_lowQbeads) bead barcode(s) for having 2+ N bases")
    sb_whitelist = [str for str in sb_whitelist if count(c -> c == 'N', str) < 2]
end
println("Total number of unique spatial barcodes: $(length(sb_whitelist))")
@assert length(sb_whitelist) < 2^32-1 "Must change type of sb_i from U32 to U64"

##### Helper methods ###########################################################

# Returns a set of all strings within a certain hamming distance of the input
function listHDneighbors(str, hd, charlist = ['A','C','G','T','N'])::Set{String}
    res = Set{String}()
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

# Returns a (string, index) tuple for all strings withing 1 HD of the input
# the index stores the 1-indexed position of the mismatched base
function listHD1neighbors(str, charlist = ['A','C','G','T','N'])::Vector{Tuple{String, UInt8}}
    res = Vector{Tuple{String, UInt8}}()
    for i in 1:length(str)
        for char in setdiff(charlist, str[i])
            s = String15(str[1:i-1]*string(char)*str[i+1:end])
            push!(res, (s,i))
        end
    end
    return(res)
end

# Given a string, return a list of all strings made by removing one character
function listDEL1neighbors(str)::Vector{String}
    return([string(str[1:i-1], str[i+1:end]) for i in 1:length(str)])
end

# Given a string with N, return a list of strings with all possible replacements
function expand_N(s::String15, charlist = ['A','C','G','T'])::Vector{String15}
    if !occursin('N', s)
        return [s]
    end
    combins = String15[]
    for nucleotide in charlist
        new_str = String15(replace(s, 'N' => nucleotide, count=1))
        append!(combins, expand_N(new_str))
    end
    return combins
end

# Create spatial barcode matching dictionary (sb1, sb2) -> (index, pos)
# index: ==0 is ambiguous, >0 is the unique index into sb_whitelist
# pos: the position of the fuzzy matched base, 0 for exact/ambig matches
function create_SBtoindex(sb_whitelist)
    SBtoindex = Dict{Tuple{String15, String7}, Tuple{UInt32, Int8}}()

    # Fuzzy matches
    for (i, sb) in enumerate(sb_whitelist)
        numN = count(c -> c == 'N', sb)
        if numN == 0
            # Add Hamming distance = 1 strings
            for (sb_f, ind) in listHD1neighbors(sb)
                sbtuple = (sb_f[1:8], sb_f[9:14])
                haskey(SBtoindex, sbtuple) ? SBtoindex[sbtuple] = (0, 0) : SBtoindex[sbtuple] = (i, ind)
            end
            # Add strings with one deletion
            sb2 = sb[9:14]
            for (ind, sb1) in enumerate(listDEL1neighbors(sb[1:8]))
                sbtuple = (sb1, sb2)
                haskey(SBtoindex, sbtuple) ? SBtoindex[sbtuple] = (0, 0) : SBtoindex[sbtuple] = (i, -ind)
            end
        elseif numN == 1
            ind = findfirst(isequal('N'), sb)
            sbtuple = (sb[1:8], sb[9:14])
            haskey(SBtoindex, sbtuple) ? SBtoindex[sbtuple] = (0, 0) : SBtoindex[sbtuple] = (i, ind)
            for sb_f in expand_N(sb)
                sbtuple = (sb_f[1:8], sb_f[9:14])
                haskey(SBtoindex, sbtuple) ? SBtoindex[sbtuple] = (0, 0) : SBtoindex[sbtuple] = (i, ind)
            end
        end
    end

    # Exact matches
    for (i, sb) in enumerate(sb_whitelist)
        if !occursin('N', sb)
            sbtuple = (sb[1:8], sb[9:14])
            SBtoindex[sbtuple] = (i, 0)
        end
    end
        
    return(SBtoindex)
end

# UMI compressing (between 0x00000000 and 0x00ffffff for a 12bp UMI)
const px = [convert(UInt32, 4^i) for i in 0:(12-1)]
function UMItoindex(UMI::StringView{SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int64}}, true}})::UInt32
    return(dot(px, (codeunits(UMI).>>1).&3))
end
# Convert compressed UMIs back into strings
bases = ['A','C','T','G'] # MUST NOT change this order
function indextoUMI(i::UInt32)::String15
    return(String15(String([bases[(i>>n)&3+1] for n in 0:2:22])))
end

##### Read the FASTQS ##########################################################

function process_fastqs(R1s, R2s, sb_whitelist)
    # Create data structures
    reads = 0
    m = Dict("exact"=>0, "GG"=>0, "none"=>0, "-1X"=>0, "1D-"=>0, "1D-1X"=>0, "-1D"=>0, "-2X"=>0, "umi_N"=>0, "umi_homopolymer"=>0)
    p = Dict("exact"=>0, "HD1"=>0, "HD1ambig"=>0, "none"=>0)
    l = Dict(i=>0 for i in collect(-8:14))
    cb_dictionary = Dict{String31, UInt32}() # cb -> cb_i
    mat = Dict{Tuple{UInt32, UInt32, UInt32}, UInt32}() # (cb_i, umi_i, sb_i) -> reads
    
    print("Creating matching dictionaries... ") ; flush(stdout)
    
    homopolymer_whitelist = Set{String15}()
    for i in 0:2
        for str in [String15(c^12) for c in ["A","C","G","T","N"]]
            union!(homopolymer_whitelist, listHDneighbors(str, i))
        end
    end
    
    UPseq = String31("TCTTCAGCGTTCCCGAGA")
    UPseqHD1 = Set{String31}(listHDneighbors(UPseq, 1))
    UPseqHD2 = Set{String31}(listHDneighbors(UPseq, 2))
    UPseqLD1 = Set{String31}(listDEL1neighbors(UPseq))
    GG_whitelist = Set{String31}(reduce(union, [listHDneighbors("G"^18, i) for i in 0:3]))

    SBtoindex = create_SBtoindex(sb_whitelist)
    
    println("done") ; flush(stdout) ; GC.gc()
    
    print("Reading FASTQs... ") ; flush(stdout)
    
    for fastqpair in zip(R1s, R2s)
        it1 = fastqpair[1] |> open |> GzipDecompressorStream |> FASTQ.Reader;
        it2 = fastqpair[2] |> open |> GzipDecompressorStream |> FASTQ.Reader;
        for record in zip(it1, it2)
            reads += 1

            # Load R1
            cb  = FASTQ.sequence(record[1], 1:16)
            umi = FASTQ.sequence(record[1], 17:28)
            if occursin('N', umi)
                m["umi_N"] += 1
                continue
            elseif in(umi, homopolymer_whitelist)
                m["umi_homopolymer"] += 1
                continue
            end

            # Load R2
            r2 = FASTQ.sequence(record[2], 1:32)
            if r2[9:26] == UPseq # exact match
                sb1=r2[1:8]; sb2=r2[27:32]; m["exact"]+=1
            elseif in(r2[9:26], UPseqHD1) # one UP base flip (could be del at last base - could check cap seq)
                sb1=r2[1:8]; sb2=r2[27:32]; m["-1X"]+=1
            elseif in(r2[9:26], GG_whitelist) # discard
                m["GG"]+=1; continue
            elseif r2[8:25]==UPseq # deletion in sb1, exact UP match
                sb1=r2[1:7]; sb2=r2[26:31]; m["1D-"]+=1
            elseif in(r2[8:25], UPseqHD1) # deletion in sb1, one UP base flip
                sb1=r2[1:7]; sb2=r2[26:31]; m["1D-1X"]+=1
            elseif in(r2[9:25], UPseqLD1) # deletion in UP
                sb1=r2[1:8]; sb2=r2[26:31]; m["-1D"]+=1
            elseif in(r2[9:26], UPseqHD2) # two UP base flips
                sb1=r2[1:8]; sb2=r2[27:32]; m["-2X"]+=1
            else # No detectable UP sequence
                m["none"]+=1; continue
            end

            # match sb -> sb_i
            sb_i, ind = get(SBtoindex, (sb1, sb2), (-1, 0))
            if sb_i > 0 && ind == 0
                p["exact"]+=1 ; l[ind]+=1
            elseif sb_i > 0
                p["HD1"]+=1 ; l[ind]+=1
            elseif sb_i == 0
                p["HD1ambig"]+=1
                continue
            else
                p["none"]+=1
                continue
            end

            # update counts
            cb_i = get!(cb_dictionary, cb, length(cb_dictionary) + 1)
            umi_i = UMItoindex(umi)
            key = (cb_i, umi_i, sb_i)
            mat[key] = get(mat, key, 0) + 1
        end
    end
    
    println("done") ; flush(stdout) ; GC.gc()

    print("Processing results... ") ; flush(stdout)
    
    # Turn matrix from dictionary into dataframe
    df = DataFrame(cb_i = UInt32[], umi_i = UInt32[], sb_i = UInt32[], reads = UInt32[])
    for key in keys(mat)
        value = pop!(mat, key)
        push!(df, (key[1], key[2], key[3], value))
    end

    # Turn cb_dictionary into cb_whitelist
    cb_whitelist = DataFrame(cb = collect(String31, keys(cb_dictionary)), cb_i = collect(UInt32, values(cb_dictionary)))
    sort!(cb_whitelist, :cb_i)
    @assert cb_whitelist.cb_i == 1:size(cb_whitelist, 1)

    println("done") ; flush(stdout) ; GC.gc()

    return(reads, m, p, l, cb_whitelist.cb, df)
end

# reads: # of total reads
# m: UMI and UP filter statistics
# p: spatial barcode matching statistics
# l: spatial barcode fuzzy matching location

reads, m, p, l, cb_whitelist, df = process_fastqs(R1s, R2s, sb_whitelist)
@assert reads == sum(values(m))
@assert m["1D-"] + m["1D-1X"] + m["exact"] + m["-1X"] + m["-1D"] + m["-2X"] == sum(values(p))
@assert p["exact"] + p["HD1"] == sum(df.reads)
@assert p["exact"] == l[0] && sum(values(l))-l[0] == p["HD1"]
l = sort(DataFrame(k = l |> keys |> collect, v = l |> values |> collect), :k)

# Create a downsampling curve
downsampling = UInt32[]
table = countmap(df.reads)
for prob in 0:0.05:1
    s = [length(unique(floor.(sample(0:k*v-1, round(Int,k*v*prob), replace=false)/k))) for (k,v) in zip(keys(table),values(table))]
    append!(downsampling, sum(s))
    GC.gc()
end

##### Save results #############################################################

print("Saving results... ") ; flush(stdout)

h5open("SBcounts.h5", "w") do file
    create_group(file, "lists")
    file["lists/cb_list", compress=9] = cb_whitelist # Vector{String31}
    file["lists/sb_list", compress=9] = sb_whitelist # Vector{String15}
    file["lists/puck_list"] = basename.(pucks)        # Vector{String}

    create_group(file, "matrix")
    file["matrix/cb_index", compress=9] = df.cb_i # Vector{UInt32}
    file["matrix/umi", compress=9] = df.umi_i     # Vector{UInt32}
    file["matrix/sb_index", compress=9] = df.sb_i # Vector{UInt32}
    file["matrix/reads", compress=9] = df.reads   # Vector{UInt32}

    create_group(file, "puck")
    file["puck/sb", compress=9] = puckdf.sb                 # Vector{String15}
    file["puck/x", compress=9] = puckdf.x                   # Vector{Float64}
    file["puck/y", compress=9] = puckdf.y                   # Vector{Float64}
    file["puck/puck_index", compress=9] = puckdf.puck_index # Vector{UInt8}
    
    create_group(file, "metadata")
    
    file["metadata/R1s"] = R1s
    file["metadata/R2s"] = R2s
    file["metadata/switch"] = convert(Int8, switch)
    file["metadata/num_reads"] = reads
    file["metadata/num_lowQbeads"] = num_lowQbeads

    create_group(file, "metadata/UP_matching")
    file["metadata/UP_matching/type"] = keys(m) |> collect
    file["metadata/UP_matching/count"] = values(m) |> collect

    create_group(file, "metadata/SB_matching")
    file["metadata/SB_matching/type"] = keys(p) |> collect
    file["metadata/SB_matching/count"] = values(p) |> collect
    file["metadata/SB_matching/position"] = l.k |> collect
    file["metadata/SB_matching/position_count"] = l.v |> collect
    
    file["metadata/downsampling"] = downsampling # Vector{UInt32}
end;

println("Done") ; flush(stdout) ; GC.gc()

##### Documentation ############################################################

# The workflow for processing the FASTQs is as follows:
#   1) UMI filter: throw out reads that have an N in the umi or a homopolymer umi
#   2) throw out reads that don't have a detectable UP site
#      - UP site does fuzzy matching - allows 1 deletion or up to 2 mismatches
#   3) throw out reads whose spatial barcode (sb) doesn't match the puck whitelist
#      - allow fuzzy matching with either 1 mismatch or 1 deletion
#   4) matrix update
#      - for reads that made it this far, the (cb, umi, sb) key of the dictionary will be incremented by 1
#      - cb_i, sb_i represent an index into the cb_whitelist, sb_whitelist vectors respectively
#      - umi is 2-bit encoded