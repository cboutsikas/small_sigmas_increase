# Input parameters
########################################################################################
# m: Int, number of rows

# n: Int number of columns (we are interested in the case where m >=n)

# num_cl: Int, number of singular values clusters
#        The script runs for arbitrary num_cl, in the paper we focus on num_cl = 2.

# k: 1-d Int array: contains the number of elements for each custer except 
#    for the last one, where it can be inferred as n - sum(k). This way ensures
#    that sum(k) == n. Notice that length(k) == num_cl - 1.

# g: 1-d array: contains the exponents (with base 10) that control the dintance 
#    between the largest and smallest singular value within the cluster. For instance,
#    assume the clustres \Sigma_{1} \in [10^{-2},10^{3}] and \Sigma_{2} = [10^{-4}],
#    then g = [5,0].

# gi: 1-d array: contains the exponents (with base 10) that control the distance between
#     the smallest and largest singular value for two consequtive clusters. USing again 
#     the example above we would have gi=[2]. Notice that length(gi) = num_cl -1.

########################################################################################

# Output parameters
########################################################################################
# name: A string that describes the parameters of the experiment. It has the form
#       of num_cl_smax_k[1]_k[2]..._k[num_cl]_g[1]_gi[1]..._g[num_cl-1]_gi[num_cl-1]_smin
#       For example if Sigma_{1} \in [10^{-2},10^{3}] and \Sigma_{2} = [10^{-4}],
#       name = "2_3_228_28_5_2_0_-4".

# Att: The generated matrix with fixed S singular values.

# M: A matrix that has columns, M[:,1] -> fixed singular values S,
#    M[:,2] -> singular values computed in double Sjd,
#    M[:,3] -> singular values computed in single Sjs
#    M[:,4] -> singular values computed in half Sjh (if possible),

# avg_in_s: The average increase of the r = k[end] smallest singular values,
#           in single precision.

# avg_in_h: The average increase of the r = k[end] smallest singular values,
#           in half precision (if possible).

# half_op_tt: A binary flag that is true for computed singular values in half,
#             and false if we had overflow/undeflow instacnes. If is false,
#             we do not compute Sjh. 
########################################################################################

using Pkg, LinearAlgebra, Printf, SparseArrays, 
MatrixMarket, Distributions, DelimitedFiles, FileIO, JLD2

include("gen_matrices_k_clusters.jl");

function ComputeAndStore(m::Int,n::Int, s_max::Number, num_cl::Int, k::AbstractArray{Int}, g::AbstractArray, gi::AbstractArray)
    @assert(sum(k)==n,"The number of elements of k do not sum up to n!");
    # fix the name based on the input parameters
    name = string(num_cl);
    for j=1:num_cl
        name = name*"_"*string(k[j]);
    end
    name = name*"_"*string(s_max);
    # for the arrays g, gi
    for j=1:num_cl
        name = name*"_"*string(g[j]);
        if j < num_cl
            name = name*"_"*string(gi[j]);
        end
    end
    # include the sigma min exponent in the name
    s_min = s_max -sum(g) - sum(gi);
    name = name*"_"*string(s_min);

    Att, S = gen_matrices_k_clusters(m,n,num_cl, g, gi, k[1:end-1], s_max);

    # convert to lower precisions
    As = convert(Matrix{Float32},Att);
    Ah = convert(Matrix{Float16},Att);

    half_op_tt = true;
    # check for overflow/underflow incidents
    if any(isnan.(Ah))
        @printf("Nan in tall & thin half matrix! \n");
        half_op_tt = false;
    end

    if any(isinf.(Ah))
        @printf("Inf in tall & thin half matrix! \n");
        half_op_tt = false;
    end

    if length(findall(==(abs(0.0)),Ah)) == m*n
        @printf("All zeros in tall & thin half matrix!\n");
        half_op_tt = false;
    end

    l = 4;
    M = zeros(length(S),l); M[:,1] = copy(S);
    F = svd(Att); Sjd = F.S; # double
    M[:,2] = copy(Sjd);
    # single
    F = svd( convert(Matrix{Float64},As) ); Sjs = F.S;
    M[:,3] = copy(Sjs);
    if half_op_tt
        # half
        F = svd( convert(Matrix{Float64},Ah) ); Sjh = F.S;
        M[:,4] = copy(Sjh);
    else
        l = 3;
    end
    
    # save the results as Julia's variables
    # and M for python plotting
    cd("data/");
    # calculate the avererage increase for single and half precision
    r  = k[end];
    avg_in_s = mean(Sjs[n-r+1:end] .- Sjd[n-r+1:end]);
    mmwrite("M_"*name*".mtx", convert(SparseMatrixCSC, M[:,1:l]));
    # write also into .txts the avererage increase 
    # of the r smallest singular values
    writedlm("avg_in_s"*name*".txt",mean(avg_in_s));
    FileIO.save("A_"*name*".jld2","A",Att);
    if half_op_tt
        FileIO.save("Sjd_"*name*".jld2","Sjd",Sjd);
        FileIO.save("Sjs_"*name*".jld2","Sjs",Sjs);
        FileIO.save("Sjh_"*name*".jld2","Sjh",Sjh);
        FileIO.save("S_"*name*".jld2","S",S);
        avg_in_h = mean(Sjh[n-r+1:end] .- Sjd[n-r+1:end]);
        writedlm("avg_in_h"*name*".txt",mean(avg_in_h));
    else
        FileIO.save("Sjd_"*name*".jld2","Sjd",Sjd);
        FileIO.save("Sjs_"*name*".jld2","Sjs",Sjs);
        FileIO.save("S_"*name*".jld2","S",S);
        # Sjh = zeros(2,1);
    end
    # write into a .txt the smallest singular values
    # for "exact" S, double Sjd, single Sjs, and half Sjh
    writedlm("min_sigmas"*name*".txt",M[end,:]);
    cd("../");

    return name,Att,M,avg_in_s, avg_in_h,half_op_tt;
    
end

# for instance if you want to generate a matrix A \in Real^{4096 \times 256}
# with fixed singular values (S) as
# \Sigma_{1} \in [10^{-2},10^{3}], 
# \Sigma_{2} \in [10^{-4},10^{-5}]
# where \Sigma_{1} containts the 228 largest singular values
# and \Sigma_{2} contains the 28 smallest ones.

m = 4096; n = 256; s_max = 3; num_cl = 2;
g = [5,1]; gi = [2]; k = [228,28];
name,Att,M,avg_in_s,avg_in_h,half_op_tt = ComputeAndStore(m,n, s_max, num_cl, k, g, gi);

# print the smallest singular values for each case
# as well as the average increase of the smallest r
# singular values

@printf("S[n] = %.10f\n",M[end,1]);
@printf("Sjd[n] = %.10f\n",M[end,2]);
@printf("Sjs[n] = %.10f\n",M[end,3]);
@printf("avg(Sigma_{2}^{s}) = %.10f\n", avg_in_s);
if half_op_tt
    @printf("Sjh[n] = %.10f\n",M[end,4]);
    @printf("avg(Sigma_{2}^{h}) = %.10f\n", avg_in_h);
end