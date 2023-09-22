using Pkg, LinearAlgebra, Printf, SparseArrays, 
MatrixMarket, Distributions, DelimitedFiles, FileIO, JLD2

# this script aims to reprodoce the fixed cases from the paper

@printf("Give the name of the figure, e.g. 4.2 for Figure 4.2:\n");
fig = readline();

if fig == "4.1"
    name = "2_255_1_2_6_3_0_-7";
    k = [255,1];
elseif fig == "4.2"
    name = "2_255_1_1_2_2_0_-3";
    k = [255,1];
elseif fig == "4.3"
    name = "2_228_28_5_6_2_2_-5";
    k = [228,28];
elseif fig == "4.4"
    name = "2_228_28_2_2_2_2_-4";
    k = [228,28];
elseif fig =="sup_1" # figure 1 from summplementary material
    name = "2_228_28_4_7_1_3_-7"
    k = [228,28];
elseif fig =="sup_2" # figure 2 from summplementary material
    name = "2_255_1_2_4_1_0_-3"
    k = [255,1];
else
    @assert(false,"Wrong input for data!");
end

cd("data/");
# Load the corresponding matrix
# and its fixed singular values
A = FileIO.load("A_"*name*".jld2","A");
m,n = size(A);
S = FileIO.load("S_"*name*".jld2","S");

# run the SVD for each precision
# convert to lower precisions
As = convert(Matrix{Float32},A);
Ah = convert(Matrix{Float16},A);

F = svd(A); Sjd = F.S; # double
F = svd( convert(Matrix{Float64},As) ); Sjs = F.S;

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

r  = k[end];
avg_in_s = mean(Sjs[n-r+1:end] .- Sjd[n-r+1:end]);

# print the smallest singular values for each case
# as well as the average increase of the smallest r
# singular values 
@printf("S[n] = %.10f\n",S[end]);
@printf("Sjd[n] = %.10f\n",Sjd[end]);
@printf("Sjs[n] = %.10f\n",Sjs[end]);
@printf("avg(Sigma_{2}^{s}) = %.10f\n", avg_in_s);
if half_op_tt
    F = svd( convert(Matrix{Float64},Ah) ); Sjh = F.S;
    avg_in_h = mean(Sjh[n-r+1:end] .- Sjd[n-r+1:end]);
    @printf("Sjh[n] = %.10f\n",Sjh[end]);
    @printf("avg(Sigma_{2}^{h}) = %.10f\n", avg_in_h);
end

cd("../");