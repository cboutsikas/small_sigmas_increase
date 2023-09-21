
include("create_interval_new.jl");
function gen_matrices_k_clusters(m::Int,n::Int,i::Int, g::AbstractArray, gi::AbstractArray, k::AbstractArray, s1::Int)
    #=
    Input arguments:
    - m: Number of rows
    - n: Number of columns (n<=m)
    - i: Number of clusters (i>0)
    - g: Array with gaps within the clusters (g >= 0)
    - gi: Array with gaps between the clusters (gi > 0)
    - k: Array with number of elements for each cluster..
    - k: .. The last one is ommited
    - s1: The exponent of the dominant singular value
    =#

    @assert(i>0,"Number of clusters should be positive integer!");
    @assert(length(g[g.<0])==0,"Gap within a cluster cannot be negative!");
    @assert(length(g)==i,"Wrong length of array within the clusters!");
    @assert(length(gi[gi.<=0])==0,"Gap between clusters cannot be zero or negative!");
    @assert(length(gi)== i -1, "Wrong length of array between the clusters!");
    @assert(length(k[k.<=0])==0, "Number of elements for each cluster should be positive!");
    @assert(length(k) == i-1, "Wrong length of array with number of elements for each cluster!");
    @assert(sum(k) < n, "Sum of k should be less than n!");

    Att = randn(m,n); Ftt = svd(Att); S = Ftt.S;
    # change the k as desired
    k_new = Vector{Int64}(zeros(length(k)+1)); k_new[1:length(k)] = copy(k);
    k_new[end] = n - sum(k);
    f1 = findall(x->(x==1),k_new); f2 = findall(x->(x==0),vec(g));
    for j in eachindex(f1)
        @assert(in(f1[j],f2),"Some cluster has size of 1 and the corresponding gap is not zero!");
    end
    # for each cluster
    # initiliaze where counting starts
    t1 = 1;
    # initialize for gaps
    c = s1;
    for j=1:i
        c -= g[j]; 
        S[t1:sum(k_new[1:j])] .= create_interval_new(k_new[j],c,c+g[j]);
        t1 += k_new[j]; 
        if j != i
            c -= gi[j];
        end 
    end
    Att = Ftt.U*diagm(S)*Ftt.Vt; 

    return Att, S;
end