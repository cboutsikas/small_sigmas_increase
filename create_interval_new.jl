function create_interval_new(l::Int,s1::Int,s2::Int)
    @assert(l>0, "Number of elements should be positive!");

    if l == 1 # in this case only s1 matters
        return 10.0^s1;
    end
    s = zeros(l); s[1] = 10.0^s2; s[l] = 10.0^(s1);
    left_lim = 0.1; right_lim = 1; d = Uniform();
    if s1 == s2
        s1 -= 1;
        left_lim = 1;
        d = Normal();
    end

    y1 = rand(s1+1:s2,l-2); y1 = broadcast(^,10.0,y1);
    y2 = rand(Truncated(d,left_lim,right_lim),l-2);
    s[2:l-1] = sort((y1.*y2),rev=true);

    if norm(sort(s,rev=true) - s) != 0
        @printf("Something is messed up with the limits!\n");
    end
    
    return s
end