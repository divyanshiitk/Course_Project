function [K,M] = assembly(k,mass,m,Ne,N,C)

M = zeros(N);
K = zeros(N);

for e = 1:Ne
    for p = 1:m
        for q = 1:m
            r = C(e,p);
            s = C(e,q);
         K(r,s) = K(r,s) + k(p,q);
         M(r,s) = M(r,s) + mass(p,q);
        end 
    end
end
end