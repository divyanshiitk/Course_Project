function C = connectivity_matrix(Ne,m)
C = zeros(Ne,m);
o = 1:m;
    for i = 1:Ne
    C(i,o) = 1 + (i-1)*(m-1):m+(i-1)*(m-1); %connectivity matrix
    end 
end