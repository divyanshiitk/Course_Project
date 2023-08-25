clc 
clear
%% Pre Processor

Ne = 15;
m = input("Enter no of nodes per element "); % number of nodes per element
Nn = Ne*(m-1) + 1; % total no of nodes
Nd = 1; %no of dof per node
N = Nd*Nn; %total no of dof
C = connectivity_matrix(Ne,m); %defining connectivity matrix
L = 20; %length of each element    
Le = L./Ne;
Nu = 1; % Number of nodes at which u is specified
Np = 0; % Number of nodes at which p is specified
X = linspace(0,L,Nn); %Global Coordinate Vector
x= linspace(0,Le,m); %Local Coordinate Vector
Bu = [1 0]; 
Bp = [0 0];
A = @(y) y./L + 1;
E = 1;
P = @(y) y./L + 1;

%% Processor

%define Gaussian Points, zhi, weights

[Ng,Ng1,zhi,zhi1,weights,w1] = Gaussian(m);

% Shape Function

[Ni,Bi] = Shape_Functions(m,zhi);

% defining {k} and {m} for a single element
k = zeros(m,m);
if m == 2
    for i = 1:Ng
    k = k + weights(i).*((A(((x(1)+x(end))/2)+((Le/2).*zhi(i))))./Le).*((Bi(:,1))*((Bi(:,1)).'));
    end
else
    for i = 1:Ng
        k = k + weights(i).*((A(((x(1)+x(end))/2)+((Le/2).*zhi(i))))./Le).*((Bi(:,i))*((Bi(:,i)).'));
    end
end

[Ni1,Bi1] = Shape_Functions(m,zhi1);
mass = zeros(m,m);
for i = 1:Ng1
    mass = mass + w1(i).*(A(((x(1)+x(end))/2)+((Le/2).*zhi1(i)))).*(P((x(1)+x(end))/2)+((Le/2).*zhi1(i))).*(Le/2).*((Ni1(:,i))*((Ni1(:,i)).'));
end


% Assemble of [K] and [M] for the domain 
[K,M] = assembly(k,mass,m,Ne,N,C);  

%% Post Processer

% Applying Boundary Conditions

K(1,:) = [];
K(:,1) = [];
M(1,:) = [];
M(:,1) = [];
X(1) = [];

[V,D] = eig(K,M);

Omega = real(sqrt(D));


Answer = table(X',Omega);
Answer.Properties.VariableNames = ["Global coordinate Vector","Omega"];

for i = 1:4
hold on
plot(X,V(:,i))
hold off
end



