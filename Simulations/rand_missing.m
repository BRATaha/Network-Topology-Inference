clc
clear all
close all
graph_v = 0;
N = 80; % Number of nodes
N_av = 5:5:N/2; % Number of anchors
R = 100; % Range
Rn = 25; % Communication range
Nsim = 100; % Number of MC simulations
pc1=1; pc2=5;
per_err = [0.05 0.1 0.3];
bet = 0.5;
tau_v = 20:20:100;
epsln = 0.05;

% Generating a connected graph
connected_var = 0;
while connected_var == 0
Nxy = randi(R,N,2);
A = [(1:N)' Nxy];
E_rand = zeros(N,N);
for i = 1:N
    for j=1:N
        if (i~=j) && (norm(A(i,2:3) - A(j,2:3),2) < Rn)  
            E_rand(i,j) = 1;
        end
    end
end
E_rand_ = E_rand; 
E_rand_(1:N+1:end) = 1;
[~,~,r,~] = dmperm(E_rand_);
if length(r) == 2 
    connected_var = 1;
end
end
Eu = triu(E_rand);
[xe,ye,~] = find(Eu);
M = length(xe);
E = [(1:M)' xe ye];
DG = sparse(xe',ye',ones(1,M),N,N); 
UG = tril(DG + DG');

% Testing accuracy of different inference methods
c_test = 0;
for na_idx = 1:length(N_av)   
N_a = N_av(na_idx);
for sim_idx = 1:Nsim
[N_a sim_idx]    

a = randperm(N,N_a);
d = zeros(N,length(a));
for j = 1:length(a)
for i = 1:N
    [d(i,j),~] = graphshortestpath(UG,i, a(j),'directed',false);
end
end
Gm = 2*ones(N,N) - 2*eye(N);
Gm = triu(Gm);
for i = 1:length(a)
    a1 = a(i);
    d1 = d(:,i);
    Gm = Gen_G0(a1,d1,Gm);
end
% PCA
[iu1,ju1] = find(Gm == 1);
[iu2,ju2] = find(Gm == 2);
iu = [iu1;iu2];
ju = [ju1;ju2];
M0 = length(iu);
kn = []; % known nodes
ukn = []; % unknown nodes
for i = 1:N
    if (sum(ismember([Gm(1:i,i)' Gm(i,i+1:N)],2)) == 0)
        kn = [kn i];
    else
        ukn = [ukn i];
    end
end
Euk = E;
iv = [];
for i = 1:M
    if any(kn == E(i,2)) || any(kn == E(i,3))
        iv = [iv i];
    end
end
Euk(iv,:) = [];
Muk = size(Euk,1);
if (Muk == 0)
    c_test = c_test+1;
    e4(na_idx,sim_idx) = length(find(Gm~=Eu));
    break
end
d1 = d;
d1(kn,:) = [];
[U1,S1,V1] = svd(d1);
P_SVD1 = d1*V1;
PS1 = P_SVD1(:,pc1:pc2);
d_v1 = zeros(M0,1);
for m = 1:M0
    if (sum(ismember(ukn,iu(m))) == 1) && (sum(ismember(ukn,ju(m))) == 1)
        i1 = find(ukn == iu(m)); i2 = find(ukn == ju(m));
        d_v1(m) = norm(PS1(i1,:) - PS1(i2,:),2);
    end
end
dm1 = max(d_v1);
E4 = Gm;
for i = 1:length(ukn)-1
    for j = i+1:length(ukn)
        dij2 = norm(PS1(i,:) - PS1(j,:),2);
        if (dij2 <= dm1) && (E4(ukn(i),ukn(j)) == 2)
            E4(ukn(i),ukn(j)) = 1;
        elseif (E4(ukn(i),ukn(j)) == 2)
            E4(ukn(i),ukn(j)) = 0;
        end
    end
end
e4(na_idx,sim_idx) = length(find(E4~=Eu));
% -------------------------
for i_err = 1:length(per_err)
N_err = randperm(N*N_a,floor(per_err(i_err)*(N*N_a)));
W = ones(N,N_a);
W(N_err) = zeros(1,length(N_err));
d_err = d;
d_err(N_err) = -ones(1,length(N_err));
Gr = 2*ones(N,N) - 2*eye(N);
Gr = triu(Gr);
for i = 1:length(a)
    a1 = a(i);
    d1 = d_err(:,i);
    Gr = Gen_G0(a1,d1,Gr);
end

% PCA
[iu1,ju1] = find(Gr == 1);
[iu2,ju2] = find(Gr == 2);
iu = [iu1;iu2];
ju = [ju1;ju2];
M0 = length(iu);

d_m2 = d;
mu_N = sum(W.*d_m2,1)./(sum(W,1));
for i = 1:length(a)
    for j = 1:length(a)
        sig_n(i,j) = sum(W(:,i).*W(:,j).*(d_m2(:,i)-mu_N(i)*ones(N,1)).*(d_m2(:,j)-mu_N(j)*ones(N,1))/sum(W(:,i).*(sum(W(:,j)))));
    end
end

[U1,S1,V1] = svd(sig_n);
P_SVD1 = (d_m2.*W)*V1;
PS1 = P_SVD1(:,pc1:pc2);
d_v1 = zeros(M0,1);
for m = 1:M
        d_v1(m) = norm(PS1(iu(m),:) - PS1(ju(m),:),2);
end
dm1 = max(d_v1);
E44 = Gr;
for i = 1:N-1
    for j = i+1:N
        dij2 = norm(PS1(i,:) - PS1(j,:),2);
        if (dij2 <= dm1) && (E44(i,j) == 2)
            E44(i,j) = 1;
        elseif (dij2 > dm1) && (E44(i,j) == 2)
            E44(i,j) = 0;
        end
    end
end
e44(na_idx,sim_idx,i_err) = length(find(E44~=Eu));

% matrix comletion -------------------------------------- 

Omeg = ones(N,N_a);
Omeg(N_err) = zeros(1,length(N_err));
X = d.*Omeg;
for id_tau = 1:length(tau_v)
tau = tau_v(id_tau);
Z0 = zeros(N,N_a);
Z = Z0;
[Us,Ss,Vs] = svd(Omeg.*Z);
d_c = Us*(sign(Ss).*max(abs(Ss)-tau*ones(size(Ss)),zeros(size(Ss))))*Vs';
Z = Z + bet*(Omeg.*X - Omeg.*d_c);
while norm(Z-Z0) > epsln
    [tau_v(id_tau) norm(Z-Z0)];
    Z0 = Z;
    [Us,Ss,Vs] = svd(Omeg.*Z);
    d_c = Us*(sign(Ss).*max(abs(Ss)-tau*ones(size(Ss)),zeros(size(Ss))))*Vs';
    Z = Z + bet*(Omeg.*X - Omeg.*d_c);
end
end
% PCA
[U,S,V] = svd(d_c);
P_SVD = d_c*V;
PS = P_SVD(:,pc1:pc2);

[iu1,ju1] = find(Gr == 1);
[iu2,ju2] = find(Gr == 2);
iu = [iu1;iu2];
ju = [ju1;ju2];
M0 = length(iu);
d_v0 = zeros(M0,1);
for m = 1:M0
    d_v0(m) = norm(PS(iu(m),:) - PS(ju(m),:),2); 
end
dm0 = max(d_v0);
% PCA modified
kn = []; % known nodes
ukn = []; % unknown nodes
for i = 1:N
    if (sum(ismember([Gr(1:i,i)' Gr(i,i+1:N)],2)) == 0)
        kn = [kn i];
    else
        ukn = [ukn i];
    end
end
Euk = E;
iv = [];
for i = 1:M
    if any(kn == E(i,2)) || any(kn == E(i,3))
        iv = [iv i];
    end
end
Euk(iv,:) = [];
Muk = size(Euk,1);

if (Muk == 0)
    c_test = c_test+1;
    e45(na_idx,sim_idx,i_err) = length(find(Gr~=Eu));
    break
end
    
d1 = d_c;
d1(kn,:) = [];
[U1,S1,V1] = svd(d1);
P_SVD1 = d1*V1;
PS1 = P_SVD1(:,pc1:pc2);

d_v1 = zeros(M0,1);
for m = 1:M0
    if (sum(ismember(ukn,iu(m))) == 1) && (sum(ismember(ukn,ju(m))) == 1)
        i1 = find(ukn == iu(m)); i2 = find(ukn == ju(m));
        d_v1(m) = norm(PS1(i1,:) - PS1(i2,:),2);
    end
end
dm1 = max(d_v1);
E45 = Gr;
for i = 1:length(ukn)-1
    for j = i+1:length(ukn)
        dij2 = norm(PS1(i,:) - PS1(j,:),2);    
        if (dij2 <= dm1) && (E45(ukn(i),ukn(j)) == 2)
            E45(ukn(i),ukn(j)) = 1;
        else
            E45(ukn(i),ukn(j)) = 0;
        end
    end
end
e45(na_idx,sim_idx,i_err) = length(find(E45~=Eu));

end
end
end
e4 = 100*mean(e4,2)/(N*(N-1)/2);
e44 = 100*mean(e44,2)/(N*(N-1)/2);
e45 = 100*mean(e45,2)/(N*(N-1)/2);

figure 
plot(N_av,e4,'m-*','linewidth',2)
hold on
plot(N_av,e44(:,:,1),'r-*','linewidth',1.5)
plot(N_av,e44(:,:,2),'k-*','linewidth',1.5)
plot(N_av,e44(:,:,3),'b-*','linewidth',1.5)
plot(N_av,e45(:,:,1),'r-s','linewidth',1.5)
plot(N_av,e45(:,:,2),'k-s','linewidth',1.5)
plot(N_av,e45(:,:,3),'b-s','linewidth',1.5)

legend('PCA','IPCA (5% missing)','IPCA (10% missing)','IPCA (30% missing)',...
    'Matrix completion (5% missing)','Matrix completion (10% missing)','Matrix completion (30% missing)')

xlabel('Number of anchors')
ylabel('Percentage of edge erorr (%)')
grid on





