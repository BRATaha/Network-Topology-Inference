clc
clear all
close all
% Generating a network from a grid
N = 0;
for i = 0:10
    N = N+1; A(N,:) = [N i 0];
    N = N+1; A(N,:) = [N i 10];
    if ((i>=1) && (i<=9)) N = N+1; A(N,:) = [N 0 i];end
    if ((i>=1) && (i<=9)) N = N+1; A(N,:) = [N 10 i];end
    if ((i>=1) && (i<=9)) N = N+1; A(N,:) = [N 2 i];end
    if ((i>=1) && (i<=5)) N = N+1; A(N,:) = [N 6 i];end
    if ((i>=7) && (i<=9)) N = N+1; A(N,:) = [N i 5];end
end

% Generating a set of anchors
% Ba =  [5 0;
%        6 1;
%        7 5;
%        6 3;
%        8 0;
%        10 3;
%        10 0];

Ba =  [5 0;
       6 1;
       7 5;
       6 3;
       6 5; %8 0;
       10 4;
       10 6];      
    

% Ba =  [0 0;
%        0 10;
%        10 10;
%        10 0;
%        0 5;
%        5 0;
%        10 5;
%        5 10];
   
N_anch = size(Ba,1);   

for n = 1:N
    for m = 1:N_anch
    if A(n,[2 3]) == Ba(m,:)
        B(m,:) = [n Ba(m,:)];
    end
    end
end

% Generating network edges
M = 0;
for i = 1:N
    for j = 1:N
        if norm(A(i,2:3)-A(j,2:3),2) == 1
            M = M+1;
            E(M,:) = [M i j];
        end
    end
end



pc1 = 2;
pc2 = 3;

% Distances between anchors and nodes  
for n = 1:N
    for m = 1:N_anch
        [P(n,m), ~] = dijkstra(A, E, n, B(m,1));
    end
end

Q = (eye(N)-(ones(N,1)*ones(N,1)'/N));
Pz = Q*P;
P0 = (ones(N,1)*ones(N,1)'/N)*P;

% Principal Component Analysis
[U,S,V] = svd(P);
P_SVD = P*V;

x_c = P_SVD(:,2);
y_c = P_SVD(:,3);

x_c2 = P_SVD(:,1);
y_c2 = P_SVD(:,2);

[Uz,Sz,Vz] = svd(Pz);
P_SVDz = Pz*Vz;

[U0,S0,V0] = svd(P0);

x_cz = P_SVDz(:,2);
y_cz = P_SVDz(:,3);

x_c2z = P_SVDz(:,1);
y_c2z = P_SVDz(:,2);


figure
subplot(2,2,1);
plot(A(:,2), A(:,3),'k.','MarkerSize',15);
hold on;
ylim([-1 11]); xlim([-1 11]);
plot(B(:,2), B(:,3),'go','MarkerSize',7,'linewidth',3);
for s = 1:M
    plot(A(E(s,2:3)',2),A(E(s,2:3)',3),'k','linewidth',1.5);
end
xlabel('(a) Anchored network','Interpreter','latex')
subplot(2,2,3);
scatter(x_c,y_c,20,'filled','MarkerFaceColor','r')
xlabel('(c) $2^{nd}$ and $3^{rd}$ PC of $\bf{H}$','Interpreter','latex')

subplot(2,2,4);
scatter(x_c2z,-y_c2z,20,'filled','MarkerFaceColor','r')
xlabel('(d) $1^{st}$ and $2^{nd}$ PC of $\tilde{\bf{H}}$','Interpreter','latex' )

subplot(2,2,2);
DS = diag(S);
DSz = diag(Sz);
% DSz2 = diag(Sz2);
plot(DS,'-*b','linewidth',2)
hold on
plot([DSz(end);DSz(1:end-1)],'-.*r','linewidth',2)
l=legend('Singular values of $\bf{H}$','Singular values of $\tilde{\bf{H}}$');
set(l, 'Interpreter', 'latex')
xlabel('(b) Singular values','Interpreter','latex')






