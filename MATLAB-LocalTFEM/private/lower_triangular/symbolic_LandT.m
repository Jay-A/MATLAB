function [L, T, M] = symbolic_LandT(p);

if (p ~= 3 && p~=4 && p~=5)
    fprintf('Lower-triangular pseudo-mass methods only work for p=3,4,and 5.\n\n');
    return;
end

M = Dubiner_Mass(p);    % Dubiner Mass Matrix for p

dim = (p+1)*(p+2)/2;    % total number of basis functions for T(p)
vm = 1:3;   % vertex mode indices 
sm1 = 1+3.*[1:(p-1)]; sm2 = sm1+1; sm3 = sm2+1; % edge mode indices
im = 3*p+1:dim; % interior mode indices
nim = (p-1)*(p-2)/2; % number of interior modes
sm1m1 = sm1(1:end-1);   % side 1 functions in T(p-1)
sm2m1 = sm2(1:end-1);   % side 2 functions in T(p-1) 
sm3m1 = sm3(1:end-1);   % side 3 functions in T(p-1)
imm1 = im(1:(p-2)*(p-3)/2);   % interior functions in T(p-1)
mu = [vm, sm1m1, sm2m1, sm3m1, imm1]; % all indices in T(p-1)

DofsT = cell((p-1),3); 
DofsL = cell((p-1),3); 

%% Vertex Rows of L and T
% mu is broken into zero constraints to define T & degrees of freedom in L
T0 = sym(eye(dim,dim));
L0 = T0;

muL_A = [1]; muL_B = [2]; muL_C = [3]; % vertex degrees of freedom in L
mu0_A = [2,3, sm1m1, sm2m1, sm3m1, imm1]; % Zero constraints for T ... 
mu0_B = [1,3, sm1m1, sm2m1, sm3m1, imm1]; 
mu0_C = [1,2, sm1m1, sm2m1, sm3m1, imm1]; 
Dofs_A = [sm2, sm3, im]; Dofs_B = [sm1, sm3, im]; Dofs_C = [sm1, sm2, im];

T0(vm(1),Dofs_A) = -(M(mu0_A, Dofs_A))\(M(mu0_A, vm(1)));
T0(vm(2),Dofs_B) = -(M(mu0_B, Dofs_B))\(M(mu0_B, vm(2)));
T0(vm(3),Dofs_C) = -(M(mu0_C, Dofs_C))\(M(mu0_C, vm(3)));

L0(1:3, mu) = T0(1:3, :)*M(:, mu);

T = T0; 
L = L0;

%% Edge Rows of L and T
% mu is broken into zero constraints to define T & degrees of freedom in L
for side_count = 1:(p-1)
    mu0_Ai = [sm1m1(1:side_count-1), sm1m1(side_count+1:end), sm2m1(side_count:end), sm3m1(side_count:end), imm1]; % Zero constraints for T ... 
    mu0_Bi = [sm1m1(side_count:end), sm2m1(1:side_count-1),  sm2m1(side_count+1:end), sm3m1(side_count:end), imm1]; 
    mu0_Ci = [sm1m1(side_count:end), sm2m1(side_count:end), sm3m1(1:side_count-1),  sm3m1(side_count+1:end), imm1];
    
    Dofs_A = [sm1(side_count+1:end), im]; 
    Dofs_B = [sm2(side_count+1:end), im]; 
    Dofs_C = [sm3(side_count+1:end), im];
    
    T(sm1(side_count),Dofs_A) = -(M(mu0_Ai, Dofs_A))\(M(mu0_Ai, sm1(side_count)));
    T(sm2(side_count),Dofs_B) = -(M(mu0_Bi, Dofs_B))\(M(mu0_Bi, sm2(side_count)));
    T(sm3(side_count),Dofs_C) = -(M(mu0_Ci, Dofs_C))\(M(mu0_Ci, sm3(side_count)));
    
    L(sm1(side_count), mu) = T(sm1(side_count),:)*M(:,mu);
    L(sm2(side_count), mu) = T(sm2(side_count),:)*M(:,mu);
    L(sm3(side_count), mu) = T(sm3(side_count),:)*M(:,mu);
    
    if side_count == (p-1)
        muPA = [vm, sm1, sm2m1, sm3m1, imm1];
        muPB = [vm, sm1m1, sm2, sm3m1, imm1];
        muPC = [vm, sm1m1, sm2m1, sm3, imm1];
        L(sm1(side_count), muPA) = T(sm1(side_count),:)*M(:,muPA);
        L(sm2(side_count), muPB) = T(sm2(side_count),:)*M(:,muPB);
        L(sm3(side_count), muPC) = T(sm3(side_count),:)*M(:,muPC);
        clear muPA muPB muPC;
    end
end

L(im, :) = T(im, :)*M;

