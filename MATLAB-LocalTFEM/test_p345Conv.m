clc; close all; clear all;

addpath('./private/dubiner/', './private/quadrature/', './private/mesh/', ...
    './private/assembly/', './private/local/');

% set function to approximate F0 = @(x,y) cos(2.*pi.*y) + sin(2.*pi.*y);
F0 = @(x,y) (x-1).^2.*(y-1).^2.*(x+y).^2.*(x+1).^2.*(y+1).^2; %cos(2.*pi.*y) + sin(2.*pi.*y);

% Define p for basis degree
p = 5;

% Define L and T for p
% [Lj,Tj] = Make_LT_Nov19(p);
eval(sprintf("load('./private/libraries/lib_LT_Ts/L%d.mat');", p));
eval(sprintf("load('./private/libraries/lib_LT_Ts/T%d.mat');", p));
% T = Tg; L = Lg;
Tj = T; Mj = L;

% set local  modified dubiner functions and mass
[Phi, EvalPhi, Modes_Xi, Modes_Eta] = Dubiner_Modes(p);
% [M] = Dubiner_Mass(p);
eval(sprintf("load('./private/libraries/lib_Dub_Mass/D%d.mat');", p));

% Number of mesh refinements to loop through
m_first = 0;  % initial mesh coarseness
m_final = 2;  % final mesh coarseness
m_range = m_final - m_first + 1; % number of mesh values

% initialize the L2_Error Arrays
L2_Errors = zeros(m_range+1, 3);
L2_Errors(1,1) = 999;
L2_Errors(1,3) = 999;
L2_Errors(2:end, 1) = m_first:m_final;
L2_Errors(1, 2) = p;


syms xi eta;
[lgl_nodes,lgl_weights] = lglnodes(15);
[XI, ETA] = meshgrid(lgl_nodes, lgl_nodes);
WTS = lgl_weights*lgl_weights';
S = ETA; 
R = (1/2).*(XI+1).*(1-ETA)-1;

i_cnt = 0;
for n = m_first:m_final
    i_cnt = i_cnt+1;
    
    [Nodes, Elements, Areas, Edges, Boundary, Num_Elements, Num_Edges, Num_Nodes] = structured_mesh(n);
    L2_Errors(i_cnt+1,3) = Nodes(2,1)-Nodes(1,1);
    ErrorsD_k = zeros(2^(2*n+1),1);

    [map, sgn] = Create_Map(p, Edges, Nodes, Elements);
    Num_Global_DoFs = Num_Nodes + (p-1)*Num_Edges + (1/2)*(p-1)*(p-2)*Num_Elements;

    %%%  Create Global Mass Matrix and Right Hand Side
    GMassLT = Create_Global_Mass(p, L, Num_Global_DoFs, Num_Elements, Nodes, Elements, Areas, sgn, map);
    GRHSLT = Create_Global_RHSLT(p, Num_Global_DoFs, Num_Elements, Nodes, Elements, Areas, sgn, map, F0, Phi, T);

    ULT = GMassLT\GRHSLT;
	%%% Create Psi_Evals (size(XI), dim(T(p)))
    [tempdim,dmp] = size(XI);
    Psi_Evals = zeros(tempdim, tempdim, (1/2)*(p+1)*(p+2));
    clear dmp tempdim;
    for cnt_psi = 1:(1/2)*(p+1)*(p+2)
       for cnt_T = 1:(1/2)*(p+1)*(p+2)
            Psi_Evals(:,:,cnt_psi) = Psi_Evals(:,:,cnt_psi) + T(cnt_psi, cnt_T).*Phi{cnt_T}(XI, ETA);
       end
       clear cnt_T;
    end
    clear cnt_psi;

    for cnt_Elements = 1:Num_Elements
        
            Verts = Nodes(Elements(cnt_Elements, 1:3), :);
            B1 = (1/2)*(Verts(1,1)+Verts(3,1));  B2 = (1/2)*(Verts(1,2)+Verts(3,2));
            A11 = (1/2)*(Verts(3,1)-Verts(2,1)); A12 =  (1/2)*(Verts(1,1)-Verts(2,1));
            A21 = (1/2)*(Verts(3,2)-Verts(2,2)); A22 = (1/2)*(Verts(1,2)-Verts(2,2));

            Xmap = @(r, s) A11.*r + A12.*s + B1; Ymap = @(r, s) A21.*r + A22.*s + B2; 
            clear   A11 A12 B1 A21 A22 B2;
            X = Xmap(R,S); Y = Ymap(R,S);

            F0_k = F0(X,Y);
            ULT_k = zeros(size(XI));
            for cnt_Modes = 1:(1/2)*(p+1)*(p+2)
                ULT_k = ULT_k + sgn(cnt_Elements, cnt_Modes).*ULT(map(cnt_Elements, cnt_Modes)).*Psi_Evals(:,:,cnt_Modes);
            end
            clear X Y Xmap Ymap cnt_Modes;

            ErrorsD_k(cnt_Elements,1) = (1/4).*sum(sum((F0_k-ULT_k).^2.*WTS.*Areas(cnt_Elements).*(1-ETA)));
    end
    clear cnt_Elements;

    ErrorLT =  sqrt(sum(ErrorsD_k));
    L2_Errors(i_cnt+1,2) = ErrorLT;
end

clear m_first m_final m_range n i_cnt R S Phi Psi_Evals Modes_Xi Modes_Eta 
clear XI ETA ErrorLT ErrorsD_k GRHSLT Mj Tj T F0_k Dub_loc EvalPhi GMassLT
clear T ULT_k WTS xi eta lgl_nodes lgl_weights L Num_Edges Num_Nodes 
clear Num_Elements sgn Num_Global_DoFs

L2_Errors

