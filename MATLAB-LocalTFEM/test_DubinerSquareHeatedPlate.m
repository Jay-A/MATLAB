clc; close all; clear all;

p=3;

% add necessary paths
addpath('./private/dubiner/');
addpath('./private/quadrature/');
addpath('./private/mesh/');
addpath('./private/assembly/');
addpath('./private/local/');
addpath('./private/visualization/');

% set problem specifics
init_heat = @(x,y) 2.*(x+1).*(x-1).*(y+1).*(y-1);% cos((pi/2).*x) + sin((pi/2).*y);
% as an example only, set left boundary (x = -1), to homogeneous Dirichlet
D_Bndry = 0; % imposed on U
mesh_choice = 2; % from 0 to 8 (roughest to finest) Note CFL condition
t0 = 0; tFin = 3; % t0 = 0, tFin = final time
delT = .0001; % time step must decrease with mesh size to meet cfl condition
Time = t0:delT:tFin; 
Alpha = 1/2;  % thermal diffusivity
how_rough=2; % quality of frame renders 1-4 for video;

% set some indexing constants
tm = (1/2)*(p+1)*(p+2); % number of total local modes
vm = 3; % number of local vertex modes
sm1 = vm + (1:3:3*(p-1)); sm2 = sm1+1; sm3 = sm1+2; % local side modes 
im = sm3(end)+1:tm; % local interior modes

% set local twisted structures and modified dubiner functions
[M] = Dubiner_Mass(p); % the function call says it all
[Phi, EvalPhi, Modes_Xi, Modes_Eta] = Dubiner_Modes(p);
GradPhi = GradPhi(p);

% Some things needed to symbolically integrate Ke for the time being
syms r s; 
PhiRef = Phi;
GradPhiRef = PhiRef;
Xi = @(r,s) -1+2.*(1+r)./(1-s);
Eta = @(r,s) s;
syms x y
for fcount = 1:tm
    PhiRef{fcount} =@(r,s) Phi{fcount}(Xi(r,s), Eta(r,s));
end
for fcount = 1:tm
    temp1 = matlabFunction(diff(PhiRef{fcount}(r,s),r));
    temp2 = matlabFunction(diff(PhiRef{fcount}(r,s),s));
    
    addr1 = find(char(temp1) == ')'); addr1 = addr1(1)+1;
    temp1 = char(temp1); temp1 = temp1(addr1:end);
    GradPhiRef{fcount,1} = str2func("@(r,s) 0.*r+(" + temp1 + ")");
    
    addr2 = find(char(temp2) == ')'); addr2 = addr2(1)+1;
    temp2 = char(temp2); temp2 = temp2(addr2:end);
    GradPhiRef{fcount,2} = str2func("@(r,s) 0.*s+(" + temp2 + ")");
end

%% Load Mesh Objects
eval(sprintf("load ./private/libraries/lib_Structured_Meshes/Areas%d.mat;",mesh_choice));
eval(sprintf("load ./private/libraries/lib_Structured_Meshes/Boundary%d.mat;",mesh_choice));
eval(sprintf("load ./private/libraries/lib_Structured_Meshes/Edges%d.mat;",mesh_choice));
eval(sprintf("load ./private/libraries/lib_Structured_Meshes/Elements%d.mat;",mesh_choice));
eval(sprintf("load ./private/libraries/lib_Structured_Meshes/Nodes%d.mat;",mesh_choice));
eval(sprintf("load ./private/libraries/lib_Structured_Meshes/Num_Edges%d.mat;",mesh_choice));
eval(sprintf("load ./private/libraries/lib_Structured_Meshes/Num_Elements%d.mat;",mesh_choice));
eval(sprintf("load ./private/libraries/lib_Structured_Meshes/Num_Nodes%d.mat;",mesh_choice));

% Nodes =  [-1 1; 1 -1; 1 1 ]; % Nodes(1:3,:);
% Areas = Areas(1); Edges = Edges(1:3,:); Edges(3,4) = 0; Boundary = [1;2;3];
% Elements = Elements(1,:);  Num_Edges = 3;
% Num_Elements = 1; Num_Nodes = 3;

%% Create some Global Objects
[map, sgn] = Create_Map(p, Edges, Nodes, Elements);
Num_Global_DoFs = Num_Nodes + (p-1)*Num_Edges + (1/2)*(p-1)*(p-2)*Num_Elements;
[Global_DoF_Addrs]=Create_Global_Addresses(Nodes,Edges, Boundary,Elements,Num_Global_DoFs,p);
U = zeros(Num_Global_DoFs, numel(Time));
Ke = zeros(tm, tm, Num_Elements);

GMass = Create_Global_Mass(p, M, Num_Global_DoFs, Num_Elements, Nodes, Elements, Areas, sgn, map);
    
GK = zeros(size(GMass,1), size(GMass,2));
GRHS = zeros(Num_Global_DoFs, 1);
 
% Find U(:, 1), the initial temperature distribution
%   and build GK
for elem_cnt = 1:Num_Elements
    tic;
    
    pause(0.0001);
    % Create local RHS
    Elem_Verts = Nodes(Elements(elem_cnt, 1:3),:);     
    [Loc_RHS] = Eval_RHS(@(x,y) init_heat(x,y), Elem_Verts, Areas(elem_cnt), p, Phi);
%     GRHS(map(elem_cnt,:)) = GRHS(map(elem_cnt,:))+sgn(elem_cnt,:)'.*Loc_RHS;
    for j_cnt = 1:(1/2)*(p+1)*(p+2)
      GRHS(map(elem_cnt,j_cnt)) = GRHS(map(elem_cnt,j_cnt))...
         + (sgn(elem_cnt, j_cnt)).*Loc_RHS(j_cnt);
    end  
    
    Ke_temp = zeros(tm, tm);
    Verts = Nodes(Elements(elem_cnt, 1:3), :);
    B1 = (1/2)*(Verts(1,1)+Verts(3,1));  B2 = (1/2)*(Verts(1,2)+Verts(3,2));
    A11 = (1/2)*(Verts(3,1)-Verts(2,1)); A12 =  (1/2)*(Verts(1,1)-Verts(2,1));
    A21 = (1/2)*(Verts(3,2)-Verts(2,2)); A22 = (1/2)*(Verts(1,2)-Verts(2,2));
    Jk = [A11, A12; A21, A22]; B = [B1;B2];   % Jk is Jacobian for ref to physical
    invJk = inv(Jk); BBack = Jk\(-B);  % invJk is Jacobian for physical to ref
    %%%% Note that [Gradx(.); Grady(.)] = (JT)^{-1} [Gradr(.); Grads(.)]
    Area = (1/2).*( Verts(1,1)*(Verts(2,2)-Verts(3,2)) + Verts(2,1)*(Verts(3,2)-Verts(1,2)) + Verts(3,1)*(Verts(1,2)-Verts(2,2)));   
    H = invJk*invJk';
    
       
    for icount = 1:tm
        for jcount = 1:tm
%           % Symbolic integration takes a long time
%             Ke_temp(icount, jcount) = (Area/2).* (int(int( ...
%                 GradPhiRef{jcount,1}(r,s).* (H(1,1).*GradPhiRef{icount,1}(r,s) + H(1,2).*GradPhiRef{icount,2}(r,s)) ...
%                 + GradPhiRef{jcount,2}(r,s).*( H(2,1).*(GradPhiRef{icount,1}(r,s)) + H(2,2).*GradPhiRef{icount,2}(r,s))...
%                 ,r,-1,-s ),s,-1,1 )); 
            % Numerical integration can be exact to round off and is faster
            Ke_temp(icount, jcount) = (Area/2).* (RefTri_Quad( @(r,s) ...
                GradPhiRef{jcount,1}(r,s).* (H(1,1).*GradPhiRef{icount,1}(r,s)...
                + H(1,2).*GradPhiRef{icount,2}(r,s)) + GradPhiRef{jcount,2}(r,s)...
                .*( H(2,1).*(GradPhiRef{icount,1}(r,s)) + H(2,2).*GradPhiRef{icount,2}(r,s)), p)); 
            GK(map(elem_cnt,icount), map(elem_cnt, jcount)) = GK(map(elem_cnt,icount), map(elem_cnt, jcount))+sgn(elem_cnt, icount)*sgn(elem_cnt, jcount)*Ke_temp(icount, jcount);
       end
    end
    
    Ke(:,:,elem_cnt) = Ke_temp;   
end
% % to impose Dirichlet boundary on x = -1
    DBndry_Addrs = find(Global_DoF_Addrs(:,1) == -1);

U(:,1) = GMass\GRHS;
% enforce dirichlet boundary condition 
    U(DBndry_Addrs,1) = D_Bndry;
    
for t_cnt = 2:numel(Time)
    U(:,t_cnt) = U(:,t_cnt-1) - Alpha.*delT.*((GMass\GK)*U(:,t_cnt-1));
    % enforce dirichlet boundary condition 
        U(DBndry_Addrs,t_cnt) = D_Bndry;
    if (mod(t_cnt,100) == 0)
        fprintf('time step %d of %d \n', t_cnt, numel(Time));
        pause(0.0001);
    end
end

Make_Video;