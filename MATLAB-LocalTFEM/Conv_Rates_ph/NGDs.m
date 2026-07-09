%%%%% Main
%   The intent of this script is to create an array of L2_Errors with
%   dimension p_range x mesh_range.  

clc; close all; clear all;

% set function to approximate 
F0 = @(x,y) cos((2*pi).*x) + sin((2*pi).*y);

% add necessary paths
addpath('../SYMs_L_n_T');
addpath('../');
addpath('../A_Twisted_LT_Script')
addpath('../Quadrature/');

% polynomial degree values to loop through
p_first = 3;  % initial p value
p_final = 10;  % final p value
p_range = p_final - p_first + 1; % this defines the number of values
% Number of mesh refinements to loop through
m_first = 0;  % initial mesh coarseness
m_final = 6;  % final mesh coarseness
Ms = m_first:m_final;
m_range = m_final - m_first + 1; % number of mesh values

% initialize the L2_Error Arrays
L2_Errors = zeros(m_range+1, p_range+1);
L2_Errors(1,1) = 999;
L2_Errors(2:end, 1) = m_first:m_final;
L2_Errors(1, 2:end) = p_first:p_final;

% initialize array for observed spatial convergence rates
Os = zeros(m_range+1, p_range+1);
Os(1,1) = 999;
Os(1, 2:end) = p_first:p_final;
Os(2:end, 1) = m_first:m_final;

% initialize array for recording processing time
Times = Os;

% initialize array for Num_Global_DoFs
NGD = Os;

% Create a structure of meshes
for m = 1:numel(Ms)
    [Mesh(m).Nodes, Mesh(m).Elements, Mesh(m).Areas, Mesh(m).Edges, ...
        Mesh(m).Boundary, Mesh(m).Num_Elements, Mesh(m).Num_Edges, ...
        Mesh(m).Num_Nodes] = Create_Structured_Mesh(Ms(m));
end

p_ind = 1;
for p_cnt = p_first:p_final
    p = p_cnt;
    fprintf('p = %d of %d \n', p, p_final);
    pause(0.000001);
    p_ind = p_ind + 1;
    
    % set some indexing constants
    tm = (1/2)*(p+1)*(p+2);
    vm = 3;
    sm1 = vm + (1:3:3*(p-1)); sm2 = sm1+1; sm3 = sm1+2;
    im = sm3(end)+1:tm;
    
%     % set local twisted structures and modified dubiner functions
%     [Tt, Mt, M,T, Mtapp] = Local_Twisted_Structs(p); % the function call says it all
%     [Phi, EvalPhi, Modes_Xi, Modes_Eta] = Dubiner_Modes(p);
    
    % loop through mesh refinements
    m_ind = 1;
    for m_cnt = 1:numel(Ms)
        fprintf('p[%d of %d] with m[%d of %d] \n', (p_ind-1), p_range, m_ind, m_range);
        pause(0.000001);
        m_ind = m_ind + 1;
        % load mesh objects
        Nodes = Mesh(m_cnt).Nodes;
        Elements = Mesh(m_cnt).Elements;
        Areas = Mesh(m_cnt).Areas;
        Edges = Mesh(m_cnt).Edges;
        Boundary = Mesh(m_cnt).Boundary;
        Num_Elements = Mesh(m_cnt).Num_Elements;
        Num_Edges = Mesh(m_cnt).Num_Edges;
        Num_Nodes = Mesh(m_cnt).Num_Nodes;
%         [Nodes, Elements, Areas, Edges, Boundary, Num_Elements, Num_Edges, Num_Nodes] = Create_Structured_Mesh(m_cnt);
        
        %% Create some Global Objects
        [map, sgn] = Create_Map(p, Edges, Nodes, Elements);
        Num_Global_DoFs = Num_Nodes + (p-1)*Num_Edges + (1/2)*(p-1)*(p-2)*Num_Elements;
        [Global_DoF_Addrs]=Create_Global_Addresses(Nodes,Edges, Boundary,Elements,Num_Global_DoFs,p);
        NGD(m_ind, p_ind) = Num_Global_DoFs;
    end
end
