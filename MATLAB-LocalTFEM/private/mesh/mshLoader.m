function [Nodes, Elements, Areas] = mshLoader(fname)
    %------------------------------------------------------------------------%
    %------ Gmsh to Matlab function: Import mesh to matlab---------------------%
    %------------------------------------------------------------------------%

    %-----------------------------------------------------------------------%
    % dlmread(filename,delimiter,[R1 C1 R2 C2]) reads only the range 
    % bounded by row offsets R1 and R2 and column offsets C1 and C2.
    %-----------------------------------------------------------------------%

    % no of nodes is mentioned in 5th row and first column

    N_n      = dlmread(fname,'',[5-1 1-1 5-1 1-1]);
    N_e      = dlmread(fname,'',[7+N_n 0 7+N_n 0]);

    node_id     = dlmread(fname,'',[5 0 4+N_n 0]);
    Nodes       = dlmread(fname,'',[5 1 4+N_n 3]);
    elements    = dlmread(fname,'',[8+N_n 0 7+N_n+N_e 7]);

    %------- 2D Geometry

    two_d_nodes = Nodes(:,1:2);
    elem_type   = elements(:,2);

    %--- find the starting indices of 2D elements
    two_ind = 1;
    for i = 1:N_e
        if(elem_type(i) ~= 2)
            two_ind = two_ind+1;
        end
    end
    %----------------------------------------------

    two_d_elements(1:N_e-two_ind,1:3) = 0;
    k = 1;
    for i = two_ind:N_e
       two_d_elements(k,1:3) = elements(i,6:8);

       k = k+1;
    end
    Elements = two_d_elements;
    
    [N_e, dmp] = size(Elements);
    Areas = zeros(N_e, 1);
    for e_cnt = 1:N_e
       A = [Nodes(two_d_elements(e_cnt,1), 1), Nodes(two_d_elements(e_cnt,1), 2)];
       B = [Nodes(two_d_elements(e_cnt,2), 1), Nodes(two_d_elements(e_cnt,2), 2)];
       C = [Nodes(two_d_elements(e_cnt,3), 1), Nodes(two_d_elements(e_cnt,3), 2)];
       Areas(e_cnt) = (1/2)*det([A(1), A(2), 1; B(1), B(2), 1; C(1), C(2), 1 ]);    
    end
    OrderSwaps = find(Areas<0)
    Elements(OrderSwaps, :) = Elements(OrderSwaps, :)*[0,1,0;1,0,0;0,0,1];
    Areas = abs(Areas);
    
    Nodes = Nodes(:,1:2);

    %---- visualize in matlab ---------------------

%    figure(1);
%    triplot(two_d_elements,two_d_nodes(:,1),two_d_nodes(:,2))
%    xlabel('X','fontsize',14)
%    ylabel('Y','fontsize',14)
%    title('GMsh to MATLAB import','fontsize',14)
%    fh = figure(1);
%    set(fh, 'color', 'white'); 

    %-------------------------------------------------------------------------