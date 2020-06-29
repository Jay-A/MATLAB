function [Global_DoF_Addrs]=Create_Global_Addresses(Nodes,Edges, Boundary,Elements,Num_Global_DoFs,p)
  % First we initialize the address list 
  Global_DoF_Addrs = zeros(Num_Global_DoFs, 3);

  % Define some constants that will be helpful
  [Num_Nodes, dmp] = size(Nodes);
  [Num_Elements, dmp] = size(Elements);
  Num_Modes = (1/2)*(p+1)*(p+2);
  [Num_Edges, dmp] = size(Edges);
  
  % The first Num_Nodes rows of Global_DoF_Addrs are just the Nodes 
  Global_DoF_Addrs(1:Num_Nodes, 1:2) = Nodes(:,1:2);
  Global_DoF_Addrs(Edges(Boundary,1),3) = 1;
  Global_DoF_Addrs(Edges(Boundary,2),3) = 1;
  
  % We will scatter (p-1) symmetric degrees of freedom along each edge
  [gl_pts, dmp] = lgwt(p-1,-1,1);
  gl_pts = [-1:2/p:1];gl_pts(1)=[]; gl_pts(end) =[];%fliplr(gl_pts');
  
  % Add Edge DoFs
  for edge_cnt = 1:Num_Edges
    isBoundary = ~isempty(find(Boundary==edge_cnt));
    
    param_x = @(t) (1-(t+1)/2).*Nodes(Edges(edge_cnt,1),1)+((t+1)/2).*Nodes(Edges(edge_cnt,2),1);
    param_y = @(t) (1-(t+1)/2).*Nodes(Edges(edge_cnt,1),2)+((t+1)/2).*Nodes(Edges(edge_cnt,2),2);
        
    for order_cnt = 1:(p-1)    
      temp_Edge_addr = Num_Nodes+Num_Edges*(order_cnt-1)+edge_cnt;
      Global_DoF_Addrs(temp_Edge_addr, 1) = param_x(gl_pts(order_cnt));
      Global_DoF_Addrs(temp_Edge_addr, 2) = param_y(gl_pts(order_cnt));
      Global_DoF_Addrs(temp_Edge_addr, 3) = isBoundary;
    end
  end
  
  % Add Int DoFs
  Ref_Ints = zeros((1/2)*(p-1)*(p-2), 2);
  run_addr = 1;
  l_cntr = 0;
  for ints_cnt = (p-2):-1:1
    how_many = ints_cnt;
    l_cntr=l_cntr+1;
    Ref_Ints(run_addr:run_addr+how_many-1, 1:2) = [ones(how_many,1).*gl_pts(l_cntr), gl_pts(1:ints_cnt)'];
    run_addr = run_addr+how_many;
  end
  
  % Define address of last edge mode DoFs
  Start_Ints = Num_Nodes + (p-1)*Num_Edges;
  for elem_cnt = 1:Num_Elements
    % Create a map from reference element to global element
    A = [Nodes(Elements(elem_cnt, 1), 1), Nodes(Elements(elem_cnt, 1), 2)];
    B = [Nodes(Elements(elem_cnt, 2), 1), Nodes(Elements(elem_cnt, 2), 2)];
    C = [Nodes(Elements(elem_cnt, 3), 1), Nodes(Elements(elem_cnt, 3), 2)];
    
    Pmap = @(xref, yref) A(1).*( (1/2).*(yref+1)) + B(1).*(-(1/2).*(xref+yref)) + C(1).*((1/2).*(xref+1));
    Qmap = @(xref, yref) A(2).*( (1/2).*(yref+1)) + B(2).*(-(1/2).*(xref+yref)) + C(2).*((1/2).*(xref+1));
    
    Addrs = [Start_Ints + (1/2)*(p-1)*(p-2)*(elem_cnt-1)+1:Start_Ints + (1/2)*(p-1)*(p-2)*(elem_cnt)];
    Global_DoF_Addrs(Addrs, 1:2) = [Pmap(Ref_Ints(:,1), Ref_Ints(:,2)), Qmap(Ref_Ints(:,1),Ref_Ints(:,2))];
  end