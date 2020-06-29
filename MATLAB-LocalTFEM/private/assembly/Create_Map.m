%  Accepting as input polynomial basis order, a global edge list, 
%     a global node list, and a connectivity array for elements describing 
%     order of vertices for each element, Create_Map(p, Edges, Nodes, Conn)
%     returns two arrays, map[e][i] and sgn[e][i] describing the local to 
%     global address mapping and the sign of each mode for the assembly process.

function [map, sgn] = Create_Map(p, Edges, Nodes, Conn)
  % Initialize map and sgn
  [Num_Nodes, dmp] = size(Nodes);
  [Num_Edges, dmp] = size(Edges);
  [Num_Elems, dmp] = size(Conn);
  Loc_Dim = (1/2)*(p+1)*(p+2);
  Loc_Vert_Addrs = 1:3;
  Loc_Side_Addrs = 4:3+3*(p-1);
  Loc_Ints_Addrs = 4+3*(p-1):3+3*(p-1)+(1/2)*(p-1)*(p-2);
  
  Num_Glob_DoFs = Num_Nodes + (p-1)*Num_Edges + (1/2)*(p-1)*(p-2)*Num_Elems;
  Glob_Vert_Addrs = 1:Num_Nodes;
  Glob_Side_Addrs = Num_Nodes+1:Num_Nodes+(p-1)*Num_Edges;
  Glob_Int_Addrs = Num_Nodes+(p-1)*Num_Edges+1:Num_Glob_DoFs;
  
  map = zeros(Num_Elems, Loc_Dim);
  sgn = ones(Num_Elems, Loc_Dim);
  %--------------------------------------------------------------------------%

  % Define Local to global vertex mappings to be the connectivity array with 
  %  only positive signs
  map(:,Loc_Vert_Addrs) = Conn;
  %--------------------------------------------------------------------------%

  % Define Interior and Edges to global mappings
  for elem_cnt = 1:Num_Elems
      
    % Define Interior to global mappings
    intDim = (1/2)*(p-1)*(p-2);
    intShift = (elem_cnt-1)*intDim;
    map(elem_cnt, Loc_Ints_Addrs) = Glob_Int_Addrs(intShift+1:intShift+intDim);  
    %--------------------------------------------------------------------------%
    
    % Now define the edge mappings
    % Define global addresses of local vertices
    GlobVerts = map(elem_cnt,1:3);
    %------------------------------------------------------------------------%
    
    % Define global addresses of local edges
    tempEdges = (find(abs(Edges(:,3))==elem_cnt | abs(Edges(:,4))==elem_cnt));
    % and organize the edges for this element
    for edge_cnt = 1:3
      tempEdge(edge_cnt) = find((Edges(tempEdges,1)-GlobVerts(edge_cnt)).* ...
        (Edges(tempEdges,2)-GlobVerts(edge_cnt)));
    end

    clear edge_cnt;
    % Then the global edge addresses of the local element is defined
    GlobEdges = [tempEdges(tempEdge(1)),tempEdges(tempEdge(2)),...
        tempEdges(tempEdge(3))];
    %------------------------------------------------------------------------%    
    
    % Define edge_sign to be 1 if the edge runs in the same
    %  direction and -1 if it runs in the opposite direction
    for edge_cnt = 1:3
      edge = GlobEdges(edge_cnt);
      edge_sign(edge_cnt) = sign(Edges(edge,2+find(abs(Edges(edge,3:4))...
        ==elem_cnt)));
    end 
    %------------------------------------------------------------------------%
    
    % Define global degree of freedom and sgn for side modes
    for mode_cnt = 1:(p-1)
      % Is the current mode odd or even
      isodd = mod(mode_cnt-1,2);
      % define the map indeces of egde modes 1,2, and 3 order mode_cnt
      col_addrs = 4+(mode_cnt-1)*3:3+mode_cnt*3;
      mode_shift = (mode_cnt-1)*Num_Edges;
      A123m_Addrs = Num_Nodes+mode_shift+GlobEdges;
      % Now state the appropriate elements of map and sgn
      map(elem_cnt, col_addrs) = A123m_Addrs;
      sgn(elem_cnt, col_addrs) = edge_sign.^isodd;     
    end
  end  
  %--------------------------------------------------------------------------%