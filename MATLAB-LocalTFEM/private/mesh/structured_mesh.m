function [Nodes, Elements, Areas, Edges, Boundary, numElements, Num_Edges, numNodes] = structured_mesh(n)
% still need to add edges(node 1, node 2, +/- elem1, +/- elem2)
%   need to add boundary (edges on boundary)
%   Num_Elements, Num_Edges, and Num_Nodes
%CREATE_STRUCTURED_MESH Summary of this function goes here
%   Detailed explanation goes here

xStep = (1/2)^(n); yStep = xStep;
numNodesX = 2^(n+1)+1; numNodesY = numNodesX;
numNodes = numNodesX*numNodesY;

xAddrs = -1:xStep:1; yAddrs = xAddrs;
Nodes = zeros(numNodes, 2);

rowCnt = 0;
for cntY = 1:numNodesY
    for cntX = 1:numNodesX
        rowCnt = rowCnt + 1;
        Nodes(rowCnt,:) = [xAddrs(cntX), yAddrs(cntY)];
    end
end

% Nodes
numElements = 8*4^(n);

Elements = zeros(numElements, 3);
elemCnt = 1;
for cntY = 2:numNodesY
    for cntX = 1:numNodesX-1
        addr11 = (cntY-1)*numNodesX+cntX; 
        addr12 = addr11;
        addr21 = (cntY-2)*numNodesX+cntX;
        addr31 = addr21+1;
        addr22 = addr31; addr32 = addr11+1;
        Elements(elemCnt:elemCnt+1, :) = [addr11, addr21, addr31; addr12, addr22, addr32];
        elemCnt = elemCnt+2;
    end
end

Areas = (1/2.*xStep.^2).*ones(numElements,1);

Edges = zeros(3*numElements, 4);
 %Build Edges by looping through elements
   Edges_Row_Cntr = 0;
   for elem_cnt = 1:numElements
    % create naive edge list for Element number elem_cnt
    temp_edge = zeros(3,2);
    temp_edge(1,:) = [Elements(elem_cnt, 1), Elements(elem_cnt,2)];
    temp_edge(2,:) = [Elements(elem_cnt, 2), Elements(elem_cnt,3)];
    temp_edge(3,:) = [Elements(elem_cnt, 3), Elements(elem_cnt,1)];
    
    % these edges will be listed increasing, determine if edge runs this way
    %    for example temp_edge_1 = [3,1], sign_1 = -1 as temp_edge_1 will be 
    %    stored as temp_edge_1 = [1,3].
    signs = zeros(3,1);
    signs(1) = all(sort(temp_edge(1,:))==temp_edge(1,:)); 
    signs(1) = (-1)^(signs(1)+1);
    temp_edge(1,:) = sort(temp_edge(1,:));
    signs(2) = all(sort(temp_edge(2,:))==temp_edge(2,:)); 
    signs(2) = (-1)^(signs(2)+1);
    temp_edge(2,:) = sort(temp_edge(2,:));
    signs(3) = all(sort(temp_edge(3,:))==temp_edge(3,:)); 
    signs(3) = (-1)^(signs(3)+1);
    temp_edge(3,:) = sort(temp_edge(3,:));
    
    %loop through temp edges
    for k = 1:3
      if isempty(find(Edges(:,1)==temp_edge(k,1) & Edges(:,2)==temp_edge(k,2)))
        Edges_Row_Cntr = Edges_Row_Cntr + 1;
        Edges(Edges_Row_Cntr, 1:3) = [temp_edge(k,:), signs(k)*elem_cnt];
        else 
          Edges(find(Edges(:,1)==temp_edge(k,1) & Edges(:,2)==temp_edge(k,2)), ...
                4) = signs(k)*elem_cnt;
      end
    end
    
    %remove unnecessary rows of Edges
    Edges(Edges_Row_Cntr+1:end, :) = [];
   end
   
   % Create Boundary list. This lists the addresses of Edges that are boundaries
   Boundary = find(Edges(:,4) == 0);
   
   
[Num_Edges, dmp] = size(Edges);
clear dmp
end


