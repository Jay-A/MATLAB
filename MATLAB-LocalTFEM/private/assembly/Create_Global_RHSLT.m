function [GRHS] = Create_Global_RHSLT(p, Num_Global_DoFs, Num_Elems, ...
    Nodes, Elements, Areas, sgn, map, F0, Phi, T)
%CREATE_GLOBAL_RHS Summary of this function goes here

% First initialize GRHS
GRHS = zeros(Num_Global_DoFs, 1);
% Loop through all elements
for elem_cnt = 1:Num_Elems
    % Create local RHS
    Elem_Verts = Nodes(Elements(elem_cnt, 1:3),:);     
    [Loc_RHS] = Eval_RHSLT(@(x,y) F0(x,y), Elem_Verts, Areas(elem_cnt), p, Phi, T);
    %Do Augment sign of appropriate rows of RHS 
    % Add the appropriately signed local systems to the global via map
    for j_cnt = 1:(1/2)*(p+1)*(p+2)
      GRHS(map(elem_cnt,j_cnt)) = GRHS(map(elem_cnt,j_cnt))...
         + (sgn(elem_cnt, j_cnt)).*Loc_RHS(j_cnt);
    end      
end

end

