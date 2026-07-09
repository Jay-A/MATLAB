function [GMass] = Create_Global_Mass(p, L_Mass, Num_Global_DoFs, ...
    Num_Elems, Nodes, Elements, Areas, sgn, map)

L_Mass = double(L_Mass);

GMass = sparse(Num_Global_DoFs, Num_Global_DoFs);
% GRHS = zeros(Num_Global_DoFs, 1);
% Loop through all elements
for elem_cnt = 1:Num_Elems
    % Create local RHS
    Elem_Verts = Nodes(Elements(elem_cnt, 1:3),:);     
%         [Loc_RHS] = Eval_RHS(F0, Elem_Verts, Areas(elem_cnt), p, Psi);
    %Do Augment sign of appropriate rows of mass matrix and RHS 
    % Add the appropriately signed local systems to the global via map
    for j_cnt = 1:(1/2)*(p+1)*(p+2)
        for k_cnt = 1:(1/2)*(p+1)*(p+2)
            GMass(map(elem_cnt,j_cnt), map(elem_cnt,k_cnt)) = ... 
            GMass(map(elem_cnt,j_cnt),map(elem_cnt,k_cnt)) + ...
            (Areas(elem_cnt)/2).*sgn(elem_cnt, j_cnt)...
            .*sgn(elem_cnt, k_cnt).*L_Mass(j_cnt, k_cnt);
        end
%       GRHS(map(elem_cnt,j_cnt)) = GRHS(map(elem_cnt,j_cnt))...
%          + (sgn(elem_cnt, j_cnt)).*Loc_RHS(j_cnt);
    end      
end

end

