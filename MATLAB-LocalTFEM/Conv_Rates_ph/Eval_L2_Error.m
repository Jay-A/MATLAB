function [error] = Eval_L2_Error(p, U, F0, Phi, map, sgn, nElems, Nodes, Elements, Areas)
    [lgl_nodes,lgl_weights] = lglnodes(2*(p+1));
    [XI, ETA] = meshgrid(lgl_nodes, lgl_nodes);
    WTS = lgl_weights*lgl_weights';

    S = ETA; 
    R = (1/2).*(XI+1).*(1-ETA)-1;
    [tempdim,dmp] = size(XI);
    Phi_Evals = zeros(tempdim, tempdim, (1/2)*(p+1)*(p+2));
    
    errors_k = zeros(nElems, 1);
    for mode_cnt = 1:(1/2)*(p+1)*(p+2)
        Phi_Evals(:,:,mode_cnt) = Phi{mode_cnt}(XI, ETA);
    end
    
    for elem_cnt = 1:nElems
       
        Verts = Nodes(Elements(elem_cnt, 1:3), :);
        B1 = (1/2)*(Verts(1,1)+Verts(3,1));  B2 = (1/2)*(Verts(1,2)+Verts(3,2));
        A11 = (1/2)*(Verts(3,1)-Verts(2,1)); A12 =  (1/2)*(Verts(1,1)-Verts(2,1));
        A21 = (1/2)*(Verts(3,2)-Verts(2,2)); A22 = (1/2)*(Verts(1,2)-Verts(2,2));

        Xmap = @(r, s) A11.*r + A12.*s + B1; Ymap = @(r, s) A21.*r + A22.*s + B2;  
        X = Xmap(R,S); Y = Ymap(R,S);
        
        F0_k = F0(X,Y);
        U_k = zeros(size(XI));
        for mode_cnt = 1:(1/2)*(p+1)*(p+2)
            U_k = U_k + sgn(elem_cnt, mode_cnt).*U(map(elem_cnt, mode_cnt)).*Phi_Evals(:,:,mode_cnt);
        end
        
        errors_k(elem_cnt,1) = (1/4).*sum(sum((F0_k-U_k).^2.*WTS.*Areas(elem_cnt).*(1-ETA)));
    end
    error = sqrt(sum(errors_k));
end