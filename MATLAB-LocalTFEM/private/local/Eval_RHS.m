function [RHS] = Eval_RHS(func, Verts, Area, p, Phi)

  % Lets define Xmap and Ymap to bring (r,s) -> (x,y)
  Xmap = @(r,s) (1/2)*(Verts(3,1)-Verts(2,1)).*r ...
          + (1/2)*(Verts(1,1)-Verts(2,1)).*s + (1/2)*(Verts(1,1)+Verts(3,1));
  Ymap = @(r,s) (1/2)*(Verts(3,2)-Verts(2,2)).*r ...
          + (1/2)*(Verts(1,2)-Verts(2,2)).*s + (1/2)*(Verts(1,2)+Verts(3,2));
  

    [lgl_nodes,lgl_weights] = lglnodes(15);
    [XI, ETA] = meshgrid(lgl_nodes, lgl_nodes);
    WTS = lgl_weights*lgl_weights';

    S = ETA; 
    R = (1/2).*(XI+1).*(1-ETA)-1;
    
    X = Xmap(R,S); Y = Ymap(R,S);
    
    F0 = func(X,Y);
      
  RHS = zeros((1/2)*(p+1)*(p+2), 1);
  % Define RHS_i to be int func(xmap(rmap(xi,eta.....).*Psi{i}(xi,eta)

  for mode_cnt = 1:(1/2)*(p+1)*(p+2)
%     This_func = @(xi,eta) func(Xmap(Rmap(xi,eta),Smap(xi,eta)),...
%                     Ymap(Rmap(xi,eta),Smap(xi,eta)))...
%                     .*Phi{mode_cnt}(xi,eta).*(Area/4).*(1-eta);
    RHS(mode_cnt) = (Area/4).*sum(sum( F0.*Phi{mode_cnt}(XI,ETA).*WTS.*(1-ETA) ));
%     RHS(mode_cnt) = (int( int(func(Xmap(Rmap(x,y),y), Ymap(Rmap(x,y),y)).*Phi{mode_cnt}(x, y).*w(x,y).*(Area/2),x,-1,1 ),y,-1,1 ));
  end