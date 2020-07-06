function [] = View_Soln(U, Elements, Nodes, map, sgn, p, Phi, how_rough)%, myfunc)

% How many elements are there?
  [Num_Elems, dmp] = size(Elements);

% load reference triangle to plot solution on
 if ~exist('how_rough', 'var')
  how_rough = 2;
 end
 
 switch(how_rough)
  case 1
    [Ref_Pts, Ref_Tri] = mshLoader('./private/libraries/lib_RefTri_Meshes/Ref_Tri_Rough.msh');
  case 2
    [Ref_Pts, Ref_Tri] = mshLoader('./private/libraries/lib_RefTri_Meshes/Ref_Tri_Med.msh');
  case 3
    [Ref_Pts, Ref_Tri] = mshLoader('./private/libraries/lib_RefTri_Meshes/Ref_Tri_Fine.msh');
  case 4
    [Ref_Pts, Ref_Tri] = mshLoader('./private/libraries/lib_RefTri_Meshes/Ref_Tri_FFine.msh');
 end

 % and define r_pts and s_pts to be the evaluation points of the local basis
  r_pts = Ref_Pts(:,1);
  s_pts = Ref_Pts(:,2);
 
% define local Dubiner basis modes for p
% [Phi, EvalPhis] = Dubiner_Modes(p); 
 
% define mapping from (xi, eta) -> (r, s) for Dub. basis evaluations 
 Ximap = @(r,s) -1+2.*(1+r)./(1-s);
 Etamap = @(r,s) 1.*s;
 % May as well preload xi_pts and eta_pts to avoid calling this mapping a lot
  xi_pts = Ximap(r_pts, s_pts);
  xi_pts(isnan(xi_pts))=1;
  eta_pts = Etamap(r_pts, s_pts);
  
% Loop through elements to plot global solution
 fig = figure('Visible','Off');
  hold on;
  for elem_cnt = 1:Num_Elems
  % reinitialize z_pts to be zeros
    z_pts = zeros(size(xi_pts));
    % Loop through modes to build local visualization
      for mode = 1:(1/2)*(p+1)*(p+2)
        z_pts = z_pts + sgn(elem_cnt,mode).*U(map(elem_cnt,mode)) ...
          .*Phi{mode}(xi_pts, eta_pts);
      end
    % Define mapping to physical elements
      Verts = Nodes(Elements(elem_cnt, 1:3), :);
      B1 = (1/2)*(Verts(1,1)+Verts(3,1));  B2 = (1/2)*(Verts(1,2)+Verts(3,2));
      A11 = (1/2)*(Verts(3,1)-Verts(2,1)); A12 =  (1/2)*(Verts(1,1)-Verts(2,1));
      A21 = (1/2)*(Verts(3,2)-Verts(2,2)); A22 = (1/2)*(Verts(1,2)-Verts(2,2));
     
      Xmap = @(r, s) A11.*r + A12.*s + B1; Ymap = @(r, s) A21.*r + A22.*s + B2; 
%      z_pts = abs(z_pts - myfunc(Xmap(r_pts,s_pts), Ymap(r_pts,s_pts)));
    % map the local element to the global solution
      trisurf(Ref_Tri, Xmap(r_pts,s_pts), Ymap(r_pts,s_pts), z_pts); shading interp
  end
  hold off;
  set(fig, 'visible', 'on')
 
 