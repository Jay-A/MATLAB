% This function defines quadrature over the Reference Triangle 
%              {(r,s) | -1 <= r <= 1 and -1 <= s <= -r\{1}
% To avoid the singularity at s=1 we use gauss-legendre-radau in the s direction
% To allow for ease in continuity constraints we use gauss-legendre-lobatto
%   quadrature in the r direction
function [eval] = RefTri_Quad(myfunc, p);
  
  % Define map from (xi,eta) to (r,s) as per Helenbrook2009
  Rmap = @(xi, eta) (1/2).*(xi+1).*(1-eta)-1;
  Smap = @(xi, eta) 0.*xi + 1.*eta;
  %---------------------------------------------------------------------------%
    
  % Define Quadrature points grids for (r,s)
  [lgl_pts, lgl_wts] = lglnodes(3*p);
  [lgr_pts, lgr_wts, dmp] = lgrnodes(3*p);
  lgr_pts = fliplr(lgr_pts')';
  lgr_wts = fliplr(lgr_wts')';

  [tempxi, tempeta] = meshgrid(lgl_pts, lgr_pts);
  [xi_wts, eta_wts] = meshgrid(lgl_wts, lgr_wts);

  R_quad_pts = Rmap(tempxi, tempeta);
  S_quad_pts = Smap(tempxi, tempeta);
  %---------------------------------------------------------------------------%
  
  % Sum over product of pt_evals and weights
  eval = sum(sum(myfunc(R_quad_pts, S_quad_pts).*xi_wts.*eta_wts...
            .*polyval([-1/2, 1/2], tempeta)));