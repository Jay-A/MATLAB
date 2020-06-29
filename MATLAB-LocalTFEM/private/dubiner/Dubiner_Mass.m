% This function creates the Dubiner Mass matrix of order p following 
%  Helenbrook2009 ordering
function [myD] = Dubiner_Mass(p, alpha, beta)

   if ~exist('alpha', 'var')
    alpha = 2; beta = 2;
   end

  dim = (1/2)*(p+1)*(p+2);
  syms xi eta;
  
%  Ximap = @(r,s) -1+2.*(1+r)./(1-s);
%  Etamap = @(r,s) 1.*s;
  
  % Define basis functions
  [Phi, EvalPhis] = Dubiner_Modes(p, alpha, beta);
  % and the weighting function
  w = @(xi, eta) (1-eta)./2;
  %--------------------------------------------------------------------------%
  
  % Define Dubiner mass matrix by looping through modes and calling RefTri_Quad
  myD = sym(zeros(dim, dim));
  for m = 1:dim
    for n = 1:dim
      myD(m,n) = sym(int(int(Phi{m}(xi,eta).*Phi{n}(xi,eta).*w(xi,eta), xi, -1,1), eta, -1,1));
    end
  end
  
%  restoredefaultpath;