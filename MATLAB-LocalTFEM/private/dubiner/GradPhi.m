function [DPhi] = GradPhi(p)

if ~exist('alpha', 'var')
  alpha = 2; beta = 2;
end

% Initialize a cell array structure of Dubiner modes ordered via Helenbrook's 
% ordering of:
%   A,B,C,  A0,B0,C0,  A1,B1,C1, ... A(p-2),B(p-2),C(p-2),  interiors by j(m,n)
DPhi = cell((1/2)*(p+1)*(p+2),2);
j_ind = @(m,n) (1/2)*(m+n)*(m+n+1)+n+1;

% % PolyModes provides a storage structure of the basis functions that allows usejava
% %  of polyder to compute partial derivatives.
% Modes_Xi = cell((1/2)*(p+1)*(p+2),1); Modes_Eta = cell((1/2)*(p+1)*(p+2),1);
% %-----------------------------------------------------------------------------%

% Define vertex modes 
DPhi{1,1} = @(xi, eta) 0; 
DPhi{1,2} = @(xi, eta) 1/2;

DPhi{2,1} = @(xi, eta) ((-1+eta)./4);
DPhi{2,2} = @(xi, eta) ((-1+xi)./4);

DPhi{3,1} = @(xi, eta) (1/4).*(1-eta);
DPhi{3,2} = @(xi, eta) (-1/4).*(1+xi);

%-----------------------------------------------------------------------------%

% Define the side modes Am,Bm,Cm for 1 <= m <= p-1 as per Appleton2019
% Special consideration is taken when m = 1, because of the P_{m-2}^{a,b}
% term
for m = 1:(p-1)

  f_mab = Create_Poly(Jacobi_Coeffs(m-1, alpha, beta));
  df_mab = @(xi) 0;
%   if m~=1
%       df_mab = Create_Poly(Jacobi_Coeffs(m-2, alpha+1, beta+1));
%   else
%       df_mab = @(xi) 0;
%   end 
  
  % partial phi_Am, partial xi then partial eta
  DPhi{3*(m)+1,1} = @(xi, eta) (1/2).*((1-eta)./2).^((m)+1).* ...
      (-xi.*f_mab(xi)+((1-xi.^2)./4).*(4+m).*df_mab(xi));
  DPhi{3*(m)+1,2} = @(xi, eta) -((m+1)./8).*(1-xi.^2).* ...
      f_mab(xi).*((1-eta)./2).^m;
  
  % partial phi_Bm, partial xi then partial eta
  DPhi{3*(m)+2,1} = @(xi, eta) ((1-eta.^2)./8).*f_mab(eta);
  DPhi{3*(m)+2,2} = @(xi, eta) ((1+xi)./4).* ...
      (-eta.*f_mab(eta)+((1-eta.^2)./4).*(4+m).*df_mab(eta));
  
  % partial phi_Cm, partial xi then partial eta      
  DPhi{3*(m)+3,1} = @(xi, eta) (((-1).^m)./8).*(1-eta.^2).*f_mab(eta);
  DPhi{3*(m)+3,2} = @(xi, eta) ((-1).^(m-1)).*((1-xi)./4).* ...
      (-eta.*f_mab(eta)+(4+m).*((1-eta.^2)./4).*df_mab(eta));
end
clear m;
%-----------------------------------------------------------------------------%

% Define the interior modes organized by j_ind(m,n) for 0 <= m <= p-3 
%  ...                                              and 0 <= n <= p-3-m
for m = 0:(p-3)
  for n = 0:(p-3)-m
    f_mab = Create_Poly(Jacobi_Coeffs(m, alpha, beta));
    f_mab2 = Create_Poly(Jacobi_Coeffs(n, 2*m+4*(alpha-1)+1, beta));
    
      if m~=0
        df_mab = Create_Poly(Jacobi_Coeffs(m-1, 3,3));
      else
        df_mab = @(xi) 0;
      end
      
      if n~=0
        df_mab2 = Create_Poly(Jacobi_Coeffs(n-1, 2*m+6,3));
      else
        df_mab2 = @(eta) 0;
      end
      
    DphiImnXidxi = @(xi) -(xi./2).*f_mab(xi)+(5+m).* ...
      ((1-xi.^2)./8).*df_mab(xi);
    DphiImnEtadeta =@(eta) (((1-eta)./4)-(m+2)).*((1-eta)./2).^(m+1).*...
      f_mab2(eta) + (8+2*m+n).*((1+eta)./4).* ...
      ((1-eta)./2).^(m+2).*df_mab2(eta);
  
    DPhi{3*p+j_ind(m,n),1} = @(xi, eta) (-(xi./2).*f_mab(xi)+(5+m).* ...
      ((1-xi.^2)./8).*df_mab(xi)) .* ...
      ((1-eta)./2).^(m+2).*((1+eta)./2).*f_mab2(eta);
    DPhi{3*p+j_ind(m,n),2} = @(xi, eta) (1./4).*(1-xi.^2).*f_mab(xi).* ...
      ((((1-eta)./4)-(m+2)).*((1-eta)./2).^(m+1).*...
      f_mab2(eta) + (8+2*m+n).*((1+eta)./4).* ...
      ((1-eta)./2).^(m+2).*df_mab2(eta));      
  end
end
% %-----------------------------------------------------------------------------%