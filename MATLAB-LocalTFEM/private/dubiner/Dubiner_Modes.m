% a file that creates a cell array of Dubiner modes for p
% Defined over the reference triangle with a singularity at vertex 1
% therefore Phi{i}:(r,s)->R with  -1 <= s < 1 and -1 <= r <= s U r=1
function [Phi, EvalPhi, Modes_Xi, Modes_Eta] = Dubiner_Modes(p, alpha, beta)

if ~exist('alpha', 'var')
  alpha = 2; beta = 2;
end

% Initialize a cell array structure of Dubiner modes ordered via Helenbrook's 
% ordering of:
%   A,B,C,  A0,B0,C0,  A1,B1,C1, ... A(p-2),B(p-2),C(p-2),  interiors by j(m,n)
Phi = cell((1/2)*(p+1)*(p+2),1);
j_ind = @(m,n) (1/2)*(m+n)*(m+n+1)+n+1;

% % PolyModes provides a storage structure of the basis functions that allows usejava
% %  of polyder to compute partial derivatives.
% Modes_Xi = cell((1/2)*(p+1)*(p+2),1); Modes_Eta = cell((1/2)*(p+1)*(p+2),1);
% %-----------------------------------------------------------------------------%

% Define vertex modes 
Phi{1} = @(xi, eta) (1/2).*(1+eta); 
  Modes_Xi{1} = [0,1]; Modes_Eta{1} = (1/2).*[1, 1];
Phi{2} = @(xi, eta) (1/4).*(1-xi).*(1-eta);
  Modes_Xi{2} = (1/2).*[-1,1]; Modes_Eta{2} = (1/2).*[-1,1];
Phi{3} = @(xi, eta) (1/4).*(1+xi).*(1-eta);
  Modes_Xi{3} = (1/2).*[1,1]; Modes_Eta{3} = (1/2).*[-1,1];
%-----------------------------------------------------------------------------%

% Define the side modes Am,Bm,Cm for 0 <= m <= p-2 as per Helenbrook2009
for m = 0:(p-2)
% Evaluate (1-eta)^(m+2)
  temp_eta_term = [0,1];
  for pwr_cnt = 1:m+2
    temp_eta_term = conv(temp_eta_term, (1/2).*[-1,1]);
  end

  f_mab = Create_Poly(Jacobi_Coeffs(m, alpha, beta));
  
  Phi{3*(m+1)+1} = @(xi, eta) ((1/2)^(m+4)).*(1+xi).*(1-xi).*(1-eta).^(m+2)... 
                    .* f_mab(xi);               
    Modes_Xi{3*(m+1)+1} =  conv(conv((1/2).*[1,1],(1/2).*[-1,1]), ...
        Jacobi_Coeffs(m,alpha,beta));
    Modes_Eta{3*(m+1)+1} = temp_eta_term;
    
  Phi{3*(m+1)+2} = @(xi, eta) ((1/2)^3).*(1+xi).*(1-eta).*(1+eta)...
                    .*f_mab(eta);
    Modes_Xi{3*(m+1)+2} =  (1/2).*[1,1];
    Modes_Eta{3*(m+1)+2} = conv(conv((1/2).*[1,1], (1/2).*[-1,1]), ...
        Jacobi_Coeffs(m,alpha,beta));
        
  Phi{3*(m+1)+3} = @(xi, eta) ((-1)^m).*((1/2)^3).*(1-xi).*(1-eta).*(1+eta)...
                    .*f_mab(eta);
    Modes_Xi{3*(m+1)+3} = (-1).^m.*(1/2).*[-1,1];
    Modes_Eta{3*(m+1)+3} = conv(conv((1/2).*[1,1], (1/2).*[-1,1]), ...
        Jacobi_Coeffs(m,alpha,beta));
end
clear m;
%-----------------------------------------------------------------------------%

% Define the interior modes organized by j_ind(m,n) for 0 <= m <= p-3 
%  ...                                              and 0 <= n <= p-3-m
for m = 0:(p-3)

  temp_eta_term = [0,1];
    for pwr_cnt = 1:m+2
      temp_eta_term = conv(temp_eta_term, (1/2).*[-1,1]);
    end
    
  for n = 0:(p-3)-m
      f_mab = Create_Poly(Jacobi_Coeffs(m, alpha, beta));
      f_mab_2 = Create_Poly(Jacobi_Coeffs(n, 2*m+4*(alpha-1)+1, beta));
      
    Phi{3*p+j_ind(m,n)} = @(xi, eta) ((1/2)^(m+5)).*(1-xi).*(1+xi)...
          .*(1+eta).*((1-eta).^(m+2)).*f_mab(xi) .*f_mab_2(eta);
      
      Modes_Xi{3*p+j_ind(m,n)} = conv(conv((1/2).*[1,1], (1/2).*[-1,1]), ...
          Jacobi_Coeffs(m,alpha,beta));
      Modes_Eta{3*p+j_ind(m,n)} = conv(conv(temp_eta_term, (1/2).*[1,1]), ...
          Jacobi_Coeffs(n,2*m+4*(alpha-1)+1,beta));
  end
end
%-----------------------------------------------------------------------------%

% Define the vector evaluation of all Phi modes
EvalPhi = @(xi,eta) cellfun(@(f)f(xi,eta),Phi);