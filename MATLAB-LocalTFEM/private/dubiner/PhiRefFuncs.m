function [PhiRef, GradPhiRef] = PhiRefFuncs(p)

dim = (1/2)*(p+1)*(p+2);
[Phi, EvalPhi, Modes_Xi, Modes_Eta] = Dubiner_Modes(p);

syms r s; 
PhiRef = cell(dim,1);
GradPhiRef = cell(dim,2);
Xi = @(r,s) -1+2.*(1+r)./(1-s);
Eta = @(r,s) s;
syms x y
for fcount = 1:dim
    PhiRef{fcount} =@(r,s) Phi{fcount}(Xi(r,s), Eta(r,s));
end
for fcount = 1:dim
    temp1 = matlabFunction(diff(PhiRef{fcount}(r,s),r));
    temp2 = matlabFunction(diff(PhiRef{fcount}(r,s),s));
    
    addr1 = find(char(temp1) == ')'); addr1 = addr1(1)+1;
    temp1 = char(temp1); temp1 = temp1(addr1:end);
    GradPhiRef{fcount,1} = str2func("@(r,s) 0.*r+(" + temp1 + ")");
    
    addr2 = find(char(temp2) == ')'); addr2 = addr2(1)+1;
    temp2 = char(temp2); temp2 = temp2(addr2:end);
    GradPhiRef{fcount,2} = str2func("@(r,s) 0.*s+(" + temp2 + ")");
end