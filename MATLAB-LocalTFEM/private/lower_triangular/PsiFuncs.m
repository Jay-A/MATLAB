%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PsiRef{dim,1} = PsiRefFuncs(T) is a function that accepts the input
% parameter T, the change of basis, and returns a cell PsiRef of size
% dim(p) in which each element is a function of (r,s) \in RefTri.  It is
% organized by Helenbrook's cite(Helenbrook2009) organization -- first
% vertex A, B, C, then edge functions cA1,CB1, cC1, cA2, ...., and finally
% the interior functions by j(m,n) cite(Helenbrook2009).  The new basis
% functions are all defined inline to allow the eventual inline gradients
% to be defined as well as to possibly ease the visualization of any
% solutions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Psi, PsiRef, GradPsiRef] = PsiFuncs(p)

dim = (1/2)*(p+1)*(p+2);
[L, T, M] = symbolic_LandT(p);

[Phi, EvalPhi, Modes_Xi, Modes_Eta] = Dubiner_Modes(p);

syms r s xi eta; 
PhiRef = cell(dim,1);
PsiRef = cell(dim,1);
Psi = cell(dim,1);
GradPsiRef = cell(dim,2);

Xi = @(r,s) -1+2.*(1+r)./(1-s);
Eta = @(r,s) s;

for fcount = 1:dim
    PhiRef{fcount} =@(r,s) Phi{fcount}(Xi(r,s), Eta(r,s));
end

for fcount = 1:dim
    tempstr = " ";
    tempstrxi = " ";
    for mode_cnt = fcount:dim
        temp1 = char(PhiRef{mode_cnt,1}(r,s));
        tempxi = char(Phi{mode_cnt,1}(xi,eta));
        mults = find(temp1 == '*');
        multsxi = find(tempxi == '*');
        
        for k = numel(temp1(mults)):-1:1
            temp1 = [temp1(1:mults(k)-1),'.',temp1(mults(k):end)];
        end
        for k = numel(tempxi(multsxi)):-1:1
            tempxi = [tempxi(1:multsxi(k)-1),'.',tempxi(multsxi(k):end)];
        end
        
        divs = find(temp1 == '/');
        divsxi = find(tempxi == '/');
        for k = numel(temp1(divs)):-1:1
            temp1 = [temp1(1:divs(k)-1),'.',temp1(divs(k):end)];
        end
        for k = numel(tempxi(divsxi)):-1:1
            tempxi = [tempxi(1:divsxi(k)-1),'.',tempxi(divsxi(k):end)];
        end
        
        pows = find(temp1 == '^');
        powsxi = find(tempxi == '^');
        for k = numel(temp1(pows)):-1:1
            temp1 = [temp1(1:pows(k)-1),'.',temp1(pows(k):end)];
        end
        for k = numel(tempxi(powsxi)):-1:1
            tempxi = [tempxi(1:powsxi(k)-1),'.',tempxi(powsxi(k):end)];
        end
        
        temp1 = convertCharsToStrings(temp1);
        tempxi = convertCharsToStrings(tempxi);

        if T(fcount, mode_cnt) < 0
            tempstr = tempstr + string(T(fcount, mode_cnt)) + ".*(" + temp1 + ")";
            tempstrxi = tempstrxi + string(T(fcount, mode_cnt)) + ".*(" + tempxi + ")";
        end
        if T(fcount, mode_cnt) > 0
            tempstr = tempstr + "+" + string(T(fcount, mode_cnt)) + ".*(" + temp1 + ")";
            tempstrxi = tempstrxi + "+" + string(T(fcount, mode_cnt)) + ".*(" + tempxi + ")";
        end
    end
    PsiRef{fcount, 1} = str2func("@(r,s) " + tempstr  );
    Psi{fcount, 1} = str2func("@(xi,eta) " + tempstrxi  );
end

for fcount = 1:dim
    temp1 = matlabFunction(diff(PsiRef{fcount}(r,s),r));
    temp2 = matlabFunction(diff(PsiRef{fcount}(r,s),s));
    
    addr1 = find(char(temp1) == ')'); addr1 = addr1(1)+1;
    temp1 = char(temp1); temp1 = temp1(addr1:end);
    GradPsiRef{fcount,1} = str2func("@(r,s) 0.*r+(" + temp1 + ")");
    
    addr2 = find(char(temp2) == ')'); addr2 = addr2(1)+1;
    temp2 = char(temp2); temp2 = temp2(addr2:end);
    GradPsiRef{fcount,2} = str2func("@(r,s) 0.*s+(" + temp2 + ")");
end