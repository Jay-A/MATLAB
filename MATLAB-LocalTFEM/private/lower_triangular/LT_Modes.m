% a file that creates a cell array of Dubiner modes for p
% Defined over the reference triangle with a singularity at vertex 1
% therefore Phi{i}:(r,s)->R with  -1 <= s < 1 and -1 <= r <= s U r=1
function [Psi] = LT_Modes(p)

[Phi, EvalPhi, Modes_Xi, Modes_Eta] = Dubiner_Modes(p);
