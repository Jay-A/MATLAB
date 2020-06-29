function [RHS] = Eval_RHSLT(func, Verts, Area, p, Phi, T)
    
    [RHS] = Eval_RHS(func, Verts, Area, p, Phi);
  
  RHS = double(T*RHS);