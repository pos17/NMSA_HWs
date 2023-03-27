function [u_ex]=C_eval_exact_sol(femregion, Dati)
%% [u_ex]=C_eval_exact_sol(femregion, exact_sol)
%==========================================================================
% EVALUATION OF THE EXACT SOLUTION ON THE GRID POINTS
%==========================================================================
%    called in C_postprocessing.m
%
%    INPUT:
%          femregion   : (struct)  see C_create_femregion.m
%          exact_sol   : (string)  expression of exact solution
%          t           : (real)    time
%
%    OUTPUT:
%          u_ex        : (spare(ndof,1) real) exact solution vector
%

rho = Dati.rho;
omega = Dati.omega;
vel = Dati.vel;


dof = femregion.dof;
x = dof(:,1);

u_ex = eval(Dati.exact_sol);

    

