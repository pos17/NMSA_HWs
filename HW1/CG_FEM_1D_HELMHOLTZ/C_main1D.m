function [errors,solutions,femregion,Dati] = C_main1D(TestName,nRef)
%==========================================================================
% Solution of the Wave Equation with linear finite elements
% (non homogeneous Dirichlet boundary conditions)
%==========================================================================
%
%    INPUT:
%          TestName    : (string)  see C_dati.m
%          nRef        : (int)     refinement level
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) infos about finite elements
%                        discretization
%          Dati        : (struct)  see C_dati.m
%
% Usage:
%    [errors,solutions,femregion,Dati] = C_main1D('Test1',3)



addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing


%==========================================================================
% LOAD DATA FOR TEST CASE
%==========================================================================

Dati = C_dati(TestName);
Dati.nRefinement = nRef;

%==========================================================================
% MESH GENERATION
%==========================================================================

[Region] = C_create_mesh(Dati);

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[femregion] = C_create_femregion(Dati,Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES
%==========================================================================

[M_nbc,A_nbc] = C_matrix1D(Dati,femregion);

%==========================================================================
% BUILD FINITE ELEMENTS RHS a time 0
%==========================================================================
[b_nbc] = C_rhs1D(Dati,femregion);


rho = Dati.rho;
omega = Dati.omega;
vel = Dati.vel;

t = 0;
if (strcmp(Dati.bc,'NN'))
    b_nbc(1)   = b_nbc(1)   + eval(Dati.neumann1);
    b_nbc(end) = b_nbc(end) + eval(Dati.neumann2);
elseif (strcmp(Dati.bc,'ND'))
    b_nbc(end) = b_nbc(end) + eval(Dati.neumann1);    
elseif (strcmp(Dati.bc,'DN'))
    b_nbc(end) = b_nbc(end) + eval(Dati.neumann2);
end


% Solve the reduced system
K_nbc =  -A_nbc + Dati.omega^2/vel^2 * M_nbc;


%% eigenmodes(resonant frequencies) and eingevectors
if (strcmp(Dati.bc,'NN'))
    M = M_nbc;
    A = A_nbc;
else
    [M,~,~] = C_bound_cond1D(M_nbc,b_nbc,femregion,Dati);
    [A,~,~] = C_bound_cond1D(A_nbc,b_nbc,femregion,Dati);
end

% eigenvalues D -- eigenvector V
% A x = k^2 M x
[V,D] = eigs(A,M,40,0.1);

k2 = diag(D);
omega_h = sqrt(k2*Dati.vel^2);
modes = omega_h;
%% plot of 4 eigenvectors  
if (strcmp(Dati.plot_eigvct,'Y'))
    for i = 1 : 6
        v = V(:,i);
        disp(modes);
        strTitle = ['eigenvector n ', num2str(i)];
        [v] = C_snapshot_1D(femregion,v,10+i,strTitle);
    end
end

fprintf("==================================================== \n" + ...
     "                     Modes                           \n"    + ...
     "==================================================== \n")
disp(modes)

%==========================================================================
% COMPUTE BOUNDARY CONDITIONS -- MODIFICATION OF M an b
%==========================================================================

if (strcmp(Dati.bc,'NN'))
    p = K_nbc\b_nbc;
else
    [M,b,u_g] = C_bound_cond1D(K_nbc,b_nbc,femregion,Dati);
    p = M\b;
    p = p + u_g;
end


%[p] = C_snapshot_1D(femregion,p,1);

%% computation of the displacement
n = length(p);
v(2:n-1,1) = (p(3:n,1)-p(1:n-2,1))/(2*femregion.h);
v(1) = eval(Dati.neumann1); 
v(n) = eval(Dati.neumann2); 
v = 1/(rho*omega).^2 * v;

%% plot of displacement
[v] = C_snapshot_1D(femregion,v,2,'displ.');


%% impedance in x=0
Z_imp = p(1)/eval(Dati.neumann1);

%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

[solutions] = C_postprocessing(Dati,femregion,p);

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
errors = [];
if (Dati.plot_errors)
     [errors] = C_compute_errors(Dati,femregion,solutions);
end



