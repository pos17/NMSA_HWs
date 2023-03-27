function [errors,solutions,femregion,Dati,u_snp] = C_main1D(TestName,nRef)
%==========================================================================
% Solution of the Wave Equation with linear finite elements 
% coupled with the leap-frog scheme 
% 
%  u_tt - c^2 u_xx = f  in (a,b)x(0,T)
%  u_t(x,0) = v_0       in (a,b) 
%  u(x,0)   = u_0       in (a,b)
%  + boundary conditions:
%  Dirichlet:  u(s,t) = g  s=a,b
%  Neumann  : c^2du/dn(s,t) = g s=a,b
%  Periodic : u(a,t) = a(b,t)
%  Absorbing: du/dt(s,t) + cdu/dn(s,t) = 0  s=a,b 
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
gamma = sqrt(Dati.c2/Dati.domain(2));

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

% [M_nbc,A_nbc,I_nbc] = C_matrix1D(Dati,femregion);
[M_nbc, A_nbc, ~] = C_matrix1D(Dati,femregion);

if(strcmp(Dati.bc,'R'))
    A_nbc(1,1)     = A_nbc(1,1) + Dati.c2;
    A_nbc(end,end) = A_nbc(end,end) + Dati.c2;
end

%==========================================================================
% BUILD FINITE ELEMENTS RHS a time 0
%==========================================================================
Dati.t = 0;
[b_nbc] = C_rhs1D(Dati,femregion);
t = 0;

F_1 = b_nbc;

if(strcmp(Dati.bc,'N'))
    b_nbc(1)   = b_nbc(1)   + eval(Dati.neumann1);
    b_nbc(end) = b_nbc(end) + eval(Dati.neumann2);
elseif(strcmp(Dati.bc,'R'))
    b_nbc(1)   = b_nbc(1)   - Dati.c2*eval(Dati.neumann1);
    b_nbc(end) = b_nbc(end) + Dati.c2*eval(Dati.neumann2);
% elseif(strcmp(Dati.bc,'M'))
%     b_nbc(1)   = b_nbc(1)   - M_nbc(1,1)*eval(Dati.Lifting1tt) - A_nbc(1,1)*eval(Dati.Lifting1);
end
% omega_range = eval(Dati.omega);
% omega_linspace = linspace(omega_range(1),omega_range(2));
% Solve the reduced system
% K_nbc =  -A_nbc + omega_linspace.^2/Dati.c2 .* M_nbc;


%% eigenmodes(resonant frequencies) and eingevectors
if (strcmp(Dati.bc,'NN'))
    M = M_nbc;
    A = A_nbc;
elseif (strcmp(Dati.bc,'P'))
%     [M,~,~] = C_bound_cond1D(M_nbc,b_nbc,femregion,Dati)
%     [A,~,~] = C_bound_cond1D(A_nbc,b_nbc,femregion,Dati);
    mat = full(M_nbc)\full(A_nbc);
    [eig_mat,~,~] = C_bound_cond1D(mat,b_nbc,femregion,Dati);
end

[V, D] = eigs(eig_mat);

w =sort(sqrt(diag(D)));
% omega_h = sqrt(k2*Dati.c2);
% modes = omega_h
modes = w



M_mod = V' * M_nbc * V;
A_mod = V' * A_nbc * V;

omega_range = eval(Dati.omega);
omega_linsp = linspace(omega_range(1), omega_range(end), 10000);
H = zeros(1,length(omega_linsp));
for i= 1:length(modes)
% for i= 1
    H = H + 1./(-omega_linsp.^2* M_mod(i,i) + A_mod(i,i));
    disp(i)
end

%figure(13)
semilogy(omega_linsp, abs(H));
xlim(omega_range);
    %view(0,90)
    %title('Space time evolution')
    %saveas(gcf,strcat("Plots/",Dati.name,"_SpaceTimeEvolution.png"))

fprintf("======================================================== \n" + ...
        "                     Modes                               \n" + ...
        "======================================================== \n")
disp(modes)

fprintf("======================================================== \n" + ...
        "                    Modes 2                              \n" + ...
        "======================================================== \n")
disp(D)



%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

% [solutions] = C_postprocessing(Dati,femregion,u2);

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
errors = [];
if (Dati.plot_errors)
%     [errors] = C_compute_errors(Dati,femregion,solutions);
end

if(strcmp(Dati.st_plot, "Y"))
    figure(10)
    surf(u_snp, EdgeColor="none");
    view(0,90)
    title('Space time evolution')
    saveas(gcf,strcat("Plots/",Dati.name,"_SpaceTimeEvolution.png"))

    figure(11)
    surf(p_snp, EdgeColor="none");
    view(0,90)
    title('Pressure Field evolution')
    saveas(gcf,strcat("Plots/",Dati.name,"_PresFieldEvol.png"))
    
end 



