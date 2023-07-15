function [femregion, Dati] = C_main1D_Ste(TestName,nRef)

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

[femregion] = C_create_femregion(Dati, Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES
%==========================================================================

[M_nbc, A_nbc, ~] = C_matrix1D(Dati,femregion);

% A = C_bound_cond1D(A_nbc, femregion, Dati);

mat = full(M_nbc)\full(A_nbc);
mat_eig = C_bound_cond1D(mat, femregion, Dati);

% [V,D] = eig(mat_eig);
[V,D] = eig(mat_eig);

k2 = sort(diag(D));
% w_r = sqrt(Dati.c2*k2);
w_r=sqrt(k2);
fprintf("=================================================\n" + ...
        "                    omega res                    \n" + ...
        "=================================================\n");
disp(w_r);

M_m = V' * M_nbc * V;
A_m = V' * A_nbc * V;

omega_range = eval(Dati.omega);
% omega_linsp = linspace(omega_range(1)/Dati.vel, omega_range(2)/Dati.vel, 1000);
omega_linsp = linspace(0, omega_range(2), 1000);
sum = zeros(1, length(omega_linsp));
for i= 1:length(w_r)
    H = 1./(-omega_linsp.^2* M_m(i,i) + A_m(i,i));
    sum = sum+H;
%     disp(i)
end



% semilogy(omega_linsp*Dati.vel, abs(sum))
semilogy(omega_linsp, abs(sum), LineWidth=3)
xlim(omega_range)
grid minor
xlabel('\omega', FontSize=16)
ylabel('|\psi|', FontSize=16)
disp("fine")






