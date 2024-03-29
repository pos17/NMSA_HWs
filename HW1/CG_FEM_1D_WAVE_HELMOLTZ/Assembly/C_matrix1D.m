function [M,A,I]=C_matrix1D(Dati,femregion)
%% [M,A] = C_matrix1D(Dati,femregion)
%==========================================================================
% Assembly of the mass matrix M and the stiffness matrix A 
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%
%    OUTPUT:
%          M           : (sparse(ndof,ndof) real) mass matrix
%          A           : (sparse(ndof,ndof) real) stiffnes matrix


addpath FESpace
addpath Assembly

fprintf('============================================================\n')
fprintf('Assembling matrices M and A ... \n');
fprintf('============================================================\n')


% connectivity infos
ndof         = femregion.ndof; % degrees of freedom
nln          = femregion.nln;  % local degrees of freedom
ne           = femregion.ne;   % number of elements
connectivity = femregion.connectivity; % connectivity matrix


% shape functions
[basis] = C_shape_basis(Dati);

% quadrature nodes and weights for integrals
[nodes_1D, w_1D] = C_quadrature(Dati);

% evaluation of shape bases 
[dphiq,Grad] = C_evalshape(basis,nodes_1D);


% Assembly begin ...
M = sparse(ndof,ndof);  % Global Mass matrix
A = sparse(ndof,ndof);  % Global Stiffness matrix
I = sparse(ndof,ndof);  % Global Stiffness matrix
for ie = 1 : ne
     
    % Local to global map --> To be used in the assembly phase
    iglo = connectivity(1:nln,ie);
  
    
    [BJ, pphys_1D] = C_get_Jacobian(femregion.coord(iglo,:), nodes_1D);
    % BJ        = Jacobian of the elemental map 
    % pphys_2D = vertex coordinates in the physical domain 
   
    %=============================================================%
    % STIFFNESS MATRIX
    %=============================================================%
    
    x = pphys_1D;
    pphys_1D = eval(Dati.shape)
    
    % Local stiffness matrix 
    [A_loc] = C_lap_loc(Grad,w_1D,nln,BJ,pphys_1D);

    % Assembly phase for stiffness matrix
    A(iglo,iglo) = A(iglo,iglo) + A_loc; 

    %=============================================================%
    % MASS MATRIX
    %=============================================================%
    
    % Local mass matrix 
    [M_loc] = C_mass_loc(dphiq,w_1D,nln,BJ,pphys_1D);

    % Assembly phase for mass matrix
    M(iglo,iglo) = M(iglo,iglo) + M_loc;   
   

end
if(strcmp(Dati.bc,'LR'))
    x = Dati.domain(1);
    S0 = eval(Dati.shape);
    x = Dati.domain(2);
    S1 = eval(Dati.shape);
    A(end,end) = A(end,end) + Dati.c2/Dati.domain(2) * S1*(Dati.domain(2)/(0.8216*sqrt(S0*S1/pi)))
    I(end,end) = Dati.c2/Dati.domain(2) * S1*   1/(2*0.8216^2*sqrt(Dati.c2/(Dati.domain(2))))
end
