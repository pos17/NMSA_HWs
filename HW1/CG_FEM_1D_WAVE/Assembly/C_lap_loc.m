function [K_loc] = C_lap_loc(Grad,w_1D,nln,BJ,shape_1D)
%% [K_loc] = C_lap_loc(Grad,w_1D,nln,BJ)
%==========================================================================
% Build the local stiffness matrix for the term grad(u)grad(v)
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          Grad        : (array real) evaluation of the gradient on
%                        quadrature nodes
%          w_1D        : (array real) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : (array real) Jacobian of the map 
%          phys_1D     : physical coordinates of the nodes
%
%    OUTPUT:
%          K_loc       :  (array real) Local stiffness matrix


K_loc = zeros(nln,nln);

%% General implementation -- to be used with general finite element spaces
for i=1:nln
    for j=1:nln
        for k=1:length(w_1D)
            Binv = 1./BJ;    % inverse
            Jdet = BJ;       % determinant 
            %disp(x)
            K_loc(i,j) = K_loc(i,j) + (Jdet.*w_1D(k)) .* ( shape_1D(k) * (Grad(k,:,i) * Binv) * (Grad(k,:,j) * Binv )');
        end
    end
end



                                              
                                              

