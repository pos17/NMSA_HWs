function C_pointwise_sol(femregion, uh, u_ex, Dati)
%% C_pointwise_sol(femregion, uh, u_ex,Dati)
%==========================================================================
% PLOT THE EXACT SOLUTION ON THE DOFS
%==========================================================================
%    called in C_postprocessing.m
%
%    INPUT:
%          femregion   : (struct)  see C_create_femregion.m
%          uh          : (sparse(ndof,1) real) solution vector
%          u_ex        : (sparse(ndof,1) real) exact solution vector
%          Dati        : (struct) see C_Dati.m
%



dof=femregion.dof;

x1 = femregion.domain(1,1);
x2 = femregion.domain(1,2);


M=max(uh);
m=min(uh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PLOT OF SOLUTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(femregion.coord(:,1),full(uh));
title('u_h(x)'); xlabel('x-axis'); ylabel('y-axis');
axis([x1,x2,m,M]); 

if(min(u_ex) ~= 0 && max(u_ex) ~= 0)
    hold on
    stem(femregion.coord(:,1),u_ex,'-r*');
    legend('computed', 'exact')
end

