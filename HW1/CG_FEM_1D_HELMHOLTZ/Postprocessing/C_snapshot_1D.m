function [uh] = C_snapshot_1D(femregion, uh, figID, figTitle)
%% C_snapshot_1D(femregion, uh, Dati)
%==========================================================================
% PLOT THE EXACT SOLUTION ON THE DOFS
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          femregion   : (struct)  see C_create_femregion.m
%          uh          : (sparse(ndof,1) real) solution vector
%          Dati        : (struct) see C_Dati.m
%


x1 = femregion.domain(1,1);
x2 = femregion.domain(1,2);


% M =  1; % max(uh)
% m = -1; % min(uh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% PLOT OF SOLUTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figID);
plot(femregion.coord(:,1),full(uh));
title(figTitle); 
xlabel('x-axis'); ylabel('y-axis');
% axis([x1,x2,m,M]); 
