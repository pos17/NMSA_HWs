% Solution and error of the transport equation
%
%    u_t + a*u_x = 0
%
% x \in I=[a,b]  and  t \in [0,T] with the
% initial data  u(x,0) = u_0(x) with different methods
%   Euler forward for time / central for space
%   Lax-Friedrichs
%   Lax-Wendroff
%   Upwind


clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial data
u0 = @(x) sin(2*pi*x);

% velocity
a = 1;

% exact solution for computing the error
u = @(x,t) sin(2*pi*(x-a*t));


% Intervall
T = 1;
I = [0 1];

% Number of time steps
NT = 1000;
% Number of space steps
NX = 200;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% menu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag = menu('Please select a method',...
%     'Forward-Euler/Central',...
%     'Lax-Friedrichs',...
%     'Lax-Wendroff',...
%     'Upwind');
%
% disp('Please wait ...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Euler forward for time / central for space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = T/NT;
dx = (I(2)-I(1))/NX;
lambda = dt/dx;

% Initial conditions
SOL = zeros(NX+1,NT+1);
for j = 1:NX+1
    SOL(j,1) = u0(I(1) + (j-1)*dx);
end

% Finite Difference-Scheme: Euler-forward in time / central in space
%     figure;
for n = 1:NT
    for j = 2:NX
        SOL(j,n+1) = SOL(j,n)- 0.5*lambda*a*(SOL(j+1,n)-SOL(j-1,n));
    end
    SOL(1,n+1) = u(I(1),n*dt);
    %periodic condition
    SOL(end,n+1) = SOL(1,n+1);
end

% Computation of the L2-error for t=T
SOL_EX = zeros(NX+1,1);
for j = 1:NX+1
    SOL_EX(j,1) = u(I(1) + (j-1)*dx,NT*dt);
end
L2_ERR   = norm(SOL_EX-SOL(:,n+1),2);
L1_ERR   = norm(SOL_EX-SOL(:,n+1),1);
LINF_ERR = norm(SOL_EX-SOL(:,n+1),Inf);



% Visualization
time = ones(NX+1,1)*linspace(0,T,NT+1);
space = linspace(I(2),I(1),NX+1)'*ones(1,NT+1);
surf(time,space,SOL,'EdgeColor','none');
xlabel('time t','FontSize',16);
ylabel('space x','FontSize',16);






