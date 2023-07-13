% function [L2_error] = hyp_HW2(NX, NT)
clear; close all;% clc;

% problem
% u_tt - c^2 u_xx = 0
% x \in I=[0,L]  and  t \in [0,T] with the
% initial data  u(x,0) = u0(x) and v(x,0) = v0(x)
% and left Dirichlet b.c. u(0,t) = g1(t)
% and right Neumann b.c. c^2*u_x(L,t) = g2(t)

% the problem develops around the variable transformation
% w1 = u_t
% w2 = u_x
% w = [w1, w2]'
% r = T^(-1) w = [r1, r2]'


% exact solution
% u = @(x,t) sin(2*pi*x)*cos(2*pi*t);
% ut = @(x,t) -2*pi*sin(2*pi*x)*sin(2*pi*t);
% ux = @(x, t) 2*pi*cos(2*pi*x)*cos(2*pi*t);

% velocity and eigenvalues
c = 2;
a1 = c;
a2 = -c;

% initial condition for u_t and u_x
v0 = @(x) 0.*x;     % w1(x,0)
ux0 = @(x) 0.*x;    % w2(x,0)

% initial condition for r1 and r2
r1_init = @(x) sqrt(c^2 + 1)/(2*c) .* (v0(x) + c*ux0(x));
r2_init = @(x) sqrt(c^2 + 1)/(2*c) .* (-v0(x) + c*ux0(x));

% Intervals
T = 1;
I = [0, 1];

% boundary conditions
gp1 = @(t) 10*pi * sin (2 * pi * (1 - 10 * t) );
g2 = @(t) 0.*t;
% g2 = @(t) c^2 * 2*pi*cos(2*pi*t);

% Number of time steps
NT = 10000;
% Number of space steps
NX = 256;

dt = T/NT;
dx = (I(2)-I(1))/NX;
lambda = dt/dx;

SOL_r1 = zeros(NX+1, NT+1);
SOL_r2 = zeros(NX+1, NT+1);
SOL_w1 = zeros(NX+1, NT+1);
SOL_w2 = zeros(NX+1, NT+1);

% setting initial conditions on solutions
for j = 1 : NX+1
    SOL_r1(j,1)      = r1_init(I(1) + (j-1)*dx);
    SOL_r2(j,1)      = r2_init(I(1) + (j-1)*dx);
    SOL_w1(j,1)      = v0(I(1) + (j-1)*dx);
    SOL_w2(j,1)      = ux0(I(1) + (j-1)*dx);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% menu:
%       flag: 1 => LAX-WENDROFF SCHEME
%       flag: 2 => LAX-FRIEDRICHS SCHEME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag = 1;


% applying finite difference scheme
for k = 1:NT
    for j = 2:NX % takes only internal nodes, no boundaries

        if flag == 1 % LAX-WENDROFF SCHEME
            SOL_r1(j, k+1) = SOL_r1(j,k) - 0.5*lambda*a1*(SOL_r1(j+1,k)-SOL_r1(j-1,k))...
                + 0.5*lambda^2*a1^2*(SOL_r1(j+1,k)-2*SOL_r1(j,k)+SOL_r1(j-1,k));
            SOL_r2(j, k+1) = SOL_r2(j,k) - 0.5*lambda*a2*(SOL_r2(j+1,k)-SOL_r2(j-1,k))...
                + 0.5*lambda^2*a2^2*(SOL_r2(j+1,k)-2*SOL_r2(j,k)+SOL_r2(j-1,k));
        elseif flag == 2 % LAX-FRIEDRICHS SCHEME
            SOL_r1(j, k+1) = 0.5*(SOL_r1(j+1,k)+SOL_r1(j-1,k)) - 0.5*lambda*a1*(SOL_r1(j+1,k)-SOL_r1(j-1,k));
            SOL_r2(j, k+1) = 0.5*(SOL_r2(j+1,k)+SOL_r2(j-1,k)) - 0.5*lambda*a2*(SOL_r2(j+1,k)-SOL_r2(j-1,k));

        else
            disp ("flag not valid")
        end

    end

    if abs((k+1)*dt - 0.1) <= 0.05
        ggp1 = gp1((k+1)*dt);
    else
        ggp1=0;
    end

    % fixing boundary conditions
    % x=0
    SOL_r2(1, k+1) = SOL_r2(2, k+1); % constant extrapolation
%     SOL_r2(1, k+1) = 2*SOL_r2(2, k+1)-SOL_r2(3, k+1); % linear extrapolation
    SOL_r1(1, k+1) = sqrt(c^2+1)/c*ggp1 + SOL_r2(1, k+1);
    
    % x=L
    SOL_r1(NX+1, k+1) = SOL_r1(NX, k+1); % constant extrapolation
%     SOL_r1(NX+1, k+1) = 2*SOL_r1(NX, k+1) - SOL_r1(NX-1, k+1); % linear extrapolation
    SOL_r2(NX+1, k+1) = sqrt(c^2+1)/c^2 * g2((k+1)*dt) - SOL_r1(NX+1, k+1);

   
    % compute w using transformation w = T r
    SOL_w1(:, k+1) = c/sqrt(1+c^2) * (SOL_r1(:, k+1) - SOL_r2(:, k+1));
    SOL_w2(:, k+1) = 1/sqrt(1+c^2) * (SOL_r1(:, k+1) + SOL_r2(:, k+1));
    
    
end

%% VISUALIZATION 
time = linspace(0,T,NT+1); % ones(NX+1,1)*
space = linspace(I(1),I(2),NX+1);%'*ones(1,NT+1);

figure
% surf(time, space, SOL_w1,'EdgeColor','none');
surf(space, time, SOL_w1','EdgeColor','none');
title('u_t(x,t)','FontSize',16);
colorbar
caxis([-60,60])
ylabel('time t','FontSize',16);
xlabel('space x','FontSize',16);
view(2);
filename=["ut_LW_forced"];
saveas(gcf, "./Plots/"+filename+".png")

figure
% surf(time, space, -SOL_w2,'EdgeColor','none');
surf(space,time, -SOL_w2','EdgeColor','none');
colorbar
caxis([-30,30])
title('u_x(x,t)','FontSize',16);
ylabel('time t','FontSize',16);
xlabel('space x','FontSize',16);
view(2);
filename=["ux_LW_forced"];
saveas(gcf, "./Plots/"+filename+".png")

% CFL_cond = lambda*a1
% close all


%% error computing

% % Computation of errors for t=T
% SOLu_EX  = zeros(NX+1,1);
% SOLut_EX = zeros(NX+1,1);
% SOLux_EX = zeros(NX+1,1);
% 
% for j = 1:NX+1
%     SOLu_EX(j,1)  = u(I(1) + (j-1)*dx, NT*dt);
%     SOLut_EX(j,1) = ut(I(1) + (j-1)*dx, NT*dt);
%     SOLux_EX(j,1) = ux(I(1) + (j-1)*dx, NT*dt);
% end
% 
% 
% % for ut=w1
% L2ut_ERR   = norm(SOLut_EX-SOL_w1(:,k+1),2)*dx^0.5;
% % L1ut_ERR   = norm(SOLut_EX-SOL_w2(:,k+1),1)*dx;
% % LINFut_ERR = norm(SOLut_EX-SOL_w2(:,k+1),Inf);
% 
% L2ux_ERR   = norm(SOLux_EX-SOL_w2(:,k+1),2)*dx^0.5;
% % L1ux_ERR   = norm(SOLux_EX-SOL_w2(:,k+1),1)*dx;
% % LINFux_ERR = norm(SOLux_EX-SOL_w2(:,k+1),Inf);
% % 
% % Err_t = [L2ut_ERR, L1ut_ERR, LINFut_ERR]
% % Err_x = [L2ux_ERR, L1ux_ERR, LINFux_ERR]
% 
% L2_error = [L2ut_ERR, L2ux_ERR];







