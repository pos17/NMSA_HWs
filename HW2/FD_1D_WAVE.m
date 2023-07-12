% Solution and error of the wave equation by finite difference method
%
%    u_tt - c^2 u_xx = 0
%
% x \in I=[0,L]  and  t \in [0,T] with the
% initial data  u(x,0) = u0(x) and v(x,0) = v0(x)
% and Dirichlet boundary conditions u(0,t) = g1(t)
% and u(L,t) = g2(t)
%   
%   Finite difference methods:
%   Euler forward for time / central for space
%   Lax-Friedrichs
%   Lax-Wendroff
%   Upwind


clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exact solution for computing the error
u  = @(x,t) sin(2*pi*x)*sin(2*pi*t);
ut = @(x,t) 2*pi*sin(2*pi*x)*cos(2*pi*t);
ux = @(x,t) 2*pi*cos(2*pi*x)*sin(2*pi*t);


% velocity
c = 1;
a1 =  c;
a2 = -c;

% Initial data
u0  = @(x) 0.*x;

% Initial conditions for w_t = A w_x
v0  = @(x) 2*pi*sin(2*pi*x); % for system of eq
up0 = @(x) 0.*x;             % for system of eq

% Initial conditions for w~_t = D w~_x
% by using the transformation w~(x,0) = T^(-1) w(x,0)

w0 = @(x) sqrt(1+c^2)/(2*c).*( 2*pi*sin(2*pi*x) + c*0.*x);
w1 = @(x) sqrt(1+c^2)/(2*c).*(-2*pi*sin(2*pi*x) + c*0.*x);

% Intervals
T = 1;
I = [0 1];

%boundary conditions
gp1 = @(t) 0.*t;
gp2 = @(t) 0.*t;
% gp2 = @(t) c^2*ux(I(2), t);


% Number of time steps
NT = 3000;
% Number of space steps
NX = 1000;


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
flag = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = T/NT;
dx = (I(2)-I(1))/NX;
lambda = dt/dx;

% Initial conditions
% w~ = [w~_1, w~_2]'; 
% u~ = [u_1, u_2]';
% SOL_w(1:NX+1,:) -> w~_1; % SOL_w(NX+2:end,:) -> w~_2;
% the same for SOL_u
SOL_w = zeros(2*(NX+1),NT+1);
SOL_u = zeros(2*(NX+1),NT+1);


for j = 1 : NX+1
    SOL_w(j,1)      = w0(I(1) + (j-1)*dx);
    SOL_w(NX+1+j,1) = w1(I(1) + (j-1)*dx);
    SOL_u(j,1)      = v0(I(1) + (j-1)*dx);
    SOL_u(NX+1+j,1) = up0(I(1) + (j-1)*dx);

end


for n = 1:NT
    for j = 2:NX
        
        if flag == 1 % FCE scheme
            % w~_1
            SOL_w(j,n+1)      = SOL_w(j,n)      - 0.5*lambda*a1*(SOL_w(j+1,n)-SOL_w(j-1,n));
            % w~_2
            SOL_w(NX+1+j,n+1) = SOL_w(NX+1+j,n) - 0.5*lambda*a2*(SOL_w(NX+1+j+1,n)-SOL_w(NX+1+j-1,n));
        elseif flag == 2 % LF scheme
            SOL_w(j,n+1)      = 0.5*(SOL_w(j+1,n)+SOL_w(j-1,n))           - 0.5*lambda*a1*(SOL_w(j+1,n)-SOL_w(j-1,n));
            SOL_w(NX+1+j,n+1) = 0.5*(SOL_w(NX+1+j+1,n)+SOL_w(NX+1+j-1,n)) - 0.5*lambda*a2*(SOL_w(NX+1+j+1,n)-SOL_w(NX+1+j-1,n));
        elseif flag == 3 % LW scheme
            SOL_w(j,n+1) = SOL_w(j,n) - 0.5*lambda*a1*(SOL_w(j+1,n)-SOL_w(j-1,n)) ...
                              + 0.5*lambda^2*a1^2*(SOL_w(j+1,n)-2*SOL_w(j,n)+SOL_w(j-1,n));
            SOL_w(NX+1+j,n+1) = SOL_w(NX+1+j,n) - 0.5*lambda*a2*(SOL_w(NX+1+j+1,n)-SOL_w(NX+1+j-1,n)) ...
                              + 0.5*lambda^2*a2^2*(SOL_w(NX+1+j+1,n)-2*SOL_w(NX+1+j,n)+SOL_w(NX+1+j-1,n));
        elseif flag == 4 % UPWIND
            SOL_w(j,n+1) = SOL_w(j,n) - 0.5*lambda*a1*(SOL_w(j+1,n)-SOL_w(j-1,n)) ...
                              + 0.5*lambda*abs(a1)*(SOL_w(j+1,n)-2*SOL_w(j,n)+SOL_w(j-1,n));
            SOL_w(j,n+1) = SOL_w(j,n) - 0.5*lambda*a2*(SOL_w(NX+1+j+1,n)-SOL_w(NX+1+j-1,n)) ...
                              + 0.5*lambda*abs(a2)*(SOL_w(NX+1+j+1,n)-2*SOL_w(NX+1+j,n)+SOL_w(NX+1+j-1,n));
                          
        end
    end
    
    % Boundary conditions
    % can be improved by adding external end points --> SOL_w(1:2*(NX+1)+4,:)
    SOL_w(NX+2,n+1) = SOL_w(NX+3,n+1); %constant extrapolation w~2(0,t) = w~2(h,t)
    % update w~_1(0,t)
    SOL_w(1,n+1) = sqrt(1+c^2)/c*gp1((n+1)*dt) + SOL_w(NX+2,n+1);


    SOL_w(NX+1,n+1) = SOL_w(NX,n+1);  %constant extrapolation w~1(L,t) = w~1(L-h,t)
    % update w~_2(L,t)
    SOL_w(end,n+1)  = SOL_w(NX+1,n+1)-sqrt(1+c^2)/c*gp2((n+1)*dt);
%     SOL_w(end,n+1)  = sqrt(1+c^2)/c^2*gp2((n+1)*dt)-SOL_w(NX+1,n+1);

    
    % compute u by using the transformation w = T w~
    % w = [w_1, w_2] = [u_t, u_x];
    SOL_u(1:NX+1,n+1)   = c/sqrt(1+c^2)*(SOL_w(1:NX+1,n+1)-SOL_w(NX+2:end,n+1));
    SOL_u(NX+2:end,n+1) = 1/sqrt(1+c^2)*(SOL_w(1:NX+1,n+1)+SOL_w(NX+2:end,n+1));
end

% Computation of the L2-error for t=T
SOLu_EX  = zeros(NX+1,1);
SOLut_EX = zeros(NX+1,1);
SOLux_EX = zeros(NX+1,1);

for j = 1:NX+1
    SOLu_EX(j,1)  = u(I(1) + (j-1)*dx,NT*dt);
    SOLut_EX(j,1) = ut(I(1) + (j-1)*dx,NT*dt);
    SOLux_EX(j,1) = ux(I(1) + (j-1)*dx,NT*dt);

end
L2ut_ERR   = norm(SOLut_EX-SOL_u(1:NX+1,n+1),2)*dx^0.5;
L1ut_ERR   = norm(SOLut_EX-SOL_u(1:NX+1,n+1),1)*dx;
LINFut_ERR = norm(SOLut_EX-SOL_u(1:NX+1,n+1),Inf);

L2ux_ERR   = norm(SOLux_EX-SOL_u(NX+2:end,n+1),2)*dx^0.5;
L1ux_ERR   = norm(SOLux_EX-SOL_u(NX+2:end,n+1),1)*dx;
LINFux_ERR = norm(SOLux_EX-SOL_u(NX+2:end,n+1),Inf);

Err_t = [L2ut_ERR, L1ut_ERR, LINFut_ERR]
Err_x = [L2ux_ERR, L1ux_ERR, LINFux_ERR]


% Visualization
time = ones(NX+1,1)*linspace(0,T,NT+1);
space = linspace(I(2),I(1),NX+1)'*ones(1,NT+1);
surf(time, space, SOL_u(1:NX+1,:),'EdgeColor','none');
title('u_t','FontSize',16);
xlabel('time t','FontSize',16);
ylabel('space x','FontSize',16);
view(2);

figure;
surf(time, space, SOL_u(NX+2:end,:),'EdgeColor','none');
title('u_x','FontSize',16);
xlabel('time t','FontSize',16);
ylabel('space x','FontSize',16);
view(2);





