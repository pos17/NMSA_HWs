clear; close all; clc;

% the problem develops around the variable transformation
% w1 = u_t
% w2 = u_x
% w = [w1, w2]'
% r = T^(-1) w = [r1, r2]'


% exact solution
u = @(x,t) sin(2*pi*x)*cos(2*pi*t);
ut = @(x,t) -2*pi*sin(2*pi*x)*sin(2*pi*t);
ux = @(x, t) 2*pi*cos(2*pi*x)*cos(2*pi*t);

% velocity
c = 2;

% transport matrix
A = [0, c^2; 1, 0];

% eigenvalues and eigenvectors
[T_d, D] = eig(A);
T_inv = inv(T_d);
a = diag(D);
a1 = a(1);
a2 = a(2);

% initial condition for u_t and u_x
v0 = @(x) ut(x, 0);     % w1(x,0)
ux0 = @(x) ux(x, 0);    % w2(x,0)

% initial condition for r1 and r2
r1_init = @(x) sqrt(c^2 + 1)/(2*c) * (v0(x) + c*ux0(x));
r2_init = @(x) sqrt(c^2 + 1)/(2*c) * (-v0(x) + c*ux0(x));

% Intervals
T = 1;
I = [0, 1];

% boundary conditions
gp1 = @(t) 0.*t;
g2 = @(t) c^2 * ux(I(2), t);
% g2 = @(t) 0.*t;

% Number of time steps
NT = 2000;
% Number of space steps
NX = 1000;

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
            SOL_r1(j, k+1) = 1/2*(SOL_r1(j+1,k)+SOL_r1(j-1,k)) - lambda/2*a1*(SOL_r1(j+1,k)-SOL_r1(j-1,k));
            SOL_r2(j, k+1) = 1/2*(SOL_r2(j+1,k)+SOL_r2(j-1,k)) - lambda/2*a2*(SOL_r2(j+1,k)-SOL_r2(j-1,k));

        else
            disp ("flag not valid")
        end

    end

    % fixing boundary conditions
    % x=0
    SOL_r2(1, k+1) = SOL_r2(2, k+1); % constant extrapolation
    SOL_r1(1, k+1) = sqrt(c^2+1)/c*gp1((k+1)*dt) + SOL_r2(1, k+1);
    
    % x=L
    SOL_r1(NX+1, k+1) = SOL_r1(NX, k+1); % constant extrapolation
    SOL_r2(NX+1, k+1) = sqrt(c^2+1)/c^2 * g2((k+1)*dt) - SOL_r1(NX+1, k+1);
%     SOL_r2(NX+1, k+1) = SOL_r1(NX+1, k+1) - sqrt(c^2+1)/c * g2((k+1)*dt);
%     %sbagliato 
    

    % compute w using transformation w = T r
    SOL_w1(:, k+1) = c/sqrt(1+c^2) * (SOL_r1(:, k+1) - SOL_r2(:, k+1));
    SOL_w2(:, k+1) = 1/sqrt(1+c^2) * (SOL_r1(:, k+1) + SOL_r2(:, k+1));
    
    
end

%% VISUALIZATION 
time = linspace(0,T,NT+1); % ones(NX+1,1)*
space = linspace(I(2),I(1),NX+1);%'*ones(1,NT+1);

figure
surf(time, space, SOL_w1,'EdgeColor','none');
title('u_t','FontSize',16);
xlabel('time t','FontSize',16);
ylabel('space x','FontSize',16);
view(2);

figure
surf(time, space, SOL_w2,'EdgeColor','none');
title('u_x','FontSize',16);
xlabel('time t','FontSize',16);
ylabel('space x','FontSize',16);
view(2);

%%
% gg2 = zeros(1,length(time));
% gg1 = zeros(1, length(time));
% for i=1:length(time)
%     gg2(i) = g2((i-1)*dt);
%     gg1(i) = gp1((i-1)*dt);
% end
% figure
% plot(time, gg2)
% 
% figure
% plot(time, gg1)




