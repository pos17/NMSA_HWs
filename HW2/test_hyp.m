
clear; close all; clc;

NXX = [1000, 2000, 3000, 4000];
NTT = [2000, 4000, 6000, 8000];

T=1;
L=1;
dx = L./NXX;
dt=T./NTT;

err = zeros(3,2);

for ii=1:length(NXX)
    err(ii,:) = hyp_HW2(NXX(ii), NTT(ii));
end

figure
loglog(dt, err(:,1), "-o", LineWidth=2)
grid minor
xlabel('time res dt','FontSize',16);
ylabel('|| u_t-u_{t,h} ||_{\Delta,2}','FontSize',16, Interpreter='tex');

figure
loglog(dx, err(:, 2), "-o", LineWidth=2)
grid minor
xlabel('space res dx','FontSize',16);
ylabel('|| u_x-u_{x,h} ||_{\Delta,2}','FontSize',16, Interpreter='tex');