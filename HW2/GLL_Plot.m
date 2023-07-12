clear; close all; clc
path = "./Report_Assets/"
xx=linspace(-1,1,1000);
saving = 0;

% L0 = @(x) 1+0.*x;
% L1 = @(x) 1.*x;
% L2 = @(x) 3/2.*x^2 - 0.5;
% L2p = @(x) 3.*x;
% 
% L3 = @(x) 1/2.*(5.*x^3 - 3.*x);
% L3p = @(x) 1/2.*(15.*x.^2 - 3);

syms c
N = 2;
z = vpasolve(diff(legendreP(N,c)) == 0);
z = double(z)';
phi_P2 = zeros(N+1, length(xx));
nodes_2 = cat(2, -1, z, 1);
L2 = @(x) subs(legendreP(N,c),c,x);
L2p = @(x) subs(diff(legendreP(N,c)),c,x);

% n=2;
% phi_P2 = zeros(n+1, length(xx));
% nodes_2 = [-1,0,1];

for i=1:length(nodes_2)
    phi_P2(i,:) = -1/(N*(N+1)) .* (1-xx.^2)./(xx-nodes_2(i)) .* L2p(xx) ./ L2(nodes_2(i)); 
end
figure();
plot(xx, phi_P2(1,:), xx, phi_P2(2,:), xx, phi_P2(3,:), LineWidth=3);
grid minor
legend(["$\varphi_0$", "$\varphi_1$", "$\varphi_2$"], Interpreter="latex", FontSize=14, Location="northeastoutside")
xlabel("$x$", Interpreter="latex", FontSize=14);
ylabel("$\varphi$", Interpreter="latex", FontSize=14);
filename="P2_basis.png";
if (saving) saveas(gcf, [path+filename]); end
% close all


% TERZO ORDINE
N = 3;
z = vpasolve(diff(legendreP(N,c)) == 0);
z = double(z)';
phi_P3 = zeros(N+1, length(xx));
nodes_3 = cat(2, -1, z, 1);
L3 = @(x) subs(legendreP(N,c),c,x);
L3p = @(x) subs(diff(legendreP(N,c)),c,x);


for i=1:length(nodes_3)
    phi_P3(i,:) = -1/(N*(N+1)) .* (1-xx.^2)./(xx-nodes_3(i)) .* L3p(xx) ./ L3(nodes_3(i)); 
end
figure();
plot(xx, phi_P3(1,:), xx, phi_P3(2,:), xx, phi_P3(3,:), xx, phi_P3(4,:), LineWidth=3);
grid minor
legend(["$\varphi_0$", "$\varphi_1$", "$\varphi_2$", "$\varphi_3$"], Interpreter="latex", FontSize=14, Location="northeastoutside")
xlabel("$x$", Interpreter="latex", FontSize=14);
ylabel("$\varphi$", Interpreter="latex", FontSize=14);
filename="P3_basis.png";
if (saving) saveas(gcf, [path+filename]); end

% close all

% QUARTO ORDINE
N = 4;
z = vpasolve(diff(legendreP(N,c)) == 0);
z = double(z)';
L4 = @(x) subs(legendreP(N,c),c,x);
L4p = @(x) subs(diff(legendreP(N,c)),c,x);
phi_P4 = zeros(N+1, length(xx));
nodes_4 = cat(2, -1, z, 1);

figure();
for i=1:length(nodes_4)
    phi_P4(i,:) = -1/(N*(N+1)) .* (1-xx.^2)./(xx-nodes_4(i)) .* L4p(xx) ./ L4(nodes_4(i)); 
    plot(xx, phi_P4(i,:),LineWidth=3);
    hold on;
end
grid minor
legend(["$\varphi_0$", "$\varphi_1$", "$\varphi_2$", "$\varphi_3$", "$\varphi_4$"], Interpreter="latex", FontSize=14, Location="northeastoutside")
xlabel("$x$", Interpreter="latex", FontSize=14);
ylabel("$\varphi$", Interpreter="latex", FontSize=14);
filename="P4_basis.png";
if (saving) saveas(gcf, [path+filename]); end

% close all

node_basis = cat(2, phi_P2(3,:), phi_P2(1,:));
figure();
plot(linspace(0,1,1000), phi_P2(1,:), LineWidth=3)
hold on
plot(linspace(0,1,1000), phi_P2(2,:), LineWidth=3)
hold on
plot(linspace(0,2,2000), node_basis, LineWidth=3);
hold on
plot(linspace(1,2,1000), phi_P2(2,:), LineWidth=3)
hold on
plot(linspace(1,2,1000), phi_P2(3,:), LineWidth=3)
grid minor
legend(["$\Psi_0$", "$\Psi_1^{(1)}$", "$\Psi_1$", "$\Psi_1^{(2)}$", "$\Psi_2$"],Interpreter="latex", Location="northeastoutside", FontSize=14)
filename="union_basis.png";
if (saving) saveas(gcf, [path+filename]); end



%%
figure()
plot(xx, gradient(phi_P3(1,:)))
hold on
plot(xx, gradient(phi_P3(2,:)))
hold on
plot(xx, gradient(phi_P3(3,:)))
hold on
plot(xx, gradient(phi_P3(4,:)))
grid minor


% node1 277
% node2 724
grad = zeros(size(phi_P3));
grad(1,:) = gradient(phi_P3(1,:), 2/1000);
grad(2,:) = gradient(phi_P3(2,:), 2/1000);
grad(3,:) = gradient(phi_P3(3,:), 2/1000);
grad(4,:) = gradient(phi_P3(4,:), 2/1000);

figure()
plot(xx, grad(1,:))
hold on
plot(xx, grad(2,:))
hold on
plot(xx, grad(3,:))
hold on
plot(xx, grad(4,:))
grid minor
