clear; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Initial data 
u0 = @(x)   sin(pi*x);

% boudary conditions 
u0_x = @(x) 0*x;
uj_x = @(x) 0*x;

% Intervals
 T = 2.5;

I = [0 1];


% Number of time steps
NT = 2000;
% Number of space steps
NX = 1500;

%%
mu = 10


[SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_x, uj_x);
%SOL_EX = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu);

figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
space = linspace(I(1),I(2),NX+1);
index1 = round(0.1 / dt);
index2 = round(0.5 / dt);
index3 = round(2.3 / dt);

amp1 = SOL(:,index1);
amp2 = SOL(:,index2);
amp3 = SOL(:,index3);

[ampEx1] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,index1);
[ampEx2] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,index2);
[ampEx3] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,index3);


plot(space,amp1,'LineWidth',2);
hold on 
plot(space,ampEx1,"--",'LineWidth',2);
hold on
plot(space,amp2,'LineWidth',2);
hold on 
plot(space,ampEx2,"--",'LineWidth',2);
hold on
plot(space,amp3,'LineWidth',2);
hold on
plot(space,ampEx3,"--",'LineWidth',2);
% title("Numerical to exact solution comparison",'FontSize',16,"Interpreter","latex")
% xlabel('space x ','FontSize',16,"Interpreter","latex");
% ylabel('$u(x)$','FontSize',16,"Interpreter","latex");
% legend(["$T=0.1$","$T=0.1$ exact","$T=0.5$","$T=0.5$ exact","$T=2.3$","$T=2.3$ exact"],'FontSize',16,"Interpreter","latex")
title("Numerical to exact solution comparison")
xlabel('space x ');
ylabel('u(x)');
legend(["T=0.1","T=0.1 exact","T=0.5","T=0.5 exact","T=2.3","T=2.3 exact"])

saveas(gcf, ".\plots\amp_mu10.png");
%%

mu = 10

errorsMaxIndex =4;
l2_errors = zeros(errorsMaxIndex,1);
lInf_errors = zeros(errorsMaxIndex,1);
% varying dx
for i = 1:errorsMaxIndex
    h_base = 0.0005;
    NT=2000;
    NX = (I(2)-I(1))/h_base*i;
    [SOL,SOL_EX,L1_ERR,L2_ERR,LINF_ERR,dx,dt] = FD_1D_BURGER_FUN(mu,T,I,NT, NX ,u0, u0_x, uj_x,2,NT+1);
    l2_errors(i,1) = L2_ERR;
    lInf_errors(i,1) = LINF_ERR;
end 

figure(4)
plot(1:4,l2_errors);
saveas(gcf, ".\plots\L2_error.png");

figure(5)
plot(1:4,lInf_errors);
saveas(gcf, ".\plots\LINF_error.png");

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error and order of convergence 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 10

%L2_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),2)*dx^0.5;
%L1_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),1)*dx;
%LINF_ERR = 0;%norm(SOL_EX-SOL(:,NT+1),Inf);

errorsMaxIndex =4;


l2_errors_h = zeros(errorsMaxIndex,1);
lInf_errors_h = zeros(errorsMaxIndex,1);
l2_errors_t = zeros(errorsMaxIndex,1);
lInf_errors_t = zeros(errorsMaxIndex,1);
% varying dx
h_base = 0.0005;
dt_base = 0.0005;
 NT= T/dt_base;
 indexError = NT+1;
for i = 0:errorsMaxIndex-1
    
    NX = (I(2)-I(1))/(h_base*2^(i));
    
    [SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_x, uj_x);
    disp("dx="+num2str(dx));
    [SOL_EX] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,indexError);
    L2_ERR  = norm(SOL_EX-SOL(:,indexError),2)*(h_base*(2^i))^0.5;
    LINF_ERR = norm(SOL_EX-SOL(:,indexError),Inf);
    l2_errors_h(i+1,1) = L2_ERR;
    lInf_errors_h(i+1,1) = LINF_ERR;
end 

NX = (I(2)-I(1))/(h_base);

for i = 0:errorsMaxIndex-1
    
    NT= T/(dt_base*(2^i));
    indexError = NT+1;
    [SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_x, uj_x);
    disp("dt="+num2str(dt));
    [SOL_EX] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,indexError);
    L2_ERR  = norm(SOL_EX-SOL(:,indexError),2)*(dt_base*(2^i))^0.5;
    LINF_ERR = norm(SOL_EX-SOL(:,indexError),Inf);
    l2_errors_t(i+1,1) = L2_ERR;
    lInf_errors_t(i+1,1) = LINF_ERR;
end 

disp("finito");

% order of convergence 
% L2 

E_12_h = l2_errors_h(1,1);
E_22_h = l2_errors_h(end,1);

ph2 = (log(E_12_h/E_22_h))/(log(h_base/(h_base*(2^(errorsMaxIndex-1)))));
conv2_val_h = zeros(errorsMaxIndex,1);
for i = 0:errorsMaxIndex-1
    conv2_val_h(i+1,1) = (h_base*(2^i))^(ph2);
end

E_1Inf_h = lInf_errors_h(1,1);
E_2Inf_h = lInf_errors_h(end,1);
phINF = (log(E_1Inf_h/E_2Inf_h))/(log(h_base/(h_base*(2^(errorsMaxIndex-1)))));
convInf_val_h = zeros(errorsMaxIndex,1);
for i = 0:errorsMaxIndex-1
    convInf_val_h(i+1,1) = (h_base*(2^i))^(phINF);
end

E_12_t = l2_errors_t(1,1);
E_22_t = l2_errors_t(end,1);

pt2 = (log(E_12_t/E_22_t))/(log(dt_base/(dt_base*(2^(errorsMaxIndex-1)))));
conv2_val_t = zeros(errorsMaxIndex,1);
for i = 0:errorsMaxIndex-1
    conv2_val_t(i+1,1) = (dt_base*(2^i))^(pt2);
end

E_1Inf_t = lInf_errors_t(1,1);
E_2Inf_t = lInf_errors_t(end,1);
ptINF = log(E_1Inf_t/E_2Inf_t)/log(dt_base/(dt_base*(2^(errorsMaxIndex-1))));
convInf_val_t = zeros(errorsMaxIndex,1);
for i = 0:errorsMaxIndex-1
    convInf_val_t(i+1,1) = (dt_base*(2^i))^(ptINF);
end

hhtt= 0.0005*(2.^(0:3));

%figure('Renderer', 'painters', 'Position', [100 100 1000 600]);

%subplot(2,2,1)
figure()
loglog(hhtt,l2_errors_h,'-+b','LineWidth',2);
hold on 
loglog(hhtt,conv2_val_h,'-or','LineWidth',2);
% title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('h','FontSize',16,"Interpreter","latex");
% ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_2$","$h_2^{q}$"],'FontSize',16,"Interpreter","latex")
title("L^2 norm and convergence error")
xlabel('h');
ylabel('L^2 error');
legend(["||u-u_{ex}||_2","h_2^{q}"])

saveas(gcf, ".\plots\L2h_error_conv_mu10.png");

disp("q2 = "+num2str(ph2));

%subplot(2,2,2)
figure()
loglog(hhtt,lInf_errors_h,'-+b','LineWidth',2);
hold on 
loglog(hhtt,convInf_val_h,'-or','LineWidth',2);
% title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('h','FontSize',16,"Interpreter","latex");
% ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_{\infty}$","$h_{\infty}^{q}$"],'FontSize',16,"Interpreter","latex")

title("L^{\infty} norm and convergence error")
xlabel('h');
ylabel('L^{\infty} error');
legend(["||u-u_{ex}||_{\infty}","h_{\infty}^{q}"])
saveas(gcf, ".\plots\LINFh_error_conv_mu10.png");
disp("qINF = "+num2str(phINF));

%subplot(2,2,3)
figure()
loglog(hhtt,l2_errors_t,'-+b','LineWidth',2);
hold on 
loglog(hhtt,conv2_val_t,'-or','LineWidth',2);
% title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
% ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_2$","$\Delta t_2^{p}$"],'FontSize',16,"Interpreter","latex")
title("L^2 norm and convergence error")
xlabel('\Delta t');
ylabel('L^2 error');
legend(["||u-u_{ex}||_2","\Delta t_2^{p}"])
saveas(gcf, ".\plots\L2t_error_conv_mu10.png");
disp("p2 = "+num2str(pt2));

%subplot(2,2,4)
figure()
loglog(hhtt,lInf_errors_t,'-+b','LineWidth',2);
hold on 
loglog(hhtt,convInf_val_t,'-or','LineWidth',2);
% title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
% ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_{\infty}$","$\Delta t_{\infty}^{p}$"],'FontSize',16,"Interpreter","latex")
title("L^{\infty} norm and convergence error")
xlabel('\Delta t');
ylabel('L^{\infty} error');
legend(["||u-u_{ex}||_{\infty}","\Delta t_{\infty}^{p}"])
%saveas(gcf, ".\plots\_error_conv_mu10.png");
saveas(gcf, ".\plots\LINFt_error_conv_mu10.png");
disp("pINF = "+num2str(2));




%%

mu = 1

[SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_x, uj_x);

figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
space = linspace(I(1),I(2),NX+1);
index1 = round(0.1 / dt);
index2 = round(0.5 / dt);
index3 = round(2.3 / dt);

amp1 = SOL(:,index1);
amp2 = SOL(:,index2);
amp3 = SOL(:,index3);

[ampEx1] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,index1);
[ampEx2] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,index2);
[ampEx3] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,index3);

plot(space,amp1,'LineWidth',2);
hold on 
plot(space,ampEx1,"--",'LineWidth',2);
hold on
plot(space,amp2,'LineWidth',2);
hold on 
plot(space,ampEx2,"--",'LineWidth',2);
hold on
plot(space,amp3,'LineWidth',2);
hold on
plot(space,ampEx3,"--",'LineWidth',2);
% title("Numerical to exact solution comparison",'FontSize',16,"Interpreter","latex")
% xlabel('space x ','FontSize',16,"Interpreter","latex");
% ylabel('$u(x)$','FontSize',16,"Interpreter","latex");
% legend(["$T=0.1$","$T=0.1$ exact","$T=0.5$","$T=0.5$ exact","$T=2.3$","$T=2.3$ exact"],'FontSize',16,"Interpreter","latex")
title("Numerical to exact solution comparison")
xlabel('space x');
ylabel('u(x)');
legend(["T=0.1","T=0.1 exact","T=0.5","T=0.5 exact","T=2.3","T=2.3 exact"])

saveas(gcf, ".\plots\amp_mu1.png");

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error and order of convergence 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 1 

%L2_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),2)*dx^0.5;
%L1_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),1)*dx;
%LINF_ERR = 0;%norm(SOL_EX-SOL(:,NT+1),Inf);

errorsMaxIndex =4;
l2_errors_h = zeros(errorsMaxIndex,1);
lInf_errors_h = zeros(errorsMaxIndex,1);
l2_errors_t = zeros(errorsMaxIndex,1);
lInf_errors_t = zeros(errorsMaxIndex,1);
% varying dx
h_base = 0.0005;
dt_base = 0.0005;
 NT= T/dt_base;
 indexError = NT+1;
for i = 0:errorsMaxIndex-1
    
    NX = (I(2)-I(1))/(h_base*2^(i));
    
    [SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_x, uj_x);
    disp("dx="+num2str(dx));
    [SOL_EX] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,indexError);
    L2_ERR  = norm(SOL_EX-SOL(:,indexError),2)*(h_base*(2^i))^0.5;
    LINF_ERR = norm(SOL_EX-SOL(:,indexError),Inf);
    l2_errors_h(i+1,1) = L2_ERR;
    lInf_errors_h(i+1,1) = LINF_ERR;
end 

NX = (I(2)-I(1))/(h_base);

for i = 0:errorsMaxIndex-1
    
    NT= T/(dt_base*(2^i));
    indexError = NT+1;
    [SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_x, uj_x);
    disp("dt="+num2str(dt));
    [SOL_EX] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,indexError);
    L2_ERR  = norm(SOL_EX-SOL(:,indexError),2)*(dt_base*(2^i))^0.5;
    LINF_ERR = norm(SOL_EX-SOL(:,indexError),Inf);
    l2_errors_t(i+1,1) = L2_ERR;
    lInf_errors_t(i+1,1) = LINF_ERR;
end 

disp("finito");

% order of convergence 
% L2 

E_12_h = l2_errors_h(1,1);
E_22_h = l2_errors_h(end,1);

ph2 = (log(E_12_h/E_22_h))/(log(h_base/(h_base*(2^(errorsMaxIndex-1)))));
conv2_val_h = zeros(errorsMaxIndex,1);
for i = 0:errorsMaxIndex-1
    conv2_val_h(i+1,1) = (h_base*(2^i))^(ph2);
end

E_1Inf_h = lInf_errors_h(1,1);
E_2Inf_h = lInf_errors_h(end,1);
phINF = log(E_1Inf_h/E_2Inf_h)/log(h_base/(h_base*(2^(errorsMaxIndex-1))));
convInf_val_h = zeros(errorsMaxIndex,1);
for i = 0:errorsMaxIndex-1
    convInf_val_h(i+1,1) = (h_base*(2^i))^(phINF);
end

E_12_t = l2_errors_t(1,1);
E_22_t = l2_errors_t(end,1);

pt2 = (log(E_12_t/E_22_t))/(log(dt_base/(dt_base*(2^(errorsMaxIndex-1)))));
conv2_val_t = zeros(errorsMaxIndex,1);
for i = 0:errorsMaxIndex-1
    conv2_val_t(i+1,1) = (dt_base*(2^i))^(pt2);
end

E_1Inf_t = lInf_errors_t(1,1);
E_2Inf_t = lInf_errors_t(end,1);
ptINF = log(E_1Inf_t/E_2Inf_t)/log(dt_base/(dt_base*(2^(errorsMaxIndex-1))));
convInf_val_t = zeros(errorsMaxIndex,1);
for i = 0:errorsMaxIndex-1
    convInf_val_t(i+1,1) = (dt_base*(2^i))^(ptINF);
end



hhtt= 0.0005*(2.^(0:3));

%figure('Renderer', 'painters', 'Position', [100 100 1000 600]);

%subplot(2,2,1)
figure()
loglog(hhtt,l2_errors_h,'-+b','LineWidth',2);
hold on 
loglog(hhtt,conv2_val_h,'-or','LineWidth',2);
% title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('h','FontSize',16,"Interpreter","latex");
% ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_2$","$h_2^{q}$"],'FontSize',16,"Interpreter","latex")
title("L^2 norm and convergence error")
xlabel('h');
ylabel('L^2 error');
legend(["||u-u_{ex}||_2","h_2^{q}"])
saveas(gcf, ".\plots\L2h_error_conv_mu1.png");
disp("q2 = "+num2str(ph2));

%subplot(2,2,2)
figure()
loglog(hhtt,lInf_errors_h,'-+b','LineWidth',2);
hold on 
loglog(hhtt,convInf_val_h,'-or','LineWidth',2);
% title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('h','FontSize',16,"Interpreter","latex");
% ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_{\infty}$","$h_{\infty}^{q}$"],'FontSize',16,"Interpreter","latex")
title("L^{\infty} norm and convergence error")
xlabel('h');
ylabel('L^{\infty} error');
legend(["||u-u_{ex}||_{\infty}","h_{\infty}^{q}"])

saveas(gcf, ".\plots\LINFh_error_conv_mu1.png");
disp("qINF = "+num2str(phINF));

%subplot(2,2,3)
figure()
loglog(hhtt,l2_errors_t,'-+b','LineWidth',2);
hold on 
loglog(hhtt,conv2_val_t,'-or','LineWidth',2);
% title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
% ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_2$","$\Delta t_2^{p}$"],'FontSize',16,"Interpreter","latex")
title("L^2 norm and convergence error")
xlabel('\Delta t');
ylabel('L^2 error');
legend(["||u-u_{ex}||_2","\Delta t_2^{p}"])
saveas(gcf, ".\plots\L2t_error_conv_mu1.png");
disp("p2 = "+num2str(pt2));

%subplot(2,2,4)
figure()
loglog(hhtt,lInf_errors_t,'-+b','LineWidth',2);
hold on 
loglog(hhtt,convInf_val_t,'-or','LineWidth',2);
% title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
% ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_{\infty}$","$\Delta t_{\infty}^{p}$"],'FontSize',16,"Interpreter","latex")
title("$L^{\infty}$ norm and convergence error")
xlabel('$\Delta t$');
ylabel('$L^{\infty}$ error');
legend(["$||u-u_{ex}||_{\infty}$","$\Delta t_{\infty}^{p}$"])
% saveas(gcf, ".\plots\_error_conv_mu1.png");
saveas(gcf, ".\plots\LINFt_error_conv_mu1.png");
disp("pINF = "+num2str(ptINF));



%%

mu = 0.1

[SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_x, uj_x);

figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
space = linspace(I(1),I(2),NX+1);
index1 = round(0.1 / dt);
index2 = round(0.5 / dt);
index3 = round(2.3 / dt);

amp1 = SOL(:,index1);
amp2 = SOL(:,index2);
amp3 = SOL(:,index3);

[ampEx1] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,index1);
[ampEx2] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,index2);
[ampEx3] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,index3);


plot(space,amp1,'LineWidth',2);
hold on 
plot(space,ampEx1,"--",'LineWidth',2);
hold on
plot(space,amp2,'LineWidth',2);
hold on 
plot(space,ampEx2,"--",'LineWidth',2);
hold on
plot(space,amp3,'LineWidth',2);
hold on
plot(space,ampEx3,"--",'LineWidth',2);
% title("Numerical to exact solution comparison",'FontSize',16,"Interpreter","latex")
% xlabel('space x ','FontSize',16,"Interpreter","latex");
% ylabel('$u(x)$','FontSize',16,"Interpreter","latex");
% legend(["$T=0.1$","$T=0.1$ exact","$T=0.5$","$T=0.5$ exact","$T=2.3$","$T=2.3$ exact"],'FontSize',16,"Interpreter","latex")
title("Numerical to exact solution comparison")
xlabel('space x');
ylabel('u(x)');
legend(["T=0.1","T=0.1 exact","T=0.5","T=0.5 exact","T=2.3","T=2.3 exact"])

saveas(gcf, ".\plots\amp_mu01.png");

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error and order of convergence 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 0.1

%L2_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),2)*dx^0.5;
%L1_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),1)*dx;
%LINF_ERR = 0;%norm(SOL_EX-SOL(:,NT+1),Inf);

errorsMaxIndex =4;
l2_errors_h = zeros(errorsMaxIndex,1);
lInf_errors_h = zeros(errorsMaxIndex,1);
l2_errors_t = zeros(errorsMaxIndex,1);
lInf_errors_t = zeros(errorsMaxIndex,1);
% varying dx
h_base = 0.0005;
dt_base = 0.0005;
 NT= T/dt_base;
 indexError = NT+1;
for i = 0:errorsMaxIndex-1
    
    NX = (I(2)-I(1))/(h_base*2^(i));
    
    [SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_x, uj_x);
    disp("dx="+num2str(dx));
    [SOL_EX] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,indexError);
    L2_ERR  = norm(SOL_EX-SOL(:,indexError),2)*(h_base*(2^i))^0.5;
    LINF_ERR = norm(SOL_EX-SOL(:,indexError),Inf);
    l2_errors_h(i+1,1) = L2_ERR;
    lInf_errors_h(i+1,1) = LINF_ERR;
end 

NX = (I(2)-I(1))/(h_base);

for i = 0:errorsMaxIndex-1
    
    NT= T/(dt_base*(2^i));
    indexError = NT+1;
    [SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_x, uj_x);
    disp("dt="+num2str(dt));
    [SOL_EX] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,indexError);
    L2_ERR  = norm(SOL_EX-SOL(:,indexError),2)*(dt_base*(2^i))^0.5;
    LINF_ERR = norm(SOL_EX-SOL(:,indexError),Inf);
    l2_errors_t(i+1,1) = L2_ERR;
    lInf_errors_t(i+1,1) = LINF_ERR;
end 

disp("finito");

% order of convergence 
% L2 

E_12_h = l2_errors_h(1,1);
E_22_h = l2_errors_h(end,1);

ph2 = (log(E_12_h/E_22_h))/(log(h_base/(h_base*(2^(errorsMaxIndex-1)))));
conv2_val_h = zeros(errorsMaxIndex,1);
for i = 0:errorsMaxIndex-1
    conv2_val_h(i+1,1) = (h_base*(2^i))^(ph2);
end

E_1Inf_h = lInf_errors_h(1,1);
E_2Inf_h = lInf_errors_h(end,1);
phINF = log(E_1Inf_h/E_2Inf_h)/log(h_base/(h_base*(2^(errorsMaxIndex-1))));
convInf_val_h = zeros(errorsMaxIndex,1);
for i = 0:errorsMaxIndex-1
    convInf_val_h(i+1,1) = (h_base*(2^i))^(phINF);
end

E_12_t = l2_errors_t(1,1);
E_22_t = l2_errors_t(end,1);

pt2 = (log(E_12_t/E_22_t))/(log(dt_base/(dt_base*(2^(errorsMaxIndex-1)))));
conv2_val_t = zeros(errorsMaxIndex,1);
for i = 0:errorsMaxIndex-1
    conv2_val_t(i+1,1) = (dt_base*(2^i))^(pt2);
end

E_1Inf_t = lInf_errors_t(1,1);
E_2Inf_t = lInf_errors_t(end,1);
ptINF = log(E_1Inf_t/E_2Inf_t)/log(dt_base/(dt_base*(2^(errorsMaxIndex-1))));
convInf_val_t = zeros(errorsMaxIndex,1);
for i = 0:errorsMaxIndex-1
    convInf_val_t(i+1,1) = (dt_base*(2^i))^(ptINF);
end

hhtt= 0.0005*(2.^(0:3));

%figure('Renderer', 'painters', 'Position', [100 100 1000 600]);

%subplot(2,2,1)
figure()
loglog(hhtt,l2_errors_h,'-+b','LineWidth',2);
hold on 
loglog(hhtt,conv2_val_h,'-or','LineWidth',2);
% title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('h','FontSize',16,"Interpreter","latex");
% ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_2$","$h_2^{q}$"],'FontSize',16,"Interpreter","latex")
title("L^2 norm and convergence error")
xlabel('h');
ylabel('L^2 error');
legend(["||u-u_{ex}||_2","h_2^{q}"])
saveas(gcf, ".\plots\L2h_error_conv_mu01.png");
disp("q2 = "+num2str(ph2));

%subplot(2,2,2)
figure()
loglog(hhtt,lInf_errors_h,'-+b','LineWidth',2);
hold on 
loglog(hhtt,convInf_val_h,'-or','LineWidth',2);
% title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('h','FontSize',16,"Interpreter","latex");
% ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_{\infty}$","$h_{\infty}^{q}$"],'FontSize',16,"Interpreter","latex")

title("L^{\infty} norm and convergence error")
xlabel('h');
ylabel('L^{\infty} error');
legend(["||u-u_{ex}||_{\infty}","h_{\infty}^{q}"])
saveas(gcf, ".\plots\LINFh_error_conv_mu01.png");
disp("qINF = "+num2str(phINF));

%subplot(2,2,3)
figure()
loglog(hhtt,l2_errors_t,'-+b','LineWidth',2);
hold on 
loglog(hhtt,conv2_val_t,'-or','LineWidth',2);
% title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
% ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_2$","$\Delta t_2^{p}$"],'FontSize',16,"Interpreter","latex")
title("L^2 norm and convergence error")
xlabel('\Delta t');
ylabel('L^2 error');
legend(["||u-u_{ex}||_2","\Delta t_2^{p}"])
saveas(gcf, ".\plots\L2t_error_conv_mu01.png");
disp("p2 = "+num2str(pt2));

%subplot(2,2,4)
figure()
loglog(hhtt,lInf_errors_t,'-+b','LineWidth',2);
hold on 
loglog(hhtt,convInf_val_t,'-or','LineWidth',2);
% title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
% ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_{\infty}$","$\Delta t_{\infty}^{p}$"],'FontSize',16,"Interpreter","latex")
% saveas(gcf, ".\plots\_error_conv_mu01.png");
title("L^{\infty} norm and convergence error")
xlabel('\Delta t');
ylabel('L^{\infty} error');
legend(["||u-u_{ex}||_{\infty}","\Delta t_{\infty}^{p}"])
saveas(gcf, ".\plots\LINFt_error_conv_mu01.png");
disp("pINF = "+num2str(ptINF));



% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Error and order of convergence 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mu = 0.1 
% 
% %L2_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),2)*dx^0.5;
% %L1_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),1)*dx;
% %LINF_ERR = 0;%norm(SOL_EX-SOL(:,NT+1),Inf);
% 
% errorsMaxIndex =4;
% l2_errors = zeros(errorsMaxIndex,1);
% lInf_errors = zeros(errorsMaxIndex,1);
% % varying dx
% h_base = 0.0005;
% dt_base = 0.0005;
%  NT= T/dt_base;
%  indexError = NT+1;
% for i = 0:errorsMaxIndex-1
%     
%     NX = (I(2)-I(1))/(h_base*2^(i));
%     
%     [SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_x, uj_x);
%     disp("dx="+num2str(dx));
%     [SOL_EX] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,indexError);
%     L2_ERR  = norm(SOL_EX-SOL(:,indexError),2)*(h_base*2^i)^0.5;
%     LINF_ERR = norm(SOL_EX-SOL(:,indexError),Inf);
%     l2_errors(i+1,1) = L2_ERR;
%     lInf_errors(i+1,1) = LINF_ERR;
% end 
% 
% disp("finito");
% %%
% % order of convergence 
% % L2 
% E_12 = l2_errors(1,1);
% E_22 = l2_errors(end,1);
% 
% ph2 = (log(E_12/E_22))/(log(h_base/(h_base*(2^errorsMaxIndex))));
% conv2_val = zeros(errorsMaxIndex,1);
% for i = 0:errorsMaxIndex-1
%     conv2_val(i+1,1) = (h_base*(2^i))^(ph2);
% end
% 
% E_1Inf = lInf_errors(1,1);
% E_2Inf = lInf_errors(end,1);
% phINF = log(E_1Inf/E_2Inf)/log(h_base/(h_base*(2^errorsMaxIndex)));
% convInf_val = zeros(errorsMaxIndex,1);
% for i = 0:errorsMaxIndex-1
%     convInf_val(i+1,1) = (h_base*(2^i))^(phINF);
% end
% 
% 
% %% 
% figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
% plot(1:4,l2_errors,"o-");
% hold on 
% plot(1:4,conv2_val,"o-");
% title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('h','FontSize',16,"Interpreter","latex");
% ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_2$","$h_2^{q}$"],'FontSize',16,"Interpreter","latex")
% saveas(gcf, ".\plots\L2h_error_conv.png");
% disp("q = "+num2str(2));
% 
% figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
% plot(1:4,lInf_errors,"o-");
% hold on 
% plot(1:4,convInf_val,"o-");
% title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('h','FontSize',16,"Interpreter","latex");
% ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_{\infty}$","$h_{\infty}^{q}$"],'FontSize',16,"Interpreter","latex")
% saveas(gcf, ".\plots\LINFh_error_conv.png");
% disp("q = "+num2str(2));
% 
% %%
% 
% errorsMaxIndex =4;
% l2_errors = zeros(errorsMaxIndex,1);
% lInf_errors = zeros(errorsMaxIndex,1);
% % varying dx
% h_base = 0.0005;
% dt_base = 0.0005;
% NX = (I(2)-I(1))/h_base;
% for i = 1:errorsMaxIndex
%     NT= T/dt_base*2^i;
%     [SOL,SOL_EX,L1_ERR,L2_ERR,LINF_ERR,dx,dt] = FD_1D_BURGER_FUN(mu,T,I,NT, NX ,u0, u0_x, uj_x,2,NT+1)
%     L2_ERR   = norm(SOL_EX-SOL(:,NT+1),2)*dx^0.5;
%     LINF_ERR = norm(SOL_EX-SOL(:,NT+1),Inf);
%     l2_errors(i,1) = L2_ERR;
%     lInf_errors(i,1) = LINF_ERR;
% end 
% 
% %%
% % order of convergence 
% % L2 
% E_12 = l2_errors(1,1);
% E_22 = l2_errors(end,1);
% 
% ph = log(E_12/E_22)/log(h_base/(h_base*2^errorsMaxIndex));
% conv2_val = zeros(errorsMaxIndex,1);
% for i = 1:errorsMaxIndex
%     conv2_val(i,1) = (h_base*i)^(ph);
% end
% 
% E_1Inf = lInf_errors(1,1);
% E_2Inf = lInf_errors(end,1);
% ph = log(E_1Inf/E_2Inf)/log(h_base/(h_base*errorsMaxIndex));
% convInf_val = zeros(errorsMaxIndex,1);
% for i = 1:errorsMaxIndex
%     convInf_val(i,1) = (h_base*i)^(ph);
% end
% 
% 
% 
% L1_ERR   = norm(SOL_EX-SOL(:,indexError),1)*dx;
% 
% %% 
% figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
% plot(1:4,l2_errors,"o-");
% hold on 
% plot(1:4,conv2_val,"o-");
% title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
% ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_2$","$\Delta t_2^{p}$"],'FontSize',16,"Interpreter","latex")
% saveas(gcf, ".\plots\L2dt_error_conv.png");
% 
% figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
% plot(1:4,lInf_errors,"o-");
% hold on 
% plot(1:4,convInf_val,"o-");
% title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
% xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
% ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
% legend(["$||u-u_{ex}||_{\infty}$","$\Delta t_{\infty}^{p}$"],'FontSize',16,"Interpreter","latex")
% saveas(gcf, ".\plots\LINFdt_error_conv.png");