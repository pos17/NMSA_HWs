clear; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Initial data 
u0 = @(x)   1/4*cos(pi*x);

% boudary conditions 
u0_t = @(t) 0*t;
uj_t = @(t) -(1/4)*exp(-mu*t);

% Intervals
 T = 2.5;

I = [0 1];


% Number of time steps
NT = 2000;
% Number of space steps
NX = 1500;

%%
mu = 10

uj_t = @(t) -(1/4)*exp(-mu*t);

[SOL,dx,dt] = FD_1D_BURGER_FUN2_EX4(mu,T,I,NT, NX ,u0, u0_t, uj_t);

%REFERENCE SOLUTION

[SOL_EX,dxEx,dtEx] = FD_1D_BURGER_FUN2_EX4(mu,T,I,5000, 2000 ,u0, u0_t, uj_t);

%SOL_EX = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu);

figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
space = linspace(I(1),I(2),NX+1);
index1 = round(0.1 / dt);
index2 = round(0.5 / dt);
index3 = round(2.3 / dt);

amp1 = SOL(:,index1);
amp2 = SOL(:,index2);
amp3 = SOL(:,index3);

spaceEx = linspace(I(1),I(2),2000+1);
indexEx1 = round(0.1 / dtEx);
indexEx2 = round(0.5 / dtEx);
indexEx3 = round(2.3 / dtEx);

ampEx1 = SOL_EX(:,indexEx1);
ampEx2 = SOL_EX(:,indexEx2); 
ampEx3 = SOL_EX(:,indexEx3);


plot(space,amp1,'LineWidth',2);
hold on 
plot(spaceEx,ampEx1,"--",'LineWidth',2);
hold on
plot(space,amp2,'LineWidth',2);
hold on 
plot(spaceEx,ampEx2,"--",'LineWidth',2);
hold on
plot(space,amp3,'LineWidth',2);
hold on
plot(spaceEx,ampEx3,"--",'LineWidth',2);
title("Numerical to exact solution comparison",'FontSize',16,"Interpreter","latex")
xlabel('space x ','FontSize',16,"Interpreter","latex");
ylabel('$u(x)$','FontSize',16,"Interpreter","latex");
legend(["$T=0.1$","$T=0.1$ exact","$T=0.5$","$T=0.5$ exact","$T=2.3$","$T=2.3$ exact"],'FontSize',16,"Interpreter","latex")
saveas(gcf, ".\plots\amp_mu10_EX4.png");
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error and order of convergence 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 10

uj_t = @(t) -(1/4)*exp(-mu*t);

%L2_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),2)*dx^0.5;
%L1_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),1)*dx;
%LINF_ERR = 0;%norm(SOL_EX-SOL(:,NT+1),Inf);

errorsMaxIndex =3;


l2_errors_h = zeros(errorsMaxIndex,1);
lInf_errors_h = zeros(errorsMaxIndex,1);
l2_errors_t = zeros(errorsMaxIndex,1);
lInf_errors_t = zeros(errorsMaxIndex,1);
% varying dx
h_base = 0.0005;
dt_base = 0.0005;
 NT= T/dt_base;
 indexError = NT+1;
 [SOL_EX,dx,dt] = FD_1D_BURGER_FUN2_EX4(mu,T,I,5000, 2000 ,u0, u0_t, uj_t);
for i = 1:errorsMaxIndex
    
    NX = (I(2)-I(1))/(h_base*2^(i));
    
    [SOL,dx,dt] = FD_1D_BURGER_FUN2_EX4(mu,T,I,NT, NX ,u0, u0_t, uj_t);
    disp("dx="+num2str(dx));
    SOL_EXX = SOL_EX(1:2^(i):end,2001);
    L2_ERR  = norm(SOL_EXX-SOL(:,indexError),2)*(h_base*(2^i))^0.5;
    LINF_ERR = norm(SOL_EXX-SOL(:,indexError),Inf);
    l2_errors_h(i,1) = L2_ERR;
    lInf_errors_h(i,1) = LINF_ERR;
end 

NX = (I(2)-I(1))/(h_base);

for i = 1:errorsMaxIndex
    
    NT= T/(dt_base*(2^i));
    indexError = NT+1;
    [SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_t, uj_t);
    disp("dt="+num2str(dt));
    SOL_EXX = SOL_EX(:,2001);
    L2_ERR  = norm(SOL_EXX-SOL(:,indexError),2)*(dt_base*(2^i))^0.5;
    LINF_ERR = norm(SOL_EXX-SOL(:,indexError),Inf);
    l2_errors_t(i,1) = L2_ERR;
    lInf_errors_t(i,1) = LINF_ERR;
end 

disp("finito");

% order of convergence 
% L2 

E_12_h = l2_errors_h(1,1);
E_22_h = l2_errors_h(end,1);

ph2 = (log(E_12_h/E_22_h))/(log(h_base*2/(h_base*(2^(errorsMaxIndex)))));
conv2_val_h = zeros(errorsMaxIndex+1,1);
for i = 0:errorsMaxIndex
    conv2_val_h(i+1,1) = (h_base*(2^i))^(ph2);
end

E_1Inf_h = lInf_errors_h(1,1);
E_2Inf_h = lInf_errors_h(end,1);
phINF = log(E_1Inf_h/E_2Inf_h)/log(h_base*2/(h_base*(2^(errorsMaxIndex))));
convInf_val_h = zeros(errorsMaxIndex+1,1);
for i = 0:errorsMaxIndex
    convInf_val_h(i+1,1) = (h_base*(2^i))^(phINF);
end

E_12_t = l2_errors_t(1,1);
E_22_t = l2_errors_t(end,1);

pt2 = (log(E_12_t/E_22_t))/(log(dt_base*2/(dt_base*(2^(errorsMaxIndex)))));
conv2_val_t = zeros(errorsMaxIndex+1,1);
for i = 0:errorsMaxIndex
    conv2_val_t(i+1,1) = (dt_base*(2^i))^(pt2);
end

E_1Inf_t = lInf_errors_t(1,1);
E_2Inf_t = lInf_errors_t(end,1);
ptINF = log(E_1Inf_t/E_2Inf_t)/log(dt_base*2/(dt_base*(2^(errorsMaxIndex))));
convInf_val_t = zeros(errorsMaxIndex+1,1);
for i = 0:errorsMaxIndex-1
    convInf_val_t(i+1,1) = (dt_base*(2^i))^(ptINF);
end

hhtt= 0.0005*(2.^(0:3));
hhtt1= 0.0005*(2.^(1:3));
figure('Renderer', 'painters', 'Position', [100 100 1000 600]);

subplot(2,2,1)
loglog(hhtt1,l2_errors_h,"o-");
hold on 
loglog(hhtt,conv2_val_h,"o-");
title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
xlabel('h','FontSize',16,"Interpreter","latex");
ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
legend(["$||u-u_{ex}||_2$","$h_2^{q}$"],'FontSize',16,"Interpreter","latex")

%saveas(gcf, ".\plots\L2h_error_conv_mu1.png");
disp("q = "+num2str(2));
subplot(2,2,2)
loglog(hhtt1,lInf_errors_h,"o-");
hold on 
loglog(hhtt,convInf_val_h,"o-");
title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
xlabel('h','FontSize',16,"Interpreter","latex");
ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
legend(["$||u-u_{ex}||_{\infty}$","$h_{\infty}^{q}$"],'FontSize',16,"Interpreter","latex")
%saveas(gcf, ".\plots\LINFh_error_conv_mu1.png");
disp("q = "+num2str(2));

subplot(2,2,3)
loglog(hhtt1,l2_errors_t,"o-");
hold on 
loglog(hhtt,conv2_val_t,"o-");
title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
legend(["$||u-u_{ex}||_2$","$\Delta t_2^{p}$"],'FontSize',16,"Interpreter","latex")
%saveas(gcf, ".\plots\L2h_error_conv_mu1.png");
disp("q = "+num2str(2));

subplot(2,2,4)
loglog(hhtt1,lInf_errors_t,"o-");
hold on 
loglog(hhtt,convInf_val_t,"o-");
title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
legend(["$||u-u_{ex}||_{\infty}$","$\Delta t_{\infty}^{p}$"],'FontSize',16,"Interpreter","latex")
saveas(gcf, ".\plots\_error_conv_mu10_EX4.png");
disp("q = "+num2str(2));




%%
mu = 1

uj_t = @(t) -(1/4)*exp(-mu*t);

[SOL,dx,dt] = FD_1D_BURGER_FUN2_EX4(mu,T,I,NT, NX ,u0, u0_t, uj_t);

%REFERENCE SOLUTION

[SOL_EX,dxEx,dtEx] = FD_1D_BURGER_FUN2_EX4(mu,T,I,5000, 2000 ,u0, u0_t, uj_t);

%SOL_EX = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu);

figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
space = linspace(I(1),I(2),NX+1);
index1 = round(0.1 / dt);
index2 = round(0.5 / dt);
index3 = round(2.3 / dt);

amp1 = SOL(:,index1);
amp2 = SOL(:,index2);
amp3 = SOL(:,index3);

spaceEx = linspace(I(1),I(2),2000+1);
indexEx1 = round(0.1 / dtEx);
indexEx2 = round(0.5 / dtEx);
indexEx3 = round(2.3 / dtEx);

ampEx1 = SOL_EX(:,indexEx1);
ampEx2 = SOL_EX(:,indexEx2); 
ampEx3 = SOL_EX(:,indexEx3);


plot(space,amp1,'LineWidth',2);
hold on 
plot(spaceEx,ampEx1,"--",'LineWidth',2);
hold on
plot(space,amp2,'LineWidth',2);
hold on 
plot(spaceEx,ampEx2,"--",'LineWidth',2);
hold on
plot(space,amp3,'LineWidth',2);
hold on
plot(spaceEx,ampEx3,"--",'LineWidth',2);
title("Numerical to exact solution comparison",'FontSize',16,"Interpreter","latex")
xlabel('space x ','FontSize',16,"Interpreter","latex");
ylabel('$u(x)$','FontSize',16,"Interpreter","latex");
legend(["$T=0.1$","$T=0.1$ exact","$T=0.5$","$T=0.5$ exact","$T=2.3$","$T=2.3$ exact"],'FontSize',16,"Interpreter","latex")
saveas(gcf, ".\plots\amp_mu1_EX4.png");
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error and order of convergence 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 1

uj_t = @(t) -(1/4)*exp(-mu*t);

%L2_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),2)*dx^0.5;
%L1_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),1)*dx;
%LINF_ERR = 0;%norm(SOL_EX-SOL(:,NT+1),Inf);

errorsMaxIndex =3;


l2_errors_h = zeros(errorsMaxIndex,1);
lInf_errors_h = zeros(errorsMaxIndex,1);
l2_errors_t = zeros(errorsMaxIndex,1);
lInf_errors_t = zeros(errorsMaxIndex,1);
% varying dx
h_base = 0.0005;
dt_base = 0.0005;
 NT= T/dt_base;
 indexError = NT+1;
 [SOL_EX,dx,dt] = FD_1D_BURGER_FUN2_EX4(mu,T,I,5000, 2000 ,u0, u0_t, uj_t);
for i = 1:errorsMaxIndex
    
    NX = (I(2)-I(1))/(h_base*2^(i));
    
    [SOL,dx,dt] = FD_1D_BURGER_FUN2_EX4(mu,T,I,NT, NX ,u0, u0_t, uj_t);
    disp("dx="+num2str(dx));
    SOL_EXX = SOL_EX(1:2^(i):end,2001);
    L2_ERR  = norm(SOL_EXX-SOL(:,indexError),2)*(h_base*(2^i))^0.5;
    LINF_ERR = norm(SOL_EXX-SOL(:,indexError),Inf);
    l2_errors_h(i,1) = L2_ERR;
    lInf_errors_h(i,1) = LINF_ERR;
end 

NX = (I(2)-I(1))/(h_base);

for i = 1:errorsMaxIndex
    
    NT= T/(dt_base*(2^i));
    indexError = NT+1;
    [SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_t, uj_t);
    disp("dt="+num2str(dt));
    SOL_EXX = SOL_EX(:,2001);
    L2_ERR  = norm(SOL_EXX-SOL(:,indexError),2)*(dt_base*(2^i))^0.5;
    LINF_ERR = norm(SOL_EXX-SOL(:,indexError),Inf);
    l2_errors_t(i,1) = L2_ERR;
    lInf_errors_t(i,1) = LINF_ERR;
end 

disp("finito");

% order of convergence 
% L2 

E_12_h = l2_errors_h(1,1);
E_22_h = l2_errors_h(end,1);

ph2 = (log(E_12_h/E_22_h))/(log(h_base*2/(h_base*(2^(errorsMaxIndex)))));
conv2_val_h = zeros(errorsMaxIndex+1,1);
for i = 0:errorsMaxIndex
    conv2_val_h(i+1,1) = (h_base*(2^i))^(ph2);
end

E_1Inf_h = lInf_errors_h(1,1);
E_2Inf_h = lInf_errors_h(end,1);
phINF = log(E_1Inf_h/E_2Inf_h)/log(h_base*2/(h_base*(2^(errorsMaxIndex))));
convInf_val_h = zeros(errorsMaxIndex+1,1);
for i = 0:errorsMaxIndex
    convInf_val_h(i+1,1) = (h_base*(2^i))^(phINF);
end

E_12_t = l2_errors_t(1,1);
E_22_t = l2_errors_t(end,1);

pt2 = (log(E_12_t/E_22_t))/(log(dt_base*2/(dt_base*(2^(errorsMaxIndex)))));
conv2_val_t = zeros(errorsMaxIndex+1,1);
for i = 0:errorsMaxIndex
    conv2_val_t(i+1,1) = (dt_base*(2^i))^(pt2);
end

E_1Inf_t = lInf_errors_t(1,1);
E_2Inf_t = lInf_errors_t(end,1);
ptINF = log(E_1Inf_t/E_2Inf_t)/log(dt_base*2/(dt_base*(2^(errorsMaxIndex))));
convInf_val_t = zeros(errorsMaxIndex+1,1);
for i = 0:errorsMaxIndex-1
    convInf_val_t(i+1,1) = (dt_base*(2^i))^(ptINF);
end

hhtt= 0.0005*(2.^(0:3));
hhtt1= 0.0005*(2.^(1:3));
figure('Renderer', 'painters', 'Position', [100 100 1000 600]);

subplot(2,2,1)
loglog(hhtt1,l2_errors_h,"o-");
hold on 
loglog(hhtt,conv2_val_h,"o-");
title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
xlabel('h','FontSize',16,"Interpreter","latex");
ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
legend(["$||u-u_{ex}||_2$","$h_2^{q}$"],'FontSize',16,"Interpreter","latex")

%saveas(gcf, ".\plots\L2h_error_conv_mu1.png");
disp("q = "+num2str(2));
subplot(2,2,2)
loglog(hhtt1,lInf_errors_h,"o-");
hold on 
loglog(hhtt,convInf_val_h,"o-");
title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
xlabel('h','FontSize',16,"Interpreter","latex");
ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
legend(["$||u-u_{ex}||_{\infty}$","$h_{\infty}^{q}$"],'FontSize',16,"Interpreter","latex")
%saveas(gcf, ".\plots\LINFh_error_conv_mu1.png");
disp("q = "+num2str(2));

subplot(2,2,3)
loglog(hhtt1,l2_errors_t,"o-");
hold on 
loglog(hhtt,conv2_val_t,"o-");
title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
legend(["$||u-u_{ex}||_2$","$\Delta t_2^{p}$"],'FontSize',16,"Interpreter","latex")
%saveas(gcf, ".\plots\L2h_error_conv_mu1.png");
disp("q = "+num2str(2));

subplot(2,2,4)
loglog(hhtt1,lInf_errors_t,"o-");
hold on 
loglog(hhtt,convInf_val_t,"o-");
title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
legend(["$||u-u_{ex}||_{\infty}$","$\Delta t_{\infty}^{p}$"],'FontSize',16,"Interpreter","latex")
saveas(gcf, ".\plots\_error_conv_mu1_EX4.png");
disp("q = "+num2str(2));
%%
mu = 0.1

uj_t = @(t) -(1/4)*exp(-mu*t);

[SOL,dx,dt] = FD_1D_BURGER_FUN2_EX4(mu,T,I,NT, NX ,u0, u0_t, uj_t);

%REFERENCE SOLUTION

[SOL_EX,dxEx,dtEx] = FD_1D_BURGER_FUN2_EX4(mu,T,I,5000, 2000 ,u0, u0_t, uj_t);

%SOL_EX = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu);

figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
space = linspace(I(1),I(2),NX+1);
index1 = round(0.1 / dt);
index2 = round(0.5 / dt);
index3 = round(2.3 / dt);

amp1 = SOL(:,index1);
amp2 = SOL(:,index2);
amp3 = SOL(:,index3);

spaceEx = linspace(I(1),I(2),2000+1);
indexEx1 = round(0.1 / dtEx);
indexEx2 = round(0.5 / dtEx);
indexEx3 = round(2.3 / dtEx);

ampEx1 = SOL_EX(:,indexEx1);
ampEx2 = SOL_EX(:,indexEx2); 
ampEx3 = SOL_EX(:,indexEx3);


plot(space,amp1,'LineWidth',2);
hold on 
plot(spaceEx,ampEx1,"--",'LineWidth',2);
hold on
plot(space,amp2,'LineWidth',2);
hold on 
plot(spaceEx,ampEx2,"--",'LineWidth',2);
hold on
plot(space,amp3,'LineWidth',2);
hold on
plot(spaceEx,ampEx3,"--",'LineWidth',2);
title("Numerical to exact solution comparison",'FontSize',16,"Interpreter","latex")
xlabel('space x ','FontSize',16,"Interpreter","latex");
ylabel('$u(x)$','FontSize',16,"Interpreter","latex");
legend(["$T=0.1$","$T=0.1$ exact","$T=0.5$","$T=0.5$ exact","$T=2.3$","$T=2.3$ exact"],'FontSize',16,"Interpreter","latex")
saveas(gcf, ".\plots\amp_mu01_EX4.png");
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error and order of convergence 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 0.1

uj_t = @(t) -(1/4)*exp(-mu*t);

%L2_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),2)*dx^0.5;
%L1_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),1)*dx;
%LINF_ERR = 0;%norm(SOL_EX-SOL(:,NT+1),Inf);

errorsMaxIndex =3;


l2_errors_h = zeros(errorsMaxIndex,1);
lInf_errors_h = zeros(errorsMaxIndex,1);
l2_errors_t = zeros(errorsMaxIndex,1);
lInf_errors_t = zeros(errorsMaxIndex,1);
% varying dx
h_base = 0.0005;
dt_base = 0.0005;
 NT= T/dt_base;
 indexError = NT+1;
 [SOL_EX,dx,dt] = FD_1D_BURGER_FUN2_EX4(mu,T,I,5000, 2000 ,u0, u0_t, uj_t);
for i = 1:errorsMaxIndex
    
    NX = (I(2)-I(1))/(h_base*2^(i));
    
    [SOL,dx,dt] = FD_1D_BURGER_FUN2_EX4(mu,T,I,NT, NX ,u0, u0_t, uj_t);
    disp("dx="+num2str(dx));
    SOL_EXX = SOL_EX(1:2^(i):end,2001);
    L2_ERR  = norm(SOL_EXX-SOL(:,indexError),2)*(h_base*(2^i))^0.5;
    LINF_ERR = norm(SOL_EXX-SOL(:,indexError),Inf);
    l2_errors_h(i,1) = L2_ERR;
    lInf_errors_h(i,1) = LINF_ERR;
end 

NX = (I(2)-I(1))/(h_base);

for i = 1:errorsMaxIndex
    
    NT= T/(dt_base*(2^i));
    indexError = NT+1;
    [SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_t, uj_t);
    disp("dt="+num2str(dt));
    SOL_EXX = SOL_EX(:,2001);
    L2_ERR  = norm(SOL_EXX-SOL(:,indexError),2)*(dt_base*(2^i))^0.5;
    LINF_ERR = norm(SOL_EXX-SOL(:,indexError),Inf);
    l2_errors_t(i,1) = L2_ERR;
    lInf_errors_t(i,1) = LINF_ERR;
end 

disp("finito");

% order of convergence 
% L2 

E_12_h = l2_errors_h(1,1);
E_22_h = l2_errors_h(end,1);

ph2 = (log(E_12_h/E_22_h))/(log(h_base*2/(h_base*(2^(errorsMaxIndex)))));
conv2_val_h = zeros(errorsMaxIndex+1,1);
for i = 0:errorsMaxIndex
    conv2_val_h(i+1,1) = (h_base*(2^i))^(ph2);
end

E_1Inf_h = lInf_errors_h(1,1);
E_2Inf_h = lInf_errors_h(end,1);
phINF = log(E_1Inf_h/E_2Inf_h)/log(h_base*2/(h_base*(2^(errorsMaxIndex))));
convInf_val_h = zeros(errorsMaxIndex+1,1);
for i = 0:errorsMaxIndex
    convInf_val_h(i+1,1) = (h_base*(2^i))^(phINF);
end

E_12_t = l2_errors_t(1,1);
E_22_t = l2_errors_t(end,1);

pt2 = (log(E_12_t/E_22_t))/(log(dt_base*2/(dt_base*(2^(errorsMaxIndex)))));
conv2_val_t = zeros(errorsMaxIndex+1,1);
for i = 0:errorsMaxIndex
    conv2_val_t(i+1,1) = (dt_base*(2^i))^(pt2);
end

E_1Inf_t = lInf_errors_t(1,1);
E_2Inf_t = lInf_errors_t(end,1);
ptINF = log(E_1Inf_t/E_2Inf_t)/log(dt_base*2/(dt_base*(2^(errorsMaxIndex))));
convInf_val_t = zeros(errorsMaxIndex+1,1);
for i = 0:errorsMaxIndex-1
    convInf_val_t(i+1,1) = (dt_base*(2^i))^(ptINF);
end

hhtt= 0.0005*(2.^(0:3));
hhtt1= 0.0005*(2.^(1:3));
figure('Renderer', 'painters', 'Position', [100 100 1000 600]);

subplot(2,2,1)
loglog(hhtt1,l2_errors_h,"o-");
hold on 
loglog(hhtt,conv2_val_h,"o-");
title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
xlabel('h','FontSize',16,"Interpreter","latex");
ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
legend(["$||u-u_{ex}||_2$","$h_2^{q}$"],'FontSize',16,"Interpreter","latex")

%saveas(gcf, ".\plots\L2h_error_conv_mu1.png");
disp("q = "+num2str(2));
subplot(2,2,2)
loglog(hhtt1,lInf_errors_h,"o-");
hold on 
loglog(hhtt,convInf_val_h,"o-");
title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
xlabel('h','FontSize',16,"Interpreter","latex");
ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
legend(["$||u-u_{ex}||_{\infty}$","$h_{\infty}^{q}$"],'FontSize',16,"Interpreter","latex")
%saveas(gcf, ".\plots\LINFh_error_conv_mu1.png");
disp("q = "+num2str(2));

subplot(2,2,3)
loglog(hhtt1,l2_errors_t,"o-");
hold on 
loglog(hhtt,conv2_val_t,"o-");
title("$L^2$ norm and convergence error",'FontSize',16,"Interpreter","latex")
xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
ylabel('$L^2$ error','FontSize',16,"Interpreter","latex");
legend(["$||u-u_{ex}||_2$","$\Delta t_2^{p}$"],'FontSize',16,"Interpreter","latex")
%saveas(gcf, ".\plots\L2h_error_conv_mu1.png");
disp("q = "+num2str(2));

subplot(2,2,4)
loglog(hhtt1,lInf_errors_t,"o-");
hold on 
loglog(hhtt,convInf_val_t,"o-");
title("$L^{\infty}$ norm and convergence error",'FontSize',16,"Interpreter","latex")
xlabel('$\Delta t$','FontSize',16,"Interpreter","latex");
ylabel('$L^{\infty}$ error','FontSize',16,"Interpreter","latex");
legend(["$||u-u_{ex}||_{\infty}$","$\Delta t_{\infty}^{p}$"],'FontSize',16,"Interpreter","latex")
saveas(gcf, ".\plots\_error_conv_mu01_EX4.png");
disp("q = "+num2str(2));

