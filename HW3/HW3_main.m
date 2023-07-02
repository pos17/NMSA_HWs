clear; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = 10


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



[SOL,SOL_EX,L1_ERR,L2_ERR,LINF_ERR,dx,dt] = FD_1D_BURGER_FUN(mu,T,I,NT, NX ,u0, u0_x, uj_x,1,NT+1);

figure(1)
space = linspace(I(1),I(2),NX+1);
index1 = round(0.1 / dt);
index2 = round(0.5 / dt);
index3 = round(2.3 / dt);

amp1 = SOL(:,index1);
amp2 = SOL(:,index2);
amp3 = SOL(:,index3);

ampEx1 = SOL_EX(:,index1);
ampEx2 = SOL_EX(:,index2);
ampEx3 = SOL_EX(:,index3);

plot(space,amp1);
hold on 
plot(space,ampEx1,"--");
hold on
plot(space,amp2);
hold on 
plot(space,ampEx2,"--");
hold on
plot(space,amp3);
hold on
plot(space,ampEx3,"--");
xlabel('space x ','FontSize',16);
ylabel('u(x)','FontSize',16);
legend(["T=0.1","T=0.1 exact","T=0.5","T=0.5 exact","T=2.3","T=2.3 exact"])
saveas(gcf, ".\plots\amp_mu10.png");


mu = 1

[SOL,SOL_EX,L1_ERR,L2_ERR,LINF_ERR,dx,dt] = FD_1D_BURGER_FUN(mu,T,I,NT, NX ,u0, u0_x, uj_x,1,NT+1);

figure(2)
space = linspace(I(1),I(2),NX+1);
index1 = round(0.1 / dt);
index2 = round(0.5 / dt);
index3 = round(2.3 / dt);

amp1 = SOL(:,index1);
amp2 = SOL(:,index2);
amp3 = SOL(:,index3);

ampEx1 = SOL_EX(:,index1);
ampEx2 = SOL_EX(:,index2);
ampEx3 = SOL_EX(:,index3);

plot(space,amp1);
hold on 
plot(space,ampEx1,"--");
hold on
plot(space,amp2);
hold on 
plot(space,ampEx2,"--");
hold on
plot(space,amp3);
hold on
plot(space,ampEx3,"--");
xlabel('space x ','FontSize',16);
ylabel('u(x)','FontSize',16);
legend(["T=0.1","T=0.1 exact","T=0.5","T=0.5 exact","T=2.3","T=2.3 exact"])
saveas(gcf, ".\plots\amp_mu1.png");

mu = 0.1

[SOL,SOL_EX,L1_ERR,L2_ERR,LINF_ERR,dx,dt] = FD_1D_BURGER_FUN(mu,T,I,NT, NX ,u0, u0_x, uj_x,1,NT+1);

figure(3)
space = linspace(I(1),I(2),NX+1);
index1 = round(0.1 / dt);
index2 = round(0.5 / dt);
index3 = round(2.3 / dt);

amp1 = SOL(:,index1);
amp2 = SOL(:,index2);
amp3 = SOL(:,index3);

ampEx1 = SOL_EX(:,index1);
ampEx2 = SOL_EX(:,index2);
ampEx3 = SOL_EX(:,index3);

plot(space,amp1);
hold on 
plot(space,ampEx1,"--");
hold on
plot(space,amp2);
hold on 
plot(space,ampEx2,"--");
hold on
plot(space,amp3);
hold on
plot(space,ampEx3,"--");
xlabel('space x ','FontSize',16);
ylabel('u(x)','FontSize',16);
legend(["T=0.1","T=0.1 exact","T=0.5","T=0.5 exact","T=2.3","T=2.3 exact"])
saveas(gcf, ".\plots\amp_mu01.png");


L2_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),2)*dx^0.5;
L1_ERR   = 0;%norm(SOL_EX-SOL(:,NT+1),1)*dx;
LINF_ERR = 0;%norm(SOL_EX-SOL(:,NT+1),Inf);
errorsMaxIndex =4;
l2_errors = zeros(errorsMaxIndex,1);
lInf_errors = zeros(errorsMaxIndex,1);
% varying dx
for i = 1:errorsMaxIndex
    dx_base = 0.0005;
    NT=2000;
    NX = (I(2)-I(1))/dx_base*i;
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