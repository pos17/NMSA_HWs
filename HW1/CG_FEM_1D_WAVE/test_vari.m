close all; clear; clc

x = linspace(-1, 1, 3);

x0=[0,1];
hat1 = [1,0];
hat2 = [0,1,0]

figure('Renderer', 'painters', 'Position', [100 100 1200 300]);
plot(x0, hat1, "--k", LineWidth=3)
text(0, 1.05, ["\phi_{0}"], FontSize=14);
for ii=1:3
    hold on 
    plot(x+ii, hat2, LineWidth=3);
    text(ii, 1.05, ["\phi_{"+num2str(ii)+"}"], FontSize=14);
end
hold on 
plot(x+9, hat2, LineWidth=3);
text(9-0.4, 1.05, ["\phi_{N_h-1}"], FontSize=14);
hold on 
plot(x0+9, flip(hat1), LineWidth=3);
text(10-0.4, 1.05, ["\phi_{N_h}"], FontSize=14);
ylim([0,1.1])
yticks([0,1])
xticks([0,1,2,3,6,9,10])
xticklabels(["x_0", "x_1", "x_2", "x_3","...", "x_{N_h-1}", "x_{N_h}"])
grid on

saveas(gcf, "./ReportAssets/BasisFunctions.png")




%%
test = ["TestHW1E4S0";"TestHW1E4S1";"TestHW1E4S2"]
for ii=1:3
    [errors,solutions,femregion,Dati,u_snp] = C_main1D(test(ii), 7, 0.0001)
    close all
end
%%
test = ["TestHW1E5S0";"TestHW1E5S1";"TestHW1E5S2"]
for ii=1:3
    [errors,solutions,femregion,Dati,u_snp] = C_main1D(test(ii), 7, 0.0001)
    close all
end

%%
clear;close all;
t=linspace(0,1,1000);
g=t.*0;
for ii=1:length(t)
    if abs(t(ii)-0.2)<=0.1
        g(ii)=.5*(1+cos(pi*(t(ii)-0.2)/0.1));
    end
end

figure('Renderer', 'painters', 'Position', [100 100 800 250]);
plot(t, g, LineWidth=3);
xlim([0,1])
ylim([-0.2, 1.2])
xlabel("Time t", FontSize=14)
ylabel("Amplitude", FontSize=14)
title("Half Cosine Impulse", FontSize=20)
grid 

