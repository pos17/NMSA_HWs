close all;
clear;
clc;

tests=["TestHW1E6S0";"TestHW1E6S1";"TestHW1E6S2"];
figure(1)
for ii=1:3
    [femregion, Dati] = C_main1D_Ste(tests(ii),7)
    hold on
end
legend(["S_1(x)","S_2(x)","S_3(x)"], FontSize=12)
title("Frequency Response of the Acoustic Potential", FontSize=20)