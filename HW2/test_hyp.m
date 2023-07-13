clear; close all; %clc;

dx = [0.001, 0.0005, 0.001, 0.0005, 0.001, 0.0005];
dt = [0.001, 0.0005, 0.0005, 0.00025, 0.00025, 0.000125];

% err = zeros(2, length(dx));
err = struct();
for ii = 1:length(dx)
    L2err = hyp_HW2(dx(ii), dt(ii));
%     field1 = "dx="+num2str(dx(ii))+", dt="+num2str(dt(ii));
    err.dx(ii) = dx(ii);
    err.dt(ii) = dt(ii);
    err.ut(ii) = L2err(1);
    err.ux(ii) = L2err(2);
end

err
close all
%%
L2err = hyp_HW2(dx(1), dt(1));
