close all; clear; clc;
tests = ["TestHW1E3S0";"TestHW1E3S1";"TestHW1E3S2"]
dt = [0.0005,0.0001]

for tt=1:2
    for ii=1:3
        [errors_table,rates] = C_convergence_test(tests(ii), dt(tt));
    end
end