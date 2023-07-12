
function [errors_table,rates] = C_convergence_test(test_name)
%% [errors_table,rates]=C_convergence_test(test_name)
%==========================================================================
% Error analysis varying the mesh size h 
%==========================================================================
% Example of usage: [errors_table,rates] = C_mesh_test('Test1')
%
%    INPUT:
%          test_name    : (string)  test case name, see C_dati.m
%
%    OUTPUT:
%          errors_table : (struct) containing the computed errors
%          rates        : (struct) containing the computed rates

close all;
warning off;
addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing


Dati=C_dati(test_name);

refinement_vector=Dati.refinement_vector;
poly_vector = ["P2","P3", "P4", "P5"];
t_vector = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001];
deltaT = t_vector(6);
poly = poly_vector(4);
num_test=length(refinement_vector);

for k=1:num_test
    % refinement
    [errors,solutions,femregion,Dati] = C_main1D(test_name, refinement_vector(k), poly, deltaT);
    
    % polynomial
%     [errors,solutions,femregion,Dati] = C_main1D(test_name,refinement_vector(k));
    
    % delta t
%     [errors,solutions,femregion,Dati] = C_main1D(test_name,refinement_vector(k));
    Error_L2(k)=errors.Error_L2;

%     Error_SEMI_H1(k)=errors.Error_SEMI_H1;
%     Error_H1(k)=errors.Error_H1;
%     Error_inf(k)=errors.Error_inf;

    ne(k)=femregion.ne;    
    h(k)=femregion.h;
    fprintf('==========================================\n');    
    fprintf('End test %i\n',k);
    fprintf('==========================================\n');
end
p=femregion.degree;

%ERROR TABLE
errors_table=struct('ne', ne,...
                   'h', h,...
                   'Error_L2', Error_L2);
%                    'Error_SEMI_H1', Error_SEMI_H1,...
%                    'Error_H1', Error_H1,...
%                    'Error_inf',Error_inf ...
%                    );


               
%TABLE
rate_L2=log10(Error_L2(2:num_test)./Error_L2(1:num_test-1))./log10(h(2:num_test)./h(1:num_test-1));
% rate_SEMI_H1=log10( Error_SEMI_H1(2:num_test)./ Error_SEMI_H1(1:num_test-1))./log10(h(2:num_test)./h(1:num_test-1));
% rate_H1=log10( Error_H1(2:num_test)./ Error_H1(1:num_test-1))./log10(h(2:num_test)./h(1:num_test-1));
% rate_inf=log10( Error_inf(2:num_test)./ Error_inf(1:num_test-1))./log10(h(2:num_test)./h(1:num_test-1));

rates=struct('rate_L2',rate_L2);
%              'rate_SEMI_H1',rate_SEMI_H1,...
%              'rate_H1',rate_H1,...
%              'rate_inf',rate_inf ...
%              );
         
%%
% ERROR PLOTS
figure()
% hs = subplot(2,1,1);
loglog(h,h.^(p+1),'-+b','Linewidth',2);
hold on
loglog(h,Error_L2,'-or','Linewidth',2);
hold on
legend(sprintf('h^%i',p+1),'||.||_{L^2}');
% legend('||.||_{L^2}');
grid minor
title (["Errors with \Delta t = "+num2str(deltaT)+" s"]);
ylabel('L^2-error')
xlabel('h');
filename = ["Errors_t"+num2str(deltaT)+"_poly"+poly];
saveas(gcf, ["./convergence_plot/"+filename+".png"]);

% hs.FontSize = 12;

% hs = subplot(2,1,2);
% loglog(h,h.^(p),'-+b','Linewidth',2);
% hold on
% loglog(h,Error_H1,'-or','Linewidth',2);
% hold on
% legend(sprintf('h^%i',p),'||.||_{H^1}');
% ylabel('H^1-error')
% xlabel('h');
% hold off
% hs.FontSize = 12;

 