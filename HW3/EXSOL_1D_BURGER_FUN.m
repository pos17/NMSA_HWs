function [SOL_EX] = EXSOL_1D_BURGER_FUN(I,T,NX,NT,mu,indexError)
    disp("EXACT SOL START");
    dt = T/NT;
    dx = (I(2)-I(1))/NX;
    lambda = dt/dx;

    
    n_harm = 20;

    c_0_fun = @(x) exp((-1/(2*mu*pi))*(1-cos(pi*x)));
    c_0 = integral(c_0_fun,I(1),I(2));
    c_n_fun = @(n,x) exp((-1/(2*mu*pi))*(1-cos(pi*x))).*cos(n.*pi*x);
    
    
    SOL_EX = zeros(NX+1,1);
    t = (indexError-1)*dt;
    sum_den = zeros(NX+1,1);
    sum_num = zeros(NX+1,1);
    for nn = 1: n_harm
        c_n = 2 * integral(@(x) c_n_fun(nn,x),I(1),I(2)); 
        for j = 1:NX+1
            x = I(1) + (j-1)*dx;    
            sum_num(j,1) = sum_num(j,1) + (c_n*nn*exp(-nn^2*pi^2*mu*t)*sin(nn*pi*x));
            sum_den(j,1) = sum_den(j,1) +  (c_n*exp(-nn^2*pi^2*mu*t)*cos(nn*pi*x));
        end   
        
    end
    u_ex = (2* mu * pi )* (sum_num)./(c_0+(sum_den));
    SOL_EX(:,1) = u_ex;
    
    disp("EXACT SOL END");
end