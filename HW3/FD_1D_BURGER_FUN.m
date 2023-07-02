function [SOL,SOL_EX,L1_ERR,L2_ERR,LINF_ERR,dx,dt] = FD_1D_BURGER_FUN(mu,T,I,NT, NX ,u0, u0_x, uj_x, u_exf, indexError)
    %FD_1D_BURGER_FUN Summary of this function goes here
    %   Detailed explanation goes here
    
    dt = T/NT;
    dx = (I(2)-I(1))/NX;
    lambda = dt/dx;

    SOL = zeros(NX+1,NT+1);
    
    for j = 1:NX+1
         SOL(j,1) = u0(I(1) + (j-1)*dx);
    end
    
    for j = 2:NX
        SOL(j,2) = SOL(j,1) + (lambda * ((mu/dx)*(SOL(j+1,1)-2*SOL(j,1)+SOL(j-1,1)))- (SOL(j,1)/2)*(SOL(j+1,1)-SOL(j-1,1)));
    end
    
    for n = 1:NT+1
        SOL(NX+1,n)= u0_x(n);
    end

    for n = 1:NT+1
        SOL(NX+1,n)= uj_x(n);
    end

    for n = 3:NT+1
        fj = zeros(NX+1,1);
        Aj_vec = ones(NX+1,1);
        Aj = diag(Aj_vec);
        fj(1,1) = SOL(1,n);
        fj(NX+1,1) = SOL(NX+1,n);
        for j = 2:NX
            fj(j,1) = (4/3 * SOL(j,n-1))-(1/3*SOL(j,n-2));
            
            gammai = -(2/3)*lambda* ((mu/dx)+(2*SOL(j,n-1)-SOL(j,n-2))/2);
            alfai = 1+(4/3*lambda*(mu/dx));
            betai = -(2/3)*lambda* ((mu/dx)-(2*SOL(j,n-1)-SOL(j,n-2))/2);
    
            Aj(j,j-1) = gammai;
            Aj(j,j) = alfai;
            Aj(j,j+1) = betai;
                   
        end
        SOLj = Aj\fj;
        SOL(:,n) = SOLj;

    end

    if u_exf == 1
        
        n_harm = 10;
    
        c_0_fun = @(x) exp((-1/(2*mu*pi))*(1-cos(pi*x)));
        c_0 = integral(c_0_fun,I(1),I(2));
        c_n_fun = @(n,x) exp((-1/(2*mu*pi))*(1-cos(pi*x))).*cos(n.*pi*x);
        
        
        SOL_EX = zeros(NX+1,NT+1);
        for n = 1:NT+1
            t = (n-1)*dt;
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
            SOL_EX(:,n) = u_ex;
        end

        L2_ERR  = norm(SOL_EX-SOL(:,indexError),2)*dx^0.5;
        L1_ERR   = norm(SOL_EX-SOL(:,indexError),1)*dx;
        LINF_ERR = norm(SOL_EX-SOL(:,indexError),Inf);
    elseif u_exf == 2

        disp("u_exf = 2");
        n_harm = 10;
    
        c_0_fun = @(x) exp((-1/(2*mu*pi))*(1-cos(pi*x)));
        c_0 = integral(c_0_fun,I(1),I(2));
        c_n_fun = @(n,x) exp((-1/(2*mu*pi))*(1-cos(pi*x))).*cos(n.*pi*x);
        
        
        SOL_EX = zeros(NX+1,1);
        t = (indexError)*dt;
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
        SOL_EX(:,n) = u_ex;
        L2_ERR  = norm(SOL_EX-SOL(:,indexError),2)*dx^0.5;
        L1_ERR   = norm(SOL_EX-SOL(:,indexError),1)*dx;
        LINF_ERR = norm(SOL_EX-SOL(:,indexError),Inf);
    elseif u_exf == 3

        dt = 0.0005;
        dx = 0.0005;
        lambda = dt/dx;
        NX = I(2)-I(1)/dx;
        NT = T/dt;
    
        SOL_EX = zeros(NX+1,NT+1);
        
        for j = 1:NX+1
             SOL_EX(j,1) = u0(I(1) + (j-1)*dx);
        end
        
        for j = 2:NX
            SOL_EX(j,2) = SOL_EX(j,1) + (lambda * ((mu/dx)*(SOL_EX(j+1,1)-2*SOL_EX(j,1)+SOL_EX(j-1,1)))- (SOL_EX(j,1)/2)*(SOL_EX(j+1,1)-SOL_EX(j-1,1)));
        end
        
        for n = 1:NT+1
            SOL_EX(NX+1,n)= u0_x(n);
        end
    
        for n = 1:NT+1
            SOL_EX(NX+1,n)= uj_x(n);
        end
    
        for n = 3:NT+1
            fj = zeros(NX+1,1);
            Aj_vec = ones(NX+1,1);
            Aj = diag(Aj_vec);
            fj(1,1) = SOL_EX(1,n);
            fj(NX+1,1) = SOL_EX(NX+1,n);
            for j = 2:NX
                fj(j,1) = (4/3 * SOL_EX(j,n-1))-(1/3*SOL_EX(j,n-2));
                
                gammai = -(2/3)*lambda* ((mu/dx)+(2*SOL_EX(j,n-1)-SOL_EX(j,n-2))/2);
                alfai = 1+(4/3*lambda*(mu/dx));
                betai = -(2/3)*lambda* ((mu/dx)-(2*SOL_EX(j,n-1)-SOL_EX(j,n-2))/2);
        
                Aj(j,j-1) = gammai;
                Aj(j,j) = alfai;
                Aj(j,j+1) = betai;
                       
            end
            SOLj = Aj\fj;
            SOL_EX(:,n) = SOLj;
    
        end
    end
    



end

