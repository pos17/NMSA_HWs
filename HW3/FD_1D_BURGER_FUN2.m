function [SOL,dx,dt] = FD_1D_BURGER_FUN2(mu,T,I,NT, NX ,u0, u0_t, uj_t)
    disp("Burgers_fun_2_start");
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
        SOL(NX+1,n)= u0_t(n);
    end

    for n = 1:NT+1
        SOL(NX+1,n)= uj_t(n);
    end

    for n = 3:NT+1
        fj = zeros(NX+3,1);
        Aj_vec = ones(NX+3,1);
        Aj = diag(Aj_vec);
        fj(2,1) = SOL(1,n);
        fj(NX+2,1) = SOL(NX+1,n);
        for j = 3:NX+1
            fj(j,1) = (4/3 * SOL(j-1,n-1))-(1/3*SOL(j-1,n-2));
            
            gammai = -(2/3)*lambda* ((mu/dx)+(2*SOL(j-1,n-1)-SOL(j-1,n-2))/2);
            alfai = 1+(4/3*lambda*(mu/dx));
            betai = -(2/3)*lambda* ((mu/dx)-(2*SOL(j-1,n-1)-SOL(j-1,n-2))/2);
    
            Aj(j,j-1) = gammai;
            Aj(j,j) = alfai;
            Aj(j,j+1) = betai;
                   
        end
        SOLj = Aj\fj;
        SOL(:,n) = SOLj(2:NX+2,1);

    end
    disp("Burgers_fun_2_end");
end

