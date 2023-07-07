function [SOL,dx,dt] = FD_1D_BURGER_FUN2_EX4(mu,T,I,NT, NX ,u0, u0_t, uj_t)
    disp("Burgers_fun_2_start");
    %FD_1D_BURGER_FUN Summary of this function goes here
    %   Detailed explanation goes here
    
    dt = T/NT;
    dx = (I(2)-I(1))/NX;
    lambda = dt/dx;
    
    SOL = zeros(NX+3,NT+1);
    for j = 2:NX+2
         SOL(j,1) = u0(I(1) + (j-1)*dx);
    end
    SOL(1,1) = SOL(2,1);
    SOL(NX+3,1) = SOL(NX+2,1);

    SOL(2,2) = SOL(2,1) + (lambda * ((mu/dx)*(SOL(3,1)-2*SOL(2,1)+SOL(1,1))));
    
    for j = 3:NX+2
        SOL(j,2) = SOL(j,1) + (lambda * ((mu/dx)*(SOL(j+1,1)-2*SOL(j,1)+SOL(j-1,1)))- (SOL(j,1)/2)*(SOL(j+1,1)-SOL(j-1,1)));
    end
    SOL(1,2) = SOL(2,2);
    SOL(NX+3,2) = SOL(NX+2,2);
    
    %SET UP OF THE LIFTING FUNCTION
    R= zeros(NX+3,NT+1);

    for n = 1:NT+1
        R(NX+2,n)= uj_t((n-1)*dt);
        R(NX+3,n) = R(NX+2,n);
    end

    % iteration
    
    for n = 3:NT+1
        fj = zeros(NX+3,1);
        Aj_vec = ones(NX+3,1);
        Aj = diag(Aj_vec);
        fj(2,1) = ((4/3) *SOL(1,n)) - ((1/3) * SOL(1,n-1));
        
        %fj(NX+2,1) = SOL(NX+1,n);
        gamma0 = -(2/3)*dt*(mu/(dx^2));
        alfa0 = 1+ (4/3)*dt*(mu/(dx^2));
        beta0 = -(2/3)*dt*(mu/(dx^2));
        
        Aj(1,2) = -1;
        
        Aj(2,1) = gamma0;
        Aj(2,2) = alfa0;
        Aj(2,3) = beta0;
        
        for j = 3:NX+1
            gammai = -(2/3)*lambda* ((mu/dx)+(2*SOL(j-1,n-1)-SOL(j-1,n-2))/2);
            alfai = 1+(4/3*lambda*(mu/dx));
            betai = -(2/3)*lambda* ((mu/dx)-(2*SOL(j-1,n-1)-SOL(j-1,n-2))/2);
            fj(j,1) = (4/3 * SOL(j,n-1))-betai*R(j+1,n);
            %fj(j,1) = (4/3 * SOL(j,n-1))-(1/3*SOL(j,n-2));
            Aj(j,j-1) = gammai;
            Aj(j,j) = alfai;
            Aj(j,j+1) = betai;      
        end

        %dirichlet
        Aj(end-1,:) = 0;
        Aj(:,end-1) = 0;
        Aj(end-1,end-1) = 1;
        fj(end-1,1)=0;

        %Aj(NX+1,NX+2) = 0;
        %Aj(end,end-1)=-1;
        %fj(j,1) = (4/3 * SOL(j-1,n-1))-(1/3*SOL(j-1,n-2));
        %fj(end,1) =0;
        %fj(end-1,1)=0;
        SOLj = Aj\fj;
        SOLj = SOLj + R(:,n);
        SOL(:,n) = SOLj(:,1);
    end
    SOL = SOL(2:end-1,:); 
    disp("Burgers_fun_2_end");


    %%
%     SOL = zeros(NX+1,NT+1);
%     
%     for j = 1:NX+1
%          SOL(j,1) = u0(I(1) + (j-1)*dx);
%     end
% 
%     SOL(2,2) = SOL(2,1) + (lambda * ((mu/dx)*(SOL(3,1)-2*SOL(2,1)+SOL(1,1))));
%     for j = 3:NX
%         SOL(j,2) = SOL(j,1) + (lambda * ((mu/dx)*(SOL(j+1,1)-2*SOL(j,1)+SOL(j-1,1)))- (SOL(j,1)/2)*(SOL(j+1,1)-SOL(j-1,1)));
%     end
%     SOL(1,2) = SOL(2,2);
% 
% 
%     for n = 1:NT+1
%         SOL(NX+1,n)= uj_t(n*dt);
%     end
% 
%     
% 
%     for n = 3:NT+1
%         fj = zeros(NX+3,1);
%         Aj_vec = ones(NX+3,1);
%         Aj = diag(Aj_vec);
%         fj(2,1) = ((4/3) *SOL(1,n)) - ((1/3) * SOL(1,n-1));
%         fj(NX+2,1) = SOL(NX+1,n);
%         gamma0 = -(2/3)*dt*(mu/(dx^2));
%         alfa0 = 1+ (4/3)*dt*(mu/(dx^2));
%         beta0 = -(2/3)*dt*(mu/(dx^2));
%         Aj(1,2) = -1;
%         
%         Aj(2,1) = gamma0;
%         Aj(2,2) = alfa0;
%         Aj(2,3) = beta0;
%         
%         for j = 3:NX+1
%             
%             
%             
%             gammai = -(2/3)*lambda* ((mu/dx)+(2*SOL(j-1,n-1)-SOL(j-1,n-2))/2);
%             alfai = 1+(4/3*lambda*(mu/dx));
%             betai = -(2/3)*lambda* ((mu/dx)-(2*SOL(j-1,n-1)-SOL(j-1,n-2))/2);
%             
%             if(j == NX+1)
%                 fj(j,1) = (4/3 * SOL(j,n-1))-(1/3*SOL(j,n-2)) - (betai*SOL(NX+1,n));
%             else 
%                 fj(j,1) = (4/3 * SOL(j,n-1))-(1/3*SOL(j,n-2));
%             end
% 
%             Aj(j,j-1) = gammai;
%             Aj(j,j) = alfai;
%             Aj(j,j+1) = betai;
%                    
%         end
% 
%         %Aj(NX+1,NX+2) = 0;
%         %Aj(end,end-1)=-1;
%         %fj(j,1) = (4/3 * SOL(j-1,n-1))-(1/3*SOL(j-1,n-2));
%         fj(end,1) =0;
%         fj(end-1,1)=0;
%         SOLj = Aj\fj;
%         SOL(1:NX,n) = SOLj(2:NX+1,1);
% 
%     end
%     disp("Burgers_fun_2_end");
% end
% 
