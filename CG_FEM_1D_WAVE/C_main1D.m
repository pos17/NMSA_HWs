function [errors,solutions,femregion,Dati,u_snp] = C_main1D(TestName,nRef)
%==========================================================================
% Solution of the Wave Equation with linear finite elements 
% coupled with the leap-frog scheme 
% 
%  u_tt - c^2 u_xx = f  in (a,b)x(0,T)
%  u_t(x,0) = v_0       in (a,b) 
%  u(x,0)   = u_0       in (a,b)
%  + boundary conditions:
%  Dirichlet:  u(s,t) = g  s=a,b
%  Neumann  : c^2du/dn(s,t) = g s=a,b
%  Periodic : u(a,t) = a(b,t)
%  Absorbing: du/dt(s,t) + cdu/dn(s,t) = 0  s=a,b 
%==========================================================================
%
%    INPUT:
%          TestName    : (string)  see C_dati.m
%          nRef        : (int)     refinement level
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) infos about finite elements
%                        discretization
%          Dati        : (struct)  see C_dati.m
%
% Usage:
%    [errors,solutions,femregion,Dati] = C_main1D('Test1',3)



addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing


%==========================================================================
% LOAD DATA FOR TEST CASE
%==========================================================================

Dati = C_dati(TestName);
Dati.nRefinement = nRef;

%==========================================================================
% MESH GENERATION
%==========================================================================

[Region] = C_create_mesh(Dati);

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[femregion] = C_create_femregion(Dati,Region);

%==========================================================================
% BUILD FINITE ELEMENT MATRICES
%==========================================================================

[M_nbc,A_nbc] = C_matrix1D(Dati,femregion);
if(strcmp(Dati.bc,'R'))
    A_nbc(1,1)     = A_nbc(1,1) + Dati.c2;
    A_nbc(end,end) = A_nbc(end,end) + Dati.c2;
end

%==========================================================================
% BUILD FINITE ELEMENTS RHS a time 0
%==========================================================================
Dati.t = 0;
[b_nbc] = C_rhs1D(Dati,femregion);
t = 0;

F_1 = b_nbc

if(strcmp(Dati.bc,'N'))
    b_nbc(1)   = b_nbc(1)   + eval(Dati.neumann1);
    b_nbc(end) = b_nbc(end) + eval(Dati.neumann2);
elseif(strcmp(Dati.bc,'R'))
    b_nbc(1)   = b_nbc(1)   - Dati.c2*eval(Dati.neumann1);
    b_nbc(end) = b_nbc(end) + Dati.c2*eval(Dati.neumann2);

% elseif(strcmp(Dati.bc,'M'))
%     b_nbc(1)   = b_nbc(1)   - M_nbc(1,1)*eval(Dati.Lifting1tt) - A_nbc(1,1)*eval(Dati.Lifting1);
end

%==========================================================================
% BUILD INITIAL CONDITIONS
%==========================================================================
x = femregion.coord;
u0 = eval(Dati.u0)
v0 = eval(Dati.v0);


u_snp = u0;
if(strcmp(Dati.timeIntScheme,'NLF'))
    %% First step of not leapfrog
    % 1) Compute the rhs
    Dati.t = Dati.dt;
    t = Dati.dt;
    [b_nbc_1] = C_rhs1D(Dati, femregion);
    
    b_nbc = 0.5*Dati.dt^2 * (b_nbc_1) + M_nbc*u0 + Dati.dt*M_nbc*v0;

    M_nbc_nlf = M_nbc + Dati.dt^2*A_nbc/2;

    if(strcmp(Dati.bc,'D'))
        [M,b,u_g] = C_bound_cond1D(M_nbc_nlf,b_nbc,femregion,Dati);
        u1 =  M\b;
        u1 = u1 + u_g;
    end
     % 5) Plot the obtained solution --> C_snapshot_1D
    
    [u1] = C_snapshot_1D(femregion, u1, Dati);
    %% End of first step not leapfrog
    
    fprintf('============================================================\n')
    fprintf('Starting time-loop ... \n');
    fprintf('============================================================\n')
    
    u_snp(:,2) = u1;
    k = 3;
    
    for t = Dati.dt : Dati.dt : Dati.T - Dati.dt
        
        fprintf('time = %5.3e \n',t);
        
        %==========================================================================
        % BUILD FINITE ELEMENTS RHS a time t
        %==========================================================================
        Dati.t = t;
        [b_nbc] = C_rhs1D(Dati,femregion);
        
        if(strcmp(Dati.bc,'N'))
            b_nbc(1)   = b_nbc(1)   + eval(Dati.neumann1);
            b_nbc(end) = b_nbc(end) + eval(Dati.neumann2);
        elseif(strcmp(Dati.bc,'R'))
            b_nbc(1)   = b_nbc(1)   - Dati.c2*eval(Dati.neumann1);
            b_nbc(end) = b_nbc(end) + Dati.c2*eval(Dati.neumann2);
        end
        
        
        % Repeat steps 1) to 5) for the general time step
        b_nbc = Dati.dt^2 * (b_nbc) + 2 * M_nbc * u1 - M_nbc * u0;
        M_nbc_nlf = M_nbc+Dati.dt^2*A_nbc;

        if(strcmp(Dati.bc,'D'))
            [M,b,u_g] = C_bound_cond1D(M_nbc_nlf,b_nbc,femregion,Dati);
            u2 =  M\b;
            u2 = u2 + u_g;
        elseif(strcmp(Dati.bc,'N')|| strcmp(Dati.bc,'R'))
            u2 =  M_nbc\b_nbc;
        elseif(strcmp(Dati.bc,'P'))
            b_nbc(1,:) = b_nbc(1,:) + b_nbc(end,:);
            b_nbc(end) = 0;
            %solution
            u2 =  M_nbc\b_nbc;
        elseif(strcmp(Dati.bc,'A'))
            b_nbc(1)   = b_nbc(1)   - sqrt(Dati.c2)*(u1(1)-u0(1))*Dati.dt;
            b_nbc(end) = b_nbc(end) - sqrt(Dati.c2)*(u1(end)-u0(end))*Dati.dt;
            u2 =  M_nbc\b_nbc;
            
        end
        
        % plot
        [u2] = C_snapshot_1D(femregion, u2, Dati);
        
        u_snp(:,k) = u2;
        k = k +1;
        
        % Put a pause between one step and the other to see the plot
        pause(0.015);
        
        % Update the solution
        u0 = u1;
        u1 = u2;
        
        
    end



elseif(strcmp(Dati.timeIntScheme,'LF'))
    %% First step of leapfrog ...
    % 1) Compute the rhs
    
    b_nbc = 0.5*Dati.dt^2 * (b_nbc - A_nbc*u0) + M_nbc*u0 + Dati.dt*M_nbc*v0;
    
    
    if(strcmp(Dati.bc,'D'))
        [M,b,u_g] = C_bound_cond1D(M_nbc,b_nbc,femregion,Dati);
        u1 =  M\b;
        u1 = u1 + u_g;
    elseif(strcmp(Dati.bc,'N') || strcmp(Dati.bc,'R'))
        u1 =  M_nbc\b_nbc;
    elseif(strcmp(Dati.bc,'P'))
        % 1st line
        M_nbc(1,:) = M_nbc(1,:) + M_nbc(end,:);
        b_nbc(1,:) = b_nbc(1,:) + b_nbc(end,:);
        % last line
        M_nbc(end,:) = 0;
        M_nbc(end,1) = 1;  M_nbc(end,end) = -1;
        b_nbc(end) = 0;
        %solution
        u1 =  M_nbc\b_nbc;
    elseif(strcmp(Dati.bc,'A'))
        b_nbc(1)   = b_nbc(1)   - sqrt(Dati.c2)*v0(1)*Dati.dt;
        b_nbc(end) = b_nbc(end) - sqrt(Dati.c2)*v0(end)*Dati.dt;
        u1 =  M_nbc\b_nbc;
    end
    
    % 5) Plot the obtained solution --> C_snapshot_1D
    
    [u1] = C_snapshot_1D(femregion, u1, Dati);
    
    %% End of first step leapfrog
    
    fprintf('============================================================\n')
    fprintf('Starting time-loop ... \n');
    fprintf('============================================================\n')
    
    u_snp(:,2) = u1;
    k = 3;
    
    for t = Dati.dt : Dati.dt : Dati.T - Dati.dt
        
        fprintf('time = %5.3e \n',t);
        
        %==========================================================================
        % BUILD FINITE ELEMENTS RHS a time t
        %==========================================================================
        Dati.t = t;
        [b_nbc] = C_rhs1D(Dati,femregion);
        
        if(strcmp(Dati.bc,'N'))
            b_nbc(1)   = b_nbc(1)   + eval(Dati.neumann1);
            b_nbc(end) = b_nbc(end) + eval(Dati.neumann2);
        elseif(strcmp(Dati.bc,'R'))
            b_nbc(1)   = b_nbc(1)   - Dati.c2*eval(Dati.neumann1);
            b_nbc(end) = b_nbc(end) + Dati.c2*eval(Dati.neumann2);
        end
        
        
        % Repeat steps 1) to 5) for the general time step
        b_nbc = Dati.dt^2 * (b_nbc - A_nbc*u1) + 2 * M_nbc * u1 - M_nbc * u0;
        
        if(strcmp(Dati.bc,'D'))
            [M,b,u_g] = C_bound_cond1D(M_nbc,b_nbc,femregion,Dati);
            u2 =  M\b;
            u2 = u2 + u_g;
        elseif(strcmp(Dati.bc,'N')|| strcmp(Dati.bc,'R'))
            u2 =  M_nbc\b_nbc;
        elseif(strcmp(Dati.bc,'P'))
            b_nbc(1,:) = b_nbc(1,:) + b_nbc(end,:);
            b_nbc(end) = 0;
            %solution
            u2 =  M_nbc\b_nbc;
        elseif(strcmp(Dati.bc,'A'))
            b_nbc(1)   = b_nbc(1)   - sqrt(Dati.c2)*(u1(1)-u0(1))*Dati.dt;
            b_nbc(end) = b_nbc(end) - sqrt(Dati.c2)*(u1(end)-u0(end))*Dati.dt;
            u2 =  M_nbc\b_nbc;
            
        end
        
        [u2] = C_snapshot_1D(femregion, u2, Dati);
        
        u_snp(:,k) = u2;
        k = k +1;
        
        % Put a pause between one step and the other to see the plot
        pause(0.015);
        
        % Update the solution
        u0 = u1;
        u1 = u2;
        
        
    end
end
%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

[solutions] = C_postprocessing(Dati,femregion,u2);

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
errors = [];
if (Dati.plot_errors)
    [errors] = C_compute_errors(Dati,femregion,solutions);
end



