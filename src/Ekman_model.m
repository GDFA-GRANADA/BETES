%*************************************************************************%
%                                                                         %
%        Baroclinic Ekman Transport Equation Solver (BETES)               %
%                                                                         %
%        Authors:                                                         %
%        Llorente-Lázaro, Víctor Javier     ETSIAE - UPM                  %
%        Díez-Minguito, Manuel              IISTA  - UGR                  %
%                                                                         %
%        Date:                                                            %
%        September 22, 2021                                               %
%                                                                         %
%        Version 1.0:                                                     %
%        - Barotropic+atmospheric/Baroclinic pressure gradient            %
%        - Variable turbulent eddy viscosity                              %
%        - Second-order (low) and fourth-order compact (high) FD scheme   % 
%                                                                         %
%*************************************************************************%

clear variables
close all
%clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Input data                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experimental data
data = 'experimental';                                                     % experimental, no experimental, classic Ekman
if strcmp(data,'experimental')
  load datos_GoC_medios.mat
end

% Oceanographic condition
OceanCond = 'upwelling';                                                 % upwelling, downwelling

% Location of the estuary 
lat = 37 + 1/60 + 25.92/3600;                                              % Latitude 

% Water column
if strcmp(data,'experimental')
  H = abs(zcells(1));                                                      % Depth [m]
elseif strcmp(data,'classic Ekman')                                     
  H = 500;                                                                 % Depth [m]
else
  H = 20;                                                                  % Depth [m] e.g.: 20.0
end

% Longitudinal length for Rossby number                                    % (1km = 1000m)
L = 15*1000;                                                               % Length [m] 

% Wind condition
if strcmp(data,'experimental') && strcmp(OceanCond,'upwelling')
  wx = mean_wx_li_upwell_wogap;                                            % x wind component [m/s] 
  wy = mean_wy_li_upwell_wogap;                                            % y wind component [m/s]    
elseif strcmp(data,'experimental') && strcmp(OceanCond,'downwelling')
  wx = mean_wx_li_downwell_wogap;                                          % x wind component [m/s] 
  wy = mean_wy_li_downwell_wogap;                                          % y wind component [m/s]    
else
  wx = 5.0;                                                                % x wind component [m/s] e.g.: 5.0
  wy = -1.0;                                                               % y wind component [m/s] e.g.: -1.0
end

% Slip condition
if strcmp(data,'experimental') && strcmp(OceanCond,'upwelling')
  ustar = mean_residual_vparalel_li_upwell_wogap(1);                       % x bottom shear velocity component [m/s]
  vstar = mean_residual_vperpen_li_upwell_wogap(1);                        % y bottom shear velocity component [m/s] 
elseif strcmp(data,'experimental') && strcmp(OceanCond,'downwelling')
  ustar = mean_residual_vparalel_li_downwell_wogap(1);                     % x bottom shear velocity component [m/s]
  vstar = mean_residual_vperpen_li_downwell_wogap(1);                      % y bottom shear velocity component [m/s]    
elseif strcmp(data,'classic Ekman') 
  ustar = 0;                                                               % x bottom shear velocity component [m/s]
  vstar = 0;                                                               % y bottom shear velocity component [m/s]
else
  ustar = 0.001;                                                           % x bottom shear velocity component [m/s] e.g.: 0.001
  vstar = -0.0001;                                                         % y bottom shear velocity component [m/s] e.g.: -0.0001
end

% Vertical eddy viscosity
if strcmp(data,'classic Ekman')
  visc = 'eddy 0';
else
  visc = 'eddy 2';
end
switch visc
  case 'eddy 0'                                                            % Classical Ekman case, + H---->\infty
    Kz0_temp = 0.1;                                                        %   Eddy viscosity at the surface [m2/s]
    n = 0;                                                                 %   Decay parameter (Kz = Kz0 = const.) 
    zm = -0.0*H;                                                           %   Depth when Kz() has the maximum (Kz = Kz0 = const.) 
    zh_temp = -0.2*H;                                                      %   Depth where the two layers are separated. Valid except when n = 0, see profile_Kz()          
  case 'eddy 1'                                                            % Variable eddy viscosity 1 (CASE BASE)
    Kz0_temp = 0.1;                                                        %   Eddy viscosity at the surface [m2/s]
    n = 2;                                                                 %   Decay parameter 
    zm = -0.1*H;                                                           %   Depth when Kz() has the maximum
    zh_temp = -0.2*H;                                                      %   Depth where the two layers are separated. Valid except when n = 0, see profile_Kz()    
  case 'eddy 2'                                                            % Variable eddy viscosity 2
    Kz0_temp = 0.3;                                                        %   Eddy viscosity at the surface [m2/s]
    n = 3.2;                                                               %   Decay parameter 
    zm = -0.1*H;                                                          %   Depth when Kz() has the maximum
    zh_temp = -0.25*H;                                                      %   Depth where the two layers are separated. Valid except when n = 0, see profile_Kz()        
  otherwise                                                                
    warning('No viscosity')
end
linevisc = 'yes';                                                          % yes, no

% Pressures                                                                % (1 mbar = 100 Pa, 1mbar = -1cm de H20, 1km = 1000m) 
if strcmp(data,'classic Ekman')
  dpx = (0-0)*100;                                                         % Pressure difference in x direction [Pa] 
  dx = 300*1000;                                                           % Increment in x direction for pressure [m] e.g.: 300*1000
  dpy = (0-0)*100;                                                         % Pressure difference in y direction [Pa]
  dy = 300*1000;                                                           % Increment in y direction for pressure [m] e.g.: 300*1000
  drhox = (0-0);                                                           % Water density difference in x direction [kg/m^3]
  dxrho = 50*1000;                                                         % Increment in x direction for water density [m] e.g.: 50*1000
  drhoy = (0-0);                                                           % Water density difference in y direction [kg/m^3]
  dyrho = 50*1000;                                                         % Increment in y direction for water density [m] e.g.: 50*1000  
else
  dpx = (20.5-0)*100;                                                         % Pressure difference in x direction [Pa] 
  dx = 300*1000;                                                             % Increment in x direction for pressure [m] e.g.: 300*1000
  dpy = (0-23.5)*100;                                                           % Pressure difference in y direction [Pa]
  dy = 300*1000;                                                             % Increment in y direction for pressure [m] e.g.: 300*1000
  drhox = (0-5);                                                             % Water density difference in x direction [kg/m^3]
  dxrho = 5*1000;                                                           % Increment in x direction for water density [m] e.g.: 50*1000
  drhoy = (0-5);                                                             % Water density difference in y direction [kg/m^3]
  dyrho = 5*1000;                                                           % Increment in y direction for water density [m] e.g.: 50*1000
end
baroclinic = 'exact';                                                      % Integral of the baroclinic pressure: exact or numeric 

% Mesh                                                                 
N = 100;                                                                   % Number of nodes

% Discretization
order = 'high';                                                             % Order of the scheme: low, high
scheme = 'half-way';                                                       % Type of the low-order scheme: grid-point, half-way

% Constants 
g = 9.8;                                                                   % Gravity [m/s^2]
rho = 1025;                                                                % Water density [kg/m^3]
rhoa = 1.22;                                                               % Air density [kg/m^3]
T = 23*3600 + 56*60 + 4.1;                                                 % Earth period [s]
Omega = 2*pi/T;                                                            % Angular velocity [rad/s]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                             Discrete domain                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                           
z = linspace(0,-H,N)';                                                     % Discrete vertical coordinate [m], z = 0 at the surface
Dz = (H - 0)/(N - 1);                                                      % Spacing distance [m]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Wind stress at the surface                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CD = (0.8 + 0.065*(wx^2 + wy^2)^0.5)*1e-3;                                 % Air-on-water drag coefficient 
tauwx = rhoa*CD*wx*(wx^2 + wy^2)^0.5;                                      % x wind stress at the surface [kg/ms^2]
tauwy = rhoa*CD*wy*(wx^2 + wy^2)^0.5;                                      % y wind stress at the surface [kg/ms^2]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Pressure gradients                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Barotropic + atmospheric pressure gradient (includes d\eta and dpatm)
Px = 1/rho*dpx/dx;                                                         % Pressure gradient in x direction []
Py = 1/rho*dpy/dy;                                                         % Pressure gradient in y direction []

% Baroclinic pressure gradient, drho/dx(z)
eps = H/20;                                                                % e-folding depth of drho/dx [m]
if strcmp(baroclinic,'numeric')
  drdx = drhox/dxrho*exp(z/eps);                                           % Maximum baroclinic press. gradient is assumed to be at the surface
  Qx = zeros(N,1);                                                         % Definition of vector Qx
  Qx(1) = 0;
  for i=2:N
    Qx(i) = -g/rho*trapz(z(1:i),drdx(1:i));                                % Integral of the Baroclinic pressure x-gradient 
  end
else
  Qx = (g/rho)*(drhox/dxrho)*eps*(1 - exp(z/eps));                         % Integral of the Baroclinic pressure x-gradient 
end

% Baroclinic pressure gradient, drho/dy(z) 
if strcmp(baroclinic,'numeric')
  drdy = drhoy/dyrho*exp(z/eps);                                           % Maximum baroclinic press. gradient is assumed to be at the surface
  Qy = zeros(N,1);                                                         % Definition of vector Qy
  Qy(1) = 0;
  for i=2:N
    Qy(i) = -g/rho*trapz(z(1:i),drdy(1:i));                                % Integral of the Baroclinic pressure y-gradient
  end
else
  Qy = (g/rho)*(drhoy/dyrho)*eps*(1 - exp(z/eps));                         % Integral of the Baroclinic pressure y-gradient
end

% Total pressure gradient
Tx_full = -Px - Qx;                                                        % In x direction
Ty_full = -Py - Qy;                                                        % In y direction
Tx = Tx_full(2:N-1);                                                       % (Without 0 and N-1 elements)
Ty = Ty_full(2:N-1);                                                       % (Without 0 and N-1 elements)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Coriolis                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = 2*Omega*sin(lat*pi/180);                                               % Coriolis force [rad/s]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Eddy viscosity profile                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Viscosity at grid points
[Kz_full, Kz0, zh] = profile_Kz(z, Kz0_temp, n, zm, zh_temp);
Kz = Kz_full(2:N-1);                                                       % (Without 0 and N-1 elements)

% Viscosity at half way
zz = [-Dz/2:-Dz:-H+Dz/2]';
[Kzz_full, Kz0, zh] = profile_Kz(zz, Kz0_temp, n, zm, zh_temp);
Kzz = Kzz_full;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    Code                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic                                                                        % Start clock

switch order
    
  case 'low'                                                               % Low-order scheme
      
    % Sparse vectors 
    bxs   = sparse(N-2,1);
    bys   = sparse(N-2,1);
    bxBCs = sparse(N-2,1);
    byBCs = sparse(N-2,1);
    % Sparse matrices
    Ms = sparse(N-2,N-2);
    Is = speye(N-2);
    if strcmp(scheme,'grid-point')                                         % Eddy viscosity at grid points
%     1) Discrete source terms 
%     x direction 
      bxs(1    ) = -Tx(1    )*Dz^2;
      bxs(2:N-3) = -Tx(2:N-3)*Dz^2;
      bxs(  N-2) = -Tx(  N-2)*Dz^2;
%     y direction 
      bys(1    ) = -Ty(1    )*Dz^2;
      bys(2:N-3) = -Ty(2:N-3)*Dz^2;
      bys(  N-2) = -Ty(  N-2)*Dz^2;
%     2) Discrete BC source terms 
%     x direction 
      bxBCs(1    ) = -Dz/rho*Kz(1)/Kz_full(1)*tauwx;
      bxBCs(2:N-3) =  0                            ;
      bxBCs(  N-2) = -Kz_full(N)*ustar             ;
%     y direction 
      byBCs(1    ) = -Dz/rho*Kz(1)/Kz_full(1)*tauwy;
      byBCs(2:N-3) =  0                            ;
      byBCs(  N-2) = -Kz_full(N)*vstar             ;
%     3) Discrete matrix
%     Boundary points (z = 0)
      Ms(1,1) = -Kz(2);
      Ms(1,2) =  Kz(2);
%     Inner points
      for i=2:N-3
        Ms(i,i-1) =   Kz(i)           ;
        Ms(i,i  ) = -(Kz(i) + Kz(i+1));
        Ms(i,i+1) =           Kz(i+1) ;
      end
%     Boundary points (z = -H)
      Ms(N-2,N-3) =   Kz(N-2)              ;
      Ms(N-2,N-2) = -(Kz(N-2) + Kz_full(N)); 
    else                                                                   % Eddy viscosity at half-way
%     1) Discrete source terms 
%     x direction 
      bxs(1    ) = -Tx(1    )*Dz^2;
      bxs(2:N-3) = -Tx(2:N-3)*Dz^2;
      bxs(  N-2) = -Tx(  N-2)*Dz^2;
%     y direction 
      bys(1    ) = -Ty(1    )*Dz^2;
      bys(2:N-3) = -Ty(2:N-3)*Dz^2;
      bys(  N-2) = -Ty(  N-2)*Dz^2;
%     2) Discrete BC source terms 
%     x direction 
      bxBCs(1    ) = -Dz/rho*tauwx*(Kzz(1)/Kz(1));
      bxBCs(2:N-3) =  0                          ;
      bxBCs(  N-2) = -Kzz(N-1)*ustar             ;
%     y direction 
      byBCs(1    ) = -Dz/rho*tauwy*(Kzz(1)/Kz(1));
      byBCs(2:N-3) =  0                          ;
      byBCs(  N-2) = -Kzz(N-1)*vstar             ;
%     3) Discrete matrix
%     Boundary points (z = 0)
      Ms(1,1) = -Kzz(2);
      Ms(1,2) =  Kzz(2);
%     Inner points
      for i=2:N-3
        Ms(i,i-1) =   Kzz(i)            ;
        Ms(i,i  ) = -(Kzz(i) + Kzz(i+1));
        Ms(i,i+1) =            Kzz(i+1) ;
      end
%     Boundary points (z = -H)
      Ms(N-2,N-3) =              Kzz(N-2) ;
      Ms(N-2,N-2) = -(Kzz(N-1) + Kzz(N-2));      
    end
    % Coreolis matrix
    Fs = (f*Dz^2)*Is;
    % Assembling
    bs = [bxs;bys] + [bxBCs;byBCs];
    As = [Ms,Fs;-Fs,Ms];
    % Matrix solver
    Usappr = As\bs;
    Uappr = full(Usappr);                                                  % Uappr = (u1,...,uN-2,v1,...,vN-2)^T
    U = [Uappr(1  )+Dz*tauwx/(rho*Kz_full(1));Uappr(1  :    N-2  );ustar]; % Numerical u(z)
    V = [Uappr(N-1)+Dz*tauwy/(rho*Kz_full(1));Uappr(N-1:(2*(N-2)));vstar]; % Numerical v(z)
    
  case 'high'                                                              % High-order scheme
      
    % Sparse matrices
    MCDs = sparse(N,N);
    QCDs = sparse(N,N);
    Zs = sparse(N-1,N-1);
    Ks = sparse(N-1,N-1);
    % Sparse vectors
    zs = sparse(N-1,1);
%   dphi/dz-discrete matrix  (d/dz = -d/d|z|)
    MCDs(1,1) = -1;
    MCDs(1,2) = -3;
    for i=2:N-1
      MCDs(i,i-1) = -1/4;
      MCDs(i,i  ) = -1;
      MCDs(i,i+1) = -1/4;
    end
    MCDs(N,N-1) = -3;
    MCDs(N,N  ) = -1;
%   phi-discrete matrix
    QCDs(1,1) = -17/6;
    QCDs(1,2) =  3/2;
    QCDs(1,3) =  3/2;
    QCDs(1,4) = -1/6;
    for i=2:N-1
      QCDs(i,i-1) = -3/4;
      QCDs(i,i  ) =  0;
      QCDs(i,i+1) =  3/4;
    end
    QCDs(N,N-3) =  1/6;
    QCDs(N,N-2) = -3/2;
    QCDs(N,N-1) = -3/2;
    QCDs(N,N  ) =  17/6; 
%   Mq and Mu matrices
    Mqs = MCDs(1:N-1,1:N-1);
    Mus = MCDs(2:N,2:N);
%   Wq and Qu metrices
    Wqs = QCDs(1:N-1,2:N);
    Wus = QCDs(2:N,1:N-1);
%   mq and mu vectors
    mqs = -MCDs(1:N-1,N);
    mus = -MCDs(2:N,1);
%   wq and wu vectors
    wqs = QCDs(1:N-1,1);
    wus = QCDs(2:N,N);
%   Diffusion matrix
    for i=1:N-1
      Ks(i,i) = Kz_full(i+1);     
    end
    Gs = -Dz*Mus*inv(Ks);
%   Coreolis matrix
    Fs = Dz*f*Mqs;
%   pressure gradient vector
    txs = Tx_full(1:N-1);
    tys = Ty_full(1:N-1);
%   Assembling
    bs = [zs;-Dz*Mqs*txs;zs;-Dz*Mqs*tys] + ... 
         [-ustar*wus-Dz*(tauwx/(rho*Kz0))*mus;...
          -(tauwx/rho)*wqs+Dz*(Tx_full(N)+f*vstar)*mqs;...
          -vstar*wus-Dz*(tauwy/(rho*Kz0))*mus;...
          -(tauwy/rho)*wqs+Dz*(Ty_full(N)-f*ustar)*mqs];
    As = [Wus,Gs,Zs,Zs;Zs,Wqs,Fs,Zs;Zs,Zs,Wus,Gs;-Fs,Zs,Zs,Wqs];
%   Matrix solver
    Xsappr = As\bs;
    Xappr = full(Xsappr); 
    U = [Xappr(         1 :    N-1  );ustar];                              % Numerical u(z)
    V = [Xappr((2*(N-1)+1):(3*(N-1)));vstar];                              % Numerical v(z)
    
  otherwise                       
      
    warning('No scheme')
    
end

toc                                                                        % Stop clock

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Ekman and Rossby number                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximum speed
S = (U.^2 + V.^2).^0.5;
Smax = max(S);

% Maximum eddy viscosity
Kmax = max(Kz_full);

% Rossby number
Ro = Smax/(L*f);

% Ekman number
Ek = Kmax/(H^2*f);

% Warnings for Rossby
if (Ro >= 1) && (Ro < 10)
  warning('Effects of planetary rotation can be neglected')
  fprintf('Rossby number is %d\n',Ro);
elseif (Ro >= 10) && (Ro < 100)
  warning('Inertial and centrifugal forces dominate')  
  fprintf('Rossby number is %d\n',Ro);
elseif (Ro >= 100)
  warning('Geostrophic approximation fails!!!')
  fprintf('Rossby number is %d\n',Ro);
end

% Warnings for Ekman
if (Ek <= 1) && (Ek > 0.1)
  warning('Diffusion force is not balanced by Coriolis force') 
  fprintf('Ekman number is %d\n',Ek);
elseif (Ek <= 0.1) && (Ek > 0.01)
  warning('Friction effects are not capable of decaying disturbances') 
  fprintf('Ekman number is %d\n',Ek); 
elseif (Ek <= 0.01)
  warning('Ekman model fails!!!')
  fprintf('Rossby number is %d\n',Ek);    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Analytical Solution for Kz = Kz0, dp/dx = dp/dy = 0, z--> -infty       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dE = sqrt(2*Kz0/f);                                                        % Analytical estimate for Ekman layer depth for Kz0 constant
ze = linspace(0,-H,10*N)';
U0 = sqrt(2)/(rho*f*dE);                                                   % Reference speed
Ue = U0.*exp(ze/dE).*(tauwx.*cos(ze/dE - pi/4) - tauwy.*sin(ze/dE - pi/4));% Classic Ekman u(z)
Ve = U0.*exp(ze/dE).*(tauwx.*sin(ze/dE - pi/4) + tauwy.*cos(ze/dE - pi/4));% Classic Ekman v(z)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         L2-norm of the errors for the Classic Ekman problem             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(data,'classic Ekman')
% Exact solution    
  Uexact = U0.*exp(z/dE).*(tauwx.*cos(z/dE - pi/4) - tauwy.*sin(z/dE - pi/4));
  Vexact = U0.*exp(z/dE).*(tauwx.*sin(z/dE - pi/4) + tauwy.*cos(z/dE - pi/4));
% L2-norm of the error  
  error_U = (U - Uexact).^2;
  error_V = (V - Vexact).^2;
  L2norm_U = sqrt(sum(error_U)/N);
  L2norm_V = sqrt(sum(error_V)/N);
% Print values
  fprintf('The spacing distance, Dz, is %f\n',Dz);
  fprintf('The l2-norm for u is %1.5E\n',L2norm_U);
  fprintf('The l2-norm for v is %1.5E\n',L2norm_V);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Pressure gradients from the experimental datas and fitting       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(data,'experimental')
% Experimental velocity field    
  if strcmp(OceanCond,'upwelling')
    Uexpe = mean_residual_vparalel_li_upwell_wogap;
    Vexpe = mean_residual_vperpen_li_upwell_wogap;
  else
    Uexpe = mean_residual_vparalel_li_downwell_wogap;
    Vexpe = mean_residual_vperpen_li_downwell_wogap;  
  end
  UVexpe = [Uexpe';Vexpe'];
  Dzexpe = zcells(2)-zcells(1);
% Experimental pressure gradients  
%  switch order   
%    case 'low'                                                             % Low-order scheme  
      [Kz_expe, Kz0, zh] = profile_Kz(zcells, Kz0_temp, n, zm, zh_temp);
      pgx(1) = Kz_expe(1)*(2*Uexpe(1)-5*Uexpe(2)+4*Uexpe(3)-Uexpe(4))/Dz^2 + ...
               ((-3*Kz_expe(1)+4*Kz_expe(2)-Kz_expe(3))/(2*Dz))*((-3*Uexpe(1)+4*Uexpe(2)-Uexpe(3))/(2*Dz)) + ...
               f*Vexpe(1);
      pgy(1) = Kz_expe(1)*(2*Vexpe(1)-5*Vexpe(2)+4*Vexpe(3)-Vexpe(4))/Dz^2 + ...
               ((-3*Kz_expe(1)+4*Kz_expe(2)-Kz_expe(3))/(2*Dz))*((-3*Vexpe(1)+4*Vexpe(2)-Vexpe(3))/(2*Dz)) - ...
               f*Uexpe(1);     
      for i=2:ncellmax-1
        pgx(i) = Kz_expe(i)*(Uexpe(i-1)-2*Uexpe(i)+Uexpe(i+1))/Dz^2 + ...
                 ((Kz_expe(i-1)-Kz_expe(i+1))/(2*Dz))*((Uexpe(i-1)-Uexpe(i+1))/(2*Dz)) + ...
                 f*Vexpe(i);     
        pgy(i) = Kz_expe(i)*(Vexpe(i-1)-2*Vexpe(i)+Vexpe(i+1))/Dz^2 + ...
                 ((Kz_expe(i-1)-Kz_expe(i+1))/(2*Dz))*((Vexpe(i-1)-Vexpe(i+1))/(2*Dz)) - ...
                 f*Uexpe(i);      
      end
      pgx(ncellmax) = Kz_expe(ncellmax)*(-2*Uexpe(ncellmax)+5*Uexpe(ncellmax-1)-4*Uexpe(ncellmax-2)+Uexpe(ncellmax-3))/Dz^2 + ...
                      ((3*Kz_expe(ncellmax)-4*Kz_expe(ncellmax-1)+Kz_expe(ncellmax-2))/(2*Dz))*((3*Uexpe(ncellmax)-4*Uexpe(ncellmax-1)+Uexpe(ncellmax-2))/(2*Dz)) + ...
                      f*Vexpe(ncellmax);
      pgy(ncellmax) = Kz_expe(ncellmax)*(-2*Vexpe(ncellmax)+5*Vexpe(ncellmax-1)-4*Vexpe(ncellmax-2)+Vexpe(ncellmax-3))/Dz^2 + ...
                      ((3*Kz_expe(ncellmax)-4*Kz_expe(ncellmax-1)+Kz_expe(ncellmax-2))/(2*Dz))*((3*Vexpe(ncellmax)-4*Vexpe(ncellmax-1)+Vexpe(ncellmax-2))/(2*Dz)) - ...
                      f*Uexpe(ncellmax);
      Tx_expe = -pgx';   
      Ty_expe = -pgy';
%    case 'high'                                                            % High-order scheme
%     Sparse matrices
%      MCDs = sparse(ncellmax,ncellmax);
%      QCDs = sparse(ncellmax,ncellmax);
%      Zs = sparse(ncellmax,ncellmax);
%      Ks = sparse(ncellmax,ncellmax);
%     dphi/dz-discrete matrix  (d/dz = -d/d|z|)
%      MCDs(1,1) = -1;
%      MCDs(1,2) = -3;
%      for i=2:ncellmax-1
%        MCDs(i,i-1) = -1/4;
%        MCDs(i,i  ) = -1;
%        MCDs(i,i+1) = -1/4;
%      end
%      MCDs(ncellmax,ncellmax-1) = -3;
%      MCDs(ncellmax,ncellmax  ) = -1;
%     phi-discrete matrix
%      QCDs(1,1) = -17/6;
%      QCDs(1,2) =  3/2;
%      QCDs(1,3) =  3/2;
%      QCDs(1,4) = -1/6;
%      for i=2:ncellmax-1
%        QCDs(i,i-1) = -3/4;
%        QCDs(i,i  ) =  0;
%        QCDs(i,i+1) =  3/4;
%      end
%      QCDs(ncellmax,ncellmax-3) =  1/6;
%      QCDs(ncellmax,ncellmax-2) = -3/2;
%      QCDs(ncellmax,ncellmax-1) = -3/2;
%      QCDs(ncellmax,ncellmax  ) =  17/6;
%     Diffusion matrix
%      [Kz_expe, Kz0, zh] = profile_Kz(zcells, Kz0_temp, n, zm, zh_temp);
%      for i=1:ncellmax-1
%        Ks(i,i) = Kz_expe(i+1);     
%      end      
%     
%      MMCDs = [MCDs,Zs;Zs,MCDs];
%      QQCDs = [QCDs,Zs;Zs,QCDs];
%      KKs = [Ks,Zs;Zs,Ks];
%     First derivative
%      bCDs = (1/Dzexpe)*QQCDs*UVexpe;
%      dUVexpes = MMCDs\bCDs;
%     Second derivative      
%      qqexpes = KKs*dUVexpes;
%      bCDs = (1/Dzexpe)*QQCDs*qqexpes;
%      dqexpes = MMCDs\bCDs;
%     Coreolis vector
%      cfs = f*[Vexpe';-Uexpe'];
%     Experimental \/P
%      pgs = dqexpes + cfs;
%      pg = full(pgs);
%      Tx_expe = -pg(1:ncellmax);
%      Ty_expe = -pg(ncellmax+1:ncellmax+ncellmax);
%    otherwise
%      warning('No scheme')      
%  end

model = @(bx1,bx2,bx3,x) - bx1 - bx2.*bx3.*(1 - exp(x./bx3));
Tx_fit = fit(zcells',Tx_expe,model,...
             'StartPoint', [1e-7, 1e-7, 0.5],...
             'Lower', [-1e-7, -1e-7, 1]);
Ty_fit = fit(zcells',Ty_expe,model,...
             'StartPoint', [1e-7, 1e-7, 0.5],...
             'Lower', [-1e-7, -1e-7, 1]);         
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Figures and output files                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% K(z) profile
figure
plot(Kz_full, z); hold on
xlabel('Kz (m^2/s)')
ylabel('z (m)')
grid on
xlim([0 max(Kz)])
ylim([min(z) max(z)])
if strcmp(linevisc,'yes')
  x=linspace(0,max(Kz),21);
  y=zh_temp*ones(size(x));
  plot(x,y,'r-.')
  y=zm*ones(size(x));
  plot(x,y,'k-.')
end
set(gcf,'color','w')
adjustpdfpage(gcf,20)
print(gcf,'K(z)','-dpdf')

% U and V profiles vs Classic Ekman velocities
figure
plot(U,z,'ro',V,z,'bo'); hold on
plot(Ue,ze,'r-',Ve,ze,'b-'); 
if strcmp(data,'experimental') && strcmp(OceanCond,'upwelling')
  hold on  
  plot(mean_residual_vparalel_li_upwell_wogap, zcells, 'r-.'); hold on
  plot(mean_residual_vperpen_li_upwell_wogap, zcells, 'b-.');    
elseif strcmp(data,'experimental') && strcmp(OceanCond,'downwelling')
  hold on  
  plot(mean_residual_vparalel_li_downwell_wogap, zcells, 'r-.'); hold on
  plot(mean_residual_vperpen_li_downwell_wogap, zcells, 'b-.');       
end
grid on
xlabel('Current (m/s)')
ylabel('z (m)')
if strcmp(data,'experimental')
  legend('Numerical u(z)', 'Numerical v(z)', 'Classic Ekman u(z)', 'Classic Ekman v(z)', 'Experimental u(z)', 'Experimental v(z)')
else
  legend('Numerical u(z)', 'Numerical v(z)', 'Classic Ekman u(z)', 'Classic Ekman v(z)')  
end
set(gcf,'color','w')
adjustpdfpage(gcf,20)
print(gcf,'Profiles','-dpdf')

% Barotropic+atmosferic/Baroclinic pressure gradient profile
figure
plot(Tx_full,z,'ro',Ty_full,z,'bo'); hold on 
plot(Tx_expe,zcells,'r-*',Ty_expe,zcells,'b-*');
xlabel('Pressure gradient (Pa/m)')
ylabel('z (m)')
grid on
%xlim([-1e-4 1e-4])
set(gcf,'color','w')
adjustpdfpage(gcf,20)
print(gcf,'Pressure','-dpdf')

% Velocity spiral
figure
plot(U,V,'o'); hold on
plot([0 U(1)],[0 V(1)],'-')
grid on
xlabel('U (m/s)')
ylabel('V (m/s)')
legend('(U,V)', 'Surface')
set(gcf,'color','w')
adjustpdfpage(gcf,20)
print(gcf,'UvsV','-dpdf')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%