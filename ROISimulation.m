%% Define Simulation Parameters %%

% Define Simulation Grid Parameters
Nx = 2^10  ; % number of pixels
dx = 8e-6  ; % pixels size [m]
Lx = Nx*dx ; % real side length [m]
[X,Y,Rho,Phi] = sky.MakeGrid(Nx,Nx,dx,dx);

%% Define Beam %%

% General beam parameters
w0 = 1e-3 ; % Embedded Gaussian beam waist [m]
zR = 0.5*k*w0^2 ; % Rayleigh range [m]

% LG beam indices
LR = 2 ; % azimuthal index (OAM) of right-circulalry polarized component
PR = 0 ; % radial index of right-circulalry polarized component
LL = 0 ; % azimuthal index (OAM) of left-circulalry polarized component
PL = 0 ; % radial index of left-circulalry polarized component

% Define Fields
UR = sky.LG(Rho,Phi,PR,LR,1,w0) ; % right-circulalry polarized field
UL = sky.LG(Rho,Phi,PL,LL,1,w0) ; % left-circulalry polarized field
% Normalise so both beams have equal power = 1
UR = UR/sum(sum(abs(UR).^2));
UL = UL/sum(sum(abs(UL).^2));

figure(1)
subplot(2,2,1); imagesc(abs(UR).^2); axis image off; title("I_R")     
subplot(2,2,2); imagesc(abs(UL).^2); axis image off; title("I_L")
subplot(2,2,3); imagesc(angle(UR)) ; axis image off; title("arg(U_R)"); clim([-pi pi])
subplot(2,2,4); imagesc(angle(UL)) ; axis image off; title("arg(U_L)"); clim([-pi pi])
colormap jet

%% Apply ROI Mask and Compute Stoke Parameters %%

for rad_ind = 1:Nx/2

    % Define Mask
    rad = rad_ind*dx;
    mask  = zeros(Nx,Nx);
    mask(Rho<rad) = 1;
    
    % Apply Mask to Component Fields
    URmask = mask.*UR;
    ULmask = mask.*UL;
    
    % Compute Stokes Parameters
    [s0,s1,s2,s3] = sky.GetStokesRL(URmask,ULmask);
    
    % Plot Stokes
    figure(1)
    subplot(2,2,1); imagesc(s0) ; axis image off; clim([ 0 1]); title("s_0")
    subplot(2,2,2); imagesc(s1) ; axis image off; clim([-1 1]); title("s_1")
    subplot(2,2,3); imagesc(s2) ; axis image off; clim([-1 1]); title("s_2")
    subplot(2,2,4); imagesc(s3) ; axis image off; clim([-1 1]); title("s_3")
    colormap jet
        
    [Nsurf(rad_ind), sigz] = sky.SkyNumSurf(s1,s2,s3);
    
    disp(rad_ind + ": " + Nsurf(rad_ind))

end
