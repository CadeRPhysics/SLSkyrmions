%% Define Simulation Parameters %%
close all
clc

% Define Simulation Grid Parameters
Nx = 2^10  ; % number of pixels
dx = 8e-6  ; % pixels size [m]
Lx = Nx*dx ; % real side length [m]
[X,Y,Rho,Phi] = sky.MakeGrid(Nx,Nx,dx,dx);

% Define Light Properties
wvl = 633e-9 ; % wavelength [m]
k = 2*pi/wvl ; % wavenumber [1/m]

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

%% Get Stokes Paramaters %%

[s0,s1,s2,s3] = sky.GetStokesRL(UR,UL);

figure(1)
subplot(2,2,1); imagesc(s0) ; axis image off; clim([ 0 1]); title("s_0")
subplot(2,2,2); imagesc(s1) ; axis image off; clim([-1 1]); title("s_1")
subplot(2,2,3); imagesc(s2) ; axis image off; clim([-1 1]); title("s_2")
subplot(2,2,4); imagesc(s3) ; axis image off; clim([-1 1]); title("s_3")
colormap jet

%% Calculate Skyrme Number %%

% Line Integral Calculation
pixcrop = Nx/2-1        ; % normally half the simulation space, but can be adjusted
BeamSizeInPix = 2*w0/dx ; % approximate size of the beam in pixels to ensure all of it is included in calculation
thresh  = 0             ; % threshold for singularity exlcusion
[Nline, S1, S2, S3,singularities] = sky.SkyNumLine(s1,s2,s3,s0,pixcrop,BeamSizeInPix, thresh) ;
PolPhase = angle(S1+1i.*S2) ; % Polarisation phase
PlotSings = imgaussfilt(singularities,20) ; % make singularities easier to see
PlotSings = PlotSings./max(PlotSings(:));

% Surface Integral Calculation
[Nsurf, sigz] = sky.SkyNumSurf(s1,s2,s3);

figure(1)
subplot(2,3,1); imagesc(S1) ; axis image off; clim([-1 1]); title("S_1")
subplot(2,3,2); imagesc(S2) ; axis image off; clim([-1 1]); title("S_2")
subplot(2,3,3); imagesc(S3) ; axis image off; clim([-1 1]); title("S_3")
subplot(2,3,4); imagesc(PolPhase)  ; axis image off; clim([-pi pi])
subplot(2,3,5); imagesc(PlotSings) ; axis image off; clim([-1 1]); title("N_{line} = " + Nline)
subplot(2,3,6); imagesc(sigz)      ; axis image off; title("N_{surf} = " + Nsurf)
colormap jet
