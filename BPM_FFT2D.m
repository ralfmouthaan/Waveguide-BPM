% Ralf Mouthaan
% University of Adelaide
% June 2023
%
% Beam propagation method - FFT implementation, as per Okamoto,
% "Fundamentals of Optical Waveguides", Ch. 7.2
% This is valid for geometries with small refractive index changes.
% I have implemented it with two-dimensional transverse coordinates.

clc; clear variables; close all;

%% User-entered parameters

lambda = 532e-9;
k0 = 2*pi/lambda;
CoreRadius = 50e-6;
NX = 1000;
NZ = 1000;
XMAX = 100e-6;
ZMAX = 1e-3;

%% Derived parameters
x = linspace(-XMAX, XMAX, NX);
dx = x(2) - x(1);
y = x.';
z = linspace(0, ZMAX, NZ);
dz = z(2) - z(1);
radius = sqrt(x.^2 + y.^2);
drho = 1/(max(x) - min(x));
rhox = (-(NX - 1)/2:(NX-1)/2)*drho;
rhoy = rhox.';
kx = 2*pi*rhox;
ky = 2*pi*rhoy;

%% Refractive index profile

n0 = 1.45;
n = ones(NX)*n0;
n(radius < CoreRadius) = 1.48;

figure; 
imagesc(x*1e6, y*1e6, n); 
axis square;
title('Transverse Refractive Index Profile')
xlabel('\mum')
ylabel('\mum')
colorbar;

%% Starting field

w0 = 5e-6;
F = exp(-radius.^2/w0^2);

figure; 
imagesc(x*1e6, y*1e6, F); 
axis square;
title('Field: z = 0')
xlabel('\mum')
ylabel('\mum')
colorbar;


%% Calculation

% This is an adaptation of the equation from Codina.
% I cannot ge the Okamoto equation to work.
dBeta = (kx.^2 + ky.^2)./(n0*k0 + sqrt(max(0,n0^2*k0*2 - kx.^2 - ky.^2)));

figure;
for i = 1:length(z)

    Ftemp = F;

    % Free space step
    Ftemp = fftshift(fft2(fftshift(Ftemp)));
    Ftemp = Ftemp.*exp(-1i*dBeta*dz/2);
    Ftemp = fftshift(ifft2(fftshift(Ftemp)));

    % Phase propagation
    Ftemp = Ftemp.*exp(-1i*k0*(n - n0)*dz);

    % Free space step
    Ftemp = fftshift(fft2(fftshift(Ftemp)));
    Ftemp = Ftemp.*exp(-1i*dBeta*dz/2);
    Ftemp = fftshift(ifft2(fftshift(Ftemp)));

    F = Ftemp;

    imagesc(x*1e6, y*1e6, abs(F)); 
    axis square;
    title(['Field: z = ' num2str(z(i)*1e3) 'mm'])
    xlabel('\mum')
    ylabel('\mum')
    clim([0 1])
    colorbar;
    drawnow;

end



























