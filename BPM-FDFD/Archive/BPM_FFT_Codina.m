% Ralf Mouthaan
% University of Adelaide
% June 2023
%
% Beam propagation method - FFT implementation, as per Okamoto,
% "Fundamentals of Optical Waveguides", Ch. 7.2
% This is valid for geometries with small refractive index changes.
% Here I am just looking to duplicate Cordina's BPM_2Step to iron out some
% kinks

clc; clear variables; close all;

%% User-entered parameters

lambda = 1e-6;
k0 = 2*pi/lambda;
CoreRadius = 50e-6;
NX = 1000;
NZ = 1000;
XMAX = 100e-6;
ZMAX = 1e-3;

%% Derived parameters
x = linspace(-XMAX, XMAX, NX);
dx = x(2) - x(1);
z = linspace(0, ZMAX, NZ);
dz = z(2) - z(1);
drho = 1/(max(x) - min(x));
rhox = (-(NX - 1)/2:(NX-1)/2)*drho;
kx = 2*pi*rhox;

figure; plot(kx);

%% Refractive index profile

nmin = 1.499;
nmax = 1.5;
x1 = -12e-6;
x2 = -2e-6;
x3 = 2e-6; 
x4 = 12e-6;
n0 = (nmin + nmax)/2;
for j = 1:1:NX
    if (x(j) >= x1) & (x(j) <= x2)
        n(j) = nmax;
    elseif (x(j) >= x3) & (x(j) <= x4)    
        n(j) = nmax;
    else
        n(j) = nmin;
    end
end

figure; 
plot(x*1e6, n, 'LineWidth', 2); 
title('Transverse Refractive Index Profile')
xlabel('\mum')
ylabel('Refractive Index')

%% Starting field

w0 = 5e-6;
F = exp(-x.^2/w0^2);

figure; 
plot(x*1e6, F, 'LineWidth', 2); 
title('Field: z = 0')
xlabel('\mum')
ylabel('|E|')

%% Calculation

% This is what Okamoto says - it does not work:
%dBeta = sqrt(k0^2*n0^2 - kx.^2) - k0*n0;

%This is what Codina says, it does work:
dBeta = kx.^2./(n0*k0 + sqrt(max(0,n0^2*k0*2 - kx.^2)));

figure;
for i = 1:length(z)

    Ftemp = F;

    % Free space step
    Ftemp = fftshift(fft(fftshift(Ftemp)));
    Ftemp = Ftemp.*exp(-1i*dBeta*dz/2);
    Ftemp = fftshift(ifft(fftshift(Ftemp)));

    % Phase propagation
    Ftemp = Ftemp.*exp(-1i*k0*(n - n0)*dz);

    % Free space step
    Ftemp = fftshift(fft(fftshift(Ftemp)));
    Ftemp = Ftemp.*exp(-1i*dBeta*dz/2);
    Ftemp = fftshift(ifft(fftshift(Ftemp)));

    F = Ftemp;

    plot(x*1e6, abs(F).^2, 'LineWidth', 2); 
    title(['z = ' num2str(z(i)*1e3) 'mm'])
    xlabel('\mum')
    ylabel('|E|^2')
    ylim([0 1])
    drawnow;

end



























