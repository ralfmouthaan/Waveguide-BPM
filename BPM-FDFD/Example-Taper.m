% Ralf Mouthaan
% University of Adelaide
% June 2023
%
% Modelling of a fibre taper using the beam propagation method

clc; clear variables; close all;
addpath('Functions\')

%% Geometry axes

lambda = 532e-9;
NX = 1000; 
NZ = 1000;
XMAX = -100e-6;
ZMAX = 500e-6;

x = linspace(-XMAX, XMAX, NX);
z = linspace(0, ZMAX, NZ);

%% Refractive index profile

CoreRadius = 10e-6;
TipZ = 250e-6;
n0 = 1.45;
n = ones(NX, NZ)*n0;

[z_mesh, x_mesh] = meshgrid(z, x.');
n(x_mesh > CoreRadius/TipZ*z_mesh - CoreRadius & x_mesh < -CoreRadius/TipZ*z_mesh + CoreRadius) = 1.48;

%% Starting field

w0 = 10e-6;
F0 = exp(-x.^2/w0^2);

%% Solver

Fmesh = BPM_FFT1D(x, z, n, F0, lambda);

%% Plot

figure;
contourf(z*1e6, x*1e6, abs(Fmesh).^2, 'LineColor', 'none')
ylabel('x (\mum)');
xlabel('z (\mum)')
title('|E|^2')

% When plotting several contour maps, Matlab assumes they have the same
% scale. So, adjust the scale of n to simplify plotting.
ntemp = (n - min(min(n)))/max(max(n)) * max(max(abs(Fmesh).^2));

hold on;
contour(z*1e6, x*1e6, ntemp, 'k'); 
xlabel('z (\mum)')
ylabel('x (\mum)')
colorbar;
colormap parula;
set(gca, 'FontSize', 12);