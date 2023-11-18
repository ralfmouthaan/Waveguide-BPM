% Finite Difference Beam Propagation Method     3 Mayo 2007
% Edgar Guevara Codina 
% Dispositivos Optoelectronicos 
% Genera un pulso gaussiano y lo propaga 1000 um a lo largo del eje y
% Grafica superficie
close all;                          % Cierra Ventanas
clear all;                          % Limpia Variables
clc;                                % Limpia Pantalla
% ---------- Declaracion de Variables -----------------
tic                                 % Comienza Timer
x1 = -60e-6;                        % Coordenada Inicial
x2 = 60e-6;                         % Coordenada Final
num_samples = 101;                  % Numero de muestras (par)
dx = (x2-x1)/num_samples;           % Espaciado de las muestras en x
dy = 1e-6;                          % Incremento en y
x = linspace (x1, x2-dx, num_samples);     % Dominio espacial
W0 = 8e-6;                          % Radio de la cintura del pulso
n_index = ones(1,num_samples);      % Indice de refraccion del medio
nmax = 1;                           % Indice de refraccion maximo
nmin = 1;                           % Indice de refraccion minimo
n_bar = (nmax+nmin)/2;              % Indice Promedio
lambda = 0.8e-6;                    % Longitud de onda
k0 = 2*pi/lambda;                   % Numero de onda

% ---------- Generacion del pulso -----------------
modo = exp (-(x/W0).^2);          % Pulso Gaussiano

k = k0.*n_index;                     % Vector de onda
k_bar = k.*n_bar;                    % Numero de onda de referencia

% ------ Generamos la reticula para graficar ------ 
[xx,yy] = meshgrid ([x1:dx:x2-dx],[dy:dy:1000e-6]);
zz = zeros(size(xx));
h = dy;
ro = dy/(dx^2);
A = j./(2.*k_bar);
B = j.*(k.^2-k_bar.^2)./(2.*k_bar);
a = -ro.*A;    
b = 2*(1+ro.*A)-h.*B;
c = a;
d = zeros(1,num_samples);
matrix = zeros(num_samples);

% --------- Generacion de la matriz tridiagonal ---------
for m = 1:1:num_samples,
    if ((m>1) && (m<num_samples))
        matrix(m,m-1) = a(m);
        matrix(m,m) = b(m);
        matrix(m,m+1) = c(m);
    else
        matrix(1,1) = b(1);
        matrix(1,2) = c(1);
        matrix(num_samples,num_samples-1) = a(num_samples);
        matrix(num_samples,num_samples) = b(num_samples);
    end
end
    
% --------- Ciclo Principal de Propagacion ---------  
for n = 1:1:1000,
    for m = 1:1:num_samples,
        if ((m>1) && (m<num_samples))
            d(m) = (2*(1-ro*A(m))+h*B(m))*modo(m)+ro*A(m)*(modo(m-1)+modo(m+1));
        else
            d(1) = (2*(1-ro*A(1))+h*B(1))*modo(1)+ro*A(1)*(modo(2));
            d(num_samples) = (2*(1-ro*A(num_samples))+h*B(num_samples))*modo(num_samples)+ro*A(num_samples)*(modo(num_samples-1));
        end
    end
    modo = matrix\d.';
    zz(n,:) = [abs(modo.')];        % Generamos la matriz para graficar
end

% -- Graficamos la magnitud de nuestros pulsos propagados --
surf(xx,yy,zz);
shading interp;
colormap jet;
grid on;
title('Magnitud de los Pulsos Gaussianos Propagados');
toc                             % Termina Timer