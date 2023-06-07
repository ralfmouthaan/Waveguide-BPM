% Ralf Mouthaan
% University of Adelaide
% June 2023
%
% Beam propagation method - FFT implementation, as per Okamoto,
% "Fundamentals of Optical Waveguides", Ch. 7.2. One bit is taken from
% Edgar Guevara's "FFT Beam Propagation Method" on Matlab's File Exchange: 
% https://au.mathworks.com/matlabcentral/fileexchange/14795-fft-beam-propagation-method
% 
% This is only valid for geometries with small refractive index changes.

function Fmesh = BPM_FFT1D(x, z, n, F0, lambda)

    DEBUG = false;

    %% Error checks

    % Only absolute basic error checks for now. These checks avoid the
    % program crashing, but do not ensure the results are correct.

    if lambda < 100e-9 || lambda > 10e-6
        error('Check size of lambda')
    end
    if length(F0) ~= length(x)
        error('First dimension of F must match length of x');
    end

    %% Derived parameters
    
    NX = length(x);
    NZ = length(z);
    k0 = 2*pi/lambda;
    dz = z(2) - z(1);
    drho = 1/(max(x) - min(x));
    rhox = (-(NX - 1)/2:(NX-1)/2)*drho;
    kx = 2*pi*rhox;
    n0 = (min(min(n)) + max(max(n)))/2;
    
    %% Absorbing boundaries
    
    % Do I want to define this in this function, or have the ABCs passed in?
    
    alpha = zeros(1, NX);
    alpha(x > max(x)*0.9) = 1e11*abs(x(x > max(x)*0.9) - max(x)*0.9).^2;
    alpha(x < min(x)*0.9) = 1e11*abs(x(x < min(x)*0.9) + max(x)*0.9).^2;
    
    if DEBUG
        figure; 
        plot(x*1e6, alpha, 'LineWidth', 2);
        title('Absorbing boundaries');
        ylabel('Loss');
        xlabel('\mum')
    end
    
    %% Calculation
    
    % This equation is from Okamoto (Eq. 7.26). I cannot get it to work.
    % dBeta = sqrt(k0^2*n0^2 - kx.^2) - k0*n0;
    
    % This equation is from Guevara, and seems to work:
    % https://au.mathworks.com/matlabcentral/fileexchange/14795-fft-beam-propagation-method
    Phase1 = exp(-1i*dz*kx.^2./(n0*k0 + sqrt(max(0,n0^2*k0*2 - kx.^2))));
    
    Fmesh = zeros(NX, NZ);
    F = F0;
    
    for i = 1:length(z)
    
        % Free space step
        F = fftshift(fft(fftshift(F)));
        F = F.*Phase1;
        F = fftshift(ifft(fftshift(F)));
    
        % Phase propagation
        Phase2 = exp(-alpha-1i*k0*(n(:, i).' - n0)*dz);
        F = F.*Phase2;
    
        if DEBUG
            plot(x*1e6, abs(F).^2, 'LineWidth', 2); 
            title(['z = ' num2str(z(i)*1e3) 'mm'])
            xlabel('\mum')
            ylabel('|E|^2')
            ylim([0 1])
            drawnow;
        end
    
        Fmesh(:, i) = F;
    
    end

end























