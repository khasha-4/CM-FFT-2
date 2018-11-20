clear all 
%% Create the solution domain

iv = [-2,2];                        % Computational domain
fl = 8;
N  = 2^fl;                          % Number of nodes
L = (iv(2)-iv(1));                  % Domain length L
dx = L/N;                           % Meshsize
x = iv(1):dx:iv(2)-dx;              % Discrete domain


%% Generate a Hertzian pressure distribution 

inerval = find(-1 < x & x < 1);         % Pressure region
p = zeros(1,N);                         % Create a zero vector
p(inerval) = sqrt(1- x(inerval).^2);    % Scaled contact pressure eq 4.35 in compendium

%% Define the matrix for the discrete integral Kernel.
j  = 1:N;                                                           
d1 = (1-j+0.5)*dx; d2=(1-j-0.5)*dx;
k  = d1.*( log( abs( d1 ) )  - 1) - d2.*( log(abs( d2 ) ) - 1);
K  = toeplitz(-k)/(pi);                                           % Creating kernal matrix

%% Standard Method matrix, vector multiplication
% Multiplication between the Kernel and the dimensionless pressure

tic; Standard_Method = (K * p.').'; toc;
figure (1)
% Plot 
plot(x,Standard_Method-Standard_Method(1),'k-','linewidth',3);

%% Discret Convolution-Fast Fourier Method
omega  = 2*pi / L;                                   % Frequency 
k      = [k(1:N),0,fliplr(k(2:N))];                  % Creating periodicity
DCH    = real(fft(-k/pi));                           % The transfer function
pDC    = [zeros(1,N/2),p,zeros(1,N/2)];              % Extending p 
uDCFFT = real(ifft(DCH.*fft(pDC)));                  % Calc conv & IFT
uDCFFT  = uDCFFT(N/2+1:N/2+N);
hold on 
plot(x,uDCFFT-uDCFFT(1),'r--','linewidth',3);










