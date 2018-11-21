%% Assignment 2

clear all

%% Initialation. Allocating memory for variables
iv = [-2,2];                        % Computational domain
fl = 14;                             
% Here we create emty vector which gives us more flexibility to save data
% for every iteration in the loop

S_bytes =      [];
S_time =       [];
FFT_bytes =    [];
FFT_time =     [];
FFT_time_avg = [];
N =            [];                                                     


for ite = 1:fl

    N(ite)  = 2^ite;                                                        % Number of nodes 2 4 8 16 32 64 ...                          
    L = (iv(2)-iv(1));                                                      % Domain length L
    dx = L/N(ite);                                                          % Meshsize                          
    x = iv(1):dx:iv(2)-dx;                                                  % Discrete domain              
    inerval = find(-1 < x & x < 1);                                         % Pressure region        
    p = zeros(1,N(ite));                                                    % Create a zero vector                        
    p(inerval) = sqrt(1- x(inerval).^2);                                    % Scaled contact pressure eq 4.35 in compendium
    
%% Standard Method 
    j  = 1:N(ite);                                                          % Define the matrix for the discrete integral Kernel.
    d1 = (1-j+0.5)*dx; d2=(1-j-0.5)*dx;
    k  = d1.*( log( abs( d1 ) )  - 1) - d2.*( log(abs( d2 ) ) - 1);
    K  = toeplitz(-k)/(pi);                                                 % Creating kernal matrix                                          
    
    tic;
    Standard_Method = (K * p.').';                                          % Standard Method matrix, vector multiplication
    S_time(ite) = toc;
    S_SM    =  whos('K');                                                   % See how much usage of memory
    S_bytes(ite) = S_SM.bytes;                                              % Save the usage value
%% FFT - Method 
    k = [k(1:N(ite)),0,fliplr(k(2:N(ite)))];                                % Creating periodicity
        for i = 1:20
            tic;
            DCH    = real(fft(-k/pi)); % The transfer function              % Computing 20 time the same operation for making the time average                  
            FFT_time(i) = toc;
        end
    
    FFT_time_avg(ite) = mean(FFT_time);                                     % Making time average for coputation time
    pDC    = [zeros(1,N(ite)/2),p,zeros(1,N(ite)/2)];                       % Extending p                   
    uDCFFT = real(ifft(DCH.*fft(pDC)));                                     % Calc conv & IFT                                 
    uDCFFT  = uDCFFT(N(ite)/2+1:N(ite)/2+N(ite));
    S_FFT = whos('DCH');                                                    % See how much usage of memory
    FFT_bytes(ite) = S_FFT.bytes;                                           % Save the usage value
    
end

Err_CCFFT = (max(abs((Standard_Method-Standard_Method(1))-(uDCFFT-uDCFFT(1))))/...
    max(abs((Standard_Method-Standard_Method(1)))));                        % Calculating error
%% Messages to display

Message_1 = sprintf('The error is: %d',Err_CCFFT);
disp(Message_1)
Message_2 = sprintf('Maximum usage of the memory by FFT_Method is: %d',(FFT_bytes(end)));
disp(Message_2)
Message_3 = sprintf('Maximum usage of the memory by Standard_Method is: %d',(S_bytes(end)));
disp(Message_3)

%% Visualization

movegui(figure('rend','painters','pos',[10 10 1000 400]),'northwest');      % Move figure to specified location on screen

subplot(1,2,1);                                                             % Plot several plots in the figure
plot((N).*log(N),FFT_time_avg)
title('FFT Method')
xlabel('N*log (N)');ylabel('Time');
grid on;

subplot(1,2,2);                                                             
plot(N.^2,S_time)               
title('Standard Method')
xlabel('N^2');ylabel('Time');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movegui(figure('rend','painters','pos',[10 10 800 400]),'southeast');       % Move figure to specified location on screen

subplot(1,2,1);                                                             % Plot several plots in the figure
semilogy(N,FFT_bytes,N,S_bytes)
legend ('FFT bytes','S bytes')
xlabel('N (Number of nodes)');ylabel('Memory B');
grid on;

subplot(1,2,2);                                                             % Plot several plots in the figure
semilogy(N,FFT_time_avg,N,S_time)
legend ('FFT time','Standard time')
xlabel('N (Number of nodes)');ylabel('Time');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movegui(figure('rend','painters','pos',[10 10 800 400]),'northeast');       % Move figure to specified location on screen

plot(x,Standard_Method-Standard_Method(1),'k-','linewidth',3);              % Plot the deflection used in both method
hold on 
plot(x,uDCFFT-uDCFFT(1),'r--','linewidth',3);
grid on
xlabel('Dimensionless Length ');ylabel('Dimensionless Deflection');





