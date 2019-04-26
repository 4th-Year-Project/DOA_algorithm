fc = 8.68e8;                    % Operating frequency
fs = 1000;                      % Sampling frequency
c = physconst('LightSpeed');    % Propagation velocity
wl = c/fc;                      % Wavelength

d = wl/2;                       % Distance between antennas
M = 4;                          % Number of elements
L = 1;                          % Number of sources

ang1 = [-14;0];                  % [azimuth, elevation] angles
angs = [ang1];                  % Concatenate angles


Nsamp = 1024;                  % Number of snapshots
noisePwr = 0.1;                % Noise power

ula = phased.ULA('NumElements',M,'ElementSpacing',d);
pos = getElementPosition(ula)/wl;

signal = sensorsig(pos,Nsamp,angs,noisePwr);        %generate signal

ang_true = az2broadside(angs(1,:));             %broadside angle conversion

musicangle = phased.MUSICEstimator('SensorArray',ula,...
             'OperatingFrequency',fc,'ForwardBackwardAveraging',false,...
             'NumSignalsSource','Property','NumSignals',1,...
             'DOAOutputPort',true);
         
[~,mus_ang] = musicangle(signal);                     %%%%% DOA using MUSIC 


for j = 1:18000
    eta = ((j-1)*pi)/18000;
    Y(j) = -(j-9001)/100;

for x=1:M
    for k=1:L
        u(k, x) = exp((-2*1j*x*pi*d*cos(eta))/wl);
    end
end

a(:,:,j) = circshift(u, L);
a(1,1,j) = 1;
A(:,:,j) = a(:,:,j).';                   %Response through 180 degrees
    
end

Rxx = zeros(4,4);                   
snapshots = 200;

for j=1:snapshots                        % Calculate correlation matrix
    
    temp = (signal(j, 1:4)') * (signal(j, 1:4));
    Rxx = Rxx + temp;
    
end

Rxx = Rxx/snapshots;                     % R signal received, N snapshots, correlation matrix

for j = 1:18000
   a_tht(:,:,j) = A(:,:,j);
   a_thtH(:,:,j) = a_tht(:,:,j)';
   
   P_BS(j) = (a_thtH(:,:,j) * Rxx * a_tht(:,:,j));     % Calculate power spectrums
   P_MVDR(j) = 1/(a_thtH(:,:,j) * inv(Rxx) * a_tht(:,:,j));
   
end

   % Plot power spectrums
%    
% subplot(2,1,1);
% plot(Y, real(P_BS))
% title('Beamscan Spectrum')
% subplot(2,1,2);
% plot(Y, real(P_MVDR))
% title('MVDR Spectrum')
    
% Find Normalised peaks of plots

[BS_pks, BS_locs] = findpeaks(real(P_BS));
BS_DOA = [BS_pks; BS_locs];
[M,I] = max(BS_DOA (1,:));
BS_Max = BS_DOA(1,I);
BS_Norm = real(P_BS)/BS_Max;


[MVDR_pks, MVDR_locs] = findpeaks(real(P_MVDR));
MVDR_DOA = [MVDR_pks; MVDR_locs];
[M,I] = max(MVDR_DOA(1,:));
MVDR_Max = MVDR_DOA(1,I);
MVDR_Norm = real(P_MVDR)/MVDR_Max;

figure 
plot(Y, BS_Norm)
xlim([-90 90])
title('Normalised Spectrum')
xlabel('DOA Angle (Degrees)')
ylabel('Normalised Spectrum')
hold on 
plot(Y, MVDR_Norm)

legend('Beamscan','MVDR')





% fprintf(' True angle = %4.2f\n\n MATLAB MUSIC DOA angle= %4.2f\n\n Beamscan DOA angle = %4.2f\n\n MVDR DOA angle = %4.2f\n\n', ang_true, mus_ang, BS_DOA, MVDR_DOA); 

clear i j k BS_pks BS_locs MVDR_pks MVDR_locs I M eta x u