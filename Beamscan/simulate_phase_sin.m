samples = 200; % no. of samples
F_s = 500e3; % 500kHz sampling frequency
time_vect = 0 : 1/(500e3) : 200 * 1/(500e3);
phase_shift_deg = -120;
phase_shift = deg2rad(phase_shift_deg);
qps = deg2rad(90);
fc = 20e3;

theoretical_ang = rad2deg(asin(phase_shift/pi));

tpf = 2 * pi * fc;


Sig1 = sin(tpf * time_vect);
Sig2 = sin((tpf * time_vect) + phase_shift);
Sig3 = sin((tpf * time_vect) + 2*phase_shift);
Sig4 = sin((tpf * time_vect) + 3*phase_shift);

%% Noise option 
SNR = 20;

SigI(:,1) = awgn(Sig1, SNR);
SigI(:,2) = awgn(Sig2, SNR);
SigI(:,3) = awgn(Sig3, SNR);
SigI(:,4) = awgn(Sig4, SNR);

SigQ(:,1)  = awgn(sin((tpf * time_vect) + qps + phase_shift), SNR);
SigQ(:,2)  = awgn(sin((tpf * time_vect) + qps + 2*phase_shift), SNR);
SigQ(:,3)  = awgn(sin((tpf * time_vect) + qps + 3*phase_shift), SNR);
SigQ(:,4)  = awgn(sin((tpf * time_vect) + qps + 4*phase_shift), SNR);


%%


% plot(time_vect, SigQ)
% plot(time_vect, SigI)


signal = (SigI(1:200,:) + (1j*SigQ(1:200,:)));

%%%%%%%% DOA calculation with real signals
[ANGLE_BS, ANGLE_MVDR, ANGLE_QR] = DOA_calc(signal, 200, true);

method = {'Theoretical Angle', 'Beamscan', 'MVDR', 'QR '};
angle = [theoretical_ang, ANGLE_BS, ANGLE_MVDR, ANGLE_QR];

Real_data = table;
Real_data.Method = method';
Real_data.Angle = angle';
display(Real_data)


function [BS_DOA, MVDR_DOA, QR_DOA] =  DOA_calc(signal, snapshots, plot_sdoa)

fc = 8.68e8;                    % Operating frequency
fs = 100e3;                      % Sampling frequency
c = physconst('LightSpeed');    % Propagation velocity
wl = c/fc;                      % Wavelength

d = wl/2;                       % Distance between antennas
M = 4;                          % Number of elements
L = 1;                          % Number of sources

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

for j=1:snapshots                        % Calculate correlation matrix
    
    temp = (signal(j, 1:4)') * (signal(j, 1:4));
    Rxx = Rxx + temp;
    
end

Rxx = Rxx/snapshots;                     % R signal received, N snapshots, correlation matrix

[Q, R] = qr(Rxx);          % QR factorisation

Qs = Q(:,(1:L));
Qn = Q(:,(L+1):M);

QnH = Qn';

for j = 1:18000
   a_tht(:,:,j) = A(:,:,j);
   a_thtH(:,:,j) = a_tht(:,:,j)';
   
   P_BS(j) = (a_thtH(:,:,j) * Rxx * a_tht(:,:,j));     % Calculate power spectrums
   P_MVDR(j) = 1/(a_thtH(:,:,j) * inv(Rxx) * a_tht(:,:,j));
   
   P_QR(j) = (1 / (a_thtH(:,:,j) * Qn * QnH * a_tht(:,:,j)) );
end


if plot_sdoa == true
   % Plot power spectrums

figure
subplot(3,1,1);
plot(Y, real(P_BS))
title('Beamscan Spectrum')
subplot(3,1,2);
plot(Y, real(P_MVDR))
title('MVDR Spectrum')

subplot(3,1,3);
plot(Y, (real(P_QR(1,:))))
title('QR Spectrum')

end
    
%%%%%%%% Find peaks of plots

[BS_pks, BS_locs] = findpeaks(real(P_BS));
BS_DOA = [BS_pks; BS_locs];
[M,I] = max(BS_DOA (1,:));
BS_DOA = Y(BS_DOA(2,I));

[MVDR_pks, MVDR_locs] = findpeaks(real(P_MVDR));
MVDR_DOA = [MVDR_pks; MVDR_locs];
[M,I] = max(MVDR_DOA(1,:));
MVDR_DOA = Y(MVDR_DOA(2,I));

[QR_pks, QR_locs] = findpeaks(real(P_QR));
QR_DOA = [QR_pks; QR_locs];
[M,I] = max(QR_DOA(1,:));
QR_DOA = Y(QR_DOA(2,I));

method = {'Beamscan', 'MVDR', 'QR '};
angle = [BS_DOA, MVDR_DOA, QR_DOA];

% Real_data = table;
% Real_data.Method = method';
% Real_data.Angle = angle';
% display(Real_data)

end