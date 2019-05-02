%%%%%%%% signal parameters
fs = 500E3;           % calibration sample frequency
T = 1/fs;             % sample time
f0 = 9.25E3;          % calibration signal frequency
fc = 868E6;           % carrier frequency

cal_dat = readtable('cal_data2k.csv');
cal_dat = table2array(cal_dat)';

signal = readtable('0deg_test1.csv');
signal = table2array(signal)';

%%

%%%%%%%% Convert signal range 0 - 2V as measured by ADC
signal = (2*signal)/255;         

%%%%%%%% plot signals? true/false   
plot_f = false;                  

%%%%%%%% Plot normalised calibration signals   
freq = 10E3;
N_samples = 2000;
[cal_norm_i, cal_norm_q] = norm_cal(cal_dat, freq, N_samples, false);                   

%%%%%%%% Calculate phase of calibration signals   
fprintf('\n   I signal phase difference: \n')
for k = 1:3
    
x = cal_norm_i(:,1)';
y = cal_norm_i(:,(k+1))';
    
PhDiff = phdiffmeasure(x, y);
PhDiff = PhDiff*180/pi;

PhDiffstr = num2str(PhDiff);
disp(['Phase difference 1 - ' num2str(k+1) ' = ' PhDiffstr ' deg'])
Phase_from_A(k+1) = PhDiff;
end
Phase_from_A(1) = 0;
fprintf('\n')

%%%%%%%% Phase as a fraction of wavelength
Phase = (360 - Phase_from_A)/360;

%%%%%%%% Discrete phase, as sample
Sample_delay = round((Phase/f0)/T, 0);

%%%%%%%% Align calibration signals
[al_calsigi, al_calsigq] = sig_align(cal_norm_i, cal_norm_i, Sample_delay);
% plot(al_calsigi)
% title('Aligned Calibration Signals')

%%%%%%%% Normalise real signals
[ni_sig, nq_sig] = norm_cal(signal, freq, 200, false);            

%%%%%%%% Align real signals
[al_sigi, al_signq] = sig_align(ni_sig, nq_sig, Sample_delay);

%%%%%%%% DOA calculation with simulated signals
% simsig = DOA(-3, 0.01, false);
%%%%%%%% 

%%

sigsize = size(al_sigi);
sigsize(2) = [];
tt = (0:(sigsize-1))/fc;
tt = tt.';

for k = 1:4  
    MCRF(:,k) = (al_sigi(:,k).*sin(2*pi*fc*tt)) - (al_signq(:,k).*cos(2*pi*fc*tt));
    
    doa_data(:,k) = al_sigi(:,k) + (1j*al_signq(:,k));
    doa_data2(:,k) = al_calsigi(:,k) + (1j*al_calsigq(:,k));
end

%%%%%%%% DOA calculation with real signals
nsamples = size(doa_data);
[ANGLE_BS(1), ANGLE_MVDR(1), ANGLE_QR(1)] = DOA_calc(doa_data, nsamples(1), true);

% nsamples = size(MCRF);
% [ANGLE_BS(1), ANGLE_MVDR(1), ANGLE_QR(1)] = DOA_calc(MCRF, nsamples(1), true);

%%

%%%%%%%% FFT of specified signal
% FFTS = al_sigi(:,1);
% FFTS = cal_norm_i(:,1);
% FFTS = real(simsig(:,1));
% FFTS = MCRF(:,1); 

% [fft_ps fft_f fft_peak] = sig_fft(FFTS(:,1), fs);
% disp(fft_peak)


clear x y plot_f PhDiff PhDiffstr k a b sigsize nsamples N_samples

%%

function [NI, NQ] = norm_cal(Cal_dat, fs, num_samples, plot_f)

t = (0:(num_samples - 1))/fs;

%%%%% Split signals into I and Q
for k = 1:4         
    ai = round((2*k)-1);
    aq = round(2*k);
    I(:,k) = Cal_dat(:,ai);
    Q(:,k) = Cal_dat(:,round(2*k));
end

%%%%% Remove DC component
for k = 1:4
    zero_I(:,k) = (I(:,k) - mean(I(:,k)));
    zero_Q(:,k) = (Q(:,k) - mean(Q(:,k)));
end

%%%%% Find peaks
for k = 1:4
    max_I(k) = max(findpeaks((sqrt((zero_I(:,k)).^2)),fs));
    max_Q(k) = max(findpeaks((sqrt((zero_Q(:,k)).^2)),fs));
    
    norm_I(:,k) = zero_I(:,k)/max_I(k);
    norm_Q(:,k) = zero_Q(:,k)/max_Q(k);
end

limits = ([-1.1 1.1]);

if plot_f == true
    subplot(2,1,1)
    plot(t, norm_I)
    ylim(limits)
    xlabel('Time (s)')
    ylabel('Normalised Amplitude')
    title('Normalised I Calibration Signals')

    subplot(2,1,2)
    plot(t, norm_Q)
    ylim(limits)
    xlabel('Time (s)')
    ylabel('Normalised Amplitude')
    title('Normalised Q Calibration Signals')
end

NI = norm_I;
NQ = norm_Q;
end

function PhDiff = phdiffmeasure(x, y)
% function: PhDiff = phdiffmeasure(x, y)
% x - first signal in the time domain
% y - second signal in the time domain
% PhDiff - phase difference Y -> X, rad
% represent x as column-vector if it is not
if size(x, 2) > 1
    x = x';
end
% represent y as column-vector if it is not
if size(y, 2) > 1
    y = y';
end
% remove the DC component
x = x - mean(x);
y = y - mean(y);
% signals length
N = length(x);
% window preparation
win = rectwin(N);
% fft of the first signal
X = fft(x.*win);
% fft of the second signal
Y = fft(y.*win);
% phase difference calculation
[~, indx] = max(abs(X));
[~, indy] = max(abs(Y));
PhDiff = angle(Y(indy)) - angle(X(indx));
end

function [SIG_I, SIG_Q] = sig_align(norm_I_sig, norm_Q_sig, sample_delay)

SIG_I = zeros(200,4);
SIG_Q = zeros(200,4);


for k = 1:4
temp = norm_I_sig((sample_delay(k)+1):end, k);
[eend, sstart] = size(temp);

SIG_I(sstart:eend,k) = temp;

temp = norm_Q_sig((sample_delay(k)+1):end, k);
[eend, sstart] = size(temp);

SIG_Q(sstart:eend,k) = temp;
end

for k = 1:4
    a = find(SIG_I(:,k)==0, 1, 'first');
    b = find(SIG_Q(:,k)==0, 1, 'first');
    
    if a ~= 0
        z_pos_I(k) = a;
        z_pos_Q(k) = b;
    end
    
end

SIG_I(min(z_pos_I):end, :) = [];
SIG_Q(min(z_pos_Q):end, :) = [];

end

function signalout = DOA(given_angle, n_power, plot_f)

fc = 8.68e8;                    % Operating frequency
fs = 100e3;                      % Sampling frequency
c = physconst('LightSpeed');    % Propagation velocity
wl = c/fc;                      % Wavelength

d = wl/2;                       % Distance between antennas
M = 4;                          % Number of elements
L = 1;                          % Number of sources

ang1 = [given_angle;0];                  % [azimuth, elevation] angles
angs = [ang1];                  % Concatenate angles


Nsamp = 1024;                  % Number of snapshots
noisePwr = n_power;                % Noise power

ula = phased.ULA('NumElements',M,'ElementSpacing',d);
pos = getElementPosition(ula)/wl;

signal = sensorsig(pos,Nsamp,angs,noisePwr);        %generate signal
% signal = sensorsig(pos,Nsamp,angs); 

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
Simulated_data(:,:,j) = a(:,:,j).';                   %Response through 180 degrees
    
end

Rxx = zeros(4,4);                   
snapshots = 200;

for j=1:snapshots                        % Calculate correlation matrix
    
    temp = (signal(j, 1:4)') * (signal(j, 1:4));
    Rxx = Rxx + temp;
    
end

Rxx = Rxx/snapshots;                     % R signal received, N snapshots, correlation matrix

[Q, R] = qr(Rxx);          % QR factorisation

Qs = Q(:,(1:L));
Qn = Q(:,(L+1):4);

QnH = Qn';

for j = 1:18000
   a_tht(:,:,j) = Simulated_data(:,:,j);
   a_thtH(:,:,j) = a_tht(:,:,j)';
   
   P_BS(j) = (a_thtH(:,:,j) * Rxx * a_tht(:,:,j));     % Calculate power spectrums
   P_MVDR(j) = 1/(a_thtH(:,:,j) * inv(Rxx) * a_tht(:,:,j));
   
   P_QR(j) = (1 / (a_thtH(:,:,j) * Qn * QnH * a_tht(:,:,j)) );
end


if plot_f == true
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

    
   % Find peaks of plots

[BS_pks, BS_locs] = findpeaks(real(P_BS));
BS_DOA = [BS_pks; BS_locs];
[M,I] = max(BS_DOA (1,:));
BS_DOA = Y(BS_DOA(2,I));

[MVDR_pks, MVDR_locs] = findpeaks(real(P_MVDR));
MVDR_DOA = [MVDR_pks; MVDR_locs];
[M,I] = max(MVDR_DOA(1,:));
MVDR_DOA = Y(MVDR_DOA(2,I));

[QR_pks, QR_locs] = findpeaks(real(P_QR));
QR_DOA = [MVDR_pks; MVDR_locs];
[M,I] = max(QR_DOA(1,:));
QR_DOA = Y(QR_DOA(2,I));

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
QR_DOA = [MVDR_pks; MVDR_locs];
[M,I] = max(QR_DOA(1,:));
QR_DOA = Y(QR_DOA(2,I));


method = {'True angle', 'MATLAB MUSIC', 'Beamscan', 'MVDR', 'QR '};
angle = [ang_true, mus_ang, BS_DOA, MVDR_DOA, QR_DOA];

Simulated_data = table;
Simulated_data.Method = method';
Simulated_data.Angle = angle';
display(Simulated_data)

signalout = signal;

clear i j k BS_pks BS_locs MVDR_pks MVDR_locs I M eta x u
end

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

Real_data = table;
Real_data.Method = method';
Real_data.Angle = angle';
display(Real_data)

end

function [P1 f pk] = sig_fft(signal, Fs)

FFTD = signal;

Ts = 1/Fs;
L = size(FFTD);
L(2) = [];
ti = (0:(L-1))*Ts;

Y = fft(FFTD);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

figure
plot(f,P1) 

[pks locs] = findpeaks(P1);
[M I] = max(pks);
pk = f(locs(I));


end