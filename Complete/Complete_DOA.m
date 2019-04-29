%%%%%%%% signal parameters
T = 1E-7;           % sample time
fs = 1/T;           % calibration sample frequency
f0 = 10E3;          % signal frequency
fc = 868E6;         % carrier frequency
rfs = 500E3;        % real signal sample frequency   

Cal_dat = readtable('Total_data.csv');
Cal_dat = table2array(Cal_dat);

signal = readtable('data.csv');
signal = table2array(signal)';

%%%%%%%% Convert signal range 0 - 2V as measured by ADC
signal = (2*signal)/255;         

%%%%%%%% plot signals? true/false   
plot_f = false;                  

%%%%%%%% Plot normalised calibration signals   
[Cal_norm_I, Cal_norm_Q] = cal(Cal_dat, plot_f);                   

%%%%%%%% Calculate phase of calibration signals   
fprintf('\n   Q signal phase difference: \n')
for k = 1:3
    
x = Cal_norm_Q(:,1)';
y = Cal_norm_Q(:,(k+1))';
    
PhDiff = phdiffmeasure(x, y);
PhDiff = PhDiff*180/pi;

PhDiffstr = num2str(PhDiff);
disp(['Phase difference 1 - ' num2str(k+1) ' = ' PhDiffstr ' deg'])
end
fprintf('\n')

%%%%%%%% calculate phase, convert to time delay
start_delay = s_delay(Cal_norm_I, Cal_norm_Q, T);

%%%%%%%% Use time delay to calculate sample delay of real data

sample_delay = rfs * start_delay;                                   
sample_delay = round(sample_delay, 0);                              

%%%%%%%% Normalise real signals
[ni_sig, nq_sig] = sig_norm(signal);            

%%%%%%%% Align real signals
[ai_signal, aq_signal] = sig_align(ni_sig, nq_sig, sample_delay);

%%%%%%%% DOA calculation with simulated signals
DOA(-30, 0.1, false);
%%%%%%%% 

for k = 1:4  
    doa_data(:,k) = ai_signal(:,k) + (1j*aq_signal(:,k));
end

%%%%%%%% DOA calculation with real signals
[ANGLE_BS, ANGLE_MVDR, ANGLE_QR] = DOA_calc(doa_data, 170, true);


clear x y plot_f PhDiff PhDiffstr k a b

function [NI, NQ] = cal(Cal_dat, plot_f)

t = (Cal_dat(:,1)+(1E-4));       % Time data for calibration

for k = 1:4
    I(:,k) = Cal_dat(:,k+1);
    Q(:,k) = Cal_dat(:,k+5);
end

for k = 1:4
    temp1 = findpeaks(I(:,k)-1);
    aa_I(k) = mean(abs(temp1));
    
    temp2 = findpeaks(Q(:,k)-1);
    aa_Q(k) = mean(abs(temp2));
end

for k = 1:4
    norm_I(:,k) = 2 * (I(:,k)-1) * max(aa_I)/aa_I(k);
    norm_Q(:,k) = 2 * (Q(:,k)-1) * max(aa_Q)/aa_Q(k);
end

[pks,locs] = findpeaks( norm_I(:,1));
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

function start_delay = s_delay(NI, NQ, T)
    for k = 1:4
    [pks,locs] = findpeaks(NI(:,k));
    Start(k) = find(NI(:,k) == max(pks(1:110)));
    end

    for k = 1:4
    [pks,locs] = findpeaks(NQ(:,k));
    Start(k) = find(NQ(:,k) == max(pks(1:110)));
    end

start_delay = T*Start;

end

function [NI, NQ] = sig_norm(signal)
for k = 1:4
    a = round((2*k)- 1);
    b = round(2*k);

    I_signal(:,k) = signal(:,(a));
    Q_signal(:,k) = signal(:,(b));
end

for k = 1:4  
    temp1 = findpeaks(I_signal(:,k)-1);
    aa_I(k) = mean(abs(temp1));
    
    temp2 = findpeaks(Q_signal(:,k)-1);
    aa_Q(k) = mean(abs(temp2));
end

for k = 1:4  
    norm_I_sig(:,k) = (I_signal(:,k)-1) * max(aa_I)/aa_I(k);
    norm_Q_sig(:,k) = (Q_signal(:,k)-1) * max(aa_Q)/aa_Q(k);
end

	NI = norm_I_sig;
	NQ = norm_Q_sig;
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

function DOA(given_angle, n_power, plot_f)

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
QR_DOA = [MVDR_pks; MVDR_locs];
[M,I] = max(QR_DOA(1,:));
QR_DOA = Y(QR_DOA(2,I));

method = {'Beamscan', 'MVDR', 'QR '};
angle = [BS_DOA, MVDR_DOA, QR_DOA];

Real_data = table;
Real_data.Method = method';
Real_data.Angle = angle';
display(Real_data)

end