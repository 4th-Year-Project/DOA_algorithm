Cal_dat = readtable('Total_data.csv');
Cal_dat = table2array(Cal_dat);

signal = readtable('data.csv');
signal = table2array(signal)';

v_signal = (2*signal)/255;         % Signal range 0 - 2V as measured by ADC

t = (Cal_dat(:,1)+(1E-4));

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



% signal parameters
T = 1E-7;           % sample time
fs = 1/T;           % sample frequency
f0 = 10E3;          % signal frequency
fc = 868E6          % carrier frequency

fprintf('\n   I signal phase difference: \n')

for k = 1:3
    
x = norm_I(:,1)';
y = norm_I(:,(k+1))';
    
PhDiff = phdiffmeasure(x, y);
PhDiff = PhDiff*180/pi;

PhDiffstr = num2str(PhDiff);
disp(['Phase difference 1 - ' num2str(k+1) ' = ' PhDiffstr ' deg'])

I_ph_diff(k) = PhDiff;

end

fprintf('\n   Q signal phase difference: \n')

for k = 1:3
    
x = norm_Q(:,1)';
y = norm_Q(:,(k+1))';
    
PhDiff = phdiffmeasure(x, y);
PhDiff = PhDiff*180/pi;

PhDiffstr = num2str(PhDiff);
disp(['Phase difference 1 - ' num2str(k+1) ' = ' PhDiffstr ' deg'])

Q_ph_diff(k) = PhDiff;

end

I_t_diff = ((1/f0)* (I_ph_diff+360))/360;

for k=1:3
    IT_diff(k) = round(I_t_diff(k), 7);
    IT_ind(k)= find(t == IT_diff(k));
end

fprintf('\n\n')

for k = 1:4
[pks,locs] = findpeaks(norm_I(:,k));
Start(k) = find(norm_I(:,k) == max(pks(1:110)));
end

for k = 1:4
[pks,locs] = findpeaks(norm_Q(:,k));
Start(k) = find(norm_Q(:,k) == max(pks(1:110)));
end

start_delay = T*Start;

for k = 1:4
    a = round((2*k)- 1);
    b = round(2*k);
      
    temp = v_signal(:,(a)) + (1j*v_signal(:,(b)));
    
    c_signal(:,k) = temp;          % Complex signal
    
    I_signal(:,k) = v_signal(:,(a));
    Q_signal(:,k) = v_signal(:,(b));
end

for k = 1:4
    temp = findpeaks(real(c_signal(:,k)) - 1);
    average_amp(k) = mean(abs(temp));
    
    temp2 = findpeaks(I_signal(:,k)-1);
    aa_I(k) = mean(abs(temp2));
    
    temp3 = findpeaks(Q_signal(:,k)-1);
    aa_Q(k) = mean(abs(temp3));
end

for k = 1:4
    nc_signal(:,k) = (c_signal(:,k)-1) * (max(average_amp)/average_amp(k));
    
    norm_I_sig(:,k) = (I_signal(:,k)-1) * max(aa_I)/aa_I(k);
    norm_Q_sig(:,k) = (Q_signal(:,k)-1) * max(aa_Q)/aa_Q(k);
end

real_fs = 500E3;

sample_delay = real_fs * start_delay;
sample_delay = round(sample_delay, 0);

figure
hold on
for k = 1:4
plot(norm_I_sig((sample_delay(k)+1):end, k))
title('Phased aligned & Normalised I signals')
ylim(limits)
end

figure
plot(norm_I_sig)
title('Unaligned Normalised I signals')
ylim(limits)

clear temp temp1 temp2 temp3 limits k aa_I aa_Q


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