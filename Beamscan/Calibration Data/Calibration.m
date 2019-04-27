Cal_dat = readtable('Total_data.csv');
Cal_dat = table2array(Cal_dat);

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

figure
hold on
for k = 1:4
plot(norm_I(Start(k):end, k))
title('Phased aligned & Normalised I signals')
ylim(limits)
end

for k = 1:4
[pks,locs] = findpeaks(norm_Q(:,k));
Start(k) = find(norm_Q(:,k) == max(pks(1:110)));
end

figure
hold on
for k = 1:4
plot(norm_Q(Start(k):end, k))
title('Phased aligned & Normalised Q signals')
ylim(limits)
end

%%% Built in matlab signal align

% [I_A, I_B] = alignsignals(norm_I(:,1),norm_I(:,2));
% [I_A, I_C] = alignsignals(norm_I(:,1),norm_I(:,3));
% [I_A, I_D] = alignsignals(norm_I(:,1),norm_I(:,4));
% 
% figure
% hold on
% plot(I_A)
% plot(I_B)
% plot(I_C)
% plot(I_D)
% title('MATLAB Phased aligned & Normalised I signals')

%%% RF 

% carry00 = cos(2*pi*fc*t);
% carry90 = sin(2*pi*fc*t);

% for m = 1:4
% for k = 1:2000
% 
% RF(k,m) = (norm_I(k,m)*carry00(k)) - (norm_Q(k,m)*carry90(k));
% 
% end
% end
% 
% figure
% subplot(2,2,1)
% plot(RF(:,1))
% 
% subplot(2,2,2)
% plot(RF(:,2))
% 
% subplot(2,2,3)
% plot(RF(:,3))
% 
% subplot(2,2,4)
% plot(RF(:,4))


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