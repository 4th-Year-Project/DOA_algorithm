signal = readtable('data.csv');
signal = table2array(signal)';

v_signal = (2*signal)/255;         % Signal range 0 - 2V as measured by ADC

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
    
    norm_I(:,k) = (I_signal(:,k)-1) * max(aa_I)/aa_I(k);
    norm_Q(:,k) = (Q_signal(:,k)-1) * max(aa_Q)/aa_Q(k);
end

plot(norm_I)
hold on



% signal parameters
fs = 500E3;
f0 = 10E3;
T = 200/fs;

% preparation of the time vector
N = round(T*fs);
t = (0:N-1)/fs;

% plot(t, real(nc_signal))
% x = real(nc_signal(:,1))';
% y = real(nc_signal(:,4))';
% 
% 
% % phase difference calculation
% PhDiff = phdiffmeasure(x, y);
% PhDiff = PhDiff*180/pi;
% % 
% 
% % display the phase difference
% PhDiffstr = num2str(PhDiff);
% disp(['Phase difference Y->X = ' PhDiffstr ' deg'])

for k = 1:3
    
x = norm_I(:,1)';
y = norm_I(:,(k+1))';
    
PhDiff = phdiffmeasure(x, y);
PhDiff = PhDiff*180/pi;

PhDiffstr = num2str(PhDiff);
disp(['Phase difference 1 - ' num2str(k+1) ' = ' PhDiffstr ' deg'])
    
end

% figure
% plot(t, real(nc_signal))

clear temp temp2 temp3 a b k PhDiff PhDiffstr


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
