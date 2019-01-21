% function Ang = QR3_theta(n) 
Ang = [50];    % Angles of source

f = 8.68e8;     % Frequency
c = 299792458;  % Propagation velocity
wl = c/f;       % Wavelength
d = wl/2;       % Distance between antennas
M = 4;          % Number of elements
L = 1;          % Number of sources

sim('ULA_4ant_sim2a', 0.1);      % Simulate 4x1 ULA signal

for t=1:180 
    theta(t) = t*(pi)/180;  % scanning angle in radians
    A(:,:,t) = zeros(M,L);
    
    for i=1:M
%     ULA(i) = i;
        for k=1:L
        u(k, i, t) = exp((-2*i*pi*d*cos(theta(t)))/wl);
        end             % Array response vectors
    end
    
    a(:,:,t) = circshift(u(:,:,t), 1,L);
    a(1:L,1, t) = 1;
    A(:,:,t) = a(:,:,t).';      % Array response matrix     
    
end

Ryy = cov(Sim);         % Estimate Covariance
[UN, U] = lu(Ryy);       % QR factorisation

Us = U(1:L, 1:M);

Uss = Us';

Us1 = Uss(1:3, 1);  
Us2 = Uss(2:4, 1);

Om = Us2/Us1;
% Om =  (inv(Us1'*Us1))*Us1*Us2
Ei = eig(Om);
Ai = angle(Ei(3));

DOA = - acos((Ai)/(2*(pi)*d));
% ANGLE = 180*angle(DOA)/pi
DOA_Angle = angle(DOA) *180 / pi

% for i=1:3    
% D(i) = - acos((Ai(i))/(2*(pi)*d));
% DOA(i) = D(i)*pi/180;
% end

% fig = figure;
% hax = axes;
% hold on
% plot(DOA1)
% title('Direction of Arrival') 
% xlabel('Angle (Degrees)')
% ylabel('DOA amplitude')
% 
% SP = round(Ang, 0);
% line([SP SP], get(hax,'YLim'),'Color',[1 0 0])
