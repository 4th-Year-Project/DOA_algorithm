Ang = [(0.2*(pi)*57.3);(0.6*(pi)*57.3)];    % Angles of source
R_Ang = [(0.2*(pi));(0.6*(pi))];           % Angles of Source in Radians

f = 8.68e8;
c = 299792458;
wl = c/f;
d = wl/2;
M = 4;          % Number of elements
L = 2;          % Number of sources

theta = [(0.2*pi), (0.6*pi)]; %fake values of arrival

ULA = zeros([1 M]);
a= zeros([1 M]);

for i=1:M
    ULA(i) = i;
    for k=1:L
        u(k, i) = exp((-2i*pi*i*d*cos(theta(k)))/wl);

    end
end

a = circshift(u, 1,L);
a(1:L,1) = 1;
A=a.';

% Ryy = cov(Sim);         % Estimate Covariance
% 
% [Q, R] = qr(Ryy);      % QE factorisation
% 
% Qs = Q(:,1:2);
% Qn = Q(:,3:4);
% 
% QnH = Qn';
% 
% Minpeaks = norm(QnH*A);

% for Theta=0:180
%    theta = [(Theta/57.3), (Theta/57.3)] ; %fake values of arrival
% 
% ULA = zeros([1 M]);
% a= zeros([1 M]);
% 
% for i=1:M
%     ULA(i) = i;
%     for k=1:L
%         u(k, i) = exp((-2i*pi*i*d*cos(theta(k)))/wl);
% 
%     end
% end
% 
% a = circshift(u, 1,L);
% a(1:L,1) = 1;
% A=a.';
% end

