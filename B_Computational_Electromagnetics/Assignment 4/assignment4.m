clear all;
close all;
clc;
xdim = 200;
time_tot = 380;
xsource = 1;

S = 0.97;

epsilon0 = 1;
mu0 = 1;
c = 1;
delta = 1;
deltat = S*delta/c;

Ez=zeros(1,xdim);
Hy=zeros(1,xdim);

sigma = zeros(1, xdim);
mu = mu0*ones(1,xdim);
epsilon = epsilon0*ones(1,xdim);

x1 = floor((2*xdim/3)-14);
x2 = floor((2*xdim/3)+15);

epsilon(1,1:(xdim/2)) = 1.0;
epsilon(1,(xdim/2):xdim) = 1.0;
epsilon(x1: x2) = 6;
mu(x1:x2) = 2;
sigma(x1:x2) = 0.000;

gaussian = 1;
sine=0;
impulse = 0;
lambda = 20;
freq = c/lambda;

detector0 = zeros(1, time_tot);
detector1 = zeros(1, time_tot);

n=1:1:time_tot;
N_lambda = 20;
Ezsource = gaussmf(1:1:53, [7.9, 26]);
Ezsource(53:time_tot) = 0;

if sine == 1
    Ezsource = cos(2*pi*freq*n);
end

if impulse == 1
    Ezsource(1) = 1;
end

mur = (S-1)/(S+1);

for n=1:1:time_tot;
    
        if n <= 53
            Ez(xsource) = Ezsource(n);
        else
               detector0(n) = Ez(1);
               detector1(n) = Ez(xdim);
        end
        
    Hy(1:xdim-1) = Hy(1:xdim-1)+(deltat./(delta*mu(1:xdim-1))).*(Ez(2:xdim)-Ez(1:xdim-1));

    oldEz2 = Ez(2);
    oldEzxdim = Ez(xdim);
    oldEzxdim1 = Ez(xdim-1);
    oldEz1 = Ez(1);
    
    Ez(2:xdim) = Ez(2:xdim) - (deltat.*sigma(2:xdim)./epsilon(2:xdim)).*Ez(2:xdim) + (deltat./(delta*epsilon(2:xdim))).*(Hy(2:xdim)-Hy(1:xdim-1));
    
    Ez(1) = oldEz2 + mur*(Ez(2) - oldEz1);
    Ez(xdim) = oldEzxdim1 + mur*(Ez(xdim-1) - oldEzxdim);
    
    %Ez(xsource)=1;
        
    plot(1:1:xdim,Ez,'color','k','linewidth',2);
    grid on;
    axis([0, xdim, -2, 2]);
    getframe;
end

R = fft(detector0);
T = fft(detector1);
S = fft(Ezsource);

L = size(R); L = L(2);
P2 = abs(R/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = (1/deltat)*(0:(L/2))/L;

f = f/(0.333e-10);

R = abs(R/L);
R = R(1:L/2+1);
R(1:end-1) = 2*R(1:end-1);

T = abs(T/L);
T = T(1:L/2+1);
T(1:end-1) = 2*T(1:end-1);

S = abs(S/L);
S = S(1:L/2+1);
S(1:end-1) = 2*S(1:end-1);

plot(f, R);
axis([8e7, 3e9, 0, .03]);
title('Source spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

R = R(1:25); S = S(1:25); T = T(1:25);
%plot((abs(T).^2), 'r');

figure;
plot(f(1:25), (T.^2)./(S.^2), 'r');
%axis([0, 1e9, 0, 1]);
hold on;
plot(f(1:25), (R.^2)./(S.^2), 'b');
grid on;
axis([1e8, 1e9, -0.25, 1.25]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency vs Reflectance, Transmittance');
legend('Transmittance', 'Reflectance');
hold off;

figure;
X = (R.^2 + T.^2)./(S.^2);

plot(f(1:25), X, 'k');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency vs Energy Ratio');
legend('Transmittance^2 + Reflectance^2', 'Location', 'southeast');
axis([0, 1e9, 0, 1.2])
grid on;
%critical dimention : smallest critical fearure the grid should capture..
%critical timestep : minimum timestep to choose such that S is nowhere
%greater than one

%imlement mur's absorbing boundary condition