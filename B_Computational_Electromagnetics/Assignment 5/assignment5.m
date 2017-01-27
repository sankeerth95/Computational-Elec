clear all;
close all;
clc;
xdim = 300;
time_tot = 700;
xsource = 1;

S = 1;

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

x1 = floor((2*xdim/3)-29);
x2 = floor((2*xdim/3)+30);

epsilon(1,1:(xdim/2)) = 1.0;
epsilon(1,(xdim/2):xdim) = 1.0;
epsilon(x1: x2) = 12;
mu(x1:x2) = 1;
sigma(x1:x2) = 0.000;

gaussian = 0;
sine=1;
impulse = 0;
%lambda = 20;
freq = 1/25;

detector0 = zeros(1, time_tot);
detector1 = zeros(1, time_tot);

n=1:1:time_tot;
N_lambda = 20;
Ezsource = gaussmf(1:1:53, [7.9, 26]);

if sine == 1
    Ezsource = sin(2*pi*freq*n);
end

Ezsource(74:time_tot) = 0;

if impulse == 1
    Ezsource(1) = 1;
end

mur = (S-1)/(S+1);

for n=1:1:time_tot;
    
        if n <= 75
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

f = 2*f/(0.3333e-10);

R = abs(R/L);
R = R(1:L/2+1);
R(1:end-1) = 2*R(1:end-1);

T = abs(T/L);
T = T(1:L/2+1);
T(1:end-1) = 2*T(1:end-1);

S = abs(S/L);
S = S(1:L/2+1);
S(1:end-1) = 2*S(1:end-1);

plot(f, (S.^2), 'r');
grid on;

%R = R(1:25); S = S(1:25); T = T(1:25);
%plot((abs(T).^2), 'r');


figure;
plot(f, (T.^2)./(S.^2), 'r');
%axis([0, 1e9, 0, 1]);
hold on;
plot(f, (R.^2)./(S.^2), 'b');
grid on;
axis([2.2e9, 2.6e9, -0.25, 1.25]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency vs Reflectance, Transmittance');
legend('Transmittance', 'Reflectance');
hold off;

figure;
X = (R.^2 + T.^2)./(S.^2);

plot(f, X, 'k');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency vs Energy Ratio');
legend('Transmittance^2 + Reflectance^2', 'Location', 'southeast');
axis([2.2e9, 2.6e9, 0, 1.2])
grid on;
%critical dimention : smallest critical fearure the grid should capture..
%critical timestep : minimum timestep to choose such that S is nowhere
%greater than one

%imlement mur's absorbing boundary condition