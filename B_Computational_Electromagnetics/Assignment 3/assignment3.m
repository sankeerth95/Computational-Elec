clear all;
close all;
clc;
xdim=200;
time_tot=500;
xsource = 1;

S = 0.7;

epsilon0=1;
mu0 = 1;
c = 1;
delta = 1;
deltat=S*delta/c;
 
Ez=zeros(1,xdim);
Hy=zeros(1,xdim);

sigma = zeros(1, xdim);
mu=mu0*ones(1,xdim);
epsilon=epsilon0*ones(1,xdim);

x1 = floor(xdim/3);
x2 = floor(2*xdim/3);

epsilon(1,1:(xdim/2))=1.0;
epsilon(1,(xdim/2):xdim)=1.0;
sigma(x1:x2) = 0.03;


gaussian=1;
sine=0;
impulse=0;
lambda = 20;
freq = c/lambda;

detector0 = zeros(1, time_tot);
detector1 = zeros(1, time_tot);

n=1:1:time_tot;
tstart = 1;
N_lambda = 20;
Ezsource = (10 - 15*cos(n*pi/20) + 6*cos(2*n*pi/20) - cos(3*n*pi/20))/32;
Ezsource(42:time_tot) = 0;

if sine == 1
    Ezsource=cos(2*pi*freq*n);
end

if impulse == 1
    Ezsource(1) = 1;
end

mur = (S-1)/(S+1);

for n=1:1:time_tot;
    
        if n <= 43
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

shift = floor((time_tot+1)/2);
R = circshift(R', shift)';
S = circshift(S', shift)';
T = circshift(T', shift)';

n = 1:1:time_tot;
figure;
plot(n, abs(T).^2, 'r');
hold on;
plot(n, abs(S).^2, 'b');

plot(n, abs(R).^2, 'k');
hold off;

figure;
X = (abs(R).^2 + abs(T).^2)./abs(S).^2;
%X = abs(S).^2;
plot(n(220:280), X(220:280), 'r');
axis([220,280, 0,2]);

%critical dimention : smallest critical fearure the grid should capture..
%critical timestep : minimum timestep to choose such that S is nowhere
%greater than one

%imlement mur's absorbing boundary condition