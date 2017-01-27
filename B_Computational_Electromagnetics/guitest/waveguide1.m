clear all;
clc;

time_tot = 1500;
maxepsr = 11.56;
maxsigma = 5000;

S = 1/sqrt(2);

epsilon0 = (1e-9)/(36*pi);
mu0 = 4*pi*1e-7;
c = 3e+8;

delta = 15e-9;
deltat = S*delta/c;

I = imread('imageinput.png');
I = im2double(I);

[xdim, ydim] = size(I(:, :, 1)');

detector = zeros(time_tot);
storeMovie = zeros(xdim, ydim, time_tot);

Ez = zeros(xdim ,ydim);
Hy = zeros(xdim, ydim);
Hx = zeros(xdim, ydim);

IB = I(:, :, 3)';
IG = I(:, :, 2)';
IR = I(:, :, 1)';
epsilonr = 1 + maxepsr*IB;
sigma = maxsigma*IG;
sourceCoors = (IR > 0.5);

epsilon = epsilon0*epsilonr;
mu = mu0*ones(xdim, ydim);
%sigma = 4e-4*ones(xdim, ydim);
sigma_star = 4e-4*ones(xdim, ydim);

sine = 1;
lambda0 = 1.2e-6;
frequency = 3e8/lambda0;

imagesc(epsilon');
pause(2);

for n = 1:1:time_tot

    n1 = 1;
    n2 = xdim - 1;
    n11 = 1;
    n21 = ydim - 1;

    Hy(n1:n2,n11:n21)=((mu(n1:n2,n11:n21)-0.5*deltat*sigma_star(n1:n2,n11:n21))./(mu(n1:n2,n11:n21)+0.5*deltat*sigma_star(n1:n2,n11:n21))).*Hy(n1:n2,n11:n21)+(deltat/delta)*(Ez(n1+1:n2+1,n11:n21)-Ez(n1:n2,n11:n21))./(mu(n1:n2,n11:n21)+0.5*deltat*sigma_star(n1:n2,n11:n21));
    Hx(n1:n2,n11:n21)=((mu(n1:n2,n11:n21)-0.5*deltat*sigma_star(n1:n2,n11:n21))./(mu(n1:n2,n11:n21)+0.5*deltat*sigma_star(n1:n2,n11:n21))).*Hx(n1:n2,n11:n21)-(deltat/delta)*(Ez(n1:n2,n11+1:n21+1)-Ez(n1:n2,n11:n21))./(mu(n1:n2,n11:n21)+0.5*deltat*sigma_star(n1:n2,n11:n21));
    
    
    %Vector update instead of for-loop for Ez field
    Ez(n1+1:n2+1,n11+1:n21+1)=((epsilon(n1+1:n2+1,n11+1:n21+1)-0.5*deltat*sigma(n1+1:n2+1,n11+1:n21+1))./(epsilon(n1+1:n2+1,n11+1:n21+1)+0.5*deltat*sigma(n1+1:n2+1,n11+1:n21+1))).*Ez(n1+1:n2+1,n11+1:n21+1)+(deltat/delta)*(Hy(n1+1:n2+1,n11+1:n21+1)-Hy(n1:n2,n11+1:n21+1)-Hx(n1+1:n2+1,n11+1:n21+1)+Hx(n1+1:n2+1,n11:n21))./(epsilon(n1+1:n2+1,n11+1:n21+1)+0.5*deltat*sigma(n1+1:n2+1,n11+1:n21+1));

%    Ez(xdim, 1:ydim) = Ez(xdim-1, 1:ydim);
   % Ez(1:xdim, ydim) = Ez(1:xdim,ydim-1);
    
    detector(n) = Ez(120, 1);

    if sine == 1
        tstart = 1;
        N_lambda = c/(frequency*delta);
        Ez(sourceCoors) = sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
%        Ez(xsource, ydim/2+70:ydim/2+80) = sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));

    end

 %   storeMovie(:, :, n) = Ez;
   if mod(n, 2) == 0
       imagesc(Ez',[-1, 1]);
       colorbar; %colormap(jet);
    end
   getframe();

end

X = fft(detector);

%for n = 1:1:time_tot
    
%    imagesc(storeMovie(:, :, n)',[-1, 1]);
%    colorbar; %colormap(jet);
%    getframe();
%end


