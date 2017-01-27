clear all;
clc;

epsilon0 = (1e-9)/(36*pi);
mu0 = 4*pi*1e-7;
omegap = 2000e12;
vc = 2*pi*57e12;

TE = 0; TM = 1;
sine = 0; gaussian = 1;
time_tot = 1000;
N0 = 20;
chimax = 7.0;

c = 3e8;
S = 1/sqrt(2);
delta = 1.3e-8;
deltat = (S/c)*delta;

I = imread('sourceDrude.png');
I = im2double(I);
[xdim, ydim] = size(I(:, :, 1)');
IB = I(:, :, 3)';
IG = I(:, :, 2)';
IR = I(:, :, 1)';

drude = repmat((IB > 0.5), [1, 1, N0]);
sourceCoors = (IR > 0.5);

I1 = imread('permittivity.png');
I1 = im2double(I1);
IB1 = I1(:, :, 3)';
IG1 = I1(:, :, 2)';
IR1 = I1(:, :, 1)';

epsx = epsilon0*(1 + IB1*chimax);
epsy = epsilon0*(1 + IR1*chimax);
epsz = epsilon0*(1 + IG1*chimax);

%convolution non-delta componnts
epsxx = epsilon0*zeros(xdim, ydim, N0);
epsyy = epsxx;
epszz = epsxx;

for i = 1:1:N0
  
    epsxx(drude) = epsilon0*(omegap^2)*(exp(-vc*(i)*deltat));
    epsyy(drude) = epsilon0*(omegap^2)*(exp(-vc*(i)*deltat));
    epszz(drude) = epsilon0*(omegap^2)*(exp(-vc*(i)*deltat));
end

mux = mu0*ones(xdim, ydim);
muy = mu0*ones(xdim, ydim);
muz = mu0*ones(xdim, ydim);




zsource = zeros(1, time_tot);
tc = gauspuls('cutoff',2500e12,0.2,[],-40);
t = -tc : deltat : tc;
yi = 2*gauspuls(t,2500e12,0.2);
plot(t,yi);
pause(1);
sze = size(yi);
for i = 1:1:sze(2)
    
    if sine == 1
        zsource(i) = 2*sin(2*pi*2500e12*i*deltat);
    elseif gaussian == 1
        zsource(i) = yi(i);
    end
end
zsource(160:end) = 0;




%ABC
x1 = floor(xdim/6);
y1 = floor(ydim/6);
x2 = 5*x1; y2 = 5*y1;

ABC = ones(xdim, ydim);
ABC(x1:x2, y1:y2) = 0;
ABC_region = ABC > 0.5;

sigma = zeros(xdim, ydim);
sigma0 = 1e5;

for x = 1:1:x1
    sigma(x, y1:y2) = sigma0*(x/x1 - 1)^2;
end
for x = x2:1:xdim
    sigma(x, y1:y2) = sigma0*((x - x2)^2)/((x2-xdim)^2);
end
for y = 1:1:y1
    sigma(x1:x2, y) = sigma0*(y/y1 - 1)^2;
end
for y = y2:1:ydim
    sigma(x1:x2, y) = sigma0*((y - y2)^2)/((ydim-y2)^2);
end

sigma(1:x1, 1:y1) = (repmat(sigma(1:x1, y1+1), 1, y1) + repmat(sigma(x1+1, 1:y1), x1, 1));
sigma(x2:xdim, 1:y1) = (repmat(sigma(x2:xdim, y1+1), 1, y1) + repmat(sigma(x1+1, 1:y1), xdim-x2+1, 1));
sigma(x2:xdim, y2:ydim) = (repmat(sigma(x2:xdim, y1+1), 1, ydim-y2+1) + repmat(sigma(x1+1, y2:ydim), xdim-x2+1, 1));
sigma(1:x1, y2:ydim) = (repmat(sigma(1:x1, y1+1), 1, ydim-y2+1) + repmat(sigma(x1+1, y2:ydim), x1, 1));

abc = 1.05;

epsx(ABC_region) = epsilon0*abc;
epsy(ABC_region) = epsilon0*abc;
epsz(ABC_region) = epsilon0*(1/abc);

mux(ABC_region) = mu0*abc;
muy(ABC_region) = mu0*abc;
muz(ABC_region) = mu0*(1/abc);

epsx(x2+1, y1:y2) = epsilon0;
epsx(x1:x2, y2+1) = epsilon0;
epsy(x2+1, y1:y2) = epsilon0;
epsy(x1:x2, y2+1) = epsilon0;
epsz(x2+1, y1:y2) = epsilon0;
epsz(x1:x2, y2+1) = epsilon0;




Ez = zeros(xdim, ydim);
Ex = zeros(xdim, ydim);
Ey = zeros(xdim,ydim);
Hx = zeros(xdim, ydim);
Hy = zeros(xdim, ydim);
Hz = zeros(xdim,ydim);



imagesc((sigma./sigma0)',[-1, 1]);
colorbar; %colormap(jet);\
pause(2);


imagesc((epszz(:, :, 1))',[-1, 1]);
colorbar; %colormap(jet);\
pause(2);

convStore = zeros(xdim, ydim, N0);
convStorex = zeros(xdim, ydim, N0);
convStorey = zeros(xdim, ydim,N0);

nx = 1; nxx = xdim-1;
ny = 1; nyy = ydim-1;

for n = 1:1:time_tot
    
    if TE == 1
        Hy(nx:nxx, ny:nyy) = Hy(nx:nxx,ny:nyy) + (deltat/delta)*(Ez(nx+1:nxx+1, ny:nyy)-Ez(nx:nxx,ny:nyy))./muy(nx:nxx,ny:nyy);
        Hx(nx:nxx, ny:nyy) = Hx(nx:nxx,ny:nyy) - (deltat/delta)*(Ez(nx:nxx, ny+1:nyy+1)-Ez(nx:nxx,ny:nyy))./mux(nx:nxx,ny:nyy);
    
        M = zeros(xdim-1, ydim-1);
        for i = 1:1:N0-1
            M = M + convStore(nx+1:nxx+1, ny+1:nyy+1, i).*epszz(nx+1:nxx+1, ny+1:nyy+1, i+1);    
        end
    
        Ez(nx+1:nxx+1, ny+1:nyy+1) = (1./(sigma(nx+1:nxx+1, ny+1:nyy+1)*deltat + epsz(nx+1:nxx+1, ny+1:nyy+1) + epszz(nx+1:nxx+1, ny+1:nyy+1, 1)*deltat*deltat)).*(epsz(nx+1:nxx+1, ny+1:nyy+1).*Ez(nx+1:nxx+1, ny+1:nyy+1) - deltat*deltat*M + (deltat/delta)*(Hy(nx+1:nxx+1,ny+1:nyy+1)-Hy(nx:nxx,ny+1:nyy+1)-Hx(nx+1:nxx+1,ny+1:nyy+1)+Hx(nx+1:nxx+1,ny:nyy)));
    
        convStore(:, :, 2:N0) = convStore(:, :, 1:N0-1);
        convStore(:, :, 1) = Ez;
    
        Ez(sourceCoors) = zsource(n);
    
        if mod(n, 2) == 1
            imagesc(Ez',[-1, 1]);
            colorbar; %colormap(jet);
        end
        getframe();
        
    else
       
        Mx = zeros(xdim-1, ydim-1);
        for i = 1:1:N0-1
            Mx = Mx + convStorex(nx+1:nxx+1, ny+1:nyy+1, i).*epsxx(nx+1:nxx+1, ny+1:nyy+1, i+1);    
        end
    
        My = zeros(xdim-1, ydim-1);
        for i = 1:1:N0-1
            My = My + convStorey(nx+1:nxx+1, ny+1:nyy+1, i).*epsyy(nx+1:nxx+1, ny+1:nyy+1, i+1);    
        end
    
        Ey(nx+1:nxx+1, ny+1:nyy+1) = (1./(epsy(nx+1:nxx+1, ny+1:nyy+1) + epsyy(nx+1:nxx+1, ny+1:nyy+1)*deltat*deltat + sigma(nx+1: nxx+1, ny+1:nyy+1)*deltat)).*(epsy(nx+1:nxx+1,ny+1:nyy+1).*Ey(nx+1:nxx+1, ny+1:nyy+1) - My*deltat*deltat - (deltat/delta)*(Hz(nx+1:nxx+1, ny:nyy)-Hz(nx:nxx,ny:nyy)));
        Ex(nx+1:nxx+1, ny+1:nyy+1) = (1./(epsx(nx+1:nxx+1, ny+1:nyy+1) + epsxx(nx+1:nxx+1, ny+1:nyy+1)*deltat*deltat + sigma(nx+1: nxx+1, ny+1:nyy+1)*deltat)).*(epsx(nx+1:nxx+1,ny+1:nyy+1).*Ex(nx+1:nxx+1, ny+1:nyy+1) - Mx*deltat*deltat + (deltat/delta)*(Hz(nx:nxx, ny+1:nyy+1)-Hz(nx:nxx,ny:nyy)));

        Hz(nx:nxx, ny:nyy) = Hz(nx:nxx, ny:nyy) + (deltat/delta)*(Ey(nx+1:nxx+1,ny+1:nyy+1) - Ey(nx:nxx,ny+1:nyy+1) - Ex(nx+1:nxx+1,ny+1:nyy+1) + Ex(nx+1:nxx+1,ny:nyy))./(-mux(nx:nxx, ny:nyy));
    
        convStorex(:, :, 2:N0) = convStorex(:, :, 1:N0-1);
        convStorex(:, :, 1) = Ex;
    
        convStorey(:, :, 2:N0) = convStorey(:, :, 1:N0-1);
        convStorey(:, :, 1) = Ey;

        Hz(sourceCoors) = zsource(n);

        if mod(n, 2) == 1
            imagesc(Hz',[-1, 1]);
            colorbar; %colormap(jet);
        end
        getframe();

    end
        
end
