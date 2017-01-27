clear all;
clc;

time_tot = 1000;
xdim = 750;
ydim = 750;

c = 3e8;

S = 1/sqrt(2);
delta = 0.01e-6;
%delta = 0.2e-7;
deltat = (S/c)*delta;

epsilon0 = (1e-9)/(36*pi);
epsx = epsilon0*ones(xdim,ydim);
epsy = epsx;
epsz = epsx;

%ABC
x1 = floor(xdim/6);
y1 = floor(ydim/6);
x2 = 5*x1; y2 = 5*y1;

ABC = ones(xdim, ydim);
ABC(x1:x2, y1:y2) = 0;
ABC_region = ABC > 0.5;

sigma = zeros(xdim, ydim);

%convolution non-delta componnts
N0 = 15;
epsxx = epsilon0*zeros(xdim, ydim, N0);
epsyy = epsxx;
epszz = epsxx;

omegap = 2000e12; vc = 57e12;
for i = 1:1:N0
    epsxx(x1:x2, y1:ceil((y1+y2)/2), i) = epsilon0*(omegap^2)*(exp(-vc*i*deltat));
    epsyy(x1:x2, y1:ceil((y1+y2)/2), i) = epsilon0*(omegap^2)*(exp(-vc*i*deltat));
    epszz(x1:x2, y1:ceil((y1+y2)/2), i) = epsilon0*(omegap^2)*(exp(-vc*i*deltat));

end


mu0 = 4*pi*1e-7;
mux = mu0*ones(xdim, ydim);
muy = mu0*ones(xdim, ydim);
muz = mu0*ones(xdim, ydim);


sigma0 = 2e5;
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

abc = 1.1;

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


frequency = 1000e12;
ysource = floor((y1+y2)/2) - 100 : floor((y1+y2)/2) + 100;
xsource = x1;
Ezsource = zeros(1, time_tot);
tc = gauspuls('cutoff',2500e12,0.2,[],-40);
t = -tc : deltat : tc;
yi = 2*gauspuls(t,2500e12,0.2);
plot(t,yi);
pause(1);
sze = size(yi);
for i = 1:1:sze(2)
    Ezsource(i) = 2*sin(2*pi*2000e12*i*deltat);
    %Ezsource(i) = yi(i);
end
Ezsource(160:end) = 0;

Ez = zeros(xdim, ydim);
Hx = zeros(xdim, ydim);
Hy = zeros(xdim, ydim);


imagesc((sigma./sigma0)',[-1, 1]);
colorbar; %colormap(jet);\
pause(2);


imagesc((epszz(:, :, 1))',[-1, 1]);
colorbar; %colormap(jet);\
pause(2);

convStore = zeros(xdim, ydim, N0);

for n = 1:1:time_tot

    nx = 1; nxx = xdim-1;
    ny = 1; nyy = ydim-1;
    
    Hy(nx:nxx, ny:nyy) = Hy(nx:nxx,ny:nyy) + (deltat/delta)*(Ez(nx+1:nxx+1, ny:nyy)-Ez(nx:nxx,ny:nyy))./muy(nx:nxx,ny:nyy);
    Hx(nx:nxx, ny:nyy) = Hx(nx:nxx,ny:nyy) - (deltat/delta)*(Ez(nx:nxx, ny+1:nyy+1)-Ez(nx:nxx,ny:nyy))./mux(nx:nxx,ny:nyy);
    
    M = zeros(xdim-1, ydim-1);
    for i = 1:1:N0-1
        
        M = M + convStore(nx+1:nxx+1, ny+1:nyy+1, i).*epszz(nx+1:nxx+1, ny+1:nyy+1, i+1);    
    end
    
    Ez(nx+1:nxx+1, ny+1:nyy+1) = (1./(sigma(nx+1:nxx+1, ny+1:nyy+1)*deltat + epsz(nx+1:nxx+1, ny+1:nyy+1) + epszz(nx+1:nxx+1, ny+1:nyy+1, 1)*deltat*deltat)).*(epsz(nx+1:nxx+1, ny+1:nyy+1).*Ez(nx+1:nxx+1, ny+1:nyy+1) - deltat*deltat*M + (deltat/delta)*(Hy(nx+1:nxx+1,ny+1:nyy+1)-Hy(nx:nxx,ny+1:nyy+1)-Hx(nx+1:nxx+1,ny+1:nyy+1)+Hx(nx+1:nxx+1,ny:nyy)));
    
    convStore(:, :, 2:N0) = convStore(:, :, 1:N0-1);
    convStore(:, :, 1) = Ez;
    
    Ez(xsource, ysource) = Ezsource(n);
   
    if mod(n, 2) == 1
       imagesc(Ez',[-1, 1]);
       colorbar; %colormap(jet);
    end
    getframe();
    
end
