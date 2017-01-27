clear all;
clc; close all;

xdim = 300;
ydim = 300;

time_tot = 200;

xsource = 250;
ysource = 250;

S = 1/(sqrt(2));

epsilon0 = 1;
mu0 = 1;
c = 1;

delta = 1;
deltat = S*delta/c;

Hz = zeros(xdim, ydim);
Ey = zeros(xdim, ydim);
Ex = zeros(xdim, ydim);

epsilon = epsilon0*ones(xdim, ydim);
mu = mu0*ones(xdim, ydim);

gaussian = 1;
sine = 0;
frequency = 1.5e13;
impulse = 0;

for n = 1:1:time_tot

    if gaussian == 0 && sine == 0 && n ==1
        Hz(xsource, ysource) = 1;
    end
    
    Ex(1:xdim-1, 1:ydim-1) = Ex(1:xdim-1, 1:ydim-1) + (deltat./(delta.*mu(1:xdim-1, 1:ydim-1))).*(Hz(2:xdim, 1:ydim-1) - Hz(1:xdim-1, 1:ydim-1));
    Ey(1:xdim-1, 1:ydim-1) = Ey(1:xdim-1, 1:ydim-1) - (deltat./(delta.*mu(1:xdim-1, 1:ydim-1))).*(Hz(1:xdim-1, 2:ydim) - Hz(1:xdim-1, 1:ydim-1));

    Hz(2:xdim, 2:ydim) = Hz(2:xdim, 2:ydim) + (deltat./(delta.*epsilon(2:xdim, 2:ydim))).*(Ex(2:xdim, 2:ydim) - Ex(1:xdim-1, 2:ydim) - Ey(2:xdim, 2:ydim) + Ey(2:xdim, 1:ydim-1));

    if impulse == 0
        if gaussian ==0 && sine == 0
            Hz(xsource, ysource) = 1;
        end
        if sine == 1
            tstart = 1;
            N_lambda = 20*2^0.5;
            Hz(xsource, ysource) = sin(((2*pi*(c/(delta*N_lambda))*(n-tstart)*deltat)));
        end
        
        if gaussian == 1
            if n <= 42
                Hz(xsource, ysource) = 5*(10-15*cos(n*pi/20) + 6*cos(2*n*pi/20) - cos(3*n*pi/20))/32;
            else
                Hz(xsource, ysource) = 0;
            end
        end
    else
        Hz(xsource, ysource) = 0;
    end
    
    if mod(n, 5) == 0
    imagesc(Hz', [-1, 1]); colorbar;
    
    title(['\fontsize{20}Color-scaled image plot of Hz for 2D FDTD (TE)']);
    set(gca, 'FontSize', 20);
    getframe;
    end
end


