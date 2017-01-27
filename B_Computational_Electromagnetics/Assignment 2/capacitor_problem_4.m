%initialize some variables here
clear all;
n = 200;

xdim = n;
ydim = n;

epsr = 12.0;

V_now = 5*ones(xdim, ydim);
V_prev = 5*ones(xdim, ydim);

capx1 = 75;
capx2 = 125;
capy1 = 50;
capy2 = 150;
capytop1 = 95;
capytop2 = 105;
capybot1 = 50;
capybot2 = 150;

V_now(capx1, capy1:capy2) = 10;
V_now(capx2, capy1:capy2) = 0;
error = max(max(abs(V_now - V_prev)));
iter = 0;

k1 = 2.0/(1.0+epsr);
k2 = 2.0*epsr/(1.0+epsr);

%iterate as long as the change in values obtained
%in negligible. In this case, max value change is 0.001

while(error > 0.0001)
iter = iter+1;

V_now(2:xdim-1, 2:ydim-1) = (V_now(1:xdim-2, 2:ydim-1) + V_now(3:xdim, 2:ydim-1) + V_now(2:xdim-1, 1:ydim-2) + V_now(2:xdim-1, 3:ydim))*0.25;

V_now(capx1, capytop1:capytop2) = 10;
V_now(capx2, capybot1:capybot2) = 0;

V_now(capx1+1:capx2-1, capy1) = (k2*V_prev(capx1+1:capx2-1, capy1+1)+k1*V_prev(capx1+1:capx2-1, capy1-1)+V_prev(capx1+2:capx2, capy1)+V_prev(capx1:capx2-2, capy1))*0.25;
V_now(capx1+1:capx2-1, capy2) = (k1*V_prev(capx1+1:capx2-1, capy2+1)+k2*V_prev(capx1+1:capx2-1, capy2-1)+V_prev(capx1+2:capx2, capy2)+V_prev(capx1:capx2-2, capy2))*0.25;
V_now(capx1, capy1:capy2) = (k1*V_prev(capx1+1, capy1:capy2)+k2*V_prev(capx1-1, capy1:capy2)+V_prev(capx1, capy1+1:capy2+1)+V_prev(capx1, capy2-1:capy2-1))*0.25;

V_now(capx1, capytop1:capytop2) = 10;
V_now(capx2, capybot1:capybot2) = 0;

V_now(2:xdim-1, 1) =  (2*V_now(2:ydim-1, 2) + V_now(3:xdim, 1) + V_now(1:xdim-2, 1))*0.25;
V_now(2:xdim-1, ydim) =  (2*V_now(2:ydim-1, ydim-1) + V_now(3:xdim, ydim) + V_now(1:xdim-2, ydim))*0.25;
V_now(1, 2:ydim-1) =  (2*V_now(2, 2:ydim-1) + V_now(1, 3:ydim) + V_now(1, 1:ydim-2))*0.25;
V_now(xdim, 2:ydim-1) =  (2*V_now(xdim-1, 2:ydim-1) + V_now(xdim, 3:ydim) + V_now(xdim, 1:ydim-2))*0.25;

V_now(1, 1) = (V_now(2, 1) + V_now(1, 2))*0.5;
V_now(xdim, 1) = (V_now(xdim-1, 1) + V_now(xdim, 2))*0.5;
V_now(xdim, ydim) = (V_now(xdim-1, ydim) + V_now(xdim, ydim-1))*0.5;
V_now(1, ydim) = (V_now(1, ydim-1) + V_now(2, ydim))*0.5;


error = max(max(abs(V_now-V_prev)));
V_prev = V_now;

end

imagesc(V_now);

%vector plot plotted using quiver after
%normallizing the vector values
figure;
title('Vector plot of electric field(normallised) at resolution 30 points');
[ex, ey] = gradient(V_now(1:5:xdim, 1:5:ydim), 1, 1);
sz = max(max(sqrt(ex.^2 + ey.^2)));
ex = -ex;
ey = -ey;
exx = ex./sz;
eyy = ey./sz;
quiver(exx, eyy, 'k-');