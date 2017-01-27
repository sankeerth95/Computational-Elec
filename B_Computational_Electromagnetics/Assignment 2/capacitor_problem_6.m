%initialize some variables here
clear all;
n = 20;

xdim = 2*n;
ydim = 7*n;

V_now = zeros(xdim, ydim);
V_prev = zeros(xdim, ydim);
V_now(n, 3*n:4*n) = 10;
iter = 0;
error = max(max(abs(V_now - V_prev)));

epsr = 0;

k1 = 2.0/(1.0+epsr);
k2 = 2.0*epsr/(1.0+epsr);

%iterate as long as the change in values obtained
%in negligible. In this case, max value change is 0.001
while(error > 0.001)
    iter = iter+1;
    
    V_now(2:xdim-1, 2:ydim-1) = (V_now(1:xdim-2, 2:ydim-1)+V_now(3:xdim, 2:ydim-1)+V_now(2:xdim-1, 1:ydim-2)+V_now(2:xdim-1, 3:ydim))*0.25;
    V_now(n, 3*n:4*n) = 10;
    V_now(n, 2:ydim-1) = (V_prev(n, 3:ydim)+V_prev(n, 1:ydim-2)+k1*V_prev(n+1, 2:ydim-1)+k2*V_prev(n-1, 2:ydim-1))*0.25;
    V_now(n, 3*n:4*n) = 10;
    
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