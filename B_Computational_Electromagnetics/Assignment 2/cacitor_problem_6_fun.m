function [C] = cacitor_problem_6_fun(iterMax, w, epsr, width, height)
n = 20;
xdim = height*n;
ydim = width*n;

wl = w*n;
w1 = round(ydim/2 - wl/2);
w2 = round(ydim/2 + wl/2);

V_now = zeros(xdim, ydim);
V_prev = zeros(xdim, ydim);
V_now(n, 3*n:4*n) = 10;
iter = 0;
error = max(max(abs(V_now - V_prev)));

k1 = 2.0/(1.0+epsr);
k2 = 2.0*epsr/(1.0+epsr);

%iterate as long as the change in values obtained
%in negligible. In this case, max value change is 0.001
while(iter < iterMax)
    iter = iter+1;
    
    V_now(2:xdim-1, 2:ydim-1) = (V_now(1:xdim-2, 2:ydim-1)+V_now(3:xdim, 2:ydim-1)+V_now(2:xdim-1, 1:ydim-2)+V_now(2:xdim-1, 3:ydim))*0.25;
    V_now(n, w1:w2) = 10;
    V_now(n, 2:ydim-1) = (V_prev(n, 3:ydim)+V_prev(n, 1:ydim-2)+k2*V_prev(n+1, 2:ydim-1)+k1*V_prev(n-1, 2:ydim-1))*0.25;
    V_now(n, w1:w2) = 10;
    
    error = max(max(abs(V_now-V_prev)));
    V_prev = V_now;
%   imagesc(V_now);
%   getframe();
end

imagesc(V_now);
%vector plot plotted using quiver after
%normallizing the vector values
figure;
title('Vector plot of electric field(normallised) at resolution 30 points');
[ex, ey] = gradient(V_now(1:xdim, 1:ydim), 1, 1);
sz = max(max(sqrt(ex.^2 + ey.^2)));
ex = -ex;
ey = -ey;
exx = ex./sz;
eyy = ey./sz;
quiver(exx(1:5:xdim, 1:5:ydim), eyy(1:5:xdim, 1:5:ydim), 'k-');

q = (ex(n, w1) - ex(n, w2))/k1 + sum(ey(n+1, w1:w2)*epsr - ey(n-1, w1:w2));
C = (1/n)*abs(q)/10;

end