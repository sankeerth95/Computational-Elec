%initialize some variables here
clear all;
xdim = 16;
ydim = 16;

V_now = 5*ones(4*xdim, 4*ydim);
V_prev = 5*ones(4*xdim, 4*ydim);

%iterate as long as the change in values obtained
%in negligible. In this case, max value change is 0.001

capx1 = 3*xdim/2;
capx2 = 5*xdim/2;
capy1 = ydim;
capy2 = 3*ydim;

V_now(capx1, capy1:capy2) = 10;
V_now(capx2, capy1:capy2) = 0;
error = max(max(abs(V_now - V_prev)));
iter = 0;

while(error > 0.0004)
iter = iter+1;
for i = 1:1:4*xdim
    for j = 1:1:4*ydim
        if (j == 1 || j == 4*ydim) && (i == 1 || i == 4*xdim)
            V_now(i, j) = (V_now(abs(i-2)+1, j) + V_now(i, abs(j-2)+1))*0.5;
        elseif j == 1 || j == 4*ydim
            V_now(i, j) = (V_now(i-1, j) + V_now(i+1, j) + 2*V_now(i, abs(j-2)+1))*0.25;
        elseif i == 1 || i == 4*xdim
            V_now(i, j) = (2*V_now(abs(i-2)+1, j) + V_now(i, j-1) + V_now(i, j+1))*0.25;
        elseif (i == capx1 || i == capx2) && j <=  capy2 && j >= capy1
            V_now(i, j) = 10 - 5*((i/capx1) - 1);
        else
             V_now(i, j) = (V_now(i-1, j) + V_now(i+1, j) + V_now(i, j-1) + V_now(i, j+1))*0.25;
        end
    end
end


error = max(max(abs(V_now-V_prev)));
V_prev = V_now;

end

imagesc(V_now);

%vector plot plotted using quiver after
%normallizing the vector values
figure;
title('Vector plot of electric field(normallised) at resolution 30 points');
[ex, ey] = gradient(V_now(1:5:4*xdim, 1:5:4*ydim), 1, 1);
sz = sqrt(ex.^2 + ey.^2);
ex = -ex;
ey = -ey;
exx = ex./sz;
eyy = ey./sz;
quiver(exx, eyy, 'k-');