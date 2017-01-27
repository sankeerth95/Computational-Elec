%initialize some variables here
xdim = 25;
ydim = 50;


V_now = zeros(xdim, ydim);
V_prev = zeros(xdim, ydim);
V_now(1, 1:ydim) = 10;
iter = 0;
error = max(max(abs(V_now - V_prev)));

%iterate as long as the change in values obtained
%in negligible. In this case, max value change is 0.001
while(error > 0.001)
iter = iter+1;
for i = 2:1:xdim-1
    for j = 1:1:ydim
        if j == 1
            V_now(i, j) = (V_now(i-1, j) + V_now(i+1, j) + 2*V_now(i, j+1))*0.25;
        elseif j == ydim
             V_now(i, j) = (V_now(i-1, j) + V_now(i+1, j) + 2*V_now(i, j-1))*0.25;
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
[ex, ey] = gradient(V_now(1:5:xdim, 1:5:ydim), 1, 1);
sz = max(max(sqrt(ex.^2 + ey.^2)));
ex = -ex;
ey = -ey;
exx = ex./sz;
eyy = ey./sz;
quiver(exx, eyy, 'k-');