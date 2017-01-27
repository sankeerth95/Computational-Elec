%initialize some variables here
n = 20;
time = [];
iterations = [];

%for different values of x and y dimention resolution, we're
%going to solve laplace equation
for(n = 20:5:100)

xdim = n;
ydim = n;

%time count begins here
tic;

V_now = zeros(xdim, ydim);
V_prev = zeros(xdim, ydim);
V_now(1:xdim, ydim) = 10;
iter = 0;
error = max(max(abs(V_now - V_prev)));

%iterate as long as the change in values obtained
%in negligible. In this case, max value change is 0.001
while(error > 0.001)
iter = iter+1;

V_now(2:xdim-1, 2:ydim-1) = (V_now(1:xdim-2, 2:ydim-1) + V_now(3:xdim, 2:ydim-1) + V_now(2:xdim-1, 1:ydim-2) + V_now(2:xdim-1, 3:ydim))*0.25;

error = max(max(abs(V_now-V_prev)));
V_prev = V_now;
%imagesc(V_now);
%pause(0.0001)
end

%update recorded time to solve and 
% number of iterations in an array
time = [time, toc];
iterations = [iterations, iter];

%for a convenient quiver plot, store the potential
%values at a convenient resolution which, in this case, 30
if(n == 30)
    V_quiv = V_now;
end

end

%plot convergence time vs no. of points/resolution
figure;
plot(20:5:100, time, 'r+');
title('Time for convergence versus number of resolution points(no for loops)');
xlabel('Resolution');
ylabel('Time of execution');

%plot number of iterations vs resulution
figure;
plot(20:5:100, iterations, 'r+');
title('Time for convergence versus number of iterations(no for loops)');
xlabel('Resolution');
ylabel('Number of iterations');

%visualization ofthe potential field
figure;
imagesc(V_now);
title('Potential field at resolution 100');

%vector plot plotted using quiver after
%normallizing the vector values
figure;
title('Vector plot of electric field(normallised) at resolution 30 points');
[ex, ey] = gradient(V_quiv, 1, 1);
sz = sqrt(ex.^2 + ey.^2);
ex = -ex;
ey = -ey;
exx = ex./sz;
eyy = ey./sz;
quiver(exx, eyy, 'k-');