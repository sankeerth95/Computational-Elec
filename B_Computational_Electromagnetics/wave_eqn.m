dx = 1;
dt = 1;
c = 1;
L = 100;
stopTime = 150;
r = c*dt/dx;
n = L/dx + 1;
current = zeros(1, n);
past = current;
future = current;

for t = 0:dt:stopTime
    if (t == 0) current(1, 50) = 1;
    else current(1, 50) = 0;
    end
    future(1) = 0;
    future(2:n-1) = r^2*(current(1:n-2) + current(3:n)) + 2*(1-r^2)*current(2:n-1)-past(2:n-1);
    future(n) = 0;
    past = current;
    current = future;
    plot([0:dx:L], current);
    axis([0, L, -2, 2]);
    getframe();
    pause(0.0001);
end