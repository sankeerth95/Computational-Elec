clear all;
clc;
xdim=500;
time_tot=200;
xsource=250;

S=1;

epsilon0=1;
mu0=1;
c=1;
delta=1;
deltat=S*delta/c;
 
Ez=zeros(1,xdim);
Hy=zeros(1,xdim);
epsilon=epsilon0*ones(1,xdim);
epsilon(1,1:(xdim/2))=1.0;
epsilon(1,(xdim/2):xdim)=1.0;
mu=mu0*ones(1,xdim);

% gaussian=0;
sine=1;
impulse=0;
lambda = 10;
freq = c/lambda;

for n=1:1:time_tot;
    
    if sine == 1
    Ez(xsource)=cos(2*pi*freq*n);
    end
    if impulse == 1
        if n == 1
            Ez(xsource) = 1;
        else 
            Ez(xsource) = 0;
        end
    end
        
    Hy(1:xdim-1)=Hy(1:xdim-1)+(deltat./(delta*mu(1:xdim-1))).*(Ez(2:xdim)-Ez(1:xdim-1));
       
    Ez(2:xdim)=Ez(2:xdim)+(deltat./(deltat*epsilon(2:xdim))).*(Hy(2:xdim)-Hy(1:xdim-1));
    
    %Ez(xsource)=1;
        
        plot(1:1:xdim,Ez,'color','k','linewidth',2);
        grid on;
        axis([0, xdim, -2, 2]);
        getframe;
end

%critical dimention : smallest critical fearure the grid should capture..
%critical timestep : minimum timestep to choose such that S is nowhere
%greater than one

%imlement mer's absorbing boundary condition