

function [out] = linear_convec_func(dt)
%Definig a function
L = 1;      % Length
n = 81;     % Mesh points
c = 1;      % Constant
t = 0.4;     % Time at which you want to see the plot
% Initial calculations

x = linspace(0,L,n);
dx = x(2)-x(1);

% Step function

u = ones(1,n);                  


x_start = 0.1 ; 
% The ending point of  fall
x_end = 0.3 ;
% The node at which  rises
n_start = (x_start/dx) + 1;
 % The node at which  falls
n_end = (x_end/dx) + 1;        
u(n_start:n_end) = 2;          

% Initial velocity profile

u_old = u;
u_initial = u;
 % Number of time steps
nt = 50;                      
coeff = c*dt/dx;
if t < nt*dt
    nt = t/dt;
end
    
% Time integration
for k = 1:nt
    % Spacial integration
    for i =2:n

        u(i) = u_old(i) - coeff*(u_old(i) - u_old(i-1));
    end
        u_old = u;

end

plot(x,u_initial,'r')
hold on
plot(x,u,'b')
grid on
legend('Initial profile','Final profile')
xlabel('Distance')
ylabel('Velocity magnitude')
title(['Velocity profile for time step = ' num2str(dt)])

end
