
 %Defining a function with time step as argument 
clear all;
close all;
clc
iterative_solver = 1; %Enter 1 for Jacobi, 2 for Gauss-Seidel and 3 for SOR
explicit_or_impicit = 2; %Enter 1 for explicit, 2 for implicit
SOR = 1.5;
Lx=1;
Ly=Lx;
nx = 10;
ny = nx;

x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
dx = x(2) - x(1);
dy = y(2) - y(1);

alpha = 1e-4;
%timestep
dt = 0.1;
tend = 50;
tnow = 0;
% Initialize the temperature field

T = ones(nx,ny);

% Boundary Conditions
% Temperature at the top face
T(1,:) = 900;
 % Temperature at the bottom face
T(end,:) = 600;    
% Temperature at the left face
T(:,1) = 400;  
% Temperature at the right face
T(:,end) = 800;     
% If the value is 1, then the solver used is Jacobi. If it is 2, then Gauss-seidel method and if 3, then Successive over-relaxation method

%coping T values
Told=T;
T_prev_dt=T; 

error = 1;
toler = 1e-4;
k = 2*(dx^2 + dy^2)/(dx^2*dy^2);


%Transient state solver
    k1 = alpha*dt/dx^2;
    k2 = alpha*dt/dy^2;
  
    %Explicit
    if explicit_or_impicit == 1
          while(tnow < tend)
              while(error > toler)     
                  for i = 2:nx-1
                      for j = 2:ny-1
                          T(i,j) = T_prev_dt(i,j) + k1*(T(i-1,j)-2*Told(i,j)+Told(i+1,j)) + k2*(T(i,j-1)-2*Told(i,j)+Told(i,j+1));
                      endfor
                  endfor
                  explicit_iter = explicit_iter + 1;
                  error = max(max((abs(Told-T))));
         explicit_iter = 1;
             Told = SOR*T + (1-SOR)*Told;
              end     
              error = 10;
              T_prev_dt = T;
              tnow = tnow + dt;              
                  end
              end
          end
    end
          
 
  
    %Implicit
    implicit_iter = 1;
    if explicit_or_impicit == 2
          while(tnow < tend)
              while(error > toler)     
                  for i = 2:nx-1
                      for j = 2:ny-1
                          T(i,j) = (T_prev_dt(i,j) + k1*(T(i-1,j)+Told(i+1,j)) + k2*(T(i,j-1)+Told(i,j+1)))/(1+2*k1+2*k2);
                  
                     
                  implicit_iter = implicit_iter + 1;
                  error = max(max((abs(Told-T))));
                  Told = SOR*T + (1-SOR)*Told;
              end     
              error = 10;
              T_prev_dt = T;
              tnow = tnow + dt;              
                  end
              end
          end
    endif
 
  
    end

contourf(x,y,T,'Linecolor\','none\')
colorbar 
 title(sprintf('Temprature distribution on a %d x %d  grid\n using Jacobi Method at iteration no %d', Lx, Ly, k),'FontSize',12);

 