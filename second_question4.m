function [out]=second_question4(n, ax)


%input for the linear convection problem
L=1;
nt=50;
c=1;
dt=0.01;
%time start
t=0;
 % calculation
 x=linspace(0,L,n);
 dx=x(2)-x(1);
 
 %define start and end position of the square wave
 n_start=0.1/dx+1;
 n_end=0.3/dx+1;
 
 
  %initilize arrays and apply boundry conditions
  u=ones(1,n);
   u(n_start:n_end) = 2
   u_old=u;
  u_initial=u;
  
   
 
   %time marching solution
   for k=1:nt
       for i=2:n
       u(i)= u_old(i)-(c*dt/dx)*(u_old(i)-u_old(i-1));  
      
       end
    %updating velocity
      u_old=u;
      t=t+dt;
      %ploting initial velocity graph
      plot(x,u_initial,'-o');
      hold on;
      % plotting the velocity profile
      
      plot(x,u,'r','linewidth',2);
     axis([0 2 1 2.5])
     xlabel('1D','fontsize',15);
     ylabel('velocity,m/s','fontsize',15)
      title_text=sprintf('Time= %d',t);
      title(title_text);
      pause(0.3);
   
     
     
   end

