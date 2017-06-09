%//////////////////////////////////////////////////////////////////////////
%Lee Allers, Exam 3, Forced Duffings Oscillator
%//////////////////////////////////////////////////////////////////////////
clc; close all; clear all
    
param.gamma =   1/4;         % Given Gamma
param.alpha =    -1;         % Given alpha
param.beta  =     1;         % Given beta
param.w     =     1;         % Given Omega
param.t0    = 0;             % Initial time
param.tf    = 400;           % Final time

%//////////////////////////////////////////////////////////////////////////
%      Sets the ODE
%//////////////////////////////////////////////////////////////////////////
Te    = 2*pi/param.w;                             %Cycle time
res   =  1000;                              %Resolution
tspan = param.t0:(Te/res):param.tf;                     %Converts t to a vector                                                                %Vectorizes initial conditions
op    = odeset('abstol',1e-9,'reltol',1e-9);%Tolerance limit so my computer doesn't explode.

%//////////////////////////////////////////////////////////////////////////
%      Sets vectors for given ODE conditions 
%//////////////////////////////////////////////////////////////////////////
F     = [0.325 0.365 0.400];  % Forcing parameters given
IC.x10   = 0.090;                % Intial conditions given
IC.x11   = .092;                 % Intial conditions given  
IC.x20   = 0;                    % Intial conditions given

t_out = zeros(length(tspan),length(F)); %Initializes the vectors for faster computation
x     = zeros(length(tspan),length(F)); %Initializes the vectors for faster computation
xdot  = zeros(length(tspan),length(F)); %Initializes the vectors for faster computation
n  = zeros(length(tspan),length(F)); %Initializes the vectors for faster computation
%//////////////////////////////////////////////////////////////////////////
%      Gets solutions for each F in the for loop and x(0) = 0.090 
%//////////////////////////////////////////////////////////////////////////
for      i = 1:length(F)
[t,y]      = ode45(@(t,x) duffingeq(t,x,param.gamma,param.alpha,param.beta,F(i),param.w),tspan,[IC.x10 IC.x20],op); %Computes the ODE with given parameters, a time span, initial conditions, and finally the tolerance settings
t_out(:,i) = t;
x(:,i)     = y(:,1);                                                %Converts solutions to ODE into seperate vectors for F = 0.325, F = 0.365 and F = .400
xdot(:,i)  = y(:,2);
%Converts solutions to ODE into seperate vectors for F = 0.325, F = 0.365 and F = .400
end

%//////////////////////////////////////////////////////////////////////////
%      Gets solutions for each F in the for loop and x(0) = 0.092 
%//////////////////////////////////////////////////////////////////////////
for      i = 4:3+length(F)
[t,y]      = ode45(@(t,x) duffingeq(t,x,param.gamma,param.alpha,param.beta,F(i-3),param.w),tspan,[IC.x11 IC.x20],op); %Computes the ODE with given parameters, a time span, initial conditions, and finally the tolerance settings
t_out(:,i) = t;
x(:,i)     = y(:,1);                                                %Converts solutions to ODE into seperate vectors for F = 0.325, F = 0.365 and F = .400
xdot(:,i)  = y(:,2);                                                %Converts solutions to ODE into seperate vectors for F = 0.325, F = 0.365 and F = .400
end

figure(99)
hold on
plot(x(:,1),xdot(:,1)); 
plot(x(50000:end,1),xdot(50000:end,1),'r','linewidth',5)
xlabel('$x(t)$','interpreter','latex');ylabel('$\dot{x}(t)$','interpreter','latex');
string = sprintf('Forced Duffing Oscillator(usual case)\nF =%.03f, \\gamma =%.03f, \\alpha =%.01f, \\beta =%.01f, \\omega =%.01f',F(1),param.gamma,param.alpha,param.beta,param.w);
title(string); grid
hold off

%//////////////////////////////////////////////////////////////////////////
%      Computes the poincare sections
%//////////////////////////////////////////////////////////////////////////
n1 = mod(t_out(:,1), 2*pi) == 0; n2 = mod(t_out(:,2), 2*pi) == 0;
loc1 = find(n1 > 0);loc2 = find(n2 > 0);
poin_x1 = x(loc1,1); poin_y1 = xdot(loc1,1);
poin_x2 = x(loc2,2); poin_y2 = xdot(loc2,2);

%//////////////////////////////////////////////////////////////////////////
%      Setting Up General Figures (Part a) 
%//////////////////////////////////////////////////////////////////////////
figure(1)
subplot(2,3,1)
plot(x(:,1),xdot(:,1)); 
xlabel('$x(t)$','interpreter','latex');ylabel('$\dot{x}(t)$','interpreter','latex');
string = sprintf('Forced Duffing Oscillator(usual case)\nF =%.03f, \\gamma =%.03f, \\alpha =%.01f, \\beta =%.01f, \\omega =%.01f',F(1),param.gamma,param.alpha,param.beta,param.w);
title(string); grid

subplot(2,3,2) 
plot(t_out(:,1),x(:,1));
xlabel('t');ylabel('x1');
title('Time Series')


ax1 = subplot(2,3,3); 
scatter(poin_x1,poin_y1,'filled')
title('Poincare Map')                    
xlabel('x1');ylabel('x2');grid; 


subplot(2,3,4)
plot(x(:,2),xdot(:,2));  
xlabel('$x(t)$','interpreter','latex');ylabel('$\dot{x}(t)$','interpreter','latex');
string = sprintf('Forced Duffing Oscillator (usual case)\nF =%.03f, \\gamma =%.03f, \\alpha =%.01f, \\beta =%.01f, \\omega =%.01f',F(2),param.gamma,param.alpha,param.beta,param.w);
title(string); grid

subplot(2,3,5), plot(t_out(:,2),x(:,2));
xlabel('t');ylabel('x1');
title('Time Series')

subplot(2,3,6)
scatter(poin_x2,poin_y2,'filled')
xlabel('x1');ylabel('x2');grid; 
title('Poincare Map')


%------------------------------------------------------------------------------------
%%


%//////////////////////////////////////////////////////////////////////////
%      Setting Up Specific Figures (Part b) 
%//////////////////////////////////////////////////////////////////////////

figure(2)
subplot(2,1,1)
plot(x(:,2),xdot(:,2));  %plot the variable x and y
xlabel('$x(t)$','interpreter','latex');ylabel('$\dot{x}(t)$','interpreter','latex');
string = sprintf('Forced Duffing Oscillator(usual case)\nF =%.03f, \\gamma =%.03f, \\alpha =%.01f, \\beta =%.01f, \\omega =%.01f',F(2),param.gamma,param.alpha,param.beta,param.w);
title(string); grid



subplot(2,1,2)
plot(x(:,5),xdot(:,5));  %plot the variable x and y
xlabel('$x(t)$','interpreter','latex');ylabel('$\dot{x}(t)$','interpreter','latex');
string = sprintf('Forced Duffing Oscillator (nearby case)\nF =%.03f, \\gamma =%.03f, \\alpha =%.01f, \\beta =%.01f, \\omega =%.01f',F(2),param.gamma,param.alpha,param.beta,param.w);
title(string); grid
%------------------------------------------------------------------------------------


%//////////////////////////////////////////////////////////////////////////
%     Calculations and Figures (Part c)
%//////////////////////////////////////////////////////////////////////////
param.tf   = 1000;
tspan      = param.t0:(Te/res):param.tf; 
[t,y]      = ode45(@(t,x) duffingeq(t,x,param.gamma,param.alpha,param.beta,F(3),param.w),tspan,[IC.x10 IC.x20],op); %Computes the ODE with given parameters, a time span, initial conditions, and finally the tolerance settings
t_outn = t;
xn     = y(:,1);
xdotn  = y(:,2);

n3 = mod(t_outn, 2*pi) == 0; 
loc3 = find(n3 > 0);
poin_x3 = xn(loc3); poin_y3 = xdotn(loc3);


[t,y]      = ode45(@(t,x) duffingeq(t,x,param.gamma,param.alpha,param.beta,F(3),param.w),tspan,[IC.x11 IC.x20],op); %Computes the ODE with given parameters, a time span, initial conditions, and finally the tolerance settings
t_outn2 = t;
xn2     = y(:,1);
xdotn2  = y(:,2);

n3 = mod(t_outn2, 2*pi) == 0; 
loc4 = find(n3 > 0);
poin_x4 = xn2(loc4); poin_y4 = xdotn2(loc4);


figure(3)
subplot(2,3,1)
plot(xn,xdotn);  %plot the variable x and y
xlabel('$x(t)$','interpreter','latex');ylabel('$\dot{x}(t)$','interpreter','latex');
string = sprintf('Forced Duffing Oscillator(usual case)\nF =%.03f, \\gamma =%.03f, \\alpha =%.01f, \\beta =%.01f, \\omega =%.01f',F(3),param.gamma,param.alpha,param.beta,param.w);
title(string); grid

subplot(2,3,2) 
plot(t_outn,xn);
xlabel('t');ylabel('x1');
title('Time Series')

subplot(2,3,3);
scatter(poin_x3,poin_y3,'filled')
title('Poincare Map')                    
xlabel('x1');ylabel('x2');grid;

subplot(2,3,4)
plot(xn2,xdotn2);  %plot the variable x and y
xlabel('$x(t)$','interpreter','latex');ylabel('$\dot{x}(t)$','interpreter','latex');
string = sprintf('Forced Duffing Oscillator(nearby case)\nF =%.03f, \\gamma =%.03f, \\alpha =%.01f, \\beta =%.01f, \\omega =%.01f',F(3),param.gamma,param.alpha,param.beta,param.w);
title(string); grid

subplot(2,3,5) 
plot(t_outn2,xn2);
xlabel('t');ylabel('x1');
title('Time Series')

subplot(2,3,6);
scatter(poin_x4,poin_y4,'filled')
title('Poincare Map')                    
xlabel('x1');ylabel('x2');grid;




%//////////////////////////////////////////////////////////////////////////
%      Changing figure settings
%//////////////////////////////////////////////////////////////////////////
fig1 = figure(1);
fig2 = figure(2);
fig3 = figure(3);
set(fig1,'Units','normalized',...
 'Position',[0 0 .8 .8])
set(fig2,'Units','normalized',...
 'Position',[0 0 .6 .8])
set(fig3,'Units','normalized',...
 'Position',[0 0 .8 .8])


%//////////////////////////////////////////////////////////////////////////
%      Lyaponov Coefficients part(b)
%//////////////////////////////////////////////////////////////////////////

figure('Units','normalized','Position',[0 0 .6 .8])
subplot(2,2,1)
plot(t_out(:,1),log2(abs(x(:,4)-x(:,1))));
line([88 160],[-6 -31],'Color','red','LineStyle','--')
xlabel('t');ylabel('log|\Delta\phi|')
legend('','\lambda \approx -.15')
title('Approximation of Lyapunov Coefficient')
subplot(2,2,2)
plot(t_out(:,1),log2(abs(x(:,4)-x(:,1))));
line([94 160],[-6 -31],'Color','red','LineStyle','--')
xlabel('t');ylabel('log|\Delta\phi|')
legend('','\lambda \approx -.15')
title('Zoomed into region of interest of first graph')
xlim([80 170])

subplot(2,2,3)
plot(t_outn2,log(abs(xn2-xn)));
%line([88 160],[-6 -31],'Color','red','LineStyle','--')
xlabel('t');ylabel('log|\Delta\phi|')
%legend('','\lambda \approx -.4')
title('Approximation of Lyapunov Coefficient')
legend('\lambda \approx 1')

%//////////////////////////////////////////////////////////////////////////
%      Bonus Question
%//////////////////////////////////////////////////////////////////////////
%My attempt is to iterate through reducing the values until I find a value
%for labmda < 1
param.tf   = 300;
tspan      = param.t0:(Te/res):param.tf; 
u = 0;
iter = 0;
chk = 1;
while chk >= 0
 u = u + .001;
 iter = iter + 1;
 Fnew = F(3) - u;
[t,  y]      = ode45(@(t,x) duffingeq(t,x,param.gamma,param.alpha,param.beta,Fnew,param.w),tspan,[IC.x10 IC.x20],op); %Computes the ODE with given parameters, a time span, initial conditions, and finally the tolerance settings
[t ,y2]      = ode45(@(t,x) duffingeq(t,x,param.gamma,param.alpha,param.beta,Fnew,param.w),tspan,[IC.x11 IC.x20],op); %Computes the ODE with given parameters, a time span, initial conditions, and finally the tolerance settings

x_1  = y(:,1);x_2  = y2(:,1);
G = log(abs(x_2(100:end) - x_1(100:end)));
lambda = [max(G), min(G)];
if lambda(1) && lambda(2) < 0
  chk = -1;
  t2 = t;
  iteration = iter ;
else
chk = 1;
end

end
%%
%This returns .3990
param.tf   = 400;
Fn = .3650;
tspan      = param.t0:(Te/res):param.tf; 
[t,  y]    = ode45(@(t,x) duffingeq(t,x,param.gamma,param.alpha,param.beta,Fn,param.w),tspan,[IC.x10 IC.x20],op); %Computes the ODE with given parameters, a time span, initial conditions, and finally the tolerance settings
[t ,y2]    = ode45(@(t,x) duffingeq(t,x,param.gamma,param.alpha,param.beta,Fn,param.w),tspan,[IC.x11 IC.x20],op); %Computes the ODE with given parameters, a time span, initial conditions, and finally the tolerance settings

figure(5)
plot(t,log2(abs(y2(:,1)-y(:,1))));
%line([88 160],[-6 -31],'Color','red','LineStyle','--')
xlabel('t');ylabel('log|\Delta\phi|')
%legend('','\lambda \approx -.4')
title('Approximation of Lyapunov Coefficient')
legend('\lambda \approx 1')