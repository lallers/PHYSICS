%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Physics, 2016                                      /
%     University of New Mexico                           /
%/////////////////////////////////////////////////////////

%% CAR BRIDGE
clc
clear all
%Carstorm
%Distance of River
HeightHigh = 20.9;%Bridge Height Upper
HeightSmall = 1.9; %Bridge Height Lower
River =63;

g = 9.8;
%v0 =@(L,H1,H2) L*sqrt(9.8/(2*(H1-H2)));
%v =@(L,H1,H2) sqrt((L^2*9.8)/(2*(H1-H2))+(2*9.8*(H1-H2)));
A = River*sqrt(g/(2*(HeightHigh-HeightSmall)));
B = sqrt((River^2*9.8)/(2*(HeightHigh-HeightSmall))+(2*9.8*(HeightHigh-HeightSmall)));
sprintf('For Car Side= %g, River Side= %g, River Length= %g\n---------------------------------------\nA) v0 = %f How fast should be travelling\nB) v = %f Speed of the car before landing',HeightHigh,HeightSmall,River,A,B)
%% Football Quarter Back
clc
clear all

yvel = [15.1]; %Velocity in Y
xvel = [18]; %Velocity in X

A =@(y) y/9.8;
B =@(y) y^2/(2*9.8);
C =@(y) 2*y/9.8;
D =@(t1,t2) t2/t1;
E =@(y,x) 2*x*y/9.8;

Aa = A(yvel);
Ba = B(yvel);
Ca = C(yvel);
Da = Ca/Aa;
Ea = E(yvel,xvel);
sprintf('For Y Velocity= %g, X Velocity= %g\n---------------------------------------\nA) t1 = %g How much time to highest\nB)h = %f how high is the point\nC) t2 = %g How much time after thrown.\nD) t2/t1 = %g How does this comepare\nE) x = %f how far travelled horizontally.',yvel,xvel,Aa,Ba,Ca,Da,Ea)  

%% SHELL
clc
clear all
v0 = [30]; %Initial Velocity
angle = [50]; %Initial Velocity

Aa = v0*cosd(angle);
Ab = v0*sind(angle);
B = v0*sind(angle)/9.8;
C = (v0*sind(angle))^2/(2*9.8); 
D = v0^2*sind(2*angle)/9.8;
E = 9.8;
F = v0*cosd(angle);
sprintf('For v0= %g, angle= %g\n---------------------------------------\nA) v0h,v0v = %g , %g Horizontal and vert components\nB) t = %g Seconds\nC) hmax = %g max height\nD) l = %g how far is firing point\nE)ah,av = 0 , -%g highest point vert and hor components\nF) vh,vv = %g , 0', v0,angle,Aa,Ab,B,C,D,E,F)
%% BASEBALL HITTER
clc
clear all

v0 = 28.7;
angle = 35.7;
height = 10.1;

g = 9.8;
vy = v0*sind(angle);
vx = v0*cosd(angle);
r1 = (1/2)*g;
p = [-r1 vy -height]; 
A_t = roots(p);
B_HorizontalComp = vx;
C_VerticalComp = vy - g*A_t(2);
D_HorizontalComp = vx;
E_VerticalComp = vy - g*A_t(1);
F_MagnitudeVel = v0;
G_Direction = angle;
A_t = [min(A_t), max(A_t)];
sprintf('For v0= %g, angle= %g, height= %g\n---------------------------------------\nA) t1,t2 = %f , %f two times.\nB) vx = %f Horizontal at earlier time. \nC) vy = %f Vertical at earlier time.\nD) vx = %f Hoeizontal at later time.\nE) vy = %f Vertical component at later time.\nF) |v| = %f Magnitude of velocity when returns to level.\nG) \theta = %f Direction.\n',v0,angle,height,A_t,B_HorizontalComp,C_VerticalComp,D_HorizontalComp,E_VerticalComp,F_MagnitudeVel,G_Direction)

%% Flare
clc
clear all
v0 = 103;
angle = 55.2;
g_moon = 1.67;


vy = v0*sind(angle);
vx = v0*cosd(angle);
g = 9.8;
gn = -(1/2)*g;
gn_moon = -(1/2)*g_moon;
A = vy^2./(2*g);
p = [gn vy 0];
pn = roots(p);
B = vx*pn(2);
C = v0^2*sind(angle).^2./(2*g_moon);
p = [gn_moon vy 0];
pn = roots(p);
D = vx*pn(2);
sprintf('For v0= %g, angle= %g, moon gravity= %g\n---------------------------------------\nA) hmax = %f Max Height\nB) R = %f Distance from firing point \nC) hmax_moon = %f Max height on moon\nD) R = %f Distance to landing point',v0,angle,g_moon,A,B,C,D)

%% Dog
clc
clear all

v0 = 8.5;
angle = 65;
height = 10;

vy = v0*sind(angle);
vx = v0*cosd(angle);
g = -9.8;
gn = (1/2)*g;
height = abs(height);
p = [gn vy height];
T = roots(p);
d = vx*T(1);

sprintf('For v0= %g, angle= %g, tree height= %g\n---------------------------------------\nA) v = %f How fast mus dog run\nB) d = %f How far will he run',v0,angle,height,vx,d)



