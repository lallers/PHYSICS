%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Numerical Computation, 2016                        /
%     University of New Mexico                           /
%NOTE: None of my scripts are built to be robust, they   /
%      are merely an implementation of a given set of    /
%      data or instructions!                             /
%/////////////////////////////////////////////////////////

%% Exercise 8.27
%Equations to know
% P of Rebecca after collision = mass * velocity of mag * cos or sin
% Total momentum in y is [0 - PofRebecca in Y]
% Total momentum in x is [Rebecca's inialt P (m*vi)- PofRebecca in x]
% Velocity of Daniel in y is [Total Momentum in Y/mass of Daniel]
% Velocity of Daniel in x is [Total Momentum in x/mass of Daniel]
% Daniels Magnitude Velociy is [sqrt(vy^+vx^2]
% Direction is arctan(vy/vx)
% Ki = 1/2 m_R v_R^2
% Kr after collision is 1/2 m_R v_RMAG^2 

clc
clear all

m1 =      65; %Daniel
m2 =      45; %Rebecca
vR =      13; %Rebecca into Daniel
vmag =     6; % V mag after Collision
angle = 52.1; %degrees

vDi =      0; %Daniel Inital Velocity
pox = m2*vR;
poy = 0;

prx = m2*vmag*cosd(angle);
pry = m2*vmag*sind(angle);
Tpoy = poy - pry;
Tpox = pox - prx;

voy = Tpoy/m1;
vox = Tpox/m1;

vD = sqrt(voy^2+vox^2);
theta = atand(voy/vox);
ki = (1/2)*m2*vR^2;
kr = (1/2)*m2*vmag^2;
kd = (1/2)*m1*vD^2;
kf = kr+kd-ki;

string = sprintf('Exercise 8.27 -------------------');
disp(string)
string = sprintf('A) %g m/s',vD);
disp(string)
string = sprintf('B) %g Degrees',abs(theta));
disp(string)
string = sprintf('C) %g Joules',kf);
disp(string)

%% Exercise 8.49 - Neutron
%Equations to know
%Collisions = log(Speed reduction)/log(3)
clear all
clc

m1 = 2.0; %neutrons and deutrons of mass
a = 0.33;
b = 1/9;
c = [1/1.77e5];

u = 1;
v = 1;
U = m1*(u)-u;
V = -1*(m1*(v)+v); 
V2 = abs(U/V);
KE = (1/2)*m1*V2^2;
ca = abs(log(c)/log(abs(V)));


string = sprintf('Exercise 8.49 -------------------');
disp(string)
string = sprintf('A) %g',V2);
disp(string)
string = sprintf('B) %g',KE);
disp(string)
if isempty(c)
string = sprintf('C) Make sure to put in part C!\n %g',ca);
disp(string)
else
  string = sprintf('C) %g collisions',ca);
disp(string)  
end

%% Exercise 8.40 Falcons
%Equations to know 
%Velocity of the Raven is (Mfalcon*Vfalcon - Mfalcon*((-1)*v of the bounce backFalcon))/MassRaven
%Direction is arctan(Vravenf/VravenI)

clear all
clc
m1 = [55]      ; %kg falcon
v1 = [19]      ; %Falson flying at
m2 = [1.5] *100     ; %kg raven
v2 = [9]       ; %Raven flyng at
vAc = [5]      ; %Bounced back at



v_raven = (m1*v1 + m1*vAc)/m2;
direction = atand(v_raven/v2);
speed = round(sqrt(v2^2 + v_raven^2));

string = sprintf('Exercise 8.40 -------------------');
disp(string)
string = sprintf('A) %g Degrees',direction);
disp(string)
string = sprintf('B) Make sure to check rounding. 2 sig figs\n    %g m/s', speed);
disp(string)

%% Exerceise 8.38
%Equations to know
% MassB*v / (Total system mass *sin(angle))
% Speed of car A is Total system mass*cos(angle) / MassA

clear all
clc

mA = 1900  ; %kg
mB = 1400  ; %kg
v =  16    ;%m/s
angle = 65 ;%degrees

tot = mA+mB;
V2 = mB*v/(tot*sind(angle));
vA = (tot*V2*cosd(angle))/mA;
string = sprintf('Exercise 8.38 -------------------');
disp(string)
string = sprintf('A) %g m/s',V2);
disp(string)
string = sprintf('B) %g m/s', vA);
disp(string)

%% Alternative Exercise 8.132
clear all
clc

string = sprintf('Alternative Exercise 8.132 -------------------');
disp(string)
string = sprintf('A) K_A, K_B = m_B*Q/(m_A+m_B),m_A*Q/(m_A+m_B) ');
disp(string)
string = sprintf('B) KA/Q = 20, KB/Q 80 Percent');
disp(string)

%% Exercise 10.41
%Equations ti know
%Convert Days to seconds
%w2 = ri^2*w/rf^2
clear all
clc

density = 10^14;
ri =        7e5   ;%km
rf =         18    ;%km
T  =         30    ; %Days

conv = T*24*3600; %to seconds
w = 2*pi/conv;
w2 = ri^2*w/rf^2;

string = sprintf('Exercise 10.41 -------------------');
disp(string)
string = sprintf('A) w2 = %g',w2);
disp(string)

%% Exercise 10.46
%equation to know: (2/5)M + m  = 1.20 * (2/5)M
clear all
clc
percentage = 20 ;
A = (100+percentage)/100;
AA = A*(2/5)-(2/5);
string = sprintf('Exercise 10.46 -------------------');
disp(string)
string = sprintf('A) %g M',AA);
disp(string)
