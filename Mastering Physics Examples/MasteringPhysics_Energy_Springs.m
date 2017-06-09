%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Physics, 2016                                      /
%     University of New Mexico                           /
%/////////////////////////////////////////////////////////
%% 637 Graph


mass = 8;
x1 = 1;
x2 = 12;

x0 = 0;
c = 8;
xf = 12;
y0 = 0;
yf = 10;

s1 = (yf-y0)/(c-x0);
a1 = (1/2)*s1*x1;
A = sqrt(a1*2/mass);

W = (1/2)*c*yf + (1/2)*yf*(xf-c);
B = sqrt(W*2/mass);

disp('Example 6.37----------------------------')
string = sprintf('A) v = %g m/s',A);
disp(string)
string = sprintf('B) v = %g m/s',B);
disp(string)





%% Exercise 6.47 Glider **Works**
%Part A


angle =    ; %Degrees
mass  =    ; %kg
k     =    ; %N/m
d     =    ; %meters
g     =    ;
d2    =    ;%Part B, still in contact in meters

x = sqrt((mass * g * d * sind(angle) * 2)/k);
U2 = mass * g * d2 * sind(angle);
U1 = (1/2) * k * x^2;
K = U1-U2;

if x < d2 
    B = 'No';
elseif x > d2
    B = 'Yes';
elseif x == d2 
    B = 'They are still touching';
end


disp('Example 6.47\n----------------------------')
string = sprintf('A) x = %g m',x);
disp(string)
string = sprintf('B) %s',B);
disp(string)
string = sprintf('C) K = %g J',K);
disp(string)

%% Exercise 6.116 Block of Ice **Working**


mass = 5.7;

alpha = 0.2;
beta = 2.01e-2;
t1 = 3.8; %Part A
t2 = 3.8; %Part B
t3 = 3.8; %Part C


syms t
x = alpha*t^2 + beta*t^3;
v = matlabFunction(diff(x));
a = matlabFunction(diff(x,2));
W = matlabFunction(diff(x) * diff(x,2));
W = mass * integral(W,0,t3);
A = v(t1);
B = mass*a(t2);
C = W;
disp('Exercise 6.116, Block of Ice')
disp('-------------------------------------------')
fprintf('A) v = %g m/s\n',A)
fprintf('B) F = %g N\n',B)
fprintf('C) W = %g J\n',C)

%% Exercise 7.21 A spring of negligible mass ** working **


k =  1400; %N/m
E1 = 3.2 ; %Part A potential
mass = 1.5; %Part B, book mass
h = .500 ; %Part B, book height meters


g = 9.8;
Q1 =@(a,b,c) (-b + sqrt(b^2 - 4 * a * c))/(2*a);  
Q2 =@(a,b,c) (-b - sqrt(b^2 - 4 * a * c))/(2*a); 
U = sqrt(E1 * 2 / k);
a = 1; b = -2*mass*g/k; c = -2*mass*g*h/k;
roots = [Q1(a,b,c),Q2(a,b,c)];

disp('Exercise 7.21 Spring Negligible mass')
disp('-------------------------------------------')
fprintf('A) %g m\n',U)
fprintf('B) %g m\n   or %g m\n',roots)


%% Exercise 7.30 Roofer ** Working **


angle =      ; % Angle
F     =      ; % N Nudge
start =      ; % Meters Part A question
friction =   ; % N
g = 9.8;

m = F/g;
a = (F*sind(angle) - friction)/m; %m/s^2, Need m/s
v = sqrt(2*a*start);
disp('Exercise 7.30 Roofer')
disp('-------------------------------------------')
fprintf('A) %g m/s\n',v)

%% Exercise 7.33 A small Block


mass = .04 ; %kg
x1 =   0.34; %Part A, meters
y1 =   0.64; %Part B, meters


syms x y
U = 5.8 * x^2 - 3.75 * y^3;
Ux = matlabFunction(diff(U,x));
Uy = matlabFunction(diff(U,y));
ux = Ux(x1);
uy = Uy(y1);
F = (sqrt(ux^2 + uy^2));
a = F/mass;
theta = 180 + atand(uy/ux);

disp('Exercise 7.33 Small Block')
disp('-------------------------------------------')
fprintf('A) a = %g m/s^2\n',a)
fprintf('B) theta = %g degrees\n', theta)

%% Exercise 7.37 Construction Site **Working **

mass1  =    ; %kg, Bucket of concrete
mass2 =     ; %kg, Box
us =        ; %Friction Coefficient
uk =        ; %Friction Coefficient
d =         ; %Part C

g = 9.8;
T = mass1*g;
F = us*(mass2-mass1)*g;
F2 = uk*mass2*g;
dU = -mass1*g*d;
W = -F2*d;
v = sqrt( 2*(-dU + W)/(mass2 + mass1));
A = 0;
B = T;
C = v;

disp('Exercise 7.37 Construction Site')
disp('-------------------------------------------')
fprintf('A) F = %g N\n',A)
fprintf('B) F = %g N\n', B)
fprintf('C) v = %g m/s\n',C)




