% Simply call 
%
%       >> double_pendulum_init
%
% to run the double pendulum simulation with the below parameters. This
% script calls double_pendulum.
%
%   ---------------------------------------------------------------------

phi1                = 0;
dtphi1              = 0;
phi2                = pi/2;
dtphi2              = 9.8;
g                   = 9.81; 
m1                  = 8; 
m2                  = 1; 
l1                  = 1; 
l2                  = 1;
duration            = 100;
fps                 = 10;
movie               = true;

clc; figure;

interval=[0, duration];
ivp=[phi1; dtphi1; phi2; dtphi2; g; m1; m2; l1; l2];

double_pendulum(ivp, duration, fps, movie);