function [dy] = duffingeq(t,y,gamma,alpha,beta,F,w)
x1 = y(1);    x2 = y(2);
dx1=x2;
dx2=-2*gamma*x2-alpha*x1-beta*x1^3+F*cos(w*t);

dy = [dx1; dx2];