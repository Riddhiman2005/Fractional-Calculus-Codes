% The Fractional Forward Euler's Method for a System of ODES
%This is based on the paper 'Analytic study on linear systems of fractional differential equations' by Zaid M. Odibat.
clc; clear; close all;
% Inputs 
h=1/10;t(1)=0; x(1)=1.2;y(1)=4.2; alpha=0.85; tfinal=1; t=t(1):h:tfinal; N=ceil((tfinal-t(1))/h);
% Exact Solution
xExact=(1/5).*mlf(alpha, 1, (t.^alpha))+mlf(alpha, 1, (-2*t.^alpha));
yExact=(1/5).*mlf(alpha, 1, (t.^alpha))+4*mlf(alpha,1,(-2*t.^alpha));
%Fractional-Order ODE
f1 =@(t,x,y) 2*x-y; f2 =@(t,x,y) 4*x-3*y;
%Fractional Forward Euler Method
for n = 1:N
j=1:n;
t(n+1)=t(n)+h;
t(n + 1) = t(n) + h;
x(n+1)=x(1)+((h^alpha)/(gamma(alpha+1))).*sum(((n-j+1).^(alpha)-(n-j).^(alpha)).*f1(t(j),x(j),y(j)));
y(n+1)=y(1)+((h^alpha)/(gamma(alpha+1))).*sum(((n-j+1).^(alpha)-(n-j).^(alpha)).*f2(t(j), x(j),y(j)));
end
% Absolute Errors
xErrors=abs(xExact-x); yErrors=abs(yExact-y);
xLast_Error=xErrors(end),yLast_Error=yErrors(end),
xMax_Error=max(xErrors), yMax_Error=max(yErrors),v
