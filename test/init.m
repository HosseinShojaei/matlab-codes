clc
clear
close all

L   = 0.3;                       %Pendulum length
l   = 0.2175;                    %Length of center of mass of Pendulum
r   = 0.00875;                   %Radius of pulley

Ra  = 10.5;
kt  = 0.0436;
kb  = 0.0466;

c   = 3.5;
b   = 0.00008;

Rp  = 0.004;
rho = 7850;

m1  = rho*pi*Rp^2*L;
m2  = rho*pi*4/3*0.0125^3;
m   = m1+m2;
I1G = m1*L^2/12;
I2G = 2*m2*0.0125^2/5;
I1  = I1G+m1*0.0224^2;
I2  = I2G+m2*0.1376^2;
I   = I1+I2;
Mb  = 0.02442;                   %mass of belt
Mc  = 0.2051;                    %mass of cart
Me  = 0.0553;                    %mass of encoder
Mp  = 0.0131;                    %equivalent mass of pulleys
Mo  = 0.01;                      %other masses
M   = Mb + Mc + Me + Mp + Mo
g   = 9.80665;
Ed  = 2*m*g*l;

q   = I*(M+m) + M*m*(l^2);

A   = [0    0                   1                                            0
       0    0                   0                                            1
       0    (((m*l)^2)*g)/q     -((I +m*(l^2))/q)*(c + (kb*kt)/(Ra*(r^2)))   -(b*m*l)/q
       0    (m*g*l*(M+m))/q     -((m*l)/q)*(c + (kb*kt)/(Ra*(r^2)))          -((M+m)*b)/q];

B   = [0
       0
       ((I +m*(l^2))*kt)/(q*Ra*r)
       (m*l*kt)/(q*Ra*r)];

Q   = diag([1000 1000 0 0]);
R   = 0.1;   %0.05
K   = lqr(A,B,Q,R);

disp('LQR Gain ='); disp(K);
