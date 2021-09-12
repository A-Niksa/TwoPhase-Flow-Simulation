% Effects of roughness on rarefied gas flow in a microfluidic channel with porous walls 
% Amirreza Soheili, Arsha Niksa, Vahid Bazargan

% Note that for the code to work properly, the directory of this file (the
% folder which contains this MATLAB file) should be opened as the Current
% Folder. This would allow MATLAB to export the figures.

clear all
close all
clc

%% System Variables
% boundary coefficients
a = 0.49;
b = 1.28;
c = -1.0669;
d = -0.003;

% semi-arbitary variables
h = 143/2*1e-6; % m | microchannel half-depth (arbitary)
L = 143/2*1e-5; % m | microchannel length
ro = 1.165; % kg/m^3 | gas density
mu = 1.76e-5; % Pa*s | absolute viscosity
Pout = 0.002; % Pa | outlet pressure
Pin = 0.003; % Pa | inlet pressure
lan = 5.9e-5; % m | mean free path at 100 Pa (landa)
Dout = 2*h*sqrt(2)*pi/lan; % outlet inverse Knudson number 
Din = Dout*Pin/Pout; % inlet inverse Knudson number

% coordinates
n = 50; % number of samples
x = linspace(0,L,n); % m | x axis in Cartesian coordinates

% porous thickness
delt1 = 11/3*10^(-6); % m | porous media thickness (delta) | mode 1
delt2 = 143/19*10^(-6); % m | porous media thickness (delta) | mode 2
delt3 = 429/17*10^(-6); % m | porous media thickness (delta) | mode 3
delt4 = 143/3*10^(-6); % m | porous media thickness (delta) | mode 4

% permeability
k1 = 0; % permeability | mode 1
k2 = 0.00000001; % m^2 | permeability | mode 2
k3 = 0.0000001; % m^2 | permeability | mode 3
k4 = 0.000001; % m^2 | permeability | mode 4
k5 = 0.00001; % m^2 | permeability | mode 5

% inverse Knudson number
D1 = 0.01; % general inverse Knudson number | mode 1
D2 = 0.1; % general inverse Knudson number | mode 2
D3 = 1; % general inverse Knudson number | mode 3
D4 = 10; % general inverse Knudson number | mode 4

% y axis
y1 = linspace(-h-delt1,h+delt1,n); % m | y axis in Cartesian coordinates | mode 1
y2 = linspace(-h-delt2,h+delt2,n); % m | y axis in Cartesian coordinates | mode 2
y3 = linspace(-h-delt3,h+delt3,n); % m | y axis in Cartesian coordinates | mode 3
y4 = linspace(-h-delt4,h+delt4,n); % m | y axis in Cartesian coordinates | mode 4

%% Setting the Pressure Distribution
% P & dP (D1_delt4_k1) / 1
P1 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p1
    C5 = ro*p1/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p1/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p1/Pout)^d));
    C6 = ro*p1/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k1-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p1/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D1^c*h^2)/2-k1-a*lan*D1^d*h) == C5*x(i)+C6;
    P1(i) = vpasolve(eqn,p1);
end
P1 = -P1+max(P1);
dP1 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP1(1) = 0;
for i = 2:n
    dP1(i) = (P1(i)-P1(i-1))/(x(i)-x(i-1));
end

% P & dP (D1_delt4_k2) / 2
P2 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p2
    C5 = ro*p2/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p2/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p2/Pout)^d));
    C6 = ro*p2/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k2-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p2/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D1^c*h^2)/2-k2-a*lan*D1^d*h) == C5*x(i)+C6;
    P2(i) = vpasolve(eqn,p2);
end
P2 = -P2+max(P2);
dP2 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP2(1) = 0;
for i = 2:n
    dP2(i) = (P2(i)-P2(i-1))/(x(i)-x(i-1));
end

% P & dP (D1_delt4_k3) / 3
P3 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p3
    C5 = ro*p3/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p3/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p3/Pout)^d));
    C6 = ro*p3/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k3-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p3/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D1^c*h^2)/2-k3-a*lan*D1^d*h) == C5*x(i)+C6;
    P3(i) = vpasolve(eqn,p3);
end
P3 = -P3+max(P3);
dP3 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP3(1) = 0;
for i = 2:n
    dP3(i) = (P3(i)-P3(i-1))/(x(i)-x(i-1));
end

% P & dP (D1_delt4_k4) / 4
P4 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p4
    C5 = ro*p4/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p4/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p4/Pout)^d));
    C6 = ro*p4/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k4-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p4/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D1^c*h^2)/2-k4-a*lan*D1^d*h) == C5*x(i)+C6;
    P4(i) = vpasolve(eqn,p4);
end
P4 = -P4+max(P4);
dP4 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP4(1) = 0;
for i = 2:n
    dP4(i) = (P4(i)-P4(i-1))/(x(i)-x(i-1));
end

% P & dP (D1_delt4_k5) / 5
P5 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p5
    C5 = ro*p5/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p5/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p5/Pout)^d));
    C6 = ro*p5/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k5-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p5/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D1^c*h^2)/2-k5-a*lan*D1^d*h) == C5*x(i)+C6;
    P5(i) = vpasolve(eqn,p5);
end
P5 = -P5+max(P5);
dP5 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP5(1) = 0;
for i = 2:n
    dP5(i) = (P5(i)-P5(i-1))/(x(i)-x(i-1));
end

% P & dP (D2_delt4_k1) / 6
P6 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p6
    C5 = ro*p6/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p6/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p6/Pout)^d));
    C6 = ro*p6/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k1-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p6/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D2^c*h^2)/2-k1-a*lan*D2^d*h) == C5*x(i)+C6;
    P6(i) = vpasolve(eqn,p6);
end
P6 = -P6+max(P6);
dP6 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP6(1) = 0;
for i = 2:n
    dP6(i) = (P6(i)-P6(i-1))/(x(i)-x(i-1));
end

% P & dP (D2_delt4_k2) / 7
P7 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p7
    C5 = ro*p7/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p7/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p7/Pout)^d));
    C6 = ro*p7/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k2-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p7/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D2^c*h^2)/2-k2-a*lan*D2^d*h) == C5*x(i)+C6;
    P7(i) = vpasolve(eqn,p7);
end
P7 = -P7+max(P7);
dP7 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP7(1) = 0;
for i = 2:n
    dP7(i) = (P7(i)-P7(i-1))/(x(i)-x(i-1));
end

% P & dP (D2_delt4_k3) / 8
P8 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p8
    C5 = ro*p8/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p8/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p8/Pout)^d));
    C6 = ro*p8/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k3-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p8/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D2^c*h^2)/2-k3-a*lan*D2^d*h) == C5*x(i)+C6;
    P8(i) = vpasolve(eqn,p8);
end
P8 = -P8+max(P8);
dP8 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP8(1) = 0;
for i = 2:n
    dP8(i) = (P8(i)-P8(i-1))/(x(i)-x(i-1));
end

% P & dP (D2_delt4_k4) / 9
P9 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p9
    C5 = ro*p9/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p9/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p9/Pout)^d));
    C6 = ro*p9/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k4-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p9/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D2^c*h^2)/2-k4-a*lan*D2^d*h) == C5*x(i)+C6;
    P9(i) = vpasolve(eqn,p9);
end
P9 = -P9+max(P9);
dP9 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP9(1) = 0;
for i = 2:n
    dP9(i) = (P9(i)-P9(i-1))/(x(i)-x(i-1));
end

% P & dP (D2_delt4_k5) / 10
P10 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p10
    C5 = ro*p10/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p10/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p10/Pout)^d));
    C6 = ro*p10/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k5-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p10/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D2^c*h^2)/2-k5-a*lan*D2^d*h) == C5*x(i)+C6;
    P10(i) = vpasolve(eqn,p10);
end
P10 = -P10+max(P10);
dP10 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP10(1) = 0;
for i = 2:n
    dP10(i) = (P10(i)-P10(i-1))/(x(i)-x(i-1));
end

% P & dP (D3_delt4_k1) / 11
P11 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p11
    C5 = ro*p11/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p11/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p11/Pout)^d));
    C6 = ro*p11/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k1-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p11/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D3^c*h^2)/2-k1-a*lan*D3^d*h) == C5*x(i)+C6;
    P11(i) = vpasolve(eqn,p11);
end
P11 = -P11+max(P11);
dP11 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP11(1) = 0;
for i = 2:n
    dP11(i) = (P11(i)-P11(i-1))/(x(i)-x(i-1));
end

% P & dP (D3_delt4_k2) / 12
P12 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p12
    C5 = ro*p12/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p12/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p12/Pout)^d));
    C6 = ro*p12/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k2-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p12/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D3^c*h^2)/2-k2-a*lan*D3^d*h) == C5*x(i)+C6;
    P12(i) = vpasolve(eqn,p12);
end
P12 = -P12+max(P12);
dP12 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP12(1) = 0;
for i = 2:n
    dP12(i) = (P12(i)-P12(i-1))/(x(i)-x(i-1));
end

% P & dP (D3_delt4_k3) / 13
P13 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p13
    C5 = ro*p13/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p13/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p13/Pout)^d));
    C6 = ro*p13/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k3-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p13/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D3^c*h^2)/2-k3-a*lan*D3^d*h) == C5*x(i)+C6;
    P13(i) = vpasolve(eqn,p13);
end
P13 = -P13+max(P13);
dP13 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP13(1) = 0;
for i = 2:n
    dP13(i) = (P13(i)-P13(i-1))/(x(i)-x(i-1));
end

% P & dP (D3_delt4_k4) / 14
P14 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p14
    C5 = ro*p14/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p14/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p14/Pout)^d));
    C6 = ro*p14/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k4-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p14/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D3^c*h^2)/2-k4-a*lan*D3^d*h) == C5*x(i)+C6;
    P14(i) = vpasolve(eqn,p14);
end
P14 = -P14+max(P14);
dP14 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP14(1) = 0;
for i = 2:n
    dP14(i) = (P14(i)-P14(i-1))/(x(i)-x(i-1));
end

% P & dP (D3_delt4_k5) / 15
P15 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p15
    C5 = ro*p15/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p15/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p15/Pout)^d));
    C6 = ro*p15/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k5-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p15/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D3^c*h^2)/2-k5-a*lan*D3^d*h) == C5*x(i)+C6;
    P15(i) = vpasolve(eqn,p15);
end
P15 = -P15+max(P15);
dP15 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP15(1) = 0;
for i = 2:n
    dP15(i) = (P15(i)-P15(i-1))/(x(i)-x(i-1));
end

% P & dP (D4_delt4_k1) / 16
P16 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p16
    C5 = ro*p16/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p16/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p16/Pout)^d));
    C6 = ro*p16/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k1-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p16/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D4^c*h^2)/2-k1-a*lan*D4^d*h) == C5*x(i)+C6;
    P16(i) = vpasolve(eqn,p16);
end
P16 = -P16+max(P16);
dP16 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP16(1) = 0;
for i = 2:n
    dP16(i) = (P16(i)-P16(i-1))/(x(i)-x(i-1));
end

% P & dP (D4_delt4_k2) / 17
P17 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p17
    C5 = ro*p17/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p17/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p17/Pout)^d));
    C6 = ro*p17/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k2-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p17/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D4^c*h^2)/2-k2-a*lan*D4^d*h) == C5*x(i)+C6;
    P17(i) = vpasolve(eqn,p17);
end
P17 = -P17+max(P17);
dP17 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP17(1) = 0;
for i = 2:n
    dP17(i) = (P17(i)-P17(i-1))/(x(i)-x(i-1));
end

% P & dP (D4_delt4_k3) / 18
P18= zeros(n,1); % Pa | pressure
for i = 1:n
    syms p18
    C5 = ro*p18/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p18/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p18/Pout)^d));
    C6 = ro*p18/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k3-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p18/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D4^c*h^2)/2-k3-a*lan*D4^d*h) == C5*x(i)+C6;
    P18(i) = vpasolve(eqn,p18);
end
P18 = -P18+max(P18);
dP18 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP18(1) = 0;
for i = 2:n
    dP18(i) = (P18(i)-P18(i-1))/(x(i)-x(i-1));
end

% P & dP (D4_delt4_k4) / 19
P19= zeros(n,1); % Pa | pressure
for i = 1:n
    syms p19
    C5 = ro*p19/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p19/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p19/Pout)^d));
    C6 = ro*p19/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k4-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p19/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D4^c*h^2)/2-k4-a*lan*D4^d*h) == C5*x(i)+C6;
    P19(i) = vpasolve(eqn,p19);
end
P19 = -P19+max(P19);
dP19 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP19(1) = 0;
for i = 2:n
    dP19(i) = (P19(i)-P19(i-1))/(x(i)-x(i-1));
end

% P & dP (D4_delt4_k5) / 20
P20= zeros(n,1); % Pa | pressure
for i = 1:n
    syms p20
    C5 = ro*p20/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p20/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p20/Pout)^d));
    C6 = ro*p20/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k5-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p20/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D4^c*h^2)/2-k5-a*lan*D4^d*h) == C5*x(i)+C6;
    P20(i) = vpasolve(eqn,p20);
end
P20 = -P20+max(P20);
dP20 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP20(1) = 0;
for i = 2:n
    dP20(i) = (P20(i)-P20(i-1))/(x(i)-x(i-1));
end

% P & dP (D1_delt1_k4) / 21
P21 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p21
    C5 = ro*p21/(mu*L)*(4*h^2+delt1^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p21/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p21/Pout)^d));
    C6 = ro*p21/mu*(4*h^2+delt1^2)*((y1(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k4-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p21/mu*(4*h^2+delt1^2)*((y1(i)^2+h^2-b*D4^c*h^2)/2-k4-a*lan*D4^d*h) == C5*x(i)+C6;
    P21(i) = vpasolve(eqn,p21);
end
P21 = -P21+max(P21);
dP21 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP21(1) = 0;
for i = 2:n
    dP21(i) = (P21(i)-P21(i-1))/(x(i)-x(i-1));
end

% P & dP (D1_delt2_k4) / 22
P22 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p22
    C5 = ro*p22/(mu*L)*(4*h^2+delt2^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p22/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p22/Pout)^d));
    C6 = ro*p22/mu*(4*h^2+delt2^2)*((y2(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k4-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p22/mu*(4*h^2+delt2^2)*((y2(i)^2+h^2-b*D4^c*h^2)/2-k4-a*lan*D4^d*h) == C5*x(i)+C6;
    P22(i) = vpasolve(eqn,p22);
end
P22 = -P22+max(P22);
dP22 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP22(1) = 0;
for i = 2:n
    dP22(i) = (P22(i)-P22(i-1))/(x(i)-x(i-1));
end

% P & dP (D1_delt3_k4) / 23
P23 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p23
    C5 = ro*p23/(mu*L)*(4*h^2+delt3^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p23/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p23/Pout)^d));
    C6 = ro*p23/mu*(4*h^2+delt3^2)*((y3(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k4-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p23/mu*(4*h^2+delt3^2)*((y3(i)^2+h^2-b*D4^c*h^2)/2-k4-a*lan*D4^d*h) == C5*x(i)+C6;
    P23(i) = vpasolve(eqn,p23);
end
P23 = -P23+max(P23);
dP23 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP23(1) = 0;
for i = 2:n
    dP23(i) = (P23(i)-P23(i-1))/(x(i)-x(i-1));
end

% P & dP (D1_delt4_k4) / 24
P24 = zeros(n,1); % Pa | pressure
for i = 1:n
    syms p24
    C5 = ro*p24/(mu*L)*(4*h^2+delt4^2)*(b/2*h^2*Dout^c*((Pin/Pout)^c-(p24/Pout)^c)+a*lan*h*Dout^d*((Pin/Pout)^d-(p24/Pout)^d));
    C6 = ro*p24/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*(Dout*(Pin/Pout))^c*h^2)/2-k4-a*lan*(Dout*(Pin/Pout))^d*h);
    eqn = ro*p24/mu*(4*h^2+delt4^2)*((y4(i)^2+h^2-b*D4^c*h^2)/2-k4-a*lan*D4^d*h) == C5*x(i)+C6;
    P24(i) = vpasolve(eqn,p24);
end
P24 = -P24+max(P24);
dP24 = zeros(n,1); % Pa/m | pressure with respect to the x-axis
dP24(1) = 0;
for i = 2:n
    dP24(i) = (P24(i)-P24(i-1))/(x(i)-x(i-1));
end

%% Deducing the Velocity Profiles
for i = 1:n
    u1(i) = (dP1(i)*y4(i)^2)/(2*mu)+dP1(i)*(h^2/(2*mu)-k1/mu-a*lan*D1^d*h/mu-b/(2*mu)*D1^c*h^2);
    u_tilde1(i) = dP1(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u2(i) = (dP2(i)*y4(i)^2)/(2*mu)+dP2(i)*(h^2/(2*mu)-k2/mu-a*lan*D1^d*h/mu-b/(2*mu)*D1^c*h^2);
    u_tilde2(i) = dP2(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u3(i) = (dP3(i)*y4(i)^2)/(2*mu)+dP3(i)*(h^2/(2*mu)-k3/mu-a*lan*D1^d*h/mu-b/(2*mu)*D1^c*h^2);
    u_tilde3(i) = dP3(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u4(i) = (dP4(i)*y4(i)^2)/(2*mu)+dP4(i)*(h^2/(2*mu)-k4/mu-a*lan*D1^d*h/mu-b/(2*mu)*D1^c*h^2);
    u_tilde4(i) = dP4(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u5(i) = (dP5(i)*y4(i)^2)/(2*mu)+dP5(i)*(h^2/(2*mu)-k5/mu-a*lan*D1^d*h/mu-b/(2*mu)*D1^c*h^2);
    u_tilde5(i) = dP5(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u6(i) = (dP6(i)*y4(i)^2)/(2*mu)+dP6(i)*(h^2/(2*mu)-k1/mu-a*lan*D2^d*h/mu-b/(2*mu)*D2^c*h^2);
    u_tilde6(i) = dP6(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u7(i) = (dP7(i)*y4(i)^2)/(2*mu)+dP7(i)*(h^2/(2*mu)-k2/mu-a*lan*D2^d*h/mu-b/(2*mu)*D2^c*h^2);
    u_tilde7(i) = dP7(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u8(i) = (dP8(i)*y4(i)^2)/(2*mu)+dP8(i)*(h^2/(2*mu)-k3/mu-a*lan*D2^d*h/mu-b/(2*mu)*D2^c*h^2);
    u_tilde8(i) = dP8(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u9(i) = (dP9(i)*y4(i)^2)/(2*mu)+dP9(i)*(h^2/(2*mu)-k4/mu-a*lan*D2^d*h/mu-b/(2*mu)*D2^c*h^2);
    u_tilde9(i) = dP9(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u10(i) = (dP10(i)*y4(i)^2)/(2*mu)+dP10(i)*(h^2/(2*mu)-k5/mu-a*lan*D2^d*h/mu-b/(2*mu)*D2^c*h^2);
    u_tilde10(i) = dP10(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u11(i) = (dP11(i)*y4(i)^2)/(2*mu)+dP11(i)*(h^2/(2*mu)-k1/mu-a*lan*D3^d*h/mu-b/(2*mu)*D3^c*h^2);
    u_tilde11(i) = dP11(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u12(i) = (dP12(i)*y4(i)^2)/(2*mu)+dP12(i)*(h^2/(2*mu)-k2/mu-a*lan*D3^d*h/mu-b/(2*mu)*D3^c*h^2);
    u_tilde12(i) = dP12(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u13(i) = (dP13(i)*y4(i)^2)/(2*mu)+dP13(i)*(h^2/(2*mu)-k3/mu-a*lan*D3^d*h/mu-b/(2*mu)*D3^c*h^2);
    u_tilde13(i) = dP13(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u14(i) = (dP14(i)*y4(i)^2)/(2*mu)+dP14(i)*(h^2/(2*mu)-k4/mu-a*lan*D3^d*h/mu-b/(2*mu)*D3^c*h^2);
    u_tilde14(i) = dP14(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u15(i) = (dP15(i)*y4(i)^2)/(2*mu)+dP15(i)*(h^2/(2*mu)-k5/mu-a*lan*D3^d*h/mu-b/(2*mu)*D3^c*h^2);
    u_tilde15(i) = dP15(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u16(i) = (dP16(i)*y4(i)^2)/(2*mu)+dP16(i)*(h^2/(2*mu)-k1/mu-a*lan*D4^d*h/mu-b/(2*mu)*D4^c*h^2);
    u_tilde16(i) = dP16(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u17(i) = (dP17(i)*y4(i)^2)/(2*mu)+dP17(i)*(h^2/(2*mu)-k2/mu-a*lan*D4^d*h/mu-b/(2*mu)*D4^c*h^2);
    u_tilde17(i) = dP17(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u18(i) = (dP18(i)*y4(i)^2)/(2*mu)+dP18(i)*(h^2/(2*mu)-k3/mu-a*lan*D4^d*h/mu-b/(2*mu)*D4^c*h^2);
    u_tilde18(i) = dP18(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u19(i) = (dP19(i)*y4(i)^2)/(2*mu)+dP19(i)*(h^2/(2*mu)-k4/mu-a*lan*D4^d*h/mu-b/(2*mu)*D4^c*h^2);
    u_tilde19(i) = dP19(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u20(i) = (dP20(i)*y4(i)^2)/(2*mu)+dP20(i)*(h^2/(2*mu)-k5/mu-a*lan*D4^d*h/mu-b/(2*mu)*D4^c*h^2);
    u_tilde20(i) = dP20(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
    u21(i) = (dP21(i)*y1(i)^2)/(2*mu)+dP21(i)*(h^2/(2*mu)-k4/mu-a*lan*D4^d*h/mu-b/(2*mu)*D4^c*h^2);
    u_tilde21(i) = dP21(i)*(h*y1(i)/mu-(h+delt1)*h/mu);
    u22(i) = (dP22(i)*y2(i)^2)/(2*mu)+dP22(i)*(h^2/(2*mu)-k4/mu-a*lan*D4^d*h/mu-b/(2*mu)*D4^c*h^2);
    u_tilde22(i) = dP22(i)*(h*y2(i)/mu-(h+delt2)*h/mu);
    u23(i) = (dP23(i)*y3(i)^2)/(2*mu)+dP23(i)*(h^2/(2*mu)-k4/mu-a*lan*D4^d*h/mu-b/(2*mu)*D4^c*h^2);
    u_tilde23(i) = dP23(i)*(h*y3(i)/mu-(h+delt3)*h/mu);
    u24(i) = (dP24(i)*y4(i)^2)/(2*mu)+dP24(i)*(h^2/(2*mu)-k4/mu-a*lan*D4^d*h/mu-b/(2*mu)*D4^c*h^2);
    u_tilde24(i) = dP24(i)*(h*y4(i)/mu-(h+delt4)*h/mu);
end

%% Trimming the Data
Y1 = y1./h; % dimensionless coordinates | mode 1
Y2 = y2./h; % dimensionless coordinates | mode 2
Y3 = y3./h; % dimensionless coordinates | mode 3
Y4 = y4./h; % dimensionless coordinates | mode 4

for i = 1:n
    if abs(Y4(i)) >= 1
        u1(i) = NaN;
        u2(i) = NaN;
        u3(i) = NaN;
        u4(i) = NaN;
        u5(i) = NaN;
        u6(i) = NaN;
        u7(i) = NaN;
        u8(i) = NaN;
        u9(i) = NaN;
        u10(i) = NaN;
        u11(i) = NaN;
        u12(i) = NaN;
        u13(i) = NaN;
        u14(i) = NaN;
        u15(i) = NaN;
        u16(i) = NaN;
        u17(i) = NaN;
        u18(i) = NaN;
        u19(i) = NaN;
        u20(i) = NaN;
        u24(i) = NaN;
    end
    if Y4(i) <= 1
        u_tilde1(i) = NaN;
        u_tilde2(i) = NaN;
        u_tilde3(i) = NaN;
        u_tilde4(i) = NaN;
        u_tilde5(i) = NaN;
        u_tilde6(i) = NaN;
        u_tilde7(i) = NaN;
        u_tilde8(i) = NaN;
        u_tilde9(i) = NaN;
        u_tilde10(i) = NaN;
        u_tilde11(i) = NaN;
        u_tilde12(i) = NaN;
        u_tilde13(i) = NaN;
        u_tilde14(i) = NaN;
        u_tilde15(i) = NaN;
        u_tilde16(i) = NaN;
        u_tilde17(i) = NaN;
        u_tilde18(i) = NaN;
        u_tilde19(i) = NaN;
        u_tilde20(i) = NaN;
        u_tilde24(i) = NaN;
    end
end

for i = 1:n
    if abs(Y1(i)) >= 1
        u21(i) = NaN;
    end
    if Y1(i) <= 1
        u_tilde21(i) = NaN;
    end
end
for i = 1:n
    if abs(Y2(i)) >= 1
        u22(i) = NaN;
    end
    if Y2(i) <= 1
        u_tilde22(i) = NaN;
    end
end
for i = 1:n
    if abs(Y3(i)) >= 1
        u23(i) = NaN;
    end
    if Y3(i) <= 1
        u_tilde23(i) = NaN;
    end
end

%% Velocity Profiles
% 1-5
figure
C = linspecer(5,'qualitative');
plot(u1,Y4,'color',C(1,:),'LineWidth',1.5);
hold on
h1=plot(u_tilde1,Y4,'color',C(1,:),'LineWidth',1.5);
plot(u2,Y4,'color',C(2,:),'LineWidth',1.5);
h2=plot(u_tilde2,Y4,'color',C(2,:),'LineWidth',1.5);
plot(u3,Y4,'color',C(3,:),'LineWidth',1.5);
h3=plot(u_tilde3,Y4,'color',C(3,:),'LineWidth',1.5);
plot(u4,Y4,'color',C(4,:),'LineWidth',1.5);
h4=plot(u_tilde4,Y4,'color',C(4,:),'LineWidth',1.5);
plot(u5,Y4,'color',C(5,:),'LineWidth',1.5);
h5=plot(u_tilde5,Y4,'color',C(5,:),'LineWidth',1.5);
ylim([0 max(Y4)])
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'k = 0','k = 10^{-8}','k = 10^{-7}','k = 10^{-6}','k = 10^{-5}','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 0.01, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('u (^{m}/_{s})','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('Y = ^{y}/_{h}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('Fig2a.tiff','-dtiff',['-r' num2str(600)])

% 6-10
figure
C = linspecer(5,'qualitative');
plot(u6,Y4,'color',C(1,:),'LineWidth',1.5);
hold on
h1=plot(u_tilde6,Y4,'color',C(1,:),'LineWidth',1.5);
plot(u7,Y4,'color',C(2,:),'LineWidth',1.5);
h2=plot(u_tilde7,Y4,'color',C(2,:),'LineWidth',1.5);
plot(u8,Y4,'color',C(3,:),'LineWidth',1.5);
h3=plot(u_tilde8,Y4,'color',C(3,:),'LineWidth',1.5);
plot(u9,Y4,'color',C(4,:),'LineWidth',1.5);
h4=plot(u_tilde9,Y4,'color',C(4,:),'LineWidth',1.5);
plot(u10,Y4,'color',C(5,:),'LineWidth',1.5);
h5=plot(u_tilde10,Y4,'color',C(5,:),'LineWidth',1.5);
ylim([0 max(Y4)])
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'k = 0','k = 10^{-8}','k = 10^{-7}','k = 10^{-6}','k = 10^{-5}','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 0.1, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('u (^{m}/_{s})','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('Y = ^{y}/_{h}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('Fig2b.tiff','-dtiff',['-r' num2str(600)])

% 11-15
figure
C = linspecer(5,'qualitative');
plot(u11,Y4,'color',C(1,:),'LineWidth',1.5);
hold on
h1=plot(u_tilde11,Y4,'color',C(1,:),'LineWidth',1.5);
plot(u12,Y4,'color',C(2,:),'LineWidth',1.5);
h2=plot(u_tilde12,Y4,'color',C(2,:),'LineWidth',1.5);
plot(u13,Y4,'color',C(3,:),'LineWidth',1.5);
h3=plot(u_tilde13,Y4,'color',C(3,:),'LineWidth',1.5);
plot(u14,Y4,'color',C(4,:),'LineWidth',1.5);
h4=plot(u_tilde14,Y4,'color',C(4,:),'LineWidth',1.5);
plot(u15,Y4,'color',C(5,:),'LineWidth',1.5);
h5=plot(u_tilde15,Y4,'color',C(5,:),'LineWidth',1.5);
ylim([0 max(Y4)])
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'k = 0','k = 10^{-8}','k = 10^{-7}','k = 10^{-6}','k = 10^{-5}','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 1, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('u (^{m}/_{s})','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('Y = ^{y}/_{h}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('Fig2c.tiff','-dtiff',['-r' num2str(600)])

% 16-20
figure
C = linspecer(5,'qualitative');
plot(u16,Y4,'color',C(1,:),'LineWidth',1.5);
hold on
h1=plot(u_tilde16,Y4,'color',C(1,:),'LineWidth',1.5);
plot(u17,Y4,'color',C(2,:),'LineWidth',1.5);
h2=plot(u_tilde17,Y4,'color',C(2,:),'LineWidth',1.5);
plot(u18,Y4,'color',C(3,:),'LineWidth',1.5);
h3=plot(u_tilde18,Y4,'color',C(3,:),'LineWidth',1.5);
plot(u19,Y4,'color',C(4,:),'LineWidth',1.5);
h4=plot(u_tilde19,Y4,'color',C(4,:),'LineWidth',1.5);
plot(u20,Y4,'color',C(5,:),'LineWidth',1.5);
h5=plot(u_tilde20,Y4,'color',C(5,:),'LineWidth',1.5);
ylim([0 max(Y4)])
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'k = 0','k = 10^{-8}','k = 10^{-7}','k = 10^{-6}','k = 10^{-5}','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 10, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('u (^{m}/_{s})','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('Y = ^{y}/_{h}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('Fig2d.tiff','-dtiff',['-r' num2str(600)])

% 21-24
figure
C = linspecer(4,'qualitative');
plot(u21,Y4,'color',C(1,:),'LineWidth',1.5);
hold on
h1=plot(u_tilde21,Y4,'color',C(1,:),'LineWidth',1.5);
plot(u22,Y4,'color',C(2,:),'LineWidth',1.5);
h2=plot(u_tilde22,Y4,'color',C(2,:),'LineWidth',1.5);
plot(u23,Y4,'color',C(3,:),'LineWidth',1.5);
h3=plot(u_tilde23,Y4,'color',C(3,:),'LineWidth',1.5);
plot(u24,Y4,'color',C(4,:),'LineWidth',1.5);
h4=plot(u_tilde24,Y4,'color',C(4,:),'LineWidth',1.5);
ylim([0 max(Y4)])
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4],'\delta = 3.7×10^{-6} m','\delta = 7.5×10^{-6} m','\delta = 2.5×10^{-5} m','\delta = 4.8×10^{-5} m','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 10, k = 10^{-7}','FontSize',10.5,'FontName','Helvetica')
xlabel('u (^{m}/_{s})','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('Y = ^{y}/_{h}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('Fig2e.tiff','-dtiff',['-r' num2str(600)])

%% Velocity Profiles (Appendix A)
% 1-5
figure
C = linspecer(5,'qualitative');
h1=plot(u_tilde1,Y4,'color',C(1,:),'LineWidth',1.5);
hold on
h2=plot(u_tilde2,Y4,'color',C(2,:),'LineWidth',1.5);
h3=plot(u_tilde3,Y4,'color',C(3,:),'LineWidth',1.5);
h4=plot(u_tilde4,Y4,'color',C(4,:),'LineWidth',1.5);
h5=plot(u_tilde5,Y4,'color',C(5,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'k = 0','k = 10^{-8}','k = 10^{-7}','k = 10^{-6}','k = 10^{-5}','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 10, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('u (^{m}/_{s})','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('Y = ^{y}/_{h}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('FigA1a.tiff','-dtiff',['-r' num2str(600)])

% 6-10
figure
C = linspecer(5,'qualitative');
h1=plot(u_tilde6,Y4,'color',C(1,:),'LineWidth',1.5);
hold on
h2=plot(u_tilde7,Y4,'color',C(2,:),'LineWidth',1.5);
h3=plot(u_tilde8,Y4,'color',C(3,:),'LineWidth',1.5);
h4=plot(u_tilde9,Y4,'color',C(4,:),'LineWidth',1.5);
h5=plot(u_tilde10,Y4,'color',C(5,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'k = 0','k = 10^{-8}','k = 10^{-7}','k = 10^{-6}','k = 10^{-5}','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 0.1, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('u (^{m}/_{s})','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('Y = ^{y}/_{h}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('FigA1b.tiff','-dtiff',['-r' num2str(600)])

% 11-15
figure
C = linspecer(5,'qualitative');
h1=plot(u_tilde11,Y4,'color',C(1,:),'LineWidth',1.5);
hold on
h2=plot(u_tilde12,Y4,'color',C(2,:),'LineWidth',1.5);
h3=plot(u_tilde13,Y4,'color',C(3,:),'LineWidth',1.5);
h4=plot(u_tilde14,Y4,'color',C(4,:),'LineWidth',1.5);
h5=plot(u_tilde15,Y4,'color',C(5,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'k = 0','k = 10^{-8}','k = 10^{-7}','k = 10^{-6}','k = 10^{-5}','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 1, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('u (^{m}/_{s})','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('Y = ^{y}/_{h}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('FigA1c.tiff','-dtiff',['-r' num2str(600)])

% 16-20
figure
C = linspecer(5,'qualitative');
h1=plot(u_tilde16,Y4,'color',C(1,:),'LineWidth',1.5);
hold on
h2=plot(u_tilde17,Y4,'color',C(2,:),'LineWidth',1.5);
h3=plot(u_tilde18,Y4,'color',C(3,:),'LineWidth',1.5);
h4=plot(u_tilde19,Y4,'color',C(4,:),'LineWidth',1.5);
h5=plot(u_tilde20,Y4,'color',C(5,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'k = 0','k = 10^{-8}','k = 10^{-7}','k = 10^{-6}','k = 10^{-5}','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 10, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('u (^{m}/_{s})','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('Y = ^{y}/_{h}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('FigA1d.tiff','-dtiff',['-r' num2str(600)])

% 21-24
figure
C = linspecer(4,'qualitative');
h1=plot(u_tilde21,Y4,'color',C(1,:),'LineWidth',1.5);
hold on
h2=plot(u_tilde22,Y4,'color',C(2,:),'LineWidth',1.5);
h3=plot(u_tilde23,Y4,'color',C(3,:),'LineWidth',1.5);
h4=plot(u_tilde24,Y4,'color',C(4,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4],'\delta = 3.7×10^{-6} m','\delta = 7.5×10^{-6} m','\delta = 2.5×10^{-5} m','\delta = 4.8×10^{-5} m','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 10, k = 10^{-6}','FontSize',10.5,'FontName','Helvetica')
xlabel('u (^{m}/_{s})','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('Y = ^{y}/_{h}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('FigA1e.tiff','-dtiff',['-r' num2str(600)])

%% Pressure Distribution
figure
C = linspecer(4,'qualitative');
h1=plot(x,P3,'color',C(1,:),'LineWidth',1.5);
hold on
h2=plot(x,P8,'color',C(2,:),'LineWidth',1.5);
h3=plot(x,P13,'color',C(3,:),'LineWidth',1.5);
h4=plot(x,P18,'color',C(4,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'D = 0.01','D = 0.1','D = 1','D = 10','FontSize',10.75,'FontName','Helvetica');
title(leg,'k = 10^{-7}, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('x (m)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('p (Pa)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('Fig3.tiff','-dtiff',['-r' num2str(600)])

%% Pressure Distribution (Appendix A)
% 1-5
figure
C = linspecer(5,'qualitative');
h1=plot(x,P1,'color',C(1,:),'LineWidth',1.5);
hold on
h2=plot(x,P2,'color',C(2,:),'LineWidth',1.5);
h3=plot(x,P3,'color',C(3,:),'LineWidth',1.5);
h4=plot(x,P4,'color',C(4,:),'LineWidth',1.5);
h5=plot(x,P5,'color',C(5,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'k = 0','k = 10^{-8}','k = 10^{-7}','k = 10^{-6}','k = 10^{-5}','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 0.01, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('x (m)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('p (Pa)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('FigB1a.tiff','-dtiff',['-r' num2str(600)])

% 6-10
figure
C = linspecer(5,'qualitative');
h1=plot(x,P6,'color',C(1,:),'LineWidth',1.5);
hold on
h2=plot(x,P7,'color',C(2,:),'LineWidth',1.5);
h3=plot(x,P8,'color',C(3,:),'LineWidth',1.5);
h4=plot(x,P9,'color',C(4,:),'LineWidth',1.5);
h5=plot(x,P10,'color',C(5,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'k = 0','k = 10^{-8}','k = 10^{-7}','k = 10^{-6}','k = 10^{-5}','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 0.1, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('x (m)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('p (Pa)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('FigB1b.tiff','-dtiff',['-r' num2str(600)])

% 11-15
figure
C = linspecer(5,'qualitative');
h1=plot(x,P11,'color',C(1,:),'LineWidth',1.5);
hold on
h2=plot(x,P12,'color',C(2,:),'LineWidth',1.5);
h3=plot(x,P13,'color',C(3,:),'LineWidth',1.5);
h4=plot(x,P14,'color',C(4,:),'LineWidth',1.5);
h5=plot(x,P15,'color',C(5,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'k = 0','k = 10^{-8}','k = 10^{-7}','k = 10^{-6}','k = 10^{-5}','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 1, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('x (m)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('p (Pa)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('FigB1c.tiff','-dtiff',['-r' num2str(600)])

% 16-20
figure
C = linspecer(5,'qualitative');
h1=plot(x,P16,'color',C(1,:),'LineWidth',1.5);
hold on
h2=plot(x,P17,'color',C(2,:),'LineWidth',1.5);
h3=plot(x,P18,'color',C(3,:),'LineWidth',1.5);
h4=plot(x,P19,'color',C(4,:),'LineWidth',1.5);
h5=plot(x,P20,'color',C(5,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4;h5],'k = 0','k = 10^{-8}','k = 10^{-7}','k = 10^{-6}','k = 10^{-5}','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 10, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('x (m)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('p (Pa)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('FigB1d.tiff','-dtiff',['-r' num2str(600)])

% 21-24
figure
C = linspecer(4,'qualitative');
h1=plot(x,P21,'color',C(1,:),'LineWidth',1.5);
hold on
h2=plot(x,P22,'color',C(2,:),'LineWidth',1.5);
h3=plot(x,P23,'color',C(3,:),'LineWidth',1.5);
h4=plot(x,P24,'color',C(4,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2;h3;h4],'\delta = 3.7×10^{-6} m','\delta = 7.5×10^{-6} m','\delta = 2.5×10^{-5} m','\delta = 4.8×10^{-5} m','FontSize',10.75,'FontName','Helvetica');
title(leg,'D = 10, k = 10^{-6}','FontSize',10.5,'FontName','Helvetica')
xlabel('x (m)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('p (Pa)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('FigB1e.tiff','-dtiff',['-r' num2str(600)])

%% Comparison with COMSOL Data
% velocity profile
ufile = fopen('Velocity-COMSOL.txt');
U = fscanf(ufile,'%f %f'); 
for i = 1:2:length(U)
    u_COMSOL((i+1)/2) = U(i);
    y_COMSOL((i+1)/2) = U(i+1);
end
fclose(ufile);
u_COMSOL = flip(u_COMSOL); 
y_COMSOL = flip(y_COMSOL);

% pressure
pfile = fopen('Pressure-COMSOL.txt');
P = fscanf(pfile,'%f %f');
for i = 1:2:length(P)
    x_COMSOL((i+1)/2) = P(i);
    p_COMSOL((i+1)/2) = P(i+1);
end
fclose(pfile);
x_COMSOL = flip(x_COMSOL); 
p_COMSOL = flip(p_COMSOL);

%% Visualising COMSOL Data and Comparing it with MATLAB
% velocity profile
ulength = length(u_COMSOL);
figure
C = linspecer(2,'qualitative');
h1=plot(log(u4),Y4,'color',C(1,:),'LineWidth',1.5);
hold on
plot(u_tilde4,Y4,'color',C(1,:),'LineWidth',1.5);
h2=plot(log(u_COMSOL),y_COMSOL,'color',C(2,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2],'MATLAB Simulation','COMSOL Simulation','FontSize',10.75,'FontName','Helvetica');
title(leg,'k = 10^{-6}, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('ln(u)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('Y = ^{y}/_{h}','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylim([0 max(Y4)])
grid on
print('Fig6a.tiff','-dtiff',['-r' num2str(600)])

% pressure
figure
C = linspecer(2,'qualitative');
h1=plot(x,P19,'color',C(1,:),'LineWidth',1.5);
hold on
h2=plot(x_COMSOL,p_COMSOL,'color',C(2,:),'LineWidth',1.5);
hold off
set(gca,'FontSize',13)
leg = legend([h1;h2],'MATLAB Simulation','COMSOL Simulation','FontSize',10.75,'FontName','Helvetica');
title(leg,'k = 10^{-6}, \delta = 4.8×10^{-5} m','FontSize',10.5,'FontName','Helvetica')
xlabel('x (m)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
ylabel('p (Pa)','FontName','Helvetica','FontSize',14,'FontWeight','bold')
grid on
print('Fig6b.tiff','-dtiff',['-r' num2str(600)])

% Effects of roughness on rarefied gas flow in a microfluidic channel with porous walls 
% Amirreza Soheili, Arsha Niksa, Vahid Bazargan
% Please contact Amirreza Soheili in case you have any questions.