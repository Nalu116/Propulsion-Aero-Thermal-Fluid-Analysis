clear;
close all;
clc
format long g

%{
- Quasi-1D air flow in propulsive duct, one meter long (exception case 3)
- Model duct as circular cross section with linear variation along axial direction for any area variation 
%}

%------------------------------ONE------------------------------

{%
- Numerical boundry value problem
- Isentropic Flow with no work interaction
- Differential Analysis for Pressure, Temp, Velocity, Mach, and Total Cond at exit of duct
- Computes force on internal wetted surfaces of duct
- Compares results with analytical solution
%} 

%-------Input system paramaters--------

n = 100000; %Number of iterations (Step size (aka fidelity)) 

y = 1.4; %gamma air
P0 = 101325; %Static Pressure (Pa) 
T0 = 288; %Static Termperature (K) 
M0 = 0.6; %Mach Number 
A0 = 0.1; %Entrace cross sectional area (m^2)
Ae = 0.3; %Exit cross sectional area (m^2)
R = 287;
%-------End System Input----------

rho0 = P0/(R*T0)'
a0 = sqrt(y*R*T0);
u0 = M0 * a0;
mdot = rho0*u*A0;
Pt0 = P0*(1+((y-1)/2)*M0^2)^(y/(y-1));
Tt0 = T0*(1+((y-1)/2)*M0^2);

%Cf based on Reynolds Number (and dynamic pressure + skin shear by def)
mu0 = 0.0000178866;
Re0 = (rho0*u0*1)/mu0;

%Define linearly varying area differential element 
dA = (Ae-A0)/n;

{%
          INTEGRATED QUASI-1D EQUATIONS
1.) Continuity Equation : mdot = rho1 * u1 * A1 = rho2 * u2 * A2 
2.) Momentum Equation : Fx = -(mdot*u2 + P2*A2 - mdot*u1 - P1*A1)
3.) Energy : Cp(T2-T1) + 1/2(u2^2 - u1^2) = q_1->2 + w_1->2 (per unit mass)
4.) Equation of State : R = P1/rho1T1 = P2/rho2T2
5.) Entropy : dS = s2-s1 = Cpln(T2/T1) - Rln(P2/P1) = S_heat + S_gen 
%} 

% ------ Analytical Solution---------
%"Isentropic flow with no work interactions"
%Make sure to change which value of analytical M2 taken based on
%Expected subsonic or super sonic conditions for future area variations!
Tt2 = Tt0;
Pt2 = Pt0;
syms M2g
eqn1 = Ae*(1/M0*((y-1)/2*M0^2+1)^(1/(y-1))*(1+((y-1)/2)*M2g^2)^(-1/(y-1))) * ((y+1)/(2*(y-1)))^(1/2)*((y+1)/2)^(-1/(y-1))) == A0*(1/M2g*((y+1)/(2*(y-1)))^(1/2)*((y+1)/2*(M2g^2+1))^(-1/(y-1)));
sol1 = vpasolve(eqn1,M2g);
M2 = sol1(1); %take the real value that is subsonic in this case (diverging)
Pe = Pt2 / (1+((y-1)/2)*M2^2)^(y/(y-1));
Te = Tt2 / (1+((y-1)/2)*M2^2);
rho2 = Pe / (R*Te);
ue = (u0*A0*rho0) / (Ae*rho2);
Force = (rho2*ue^2 +Pe)*Ae - (rho0*u0^2 +P0)*A0;

% ------ Numerical Iterations---------

i = 1;

%local conditions of array begin at initial position of array elements prior to steps
%{
Tmat = zeros(1, n);
Pmat = zeros(1, n);
umat = zeros(1, n);
rhomat = zeros(1, n);
%}


