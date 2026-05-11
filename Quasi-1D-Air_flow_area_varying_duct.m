clear;
close all;
clc
format long g

%{
- Quasi-1D air flow in propulsive duct, one meter long 
- ModelING duct as circular cross section with linear variation along axial direction for any area variation 
- Isentropic Flow with no work interaction
- Differential Analysis for Pressure, Temp, Velocity, Mach and total cond
  at exit of duct
- Compute force on internal wetted surfaces of the duct
- Compare results with analytical solution
%}


n = 100000; %Number of Iterations (Fidelity Step Size)
%n = 500; %Temporary step size to protect against wrinkles

y = 1.4; %gamma for air
P0 = 101325; %Static Pressure (Pa)
T0 = 288; %Static Temperature (K)
M0 = 0.6; %Mach number
A0 = 0.1; %Entrance cross sectional area (m^2)
Ae = 0.3; %Exit cross sectional area (m^2)
R = 287;
rho0 = P0/(R*T0);
a0 = sqrt(y*R*T0);
u0 = M0*a0;
mdot = rho0*u0*A0;
Pt0 = P0*(1+((y-1)/2)*M0^2)^(y/(y-1));
Tt0 = T0*(1+((y-1)/2)*M0^2);

%Cf based on reynolds number
mu0 = 0.0000178866;
Re0 = (rho0*u0*1)/mu0;
dA = (Ae-A0)/n; %Definine the linearly varying area differential
%q0 = (1/2) * (rho0)*u0^2;

%Cf based on dynamic pressure and skin shear by def

%{
                    INTEGRATED QUASI-1D EQUATIONS

1.) Continuity Equation : mdot = rho1 * u1 * A1 = rho2 * u2 * A2
2.) Momentum Equation : Fx = -(mdot*u2 + P2*A2 - mdot*u1 - P1*A1)
3.) Energy : Cp(T2-T1) + 1/2(u2^2 - u1^2) = q_1->2 + w_1->2 (per unit mass)
4.) Equation of State : R = P1/rho1T1 = P2/rho2T2
5.) Entropy dS = s2-s1 = Cpln(T2/T1) - Rln(P2/P1) = Sheat +Sgen
    
%}

%---------Analytical Solution-------
%"Isentropic flow with no work interactions"
%Make sure to change which value of analytical M2 taken based on 
%Expected subsonic or super sonic conditions for future area variations!

Tt2 = Tt0;
Pt2 = Pt0;


syms M2q
eqn1 = Ae*(1/M0*((2/(y+1))*(1+((y-1)/2)*M0^2))^((y+1)/(2*(y-1)))) == A0*(1/M2q*((2/(y+1))*(1+((y-1)/2)*M2q^2))^((y+1)/(2*(y-1))));

sol1 = vpasolve(eqn1,M2q);
M2 = sol1(1); %Take the real value that is subsonic in this case (diverging)
Pe = Pt2/(1+((y-1)/2)*M2^2)^(y/(y-1));
Te = Tt2/(1+((y-1)/2)*M2^2);
rho2 = Pe/(Te*R);
ue = (u0*A0*rho0)/(Ae*rho2);
Force = (rho2*ue^2+Pe)*Ae - (rho0*u0^2+P0)*A0;



% ---------- Numerical Itterations---------------------------------------

i = 1;    %Initialize counter and starting position of array elements

%local conditions of array begin at initial conditions prior to steps

Tmat(i) = T0;
Pmat(i) = P0;
umat(i) = u0;
rhomat(i) = rho0;
Machmat(i) = M0;
Amat(i) = A0;

h = waitbar(0,'Initializing waitbar...'); %Loading bar to track progress
L=1;
c_f = 0;
eta = 1;
dW_shaft = 0;
dx = L/n;
dq = 0;
T_w = 0;
Tw = 0;
c_p = 1004.5;
dQ_rad = 0; 
u(1) = u0;
T(1) = T0;
P(1) = P0;
rho(1) = rho0;
A = linspace(A0,Ae,n);
gamma = y;

i = 2;
while i < (n+1)

    u(i) = ((R*c_p*(u(i-1)*A(i-1)*eta*dW_shaft*mdot*rho(i-1)+u(i-1)^3*A(i-1)*rho(i-1)+u(i-1)*((A(i)-2*A(i-1))*P(i-1)-(c_f*(.5*rho(i-1)*u(i-1)^2))*(2*sqrt(pi*A(i-1)))*dx))+2*u(i-1)*A(i-1)*P(i-1)*c_f*R*c_p*(dx/((2*sqrt(pi*A(i-1)))/pi))+dQ_rad)*gamma*T(i-1)+((u(i-1)^3*A(i-1)*P(i-1)*c_f-2*u(i-1)*A(i-1)*P(i-1)*c_f*T_w*R)*c_p*(dx/((2*sqrt(pi*A(i-1)))/pi))-u(i-1)*A(i-1)*P(i-1)*R*dQ_rad+(-u(i-1)*A(i-1)*P(i-1)*dW_shaft-u(i-1)^3*A(i-1)*P(i-1))*R)*gamma-u(i-1)^3*A(i-1)*P(i-1)*c_f*c_p*(dx/((2*sqrt(pi*A(i-1)))/pi)))/(R*c_p*(u(i-1)^2*A(i-1)*rho(i-1)-A(i-1)*P(i-1))*gamma*T(i-1)-u(i-1)^2*A(i-1)*P(i-1)*R*gamma);
    rho(i) = (rho(i-1)*u(i-1)*A(i-1))/(u(i)*A(i));
    T(i) = (((2*c_p*c_f*(T_w-(T(i-1)*(1+((gamma-1)/2)*((u(i-1)/(sqrt(gamma*R*T(i-1)))))^2)))*((dx)/((2*sqrt(pi*A(i-1)))/pi))+dQ_rad)+dW_shaft-u(i-1)*(u(i)-u(i-1)))/c_p)+T(i-1);
    P(i) = rho(i)*R*T(i);
    
    M(i) = u(i)/(sqrt(gamma*R*T(i)));

percent = (i/n) * 100;
waitbar(percent/100,h,sprintf('%.2f%% of Isentropic Flow Calculation Complete',percent))

    i = i+1;

end

close(h); %Close the wait bar

Pfinal = P(n);
Tfinal = T(n);
ufinal = u(n);
rhofinal = rho(n); 
Mfinal = M(n);
Ptfinal = Pfinal*((1+((y-1)/2)*Mfinal^2)^(y/(y-1)));
Ttfinal = Tfinal*(1+((y-1)/2)*Mfinal^2);
Forcefinal = (rhofinal*ufinal^2+Pfinal)*Ae - (rho0*u0^2+P0)*A0;

fprintf('------------------------BEGIN-OUTPUT-------------------------------\n ' )
fprintf('          ISENTROPIC FLOW WITH NO WORK INTERACTIONS \n')
fprintf('Analyitical Solution \n')
fprintf('Pressure          : %f Pa \n',Pe)
fprintf('Total Pressure    : %f Pa \n',Pt2)
fprintf('Temperature       : %f K \n',Te)
fprintf('Total Temperature : %f K \n',Tt2)
fprintf('Mach #            : %f \n',M2)
fprintf('Velocity          : %f m/s \n',ue)
fprintf('Force             : %f N \n',Force)
fprintf('Density           : %f kg/m^3 \n \n',rho2)

fprintf('Numerical Solution \n')
fprintf('Pressure          : %f Pa \n',Pfinal)
fprintf('Total Pressure    : %f Pa \n',Ptfinal)
fprintf('Temperature       : %f K \n',Tfinal)
fprintf('Total Temperature : %f K \n',Ttfinal)
fprintf('Mach #            : %f \n',Mfinal)
fprintf('Velocity          : %f m/s \n',ufinal)
fprintf('Force             : %f N \n',Forcefinal)
fprintf('Density           : %f kg/m^3 \n \n',rhofinal)

fprintf('--------------------------END-OUTPUT--------------------------------\n \n' )
