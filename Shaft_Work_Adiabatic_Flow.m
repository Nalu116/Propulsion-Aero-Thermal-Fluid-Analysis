clear variables

%{
Positive shaft work interactions to adiabatic flow with no wall friction
Case A - Isentropic (Eta = 1.0)
Case B - Efficiency of 95% (Eta = 0.95)
%}

n = 100000; %Axial step size
%n = 100; %Temporary Axial Step Size to lower abrasiveness to the sands of time during iteration

W = 0.8e6; %Total work per mass (J/kg)
dw = W/n;

%Define efficiency
Eta = 1.0;  %Case A
%Eta = 0.95; %Case B

P0 = 101325;
T0 = 288;
M0 = 0.5;
A0 = 0.1;
Ae = 0.002;
L = 1;
y = 1.4;
R = 287;
Cp = 1004.5;
rho0 = P0/(R*T0); 
u0 = M0*(y*R*T0)^0.5;
mdot = rho0*u0*A0;
dA = (Ae-A0)/n;
Tt0 = T0*(1+(((y-1)/2)*M0^2));
Pt0 = P0*((1+((y-1)/2)*M0^2)^(y/(y-1)));

%---Analyitical Solution------
Tt2 = Tt0 + W/Cp;
Pt2 = Pt0*(Tt2/Tt0)^(y/(y-1));

syms M2q
eqn1 = mdot == Pt2*((y/R)^0.5)*(1/(Tt2)^0.5)*M2q*Ae*(1+(((y-1)*M2q^2)/2))^-((y+1)/(2*(y-1)));

sol1 = double(vpasolve(eqn1,M2q,[0 1]));
M2 = sol1; %Take the real value that is subsonic in this case (diverging)

Pe = Pt2/(((1+((y-1)/2)*M2^2)^(y/(y-1))));
Te = Tt2/((1+((y-1)/2)*M2^2));
rho2 = Pe/(Te*R);
ue = (u0*A0*rho0)/(Ae*rho2);
Force = (rho2*ue^2+Pe)*Ae - (rho0*u0^2+P0)*A0;

%Numerical --------
i = 1;
%local conditions of array begin at initial conditions prior to steps
Amat = linspace(A0,Ae,n);
Tmat = zeros(1, n);
Pmat = zeros(1, n);
umat = zeros(1, n);
rhomat = zeros(1, n);
Machmat = zeros(1, n);
Tmat(i) = T0;
Pmat(i) = P0;
umat(i) = u0;
rhomat(i) = rho0;
Machmat(i) = M0;

h = waitbar(0,'Initializing waitbar...');
i = 1;

while i < (n+1)
            
du = - (dA/Amat(i) + (rhomat(i)*(dw*Eta))/Pmat(i) - (dw)/(Tmat(i)*Cp))/(1/umat(i) + umat(i)/(Tmat(i)*Cp) - (rhomat(i)*umat(i))/Pmat(i));
drho =  - rhomat(i)*(dA/Amat(i) + du/umat(i)); 
dP =  -rhomat(i)*(du*umat(i) - dw*Eta);
dT = (dw - du*umat(i))/Cp;

        umat(i+1) = umat(i) + du;
        Tmat(i+1) = Tmat(i) + dT;
        Pmat(i+1) = Pmat(i) + dP; 
        rhomat(i+1) = rhomat(i) + drho;
        Machmat(i+1) = umat(i+1)/sqrt((y*R*Tmat(i+1)));

percent = (i/n) * 100;
waitbar(percent/100,h,sprintf('%.2f%% of Positive Shaft Work Calculation Complete',percent))

        i = i+1;
end

close(h); %Close waitbar

Pfinal = Pmat(n);
Tfinal = Tmat(n);
ufinal = umat(n);
rhofinal = rhomat(n); 
Mfinal = Machmat(n);
Ptfinal = Pfinal*((1+((y-1)/2)*Mfinal^2)^(y/(y-1)));
Ttfinal = Tfinal*(1+((y-1)/2)*Mfinal^2);
Forcefinal = (rhofinal*ufinal^2+Pfinal)*Ae - (rho0*u0^2+P0)*A0;

fprintf('-------------------------BEGIN-OUTPUT------------------------------- \n' )
fprintf('      POSITIVE SHAFT WORK INTERACTIONS TO ADIABATIC FLOW  \n')
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
fprintf('------------------------END-OUTPUT---------------------------\n \n' )
