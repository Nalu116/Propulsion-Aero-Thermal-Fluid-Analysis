clear variables

%{
Fanno Flow --> Adiabatic, no work interaction, no area variation 
%} 

n = 100000; %Axial step size
%n = 100; %Temporary Axial Step Size for defense against aging

Cf = 0.05;
Cp = 1004.5;
P0 = 101325;
T0 = 288;
M0 = 0.2;
A = 0.1;
L = 10;
dx = L/n;
y = 1.4;
R = 287;
u0 = M0*sqrt(R*y*T0);
Pt0 = P0*((1+((y-1)/2)*M0^2)^(y/(y-1)));
Tt0 = T0*(1+((y-1)/2)*M0^2);
rho0 = P0/(R*T0);
M1 = M0;
P1 = P0;
Tt1 = Tt0;
Tt2 = Tt1; %Adiabatic with no work interactions -> Total Temp unchanged

%ASSUMES CIRCULAR CROSS SECTION
%Area = pi * r^2
%c = 2*pi*r
%r = sqrt(A/PI) then ... 
c = 2*pi*sqrt(A/pi);    

%Analyitical Fanno Flow

fvM1 = ((y+1)/2)*log((1+((y-1)*(M1^2)/2))/(M1^2))-(1/(M1^2));

syms M2q
fvM2 = ((y+1)/2)*log((1+((y-1)*(M2q^2)/2))/(M2q^2))-(1/(M2q^2));
eqn0 = y*Cf*(c/A)*L == fvM2-fvM1;
sol4 = vpasolve(eqn0,M2q,[0 1]);
M2 = double(sol4);
Pe = P1*M1/M2*sqrt( (1+(y-1)/2*M1^2)/(1+(y-1)/2*M2^2));
Pt2 = Pe*((1+((y-1)/2)*M2^2)^(y/(y-1)));
Te = Tt2/(1+(y-1)/2*M2^2);
ue = M2 * sqrt(y*R*Te);
rho2 = Pe/(R*Te);
Force = (rho2*ue^2+Pe)*A - (rho0*u0^2+P0)*A;
%Small increase in mach #, but Cf was very small so kind of makes sense

%Numerical 
i = 1;
%local conditions of array begin at initial conditions prior to steps
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

while i < (n+1)

Tauw = 0.5*Cf*rhomat(i)*umat(i)^2;    
syms duq drhoq dTq dPq

eqn1 = dTq == -(umat(i)*duq)/Cp; 
eqn2 = dPq == - rhomat(i)*(umat(i)*duq + (c*dx*Tauw)/(A*rhomat(i)));
eqn3 = drhoq == (duq/umat(i))*-rhomat(i); 
eqn4 = dPq == Pmat(i)*((drhoq/rhomat(i))+(dTq/Tmat(i))); 
eqnmat3 = [eqn1,eqn2,eqn3,eqn4];
sol4 = vpasolve(eqnmat3);
dP = sol4.dPq;
dT = sol4.dTq;
drho = sol4.drhoq;
du = sol4.duq;

        umat(i+1) = umat(i) + du;
        Tmat(i+1) = Tmat(i) + dT;
        Pmat(i+1) = Pmat(i) + dP; 
        rhomat(i+1) = rhomat(i) + drho;
        Machmat(i+1) = umat(i+1)/sqrt((y*R*Tmat(i+1)));

percent = (i/n) * 100;
waitbar(percent/100,h,sprintf('%.2f%% of Fanno Flow Calculation Complete',percent))

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
Forcefinal = (rhofinal*ufinal^2+Pfinal)*A - (rho0*u0^2+P0)*A;

fprintf('-------------------------BEGIN-OUTPUT------------------------------- \n' )
fprintf('     FANNO FLOW - ADIABATIC NO WORK INTERACTIONS CONST AREA   \n')
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
fprintf('-----------------------END-PART-OUTPUT-----------------------------\n \n' )
