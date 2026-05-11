clear variables 

%{
    Positive Shaft work interaction with convection and skin friction
%}

n = 100000; %Axial step size
%n = 500; %Temporary Axial Step Size to save the lives of RAM cards everywhere

Eta = 0.9;

y = 1.4;
R = 287;
Cp = 1004.5;
W = 2.5e6;

dw = W/n;
L = 1;
dx = L/n;
Tw = 500;
Cf = 0.02;
P0 = 20000;
T0 = 1000;
M0 = 4.0;
A0 = 0.1;
Ae = 0.3;
u0 = M0*sqrt(y*R*T0);
rho0 = P0/(R*T0); 
dA = (Ae - A0)/n;
Tt0 = T0*(1+(((y-1)/2)*M0^2));
Pt0 = P0*((1+((y-1)/2)*M0^2)^(y/(y-1)));

%-----Numerical Solution-----
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

D = 2*sqrt(Amat(i)/pi);
c = 2*pi*sqrt(Amat(i)/pi); 

Tauw = 0.5*Cf*rhomat(i)*umat(i)^2;   
Tti = Tmat(i)*(1+((y-1)/2)*Machmat(i)^2); %Local Total Temp variable (changing)

dq = 2*Cp*Cf*(Tw-Tti)*(dx/D);
du = - (dA/Amat(i) + (rhomat(i)*(dw*Eta - (c*dx*Tauw)/(Amat(i)*rhomat(i))))/Pmat(i) - (dq + dw)/(Tmat(i)*Cp))/(1/umat(i) + umat(i)/(Tmat(i)*Cp) - (rhomat(i)*umat(i))/Pmat(i));
drho =  - rhomat(i)*(dA/Amat(i) + du/umat(i));
dP =  - rhomat(i)*(du*umat(i) - dw*Eta + (c*dx*Tauw)/(Amat(i)*rhomat(i)));
dT = (dq + dw - du*umat(i))/Cp;

        umat(i+1) = umat(i) + du;
        Tmat(i+1) = Tmat(i) + dT;
        Pmat(i+1) = Pmat(i) + dP; 
        rhomat(i+1) = rhomat(i) + drho;
        Machmat(i+1) = umat(i+1)/sqrt((y*R*Tmat(i+1)));

percent = (i/n) * 100;
waitbar(percent/100,h,sprintf('%.2f%% of Combined Calculation Complete',percent))

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
Forcefinal = -((rhofinal*ufinal^2*Ae+Pfinal*1*Ae)+(rho0*u0*-u0*A0+P0*-1*A0));
HeatTransfer = Cp*(Ttfinal-Tt0)-W;

fprintf('-------------------------BEGIN-OUTPUT------------------------------- \n' )
fprintf('POSITIVE SHAFT WORK TO FLOW WITH CONVECTION AND WALL FRICTION FLOW  \n')
fprintf('Numerical Solution \n')
fprintf('Pressure          : %f Pa \n',Pfinal)
fprintf('Total Pressure    : %f Pa \n',Ptfinal)
fprintf('Temperature       : %f K \n',Tfinal)
fprintf('Total Temperature : %f K \n',Ttfinal)
fprintf('Mach #            : %f \n',Mfinal)
fprintf('Velocity          : %f m/s \n',ufinal)
fprintf('Force             : %f N \n',Forcefinal)
fprintf('Density           : %f kg/m^3 \n',rhofinal)
fprintf('Heat Transfer     : %f J \n \n',HeatTransfer)
fprintf('-----------------------END-OUTPUT----------------------------\n \n' )
