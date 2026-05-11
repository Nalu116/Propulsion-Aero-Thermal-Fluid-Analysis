clear variables

%{
Convective heat transfer in duct with no area variation
%}

n = 100000; %Axial step size
%n = 200; %Temporary Axial Step Size to run before retiring
dA = 0;
L = 1;
dx = L/n;
y = 1.4;
R = 287;
Cp = 1004.5;
Cf = 0.08;
P0 = 101325;
T0 = 288;
M0 = 0.1;
A = 0.1;
u0 = M0*sqrt(R*T0*y);
rho0 = P0/(R*T0);
Tw = 3500;
Pt0 = P0*((1+((y-1)/2)*M0^2)^(y/(y-1)));
Tt0 = T0*(1+((y-1)/2)*M0^2);

%A = pi*r^2
%A = pi * (d/2)^2
%sqrt(A/pi) = D/2
D = 2*sqrt(A/pi);
c = 2*pi*sqrt(A/pi); 

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
i = 1;

while i < (n+1)
               
Tauw = 0.5*Cf*rhomat(i)*umat(i)^2;   
Tti = Tmat(i)*(1+((y-1)/2)*Machmat(i)^2); %Local Total Temp variable (changing)
dq = 2*Cp*Cf*(Tw-Tti)*(dx/D);

syms duq drhoq dTq dPq

eqn1 = dTq == (dq + (umat(i)*duq))/Cp; 
eqn2 = dPq == - rhomat(i)*(umat(i)*duq + (c*dx*Tauw)/(A*rhomat(i)));
eqn3 = drhoq == (duq/umat(i))*-rhomat(i); 
eqn4 = dPq == Pmat(i)*((drhoq/rhomat(i))+(dTq/Tmat(i))); 
eqnmat4 = [eqn1,eqn2,eqn3,eqn4];
sol4 = vpasolve(eqnmat4);
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
waitbar(percent/100,h,sprintf('%.2f%% of Convective Heat Transfer Calculation Complete',percent))

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
HeatTransfer = Cp*(Ttfinal-Tt0);


fprintf('-------------------------BEGIN-OUTPUT------------------------------- \n' )
fprintf('        CONVECTIVE HEAT TRANSFER CONST AREA DUCT   \n')

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
fprintf('------------------------END-OUTPUT---------------------------\n \n' )
