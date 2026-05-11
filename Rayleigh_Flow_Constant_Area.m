clear variables
%{ 
-Rayleigh Flow
-Inviscid no work interaction no area variation
%}

n = 100000;  %Number of itterations (Fidelity step size)
%n = 100; %Temp itteration so dont wait 30 years for entire code

P0 = 101325;
T0 = 288;
M0 = 0.2;
A = 0.1;
dA = 0; %Constant Area duct
L = 1;
y = 1.4;
R = 287;
a0 = sqrt(y*R*T0);
Cp = 1004.5;
rho0 = P0/(R*T0);
q = 1.2e6;
u0 = M0*a0;
mdot = rho0*u0*A;
Pt0 = P0*((1+((y-1)/2)*M0^2)^(y/(y-1)));
Tt0 = T0*(1+((y-1)/2)*M0^2);
Pt1 = Pt0;
Tt1 = Tt0;

%----Analyitical Solution---(Equation referenced from Reyleigh flow from previous propulsion lecture Godbless Dr. Riggins + Dr. Abbas)
fM1 = (M0^2)*((1+((y-1)/2)*M0^2)/(1+(y*M0^2))^2);

syms M2q

fM2 = (M2q^2)*((1+((y-1)/2)*M2q^2)/(1+(y*M2q^2))^2);
eqn1 = 1+ (q/(Cp*Tt1)) == fM2/fM1;

sol1 = vpasolve(eqn1, M2q);
M2 = double(sol1(3));
Tt2 = Tt1 + q/Cp;
Te = Tt2/(1+(((y-1)/2)*M2^2));

syms Pt2q

eqn2 = mdot == ((Pt2q*(sqrt(y/R)))/(sqrt(Tt2)))*M2*A*(1+(((y-1)*M2^2)/2))^-((y+1)/(2*(y-1)));
sol2 = vpasolve(eqn2,Pt2q);
Pt2 = double(sol2);
Pe = Pt2/((1+((y-1)/2)*M2^2)^(y/(y-1)));
rho2 = Pe/(R*Te);
ue = M2 * sqrt(y*R*Te);
Force = (rho2*ue^2+Pe)*A - (rho0*u0^2+P0)*A; %A const;

% ---------- Numerical Itterations---------------------------------------
dq = q/n;    %Rayleigh flow --> Cf -> 0 so only the additive initial term
i = 1;       %Initialize counter and starting position of array elements

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

h = waitbar(0,'Initializing waitbar...'); %Loading bar to track progress

while i < (n+1)
               
syms duq drhoq dTq dPq

eqn3 = drhoq == -rhomat(i) * (duq/umat(i)); %dA -> 0 thus dA/A term nulled out
eqn4 = dPq == Pmat(i)*((drhoq/rhomat(i))+(dTq/Tmat(i)));
eqn5 = dTq == ((dq - (umat(i)*duq))/Cp); 
eqn6 = dPq == -rhomat(i)*(umat(i)*duq);
eqnmat2 = [eqn3,eqn4,eqn5,eqn6];

sol3 = vpasolve(eqnmat2); %Pressure, density decreasing, velocity, temp up
dP = sol3.dPq;
dT = sol3.dTq;
drho = sol3.drhoq;
du = sol3.duq;

        umat(i+1) = umat(i) + du;
        Tmat(i+1) = Tmat(i) + dT;
        Pmat(i+1) = Pmat(i) + dP; 
        rhomat(i+1) = rhomat(i) + drho;
        Machmat(i+1) = umat(i+1)/sqrt((y*R*Tmat(i+1)));

percent = (i/n) * 100;
waitbar(percent/100,h,sprintf('%.2f%% of Rayleigh Flow Calculation Complete',percent))

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

fprintf('-------------------------BEGIN-OUTPUT--------------------------------- \n' )
fprintf(' RAYLEIGH FLOW - INVISCID NO WORK INTERACTIONS NO AREA VARIATION \n')
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
fprintf('-----------------------END-OUTPUT-------------------------------\n \n' )


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
