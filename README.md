The scripts found in this repo are various levels of fidelity for aero-thermal-fluid analysis (specifically with propulsive use cases in mind, though any thermo-fluid use case that utilized quasi-1D flow and area/mach relations could find this useful).

Initially I utilized MATLAB’s built in syms functionality to numerically solve various systems of equations. However as this became increasingly computationally intensive, it became more efficient (and accurate) to manually define differential elements as functions of the other parameters either by hand or using equation solvers, so some of the first files made are _slighty_ less accurate (notably the initial quasi-1D file). For best results, use the "Combined_Positive_Shaft_Work_Skin_Friction_Convection.m" script and update any initial conditions manually to capture your specific use case.

For all cases, Quasi 1-D equations ( 1 - 4 ) were manipulated to determine the local differential elements (differential velocity, differential temperature, etc). Based on the particular problem and case, certain
variables drop out of the equations. For example, an inviscid flow would have the entire fraction term go to zero in equation 2 due to the lack of shear stress in an idealized reversible fluid flow. Fanno and Rayleigh flow use nested functions defined as fM and fvM to reduce clutter. For further reference, equations manipulated in these scripts are defined below:

Quasi 1-D Equations: 

$$
\begin{align}
dT = \frac{dq + dw + udu}{C_p} \quad \text{(1)}
\end{align}
$$

$$
\begin{align}
dP = −ρ(udu − ηdw +  \frac{cdxτ_w}{Aρ} ) \quad \text{(2)}
\end{align}
$$

  
$$
\begin{align}
dP =  \frac{dρ}{ρ} + \frac{du}{u} + \frac{dA}{A} = 0 \quad \text{(3)}
\end{align}
$$

  
$$
\begin{align}
dP =  \frac{dP}{P} + \frac{dρ}{ρ} + \frac{dT}{T} = 0 \quad \text{(4)}
\end{align}
$$

  
Rayleigh flow fM functions to solve Mach Number:

$$
\begin{align}
dP =  1 + \frac{q}{C_pT_t} = \frac{fM_2}{fM_1} \quad \text{(5)}
\end{align}
$$


Where the fM function to find Mach of interest is defined as :
<p align="center">
$fM_n = M_2 * \dfrac{\frac{1+γ−1}{2}M^2}{{1 + (γM^2)^2}} \quad \text{(6)}$


Fanno flow functions for Mach Number:
<p align="center">
$γC_f (\frac{c}{A}) L = fvM_2 − fvM_1 \quad \text{(7)}$


Where the fvM function is defined as :
<p align="center">
$fvM_n = \dfrac{\frac{γ−1}{2} ln (1 + (γ − 1) \frac{M^2}{2} } {M^2 - \frac{1}{M^2} } \quad \text{(8)}$


Fluid wall shear stress:

$$
\begin{align}
τ_w = \dfrac{1}{2}C_f ρu^2 \quad \text{(9)}
\end{align}
$$

Differential heating element:

$$
\begin{align}
dq = 2C_pC_f (T_w − T_t)(\frac{dx}{D}) \quad \text{(10)}
\end{align}
$$


Mass flow rate as a function of Mach and local state variables:
<p align="center">
$\dot{m} = P_{t2} \sqrt{ (\dfrac{γ}{R})} (\dfrac{1}{T_{t2}^2}) M_2A_e (1 + \dfrac{(γ − 1)M_2^2}{2}) ^ {( \dfrac{γ+1}{2*(y+1)} )}  \quad \text{(11)}$

