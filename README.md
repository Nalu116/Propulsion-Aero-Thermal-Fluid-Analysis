The scripts found in this repo are various levels of fidelity for aero-thermal-fluid analysis (specifically with propulsive use cases in mind, though any thermo-fluid use case that utilized quasi-1D flow and area/mach relations could find this useful).

Initially I utilized MATLAB’s built in syms functionality to numerically solve various systems of equations. However as this became increasingly computationally intensive, it became more efficient (and accurate) to manually define differential elements as functions of the other parameters either by hand or using equation solvers.

For all cases, Quasi 1-D equations ( 1 - 4 ) were manipulated to determine the local differential elements (differential velocity, differential temperature, etc). Based on the particular problem and case, certain
variables drop out of the equations. For example, an inviscid flow would have the entire fraction term go to zero in equation 2 due to the lack of shear stress in an idealized reversible fluid flow.

Quasi 1-D Equations: 

$$
\begin{align}
dT = \frac{dq + dw + udu}{Cp} \quad\text{(1)}
\end{align}
$$

$$
\begin{align}
dP = −ρ(udu − ηdw +  \frac{cdxτw}{Aρ} ) \quad\text{(2)}
\end{align}
$$

$$
\begin{align}
dP =  \frac{dρ}{ρ} + \frac{du}{u} + \frac{dA}{A} = 0 \quad\text{(3)}
\end{align}
$$

$$
\begin{align}
dP =  \frac{dP}{P} + \frac{dρ}{ρ} + \frac{dT}{T} = 0 \quad\text{(4)}
\end{align}
$$

Rayleigh flow fM functions of Mach Number:

$$
\begin{align}
dP =  1 + \frac{q}{CpTt} = \frac{fM2}{fM1} \quad\text{(5)}
\end{align}
$$



