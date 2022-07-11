# gpe-yoshida
Yoshida Implementation of Coupled Gross-Pitaevskii Equations

This code solves the coupled GPEs describing dilute gas cold atom Bose Einstein condensates.
Facilities are provided for doing both real-time evolution and obtaining the ground state (via imaginary time evolution).

## Real-Time Evolution

The real-time evolution equations are
$$i\hbar\dot\psi_1 = \left[-\frac{\hbar^2}{2m}\nabla^2 + g_{11}|\psi_1|^2 + g_{12}|\psi_2|^2 + V(x) - \mu  \right]\psi_1 - \nu\psi_2$$

$$i\hbar\dot\psi_2 = \left[-\frac{\hbar^2}{2m}\nabla^2 + g_{22}|\psi_2|^2 + g_{21}|\psi_1|^2 + V(x) - \mu \right]\psi_2 - \nu\psi_1$$

for user specifiable parameters $g_{ij}$ and 
$\nu$, and a user specified external trapping potential $V(x)$. 
The constraint $g_{21}=g_{12}$ is assumed to hold.
The chemical potential $\mu$ can either be specified by hand, or computed explicitly from the initial state.

## Ground State Solution

In addition to solving the dynamical equation, it can also compute the ground state for a given user specifed $V(x)$ by solving the coupled 
nonlinear eigenvalue problem

$$\mu\dot\psi_1 = \left[-\frac{\hbar^2}{2m}\nabla^2 + g_{11}|\psi_1|^2 + g_{12}|\psi_2|^2 + V(x) \right]\psi_1 - \nu\psi_2$$

$$\mu\dot\psi_2 = \left[-\frac{\hbar^2}{2m}\nabla^2 + g_{22}|\psi_2|^2 + g_{21}|\psi_1|^2 + V(x) \right]\psi_2 - \nu\psi_1$$

for the ground state.
