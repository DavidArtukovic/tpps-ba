# Model Overview – Thermal Pumped Piston Storage (TPPS)

This document provides an overview of the physical and mathematical TPPS model.


## 1. Domains
- Water
- Solid (piston, soil, insulation)
- Penstock / inlet pipe

## 2. Governing Equations

### 2.1 Water Domain
One dimensional heat transport equation **with** residual terms: 

$$
\rho_w \, c_w \, \frac{\partial T_w(z; t)}{\partial t}
= \lambda_w \, \frac{\partial^2 T_w(z; t)}{\partial z^2}
- \rho_w \, v_{z} \, c_w \, \frac{\partial T(z; t)}{\partial z}
+ q_{\text{conv, free}}
+ \sigma_{\text{piston}}
+ \sigma_{\text{wall}}
$$

Terms:
- $-\rho_w \, v_{z} \, c_w \, \partial T / \partial z$ → forced convection (if present). Shifts temperature profile to in accordance with sign of vertical (axial) velocitcy $v_z$. Is active between inlet and outlet.
- $\lambda_w \frac{\partial^2 T_w(z,t)}{\partial z^2}$ → diffusive heat transport 
- $q_{\text{conv, free}}$ → natural convection term (to be developed in the thesis)  
- $\sigma_{\text{piston}}$, $\sigma_{\text{wall}}$ → coupling terms to piston and soil (see Häuslein, 2024)


### 2.2 Solid Domain
One dimensional heat transport equation in **solid** systems:

$$
\rho_s \, c_s \frac{\partial T_s(z; t)}{\partial t} = \lambda_s \frac{\partial^2 T_s(z; t)}{\partial z^2}
$$

Two dimensional (axial=z and radial=r) heat transport equation

$$
\rho_s \, c_s \frac{\partial T_s(z, r; t)} {\partial t}
=
\lambda_s \left(
\frac{\partial^2 T_s(z,r; t)}{\partial z^2}
+
\frac{1}{r}\,\frac{\partial}{\partial r}
\left(
r\,\frac{\partial T_s(z, r; t)}{\partial r}
\right)
\right)
$$

### 2.3 Coupling Terms (σ-terms)
Coupling terms represent heat transfer between adjacent domains.
They appear as source/sink terms in each energy balance equation.


## 3. Extension: Natural Convection (Thesis Objective)
Natural convection will be introduced into the 1D water equation.
 Approach: No fluid transport via a flow but energy entrance per volume element, per incremental time unit.

### 3.1 Energy Balance Approach for Natural Convection (Gerle, Schäfer)
Thermal energy per volume to arbitrary reference point, differentiated by time under the assumption of a an incompressible fluid:
$$
\begin{aligned}
& q_{\text{conv, free}} = \rho_w \, c_w \, T \\
\iff & \frac{\partial}{\partial t} \,q = \frac{\partial}{\partial t} \,[\rho_w \, c_w T_w] \\
\iff & \dot{q}_{\text{conv, free}} = \rho_w \, c_w \,
\frac{\partial T_w}{\partial t}
\end{aligned}
$$
which corresponds to equation (3) in Schäfer (2021). 

#### Linear Approach for the Mixed Temperature

During a mixing event (natural convection), the new temperature
distribution in the mixing zone $ z \in [z_{\text{in}}, z_{\text{mix}}] $
is assumed to be linear:

$$
T_w^*(z) = a z + b
\tag{1}
$$
Note as well that $T_w^* = T_w(z, t+dt)$ and is thus the temperature distribution after an infinitesimal timestep $dt$. 

---

#### Determination of Intercept $b$ via Boundary Condition
The boundary condition at the top of the mixing zone is given by:

$$
T_w^*(z_{mix}) = T_{w,\text{in}}
\tag{2}
$$

Using (1) in (2):

$$
a z_{\text{mix}} + b = T_{w,\text{in}}
\quad\Rightarrow\quad
b = T_{w,\text{in}} - a z_{\text{mix}}
\tag{3}
$$

Insert (3) into (1):

$$
T_w^*(z)
= T_{w,\text{in}} + a (z - z_{\text{mix}})
\tag{4}
$$

---

#### Integration of the New Temperature Profile

We integrate (4) over the mixing zone:

$$
\int_{z_{\text{in}}}^{z_{\text{mix}}} T_w^*(z)\,dz
=
\int_{z_{\text{in}}}^{z_{\text{mix}}}
\left[ T_{w,\text{in}} + a(z - z_{\text{mix}}) \right] dz
\tag{5}
$$

The integral of the linear term is:

$$
\int_{z_{\text{in}}}^{z_{\text{mix}}} (z - z_{\text{mix}})\,dz
=
-\,\frac{(z_{\text{mix}} - z_{\text{in}})^2}{2}
\tag{6}
$$

We define:

$$
I := \int_{z_{\text{in}}}^{z_{\text{mix}}}
\left( T_{w,\text{in}} - T_w(z) \right) dz
\tag{7}
$$

---

#### Energy Balance to Determine the Slope $a$

The energy added by the inflowing hot water must equal the
internal energy increase of the mixing region:

$$
\begin{aligned}
\Delta U_{w,\text{mix}}
=
& \rho_w c_w A_{\text{hws}}
\int_{z_{\text{in}}}^{z_{\text{mix}}}
\left( T_w^*(z) - T_w(z) \right) dz \\
\iff 
& \rho_w c_w A_{\text{hws}}
\int_{z_{\text{in}}}^{z_{\text{mix}}}
\left( T_{w,\text{in}} + a (z - z_{\text{mix}}) - T_w(z) \right) dz
=
\dot Q_{\text{mix}} \, dt
\tag{8}
\end{aligned}
$$
where the last formulation was arrived by plugging in equation (4). 
Now use (6)–(7) to solve the remaining integral:

$$
\rho_w c_w A_{\text{hws}}
\left[
I - a\frac{(z_{\text{mix}} - z_{\text{in}})^2}{2}
\right]
=
\dot Q_{\text{mix}}\, dt
\tag{9}
$$

Solve for $a$:

$$
a
=
\frac{2}{(z_{\text{mix}} - z_{\text{in}})^2}
\left[
I - \frac{\dot Q_{\text{mix}}\, dt}
       {\rho_w c_w A_{\text{hws}}}
\right]
\tag{10}
$$
Now both constants $[a, b]$ are defined by the boundary condition and the energy equation, respectively.

---


#### Update of Linear Approach

$$
T_w^*(z) = T_{w,\text{in}} + a(z - z_{\text{mix}}),
$$

we now insert the expression for $a$ and $b$ obtained in Equation (10) and (3):
$$
T_w^*(z)
=
T_{w,\text{in}}
+
\frac{2(z - z_{\text{mix}})}{(z_{\text{mix}} - z_{\text{in}})^2}
\left[
I - \frac{1}{\rho_w c_w A_{\text{hws}}}\,\dot Q_{\text{mix}} \, dt
\right].
\tag{11}
$$

This expresses the updated temperature field after one infinitesimal mixing step in terms of the integral $I$, the inflow energy $\dot Q_{\text{mix}} dt$, and the geometric mixing factor.

---

#### Temperature Increase and Definition of the Mixing Function

The incremental temperature increase in the mixing zone is:

$$
\begin{aligned}
\Delta T(z)
&= T_w^*(z) - T_w(z) \\
&=T_{w,\text{in}} + \frac{2(z - z_{\text{mix}})}{(z_{\text{mix}} - z_{\text{in}})^2}
\left[
I - \frac{1}{\rho_w c_w A_{\text{hws}}}\,\dot Q_{\text{mix}} \, dt
\right]
- T_w(z)
\tag{12}
\end{aligned}
$$

where equation (11) has been plugged in.

---

#### Final Mixing Source Term (Equation 8 in Schäfer 2021)

Using Equation (3):

$$
\sigma_{\text{mix}}(z) = \rho_w c_w \frac{\Delta T(z)}{dt}
\tag{13}
$$

Insert (12) and (13) to obtain:

$$
\boxed{
q_{\text{conv, free}}(z) =
\frac{\rho_w c_w}{dt}
\left(
\int_{z_{\text{in}}}^{z_{\text{mix}}}
\frac{2(z-z_{\text{mix}})}{(z_{\text{mix}}-z_{\text{in}})^2}
\left(T_{w,\text{in}} - T_w(z)\right)\,dz
\right)
+
\frac{\rho_w c_w}{dt}(T_{w,\text{in}} - T_w)
-
\dot Q_{\text{mix}}
\frac{2(z - z_{\text{mix}})}
     {A_{\text{hws}}(z_{\text{mix}} - z_{\text{in}})^2}
}
\tag{15}
$$

which matches Equation (8) from Schäfer et al. (2021). Note, the unit of the convecitve term is given by: $[q_{\text{conv, free}}(z)] = \frac{W}{m^3}$. Thus, we speak of a power input per volume.


The heat flow is given by: 
 $$
 \dot{Q}_{mix} = \dot{M_w}\,c_w\,(T_{w,in}-T_{w}(z_{in}))
 $$

### 3.2. Variable Mixed Zone
The mixed zone is in the present model of dynamic size, $[z_{in}, z_{mix}(t)]$. This means in practice, every iteration step, on has to find the z value for which the temperature is equal to the inlet temperature:
 $$
 T_{w}(z_{mix}(t), t) = T_{w, in}(t)
 $$
Thus the number of cells, affected by natural convection is variable. 


## 4. TODOs

- Derive the free convection term for TPPS geometry.  
- Review Nusselt–Rayleigh correlations for low-velocity water in vertical storage.  
- Add piston–water and wall–water σ-term details.  
- Perform parameter sensitivity analysis once the extension is implemented.

