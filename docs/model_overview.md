# Model Overview – Thermal Pumped Piston Storage (TPPS)

This document provides an overview of the physical and mathematical TPPS model.

---

## 1. Domains
- Water 
- Solid (piston, soil, insulation)
- Penstock / inlet pipe

---

## 2. Governing Equations

### 2.1 Water Domain

One dimensional heat transport equation **with** residual terms: 

$$
\rho_w \, c_w \frac{\partial T_w(z; t)}{\partial t}
= \lambda_w \frac{\partial^2 T_w(z; t)}{\partial z^2}
- \dot{m} c_w \frac{\partial T(z; t)}{\partial z}
+ q_{\text{conv, free}}
+ \sigma_{\text{piston}}
+ \sigma_{\text{wall}}
$$

Terms:
- $-\dot{m} c_w \, \partial T / \partial z$ → forced convection (if present)  
- $\lambda_w \frac{\partial^2 T_w(z,t)}{\partial z^2}$ → fiffusive heat transport 
- $q_{\text{conv, free}}$ → natural convection term (to be developed in the thesis)  
- $\sigma_{\text{piston}}$, $\sigma_{\text{wall}}$ → coupling terms to piston and soil  

---

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

---

### 2.3 Coupling Terms (σ-terms)

Coupling terms represent heat transfer between adjacent domains.
They appear as source/sink terms in each energy balance equation.

Typical form:

$$
\sigma = \frac{\lambda A}{V}\left( T_{\text{neighbor}} - T \right)
$$

Where:

- $\lambda$ = thermal conductivity  
- $A$ = contact area  
- $V$ = control volume  
- $(T_{\text{neighbor}} - T)$ = temperature difference  

---

## 3. Extension: Natural Convection (Thesis Objective)

Natural convection will be introduced into the 1D water equation.
Candidate modelling approaches:

### 3.1 Rayleigh–Nusselt Correlation (buoyancy-based)

Estimate local effective mixing / enhanced heat transfer:

- Compute Rayleigh number:
  
  $$
  \mathrm{Ra} = \frac{g \beta (T - T_\text{ref}) (L)^3}{\nu \alpha}
  $$

- Use a suitable Nusselt correlation for vertical cylindrical domains:

  $$
  \mathrm{Nu} = C \, \mathrm{Ra}^n
  $$

- Convert to an effective convective term:

  $$
  q_{\text{conv, free}} = h A (T_\text{bulk} - T)
  $$

---

### 3.2 Mixing-Term Approach (Gerle, Schäfer)

Alternatively, model buoyancy-driven recirculation as a mixing velocity:

$$
q_{\text{conv, free}} = \rho c_p \, v_{\text{mix}} \,
\frac{\partial T}{\partial x}
$$

Where $v_{\text{mix}}$ is computed from stratification or temperature gradients.

---

## 4. TODOs

- Derive the free convection term for TPPS geometry.  
- Review Nusselt–Rayleigh correlations for low-velocity water in vertical storage.  
- Add piston–water and wall–water σ-term details.  
- Perform parameter sensitivity analysis once the extension is implemented.

