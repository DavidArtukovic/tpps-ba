#  01 - Code Info 
This notebook collects specific information on code junks, calculations within the matlab files and related stuff.

## scenario_1/Init.m

### Translation of Dome Zone into Cylinder
At the upper and lower end are dead zones (1m) and domes of height one meter. This volume is treated as it would be within a cylinder, thus effectively reducing the **vertical number of gridpoints** of the dome.
```matlab
h_tot = 1;                              % upper dead-zone
V_Kupp = (2/3)*r_ST^2*h_kupp*pi;        % volume upper half dome in m³
V_PP = h_in*d_in^2/4*pi;                % small volume (5cm height) of collision palte
h_tot_2 = (V_Kupp - V_PP)/(d_ST^2*pi/4) % height of dome volume if treated as cylinder
h_1D_tot = h_tot + h_tot_2;             % Total 1D height of dead-zone.
```
Thu, instead of two meters deadzone the model assumes 1.665m.

## scenario_1/Szenario1.m

### Translation of Scenario into Information on Upper/Lower Pressure Zone

The helper function `prepare_scenario_flow` preprocesses the 15-minute scenario schedule and converts it into all hydraulic quantities required for the simulation. It computes the piston position, the resulting number of vertical cells in the lower and upper pressure zones, counts total charging and discharging durations, derives the corresponding average flow velocities, and assigns the charge/discharge flow rate for each timestep.
```matlab
function [Res_900, LC, DC, Lflow, Dflow] = prepare_scenario_flow(Res_900, t_900, Nz, H, d)
```
**Calculation of flows**:
The velocities in charge and discharge direction are calculated on the basis of overall conservation of the total water volume which can be represented by its height $H(6)$. Thus the velocities are function of the total height and the cumulative charging/dischargng hours. 
$$
v_{charge, avg} = \frac{H(6)}{LC \times 3600}, \text{ with LC as the counter for charge intervals}
$$
$$
v_{discharge, avg} = \frac{H(6)}{DC \times 3600}, \text{ with DC as the counter for discharge intervals}
$$
It is thus assumed that the total water volume is once discharged and charged! Thus the vertical velocity is reffering to the the height and volume of the storage. At the same time this is the current vertical velocity within the water storage, which is needed for the forced convection part within the PDE. 

---

### Extracting Temperature Fields and Computing Thermal/Exergy Quantities

The helper function 'compute_energy_balances' compute all thermal energy contents and the exergy of the water volume based on the current temperature fields in the TPPS model. Characteristic energy calculations are exemplary shown for water and the exergy of water. 
```matlab
function [Heat_insu, Heat_Wasser, Heat_piston, Heat_Vsoil, Heat_Rsoil, Heat_Rinsu, WEX, wEX, DTRE] = ...
          compute_energy_balances(T_V, T_W, T_Sys, T_REf, Nz, dz, A, SW, z_RE, z_W, K, Spz)
        % --- Water column energy [J] ---
        Heat_Wasser = SW(1,2) * SW(1,3) * dz * A(1) * (sum(T_W(1,2:end-1)) + 0.5 * (T_W(1,1) + T_W(1,end)));  
        % with
        % SW(1,2) = 1000;                                 % Density ρ of water [kg/m³]
        % SW(1,3) = 4200;                                 % Specific heat capacity cp of water [J/kgK]

        % --- Exergy in the water volume [J] ---
        wEX = zeros(1, length(z_W));
        for p = 1:length(z_W)
            wEX(1,p) = T_W(p) - K * log((K + T_W(p)) / K);
        end
        WEX = Spz * dz * (sum(wEX(1,2:end-1)) + 0.5 * (wEX(1,1) + wEX(1,end)));
```

In general the thermal energy of a a 1D-segment is calculated as follows:

$$
E = \rho \, c \, A \int_0^H T(z)\,\mathrm{d}z
$$

and numerically discretized via:

$$
E \approx \rho\,c\,A\,\Delta z \left(\sum_{i=2}^{N-1} T_i + \frac{1}{2}(T_1 + T_N)
\right).
$$

In analogy the energy is calculated for 
- vertical insulation
- piston
- vertical soil
---

In addition the function computes the **exergy of the water column** along the vertical storage height.

First, the **specific exergy density** at each water node is defined as

$$
w_\mathrm{ex}(T_W)
= T_W - K \,\ln\!\left(\frac{T_W + K}{K}\right),
$$

where \(T_W\) is the local water temperature (or temperature surplus) and \(K\) is a constant associated with the reference temperature (e.g., return temperature).

The **total exergy** of the water column is then obtained by integrating this quantity over the height using the trapezoidal rule:

$$
W_\mathrm{ex}
\approx S_{pz}\,\Delta z
\left(
\sum_{i=2}^{N-1} w_\mathrm{ex}(T_{W,i})+ 
\frac{1}{2}\Big[w_\mathrm{ex}(T_{W,1}) + w_\mathrm{ex}(T_{W,N})\Big]
\right),
$$

where $S_{pz}$ is a scaling parameter $\rho c_P A$ and $\Delta z$ is the vertical grid spacing.

The **radial domains (soil and insulation)** are represented as a series of concentric cylindrical rings.
For each radial ring $m$ the vertical integration of temperature is carried out analogously to the 1D case, but scaled by a ring-specific geometric factor $z_{RE}(2,m)$ representing the effective lateral area of that ring.
For a given radial ring, the differential volume contribution is
$$
\Delta V_m = \Delta z \, z_{RE}(2,m)
$$
For each volume segment the sum over all vertical elements (several thousand, ```matlab sum(T_REf(2:end-Nz(9), m))```) is builded and finally multiplied with the material paramters for soil or insulation, respectively. 
```matlab
% Radial Soil before Insulation
for m = 1:Nz(12)
    DTRE(m) = dz * z_RE(2,m) * ...
        (sum(T_REf(2:end-Nz(9), m)) + ...
        0.5*(T_REf(1,m) + T_REf(end-Nz(9)+1, m)));
end
Heat_Rsoil = SW(2,2) * SW(2,3) * sum(DTRE);

% Radial Insulation
for m = 1:Nz(12)
    DTRE(m) = dz * z_RE(2,m) * ...
        (sum(T_REf(end-Nz(9)+2:end-1, m)) + ...
        0.5*(T_REf(end-Nz(9)+1, m) + T_REf(end, m)));
end
Heat_Rinsu = SW(3,2) * SW(3,3) * sum(DTRE);
```

## scenario_1_freeConv/HeatFluidSolid.m

### Position of Upper Bypass Inlet in Non-Moving Piston State.

The upper bypass inlet is placed within the replacement volume. Within the first 10 days of scenario_1 there is no piston moving,  the position is fix. Thus the water transported via the bypass enters within the ring gap, for the upper water volume.

### Contact Temperature between Water and Insulation

The goal is to have an **continious** heat flow on the contact point. Four contact temperatures are calculated
Lines: 
- 85
- 108
- 133
- 155

### Piston–Water Coupling Terms (Axial Heat Conduction)

The following code implements the **axial thermal coupling between the piston and the adjacent water nodes**, corresponding to Eqs. (4) and (5) in Häuslein's (2024) paper. 
The conductive heat flux at the piston surfaces (top and bottom) is attributed to the first discretized water layer and enters the water energy balance as a source/sink term. Note that for this calculation the contact temperatures in the previous section are needed and used.

The coupling coefficient `alpha_pw` represents the **normalized axial heat conduction** from the piston into a water control volume,

$$
\alpha_{pw}
=
\frac{\lambda_{\text{pist}} A_{\text{pist}}}
{\rho_w c_{p,w} A_{\text{TPPS}} \Delta z^2},
$$

which matches the prefactor of the piston coupling terms in Eqs. (4) and (5).

```matlab
%%%------------------------------------------%%%
% 09 Additional piston loss terms into the water nodes
%%%------------------------------------------%%%

% Coupling coefficient for axial heat conduction between piston and water
% Corresponds to λ_pist * A_pist / (ρ_w * c_p,w * A_TPPS * Δz²)
alpha_pw = (SW(2,1) * A(2)) / (A(1) * SW(1,2) * SW(1,3) * dz^2);

% Upper water node adjacent to the piston (top contact)
% Implements Eq. (5): axial heat flux from piston top into upper water node
idx_water_top_piston = Nz(3)+Nz(8)+Nz(2)+Nz(5)-2;
dTdt(idx_water_top_piston) = dTdt(idx_water_top_piston) ...
    - alpha_pw * (T(Nz(3)) - T(Nz(3)-1));

% Lower water node adjacent to the piston (bottom contact)
% Implements Eq. (4): axial heat flux from piston base into lower water node
idx_water_bottom_piston = Nz(3)+Nz(8)+Nz(2)-1;
dTdt(idx_water_bottom_piston) = dTdt(idx_water_bottom_piston) ...
    - alpha_pw * (T(1) - T(2));
```

## scenario_1_freeConv/HeattransferSzen.m

### Role of Time dt in the Wrapper

The full ODE state acknowledges both domains: the 1D cylinder system and the 2D radial soil field, where the soil is stored as a flattened vector of length `Nz(12)*Nz(13)` because `ode45` operates on 1D state vectors.
After integration, the “mapping” section reshapes the flattened soil part back to a 2D grid and synchronizes it with the 1D boundary/interface nodes to apply or enforce the 1D↔2D coupling consistently (e.g., mixed/Robin-type interface conditions and boundary bookkeeping).
Using `t2 = [0, dt/2, dt]` (`Nt2=2`) produces intermediate output states (e.g., at 450 s) so that this reshaping/interface handling can be evaluated at sub-times within one 900 s procedure step.
Note: `Nt2` mainly increases the temporal resolution of this mapping/interface handling at the stored output times; it does not directly prescribe the internal adaptive step size of `ode45`.

### ODE45 – Time Integration of the Discrete TPPS Model

The transient TPPS model is advanced in time using MATLAB’s `ode45` solver:

```matlab
[t2, T_Sys] = ode45(@HeatFluidSolid, t2, IC_Sys, [], Nz, dz, flow, T0init, SW, A, z_RE);
```

After spatial discretization of the governing heat equations (finite differences in axial and radial direction), the full system is reduced to a large set of ordinary differential equations (ODEs) of the form
$\dot{T}(t) = f(t, T)$
where the state vector $T$ contains all temperature degrees of freedom of the model (water, piston, insulation, soil), flattened into a one-dimensional vector as required by ode45. The solver `ode45` performs pure time integration of this ODE system. It does not know anything about grids, materials, or physical domains. All spatial coupling, material interfaces, and boundary conditions are encoded explicitly in the right-hand side function HeatFluidSolid. In particular:
- Diffusion and convection appear through finite-difference expressions in dTdt.
- Thermal coupling between domains occurs only if temperatures from other domains explicitly enter the corresponding derivative.
- Algebraic interface conditions (e.g. Robin-type contact temperatures) are enforced outside the ODE by direct reconstruction and are therefore not integrated as dynamic states.

Numerically, `ode45` uses an explicit adaptive **Runge–Kutta (4,5) Dormand–Prince scheme**. Conceptually, each state is updated according to: $T^{n+1} = T^n +\Delta t \dot{T}$ with the size $\Delta t$ chosen automatically to control the local truncation error. `ode45` is well suited here because the dominant thermal dynamics are moderately stiff and smooth in time, while adaptive step size control ensures numerical stability and efficiency without manual tuning of the time step.