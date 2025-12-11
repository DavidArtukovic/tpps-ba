#  01 - Code Info 
This notebook collects specific information on code junks, calculations within the matlab files and related stuff.


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
<span style="color:blue"> It is thus assumed that the total water volume is once discharged and charged??????</span>.
At the same time this is the current vertical velocity within the water storage, which is needed for the forced convection part within the PDE. 

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


## scenario_1/HeattransferSzen.m

