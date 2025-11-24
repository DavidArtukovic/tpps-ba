#  01 - Code Info 
This notebook collects specific information on code junks.


## scenario_1/Szenario1.m

**Extracting Temperature Fields and Computing Thermal/Exergy Quantities**

The thermal energy of a a 1D-segment is calculated as follows:

$$
E = \rho \, c \, A \int_0^H T(z)\,\mathrm{d}z
$$

and numerically discretized via:

$$
E \approx \rho\,c\,A\,\Delta z \left(\sum_{i=2}^{N-1} T_i + \frac{1}{2}(T_1 + T_N)
\right).
$$

```matlab

% Vertical temperature vector (full 1D model stack)
T_V(1,:) = T_Sys(1, 1:Nz(10));

% Water temperature field extracted from the system state
T_W(1,:) = T_Sys(1, Nz(3)+Nz(8) : Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3);

% --- Thermal energy in insulation (vertical) [J]
Heat_insu = SW(3,2)*SW(3,3)*dz*A(1)*...
(sum(T_V(1,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-2:end-1))...
+0.5*(T_V(1,Nz(3)+Nz(8)+Nz(2)+Nz(5)+Nz(4)-3)+T_V(1,end))); 
```
---
The script computes the **exergy of the water column** along the vertical storage height.

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

where \(S_{pz}\) is a scaling parameter $\rho c_P A$ and \(\Delta z\) is the vertical grid spacing.

```matlab
% Thermal Energy in the radial insulation
for m = 1 : Nz(12)
    DTRE(m) = dz*z_RE(2,m)...
    *(sum(T_REf(end-Nz(9)+2:end-1,m))...
    +0.5*(T_REf(end-Nz(9)+1,m)...
    +T_REf(end,m))); %'K*m³
end

% Thermal Energy in radial direction in [J] relative ot initial timestamp
Heat_Rinsu = SW(3,2)*SW(3,3)*sum(DTRE); ;

for p = 1:length(z_W) % Voranalyse des Wasservolumens
    wEX(1,p) = T_W(p)-K*log((K+T_W(p))/(K)); % Exergie im gesamten Wasservolumen auf RÃ¼cklauftemperatur bezogen ohne aktive Nullung
end 
WEX = Spz*dz*(sum(wEX(1,2:end-1))+0.5*(wEX(1,1)+wEX(1,end)));

```
---


