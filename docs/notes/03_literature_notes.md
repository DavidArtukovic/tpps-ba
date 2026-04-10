# 03 - Literature Notes – TPPS Bachelor Thesis

This document contains concise summaries of all relevant papers and theses.
For each source, the following aspects are listed (if needed):
- Key idea / problem addressed
- Model assumptions
- Relevant equations
- Relevance for my thesis

---

## Häuslein et al. (2024) – *Dynamic Modeling of Thermal Pumped Piston Storage*
**Key idea:**  
A dynamic model of the TPPS system where the system's optimal size is determined based on demand characteristics. A full cycle simulation is conducted based on consumption data from the city of Wunsiedel. For me particularly interesting is the thermal simulation. So far no free convection is incoporated. Additionally, the penstock water is modelled via a replacement volume which is to be revised. 

**Assumptions:**  
- Water modeled as a 1D continuum.
- No natural convection included (only forced or mixed-flow approximations).

<!-- **Relevant equations:**   -->
- 1D water/soil heat transport equation 
- σ-coupling terms (piston ↔ water ↔ soil)

**Relevance for my thesis:**  
Baseline model — needs extension towards natural convection.

**Questions**
The coupling terms demand further investigation and understanding!

---

## Gerle (2019) – Modeling of a Zero-Energy Storage (NES)
**Key idea:**  
1D thermal model of a hot water storage including natural convection.

**Relevance:**  
Provides foundations for modeling free convection in TPPS.

## Schäfer et al. (2021) – *Development of a zero-energy-sauna: Simulation study of thermal energy storage*
**Key idea:**  
Modeling of a zero-energy sauna via a 1D energy conservation model. Includes natural convection through an effective mixing term. Idea is thoroughly explained in Gerle (2019)

**Relevance:**  
My approach for representing free convection in low-velocity water domains.

**Questions**
Explanation of turbulent regions demands further investigation.

## Wärmeübertragung-Grundlagen und Praxis (2019) - *Conctact temperature + radial heat conduction (coupling terms)*

**Relevance - contact temperature:**  
At pages 66-67 the book explains how one can arrive at the expression of an instantantious contact temperature for two semi-infinite bodies.
Equation 2.82 is the temperatuer explained in model_overview as well as equation (8) in the paper of Häuslein (2024).


**Key idea - radial heat conduction:**  
Analytical derivation of **stationary radial heat conduction** in hollow cylinders based on Fourier’s law in cylindrical coordinates.  
Because the heat transfer area depends on the radius $A(r)=2\pi r l$, the resulting thermal resistance is **logarithmic**, which is characteristic for pipe-wall geometries.

**Core result:**
- Radial heat flux through a hollow cylinder wall:
  $$
  \dot Q
  =
  \frac{2\pi l\,\lambda}{\ln(r_2/r_1)}
  (\vartheta_1-\vartheta_2)
  $$

**Relevance for my thesis:**  
- Provides the **theoretical basis** for the radial coupling (“pipe wall”) terms used in Häuslein et al. (2024).
- The logarithmic term $\ln(R_1/R_0)$ in the TPPS model directly originates from cylindrical heat conduction.
- The radial loss terms $\sigma(z)$ can be interpreted as this pipe-wall heat flux, normalized by the local heat capacity
  $$
  \rho c_p \pi R_0^2,
  $$
  yielding a volumetric source/sink term in a 1D axial energy balance.

**Source:**  
VDI – *Wärmeübertragung: Grundlagen und Praxis*, 2019, Chapter 2.1.3 “Wärmeleitung in einem Hohlzylinder”.


## Cesaro Oliveski et al. (2003) – *Comparison between models for the simulation of hot water storage tanks*

**Key idea:**  
This paper systematically compares **two-dimensional (2D)** and **one-dimensional (1D)** numerical models for vertical hot water storage tanks under **natural and mixed convection**. A detailed 2D finite-volume model (laminar + low-Reynolds-number turbulence modeling) is **experimentally validated** and subsequently used as a **reference solution** to assess the accuracy and limitations of commonly used 1D multinode models. One objective is to determine under which conditions **1D modeling is sufficiently accurate** for long-term system simulations.

### Modeling approaches  

**Two-dimensional model (reference):**
- Axisymmetric cylindrical tank
- Full solution of:
  - Mass conservation
  - Momentum equations (axial and radial)
  - Energy equation
- Buoyancy-driven flow modeled via density variation (Boussinesq-type approach)
- Mixed convection handled using a **Launder–Sharma low-Re κ–ε turbulence model**
- Captures:
  - Inlet jet effects
  - Local turbulent regions near the inlet
  - Large laminar regions in the bulk of the tank
- Experimentally validated for both natural and mixed convection

**One-dimensional models (multinode):**
- Tank discretized into axial segments with uniform temperature per segment
- Governing equation includes:
  - Axial heat conduction
  - External heat losses
  - Advective energy transport (if mass flow is present)
- No explicit physical modeling of internal natural convection
- Require numerical corrections (“computational artifices”) to maintain physical consistency

→ In an actual water tank, the water of the vertical layer becomes denser than its surroundings and
slips towards the bottom of the tank. One-dimensional models cannot reproduce this stratifying-by-cooling
phenomenon, since physically the movement of water happens in two dimensions. However, the relevance of such **wall-induced downward flow**, which also induce upward streams in the tank middle, thus causing overall circulation movements, reduce with increasing system size or specifically with decresing surface to volume ratio.

### One-dimensional model variants  

1. **Simple multinode model**
   - Pure energy balance formulation
   - Fails for pure natural convection (cooling without flow)
   - Produces non-physical temperature profiles with a maximum in the tank center

2. **Multinode with inversion (Franke, 1997)**
   - Enforces monotonic stratification by swapping adjacent segments if  
     $ T_{i} < T_{i+1} $
   - Artificial but numerically robust
   - Corrects temperature inversion without modifying total energy directly

3. **Multinode with mean**
   - Segments involved in temperature inversion are merged via a weighted mean temperature
   - Removes inversion but introduces additional numerical diffusion
   - Requires iterative upward and downward scanning to remove discontinuities

### Key findings  

- The **2D model** reproduces experimental temperature profiles very accurately for both natural and mixed convection and serves as a reliable reference.
- In mixed convection, **turbulence is confined to the inlet region**, while most of the tank remains laminar.
- **1D models without correction are physically inconsistent** under pure natural convection.
- When combined with inversion or mean artifices, **1D models reproduce temperature stratification reasonably well** compared to the 2D reference.
- However, corrected 1D models show:
  - **Lower mean tank temperatures**
  - **Artificially increased thermal losses**
  - A slight **violation of energy conservation** compared to the 2D model
- Time-step size is critical:
  - Large time steps amplify errors in 1D models
  - 1D models still allow time steps orders of magnitude larger than 2D models

### Relevance for this thesis  

This paper provides a **methodological justification for using 1D storage models** in long-term and system-level simulations, where computational efficiency is essential and detailed flow structures are of secondary importance. At the same time, it highlights **systematic limitations of 1D models**, particularly:
- The need for artificial stratification enforcement
- The tendency to underestimate stored thermal energy
- The sensitivity to numerical time-step size

These insights are directly relevant for the present thesis, as they motivate the use of **computationally efficient 1D models** while emphasizing the importance of carefully designed source terms or correction strategies to represent **natural convective mixing** in a physically consistent manner.

### Open questions / points for deeper understanding  

- Detailed derivation and physical interpretation of the **turbulent inlet region** and its confinement.
- Qualitative understand of **multinode-with-inversion** and **multinode-with-mean** approaches.
- Transferability of the proposed correction strategies to modified 1D formulations with explicit mixing source terms.


## Franke (1997) – *Object-oriented modeling of solar heating systems*

**Key idea:** 
Franke introduces a modular, **object-oriented framework** for simulating solar thermal heating systems. The system is decomposed into physical components (e.g. collectors, pipes, storage), each described by energy and mass balances. Stratified hot water storages are modeled in 1D, with simplified representations of heat losses and mixing.

**Assumptions of Inversion:**
- Buoyancy acts instantanious
- No inertia
- No partial mixing
- No dependence from Rayleigh, Prandtl and Geometry

**Relevance for TPES Modelling:**
Franke provides the idea of multimode with inversion approach (also adopted by Cesaro (2003)):
The algorithm is a simple pair-wise **permutation policy** which conserves energy
```
for i = bottom … top-1:
    if T[i] > T[i+1]:        # instable Stratification
        permute nodes i and i+1
        set lower = T[i]
        set T[i] = T[i+1]
        set T[i+1] = lower
repeat until profile is monotone
```

- Clear separation between physical modeling and numerical implementation
- Free Convection as infinitely fast and perfectly efficient process.


→ Relevant mainly as a modeling philosophy reference, not for specific equations.
→ However, Cesaro finds that the **projected energy losses are to high**, despite solid temperature distribution.
→ Particularly not useful in penstock induced convection.

## Spall (1998) – *Transient mixed convection in cylindrical thermal storage tanks*

**Key idea:**  
Numerical CFD study (axisymmetric, cylindrical tank) of **mixed convection during charging** of thermal storage tanks. The work systematically investigates the influence of the **Archimedes / Richardson number** and Reynolds number on **stratification stability, internal mixing, and thermocline thickness**.

**Main findings (flow regime classification):**
- **Ar ≲ 1 (Ri ≲ 1):**  
  Inertia-dominated regime. Strong internal mixing, recirculation cells and descending warm plumes occur. No stable thermocline; stratification breaks down.
- **Ar ≈ 1–2:**  
  Transitional regime. Buoyancy begins to suppress inertial structures, but intermittent mixing persists. Stratification remains sensitive to geometry and turbulence modeling.
- **Ar ≥ 2 (Ri ≥ 2):**  
  Buoyancy-dominated regime. Stable vertical stratification is established. Thermocline thickness becomes well-defined and largely independent of Reynolds number.
- **Ar ≳ 5:**  
  Fully stratified flow. Inertial effects are negligible; the velocity field is governed by buoyancy-driven structures.

**Additional insights:**
- At fixed Ar, the **Reynolds number has little influence** on stratification quality (Re = 500–3000).
- Turbulence modeling strongly affects results:  
  k–ε models significantly **overpredict thermocline thickness**, while Reynolds Stress Models predict sharper, more realistic stratification.

**Relevance for TPPS modeling:**
- Confirms that **Richardson (or Archimedes) number is the governing control parameter** for buoyancy-dominated storage behavior.
- TPPS operation with low inflow velocities and large thermal gradients is expected to satisfy **Ri ≫ 1**, justifying:
  - neglect of inlet inertia effects,
  - use of **1D stratified storage models with buoyancy-based mixing closures** (Franke / Cesaro type).
- Provides quantitative support for assuming **Re-independence** of stratification once the buoyancy regime is reached.



## Xiang et al. (2022) – *CFD-based analysis of large-scale Pit Thermal Energy Storage (PTES)*

**Key idea:**  
High-resolution CFD modeling of a full-scale PTES (Dronninglund, Denmark) to analyze stratification, turbulence, and heat losses during charging and discharging. CFD is calibrated/validated against measured temperature data of a real PTES and then used as a **reference model**, while long-term simulations are explicitly delegated to simplified models due to computational cost.

**Assumptions:**  
- Full-scale PTES geometry with water and surrounding soil.
- Transient CFD simulations for selected representative days only.
- Turbulence mainly relevant near inlet diffusers; bulk storage remains quasi-laminar.
- Stratification dominated by vertical temperature gradients.

<!-- **Relevant equations:** -->
- Energy moment definition used for stratification assessment (Eq. 14).
- MIX number quantifying deviation from ideal stratification vs. fully mixed state (Eq. 15).

**Relevance for my thesis:**  
Highly relevant due to **system layout (PTES)**, **large storage size**, and **cyclic charging/discharging operation**, which closely resemble the target application. The paper strongly supports the modeling strategy of:
- using CFD only for **calibration and methodological guidance**, not for multi-year simulations,
- focusing on **vertical resolution**, as horizontal grid size is shown to be of minor importance (argument in favor of 1D models),
- selecting grid size based on temperature gradients, with **≈0.05–0.06 m** identified as a good compromise (consistent with my chosen 0.05 m),
- treating turbulence as a **localized inlet phenomenon** rather than a dominant bulk effect.

The MIX concept is particularly useful as an evaluation metric for stratification quality in simplified (1D) models.

**Questions**
- How can the MIX concept be robustly transferred to a purely 1D modeling framework?
- How to represent inlet-induced turbulence and mixing in a reduced-order model without excessive numerical diffusion?

## Xiang et al. (2022) – *A comprehensive review on pit thermal energy storage (PTES): technical elements, numerical approaches and recent applications*

**Key idea:**  
Comprehensive review of large-scale **pit thermal energy storage (PTES)** with emphasis on **design aspects**, **heat transfer mechanisms**, and especially **numerical modeling strategies** ranging from simplified 1D system models to detailed CFD approaches. The paper synthesizes results from ~160 studies and identifies best practices and limitations for long-term storage simulations.

### Numerical modeling approaches  

- **Reduced-order / system-level models (1D, TRNSYS-type):**
  - Water storage typically modeled as **1D vertically stratified**
  - Soil modeled in 2D or 3D with effective heat conduction
  - Natural convection in the water often **not resolved explicitly**
  - Stratification maintained via:
    - increased vertical resolution,
    - effective thermal conductivity,
    - artificial mixing or inversion-prevention schemes
  - Suitable for **multi-year simulations**, but sensitive to parametrization

- **Finite Element models (e.g. COMSOL):**
  - 2D/3D coupled water–soil domains
  - Can include **groundwater flow (Darcy)** and complex boundary conditions
  - High physical fidelity, but limited scalability for long-term studies

- **CFD / Finite Volume models:**
  - Used to resolve **natural convection, inlet mixing, and wall effects**
  - Typically restricted to **short time horizons**
  - Serve as **reference models** for understanding and calibration

### Key findings  

- For large PTES systems, stratification is dominated by **vertical temperature gradients**.
- **Inlet-induced mixing and turbulence** are localized phenomena; the bulk storage remains largely laminar.
- Horizontal discretization has **minor influence** compared to vertical resolution.
- Long-term storage performance predictions are highly sensitive to:
  - boundary conditions,
  - soil properties,
  - insulation and cover modeling.
- Many simplified models reproduce annual energy balances well, but:
  - mispredict thermocline thickness,
  - underestimate or overestimate heat losses on shorter time scales.

### Relevance for TPPS / this thesis  

- Although focused on PTES, the numerical insights are **directly transferable to TPPS**:
  - large water volumes,
  - slow charge/discharge,
  - dominance of buoyancy-driven stratification.
- Strongly supports the use of **1D vertically resolved water models** for system-level TPPS simulations.
- Confirms the need for:
  - carefully designed **mixing or source-term closures** to represent unresolved natural convection,
  - CFD studies as **supporting tools**, not primary simulation frameworks.
- Provides a solid **literature-backed justification** for a hierarchical modeling approach:
  - CFD / high-fidelity models → physical understanding
  - Reduced-order models → long-term and grid-coupled simulations

### Open questions / points for deeper understanding  

- How to parameterize buoyancy-driven mixing consistently in 1D models?
- How to translate CFD-based stratification metrics into reduced-order diagnostics?
- How sensitive are TPPS results to soil and boundary-condition uncertainties compared to PTES?

## Bai et al. (2019) – *Numerical and experimental study of an underground water pit for seasonal heat storage*

**Key Idea:**
Seasonal PTES modeled with a **1D multinode water model** that **does not resolve free or forced convection explicitly**, but enforces stratification numerically.

**Treatment of convection:**  
- **Forced convection (charging/discharging):**  
  Represented implicitly via **plug-flow–like node-to-node advection**
- **Free convection:**  
  Handled by **multinode with inversion (after Franke, 1997)**:
  - After each timestep, vertical temperature profile is checked
  - If $ \partial T / \partial z > 0 $ → **nodes are instantaneously mixed**
  - Ensures static stability without resolving flow physics

**Numerical implications:**  
- Free convection is **not predicted**, but **suppressed a-posteriori**
- Buoyancy effects only appear as a **numerical correction**

**Relevance for my thesis:**  
- Confirms that PTES literature **avoids explicit free-convection modeling**
- Highlights the gap between:
- Provides a **baseline approach** to improve upon with physically motivated free-convection closures.

**Relation to Franke (1997):**  
Direct application of the **multinode with inversion** concept to guarantee stratification stability.

## Yang et al. (2021) - *Seasonal thermal energy storage: A techno-economic literature review*

**Key Finding:**
Classifies different technologies of Seasonal tehrmal energy storages (STES) by technical and economic parameters. STES as overarching category of thermal energy storages where PTES are part of it (and thus TPPS). Delivers extensive review of the related literature. Important to note that one of the first large-scale PTES was developed at the University of Stuttgart in 1984. 

**Relevance for PTES/TPPS:**
Among compared PTES have higher storage efficiency given the insulation. PTES relativley indifferen to geological conditions but idealy no ground water and obviously large scale. Also with ever higher volume the efficiency starts so sink again due to thermal destraficiation (e.g. due to free convection).

## Hahne (2000) - *The ITW solar heating system: an oldtimer fully in action*
**Key Finding:**
One of the first large scale low-temperature PTES dating from 1985, operating 15 years with high reliability up to the publication of the paper. The work was conducted by the Institut of Thermodynamics and Thermal Engineering, one of the predecessors of the current IGTE. The storage had a volume of 1050m^3, wheras water composed 960m^3. The system was equipped with an extensive apparutus for thermal charging and discharging, compromising heat-exchanger coils and water exchange mechanisms at the toop and bottom. Thermal energy was supplied by solar collectors, providing 35°C warm water. The heat balance over several years showed a storage efficiency of about 80% (seite 478). While during summer months a high loss in the soil can be reported the effect is reversed in the winter months, helping in the operation of the system and achieving high efficiencies. The low temperatures demand a heat-pump which proved to be the limiting factor back then.

## Dahash (2020) - *Toward efficient numerical modeling and analysis of large-scale thermal energy storage for renewable district heating*

**Key Idea:**
Dahash et al. present a numerical model for large-scale pit thermal energy storage (PTES), focusing on the prediction of thermal stratification, heat losses, and overall system performance. The model is based on a one-dimensional energy balance equation including advection, conduction, and loss terms, while the surrounding ground is modeled separately using a multidimensional heat conduction approach. Instead of resolving fluid flow explicitly, buoyancy-driven convection is incorporated through an enhanced thermal conductivity formulation, where the effective conductivity is increased based on Nusselt–Rayleigh correlations. This approach allows the model to capture mixing effects at low computational cost, but inherently limits the physical representation of convective transport. The results show that PTES systems can achieve high efficiencies (up to ~90%), although only a fraction of the storage volume is effectively utilized due to thermal losses and mixing. Furthermore, the study highlights the strong influence of ground properties and storage geometry on system performance, as well as significant uncertainties related to model parameters and boundary conditions. Overall, the work demonstrates that simplified 1D-based models can provide reasonable accuracy for system-level analysis, while emphasizing the need for improved representations of convection and multidimensional effects.

**Relevance to my work**
Check if my increase or distribution in convection can be reasoned herewith.

**Solving of inverse thermoclines**
In the presented model, buoyancy-driven convection is not resolved explicitly but approximated using an enhanced thermal conductivity approach based on Nusselt–Rayleigh correlations. This enhancement is only applied in the absence of volumetric flow and under unstable stratification conditions, ensuring that convective effects are not double-counted when advection is present. The Nusselt number is formulated as a power-law function of the Rayleigh number, allowing convective heat transfer to be represented as an increase in effective thermal conductivity, albeit relying on empirical parameters and without explicitly resolving flow dynamics.

## Xiang et. al (2022) - *Long-term thermal performance analysis of a large-scale water pit thermal energy storage*

**Key Idea:**
Finding optimal grid size for PTES modelling. A full CFD model is validated against Dronninglund PTES (Denmark) test facility. 
 modeled with a **1D multinode water model** that **does not resolve free or forced convection explicitly**, but enforces stratification numerically.

## Dahash (2019) - *Numerical heat transfer modeling of large-scale hot water tanks and pits*
 **Key Idea**
 Numerical modelling of FK