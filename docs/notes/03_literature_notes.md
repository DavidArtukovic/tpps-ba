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
- Baseline model — needs extension towards natural convection.
- Reasons why in a system of 36m of water the piston has an optimal height of 18m (p. 8)


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

## Dahash (2017) - *Toward efficient numerical modeling and analysis of large-scale thermal energy storage for renewable district heating*

**Key Findings**
Develops a multiphysics numerical model for large-scale seasonal thermal energy storage systems (PTES and TTES) based on finite element methods. The model is validated against operational data from the Dronninglund PTES and shows good agreement in terms of temperature profiles and energy balances (deviations < 3%). The study demonstrates that although the storage achieves a high energy efficiency (~90%), only about 76% of the theoretical storage capacity is effectively utilized due to thermal mixing and losses.

**Stratification and MIX Number**
Thermal stratification is identified as a key performance driver. The study introduces the MIX number as a robust indicator based on the vertical energy distribution, showing moderate mixing behavior in real systems (typically MIX ≈ 0.2–0.55). The results highlight that operational strategies (charging/discharging schemes) significantly influence stratification quality, often more than geometry.

**Treatment of Convection**
The water domain is modeled using a one-dimensional axial discretization. Buoyancy-driven natural convection is not explicitly resolved but approximated via an enhanced thermal conductivity using Nusselt–Rayleigh correlations. This effectively introduces a mixing term, allowing the model to capture stratification degradation without resolving full fluid dynamics. Natural convection is not explicitly resolved but modeled via an enhanced thermal conductivity 
$\lambda_{\mathrm{eff}} = \lambda \cdot Nu$, with $Nu = C \cdot Ra^k$ and 
$Ra = \frac{g \beta \Delta T z^3}{\nu \alpha}$. 
The parameters $C$ and $k$ are not derived from first principles but empirically chosen or calibrated, 
so that buoyancy-driven mixing is effectively represented as a diffusion-like process rather than resolved flow.


## Xiang et. al (2022) - *Long-term thermal performance analysis of a large-scale water pit thermal energy storage*

**Key Idea:**
Finding optimal grid size for PTES modelling. A full CFD model is validated against Dronninglund PTES (Denmark) test facility. 
 modeled with a **1D multinode water model** that **does not resolve free or forced convection explicitly**, but enforces stratification numerically.

## Dahash (2019) - *Numerical heat transfer modeling of large-scale hot water tanks and pits*

 **Key Idea**
 Numerical modelling wallpaper.

## Dahash (2020) - *Techno-economic planning and construction of cost-effective largescale hot water thermal energy storage for Renewable District heating systems*

 **Key Findings**
 Investigates the techno-economic performance of large-scale seasonal TES (tank and pit) using a multiphysics model over a 10-year simulation horizon with daily time steps. The results show that storage performance improves over time due to ground preheating effects, reducing thermal losses during operation. Larger storage volumes exhibit higher efficiencies due to lower surface-to-volume ratios, confirming the importance of scaling effects.

 **Design Implications**
 Thermal insulation has a limited impact on efficiency for very large systems but remains essential to prevent excessive heating of the surrounding ground. Fully buried configurations show comparable or slightly improved performance due to reduced interaction with ambient conditions. Overall, tank TES provide higher efficiencies, while pit TES offer lower investment costs, highlighting a trade-off between performance and economic feasibility.
  

## Caputo (2021) - *District thermal systems: State of the art and promising evolutive scenarios. A focus on Italy and Switzerland*

**Key Findings**
Reviews different classification of DTS with particular focus on Italy and Switzerland.

**Relevance**
  For introduction this paper can serve as basic review for district heating systems with key findig of ever increasing number of DTS so far and in the near future.


## Dahash et al. (2021) - *Understanding the interaction between groundwater and large-scale underground hot-water tanks and pits*

**Key Findings**
Investigates TES with respect to the hydrological conditions, with groundwater flow anticipated. Darcy-Flow plays is significant for the magnitude of thermal losses. Consequently, TES do have impact on the nearby groundwater and vice-versa. In a realistic range of groundwater flow speed a representative model with 100.000m^3 tank volume loses up to 30 percentage points in thermal efficiency compared to no groundwater. The effect drastically reduces with rising system size, however does system are signficantly above 100meter in diameter. Additionaly, the paper proves the difficulty in maintainig the groundwater temperature below 25°C which is the recomendation in Germany. 

**Comsol vs FEFLOW**
COMSOL Multiphysics is used as the primary simulation tool to model the coupled heat transfer and groundwater flow processes in the TES system. In addition, the results are validated against simulations performed with FEFLOW, a specialized groundwater modeling software. Both models are applied independently to the same physical setup, and their results are compared to ensure the reliability of the COMSOL-based multiphysics approach.

**Relevance**
Particularly in the dicussion this paper can give a hint towards the groundwater topic. Especially if the ground is of high porosity the topic is increasingly dominant. Particulalrly important in face of groundwater regulations in many central European countries like Germany. Groundwater flow increases thermal losses and can push groundwater temperatures beyond regulatory limits (e.g. 20 °C), requiring mitigation measures such as cut-off walls and insulation. Larger systems are more efficient, but also more affected by groundwater interaction due to increased absolute losses. The paper highlights that purely conductive models underestimate heat transport in the presence of fluid motion. This supports your observation that enhanced transport mechanisms (e.g. free convection / flow) can act as a “thermal shortcut”, accelerating vertical heat redistribution compared to 1D diffusion-only models.

## VDI-Fachbereich Energie- und Umwelttechnik (2019) - *Verband Deutcher Ingineure*

**Relevance**
Recomendation for groundwater use up until 25°C.

## Nield (2017) - *Mechanics of Fluid Flow Through a Porous Medium*

**Relevance**
Intoductury book into Darcy flow as with the alternative for the modelling of Navier-Stokes in porous media. Darcy's law is experimentally derived.


## Borri (2021) - *Recent developments of thermal energy storage applications in the built environment: A bibliometric analysis and systematic review*

**Relevance**
Provides a structured, data-driven overview of recent research trends in thermal energy storage (TES), highlighting the increasing importance of TES in the context of decarbonization and the integration into district heating systems (DHS). The work is particularly relevant for the introduction, as it quantitatively underpins the growing research interest and application potential of TES technologies.

Furthermore, the paper offers a clear classification of TES technologies (sensible, latent, thermochemical), with a strong emphasis on sensible heat storage, which is directly relevant for TPPS. This classification can be adopted as a conceptual foundation for positioning the present work within the broader TES landscape.


## Tosato (2023) - *Simulation-based performance evaluation of large-scale thermal energy storage coupled with heat pump in district heating systems*

**Key Findings**
Investigates the coupling of large-scale TES with heat pumps in district heating systems using a combined 1D TES and multiphysics ground model. The integration of a heat pump increases the TES efficiency by about +6% for tank TES and +16% for shallow pit TES. More importantly, the usable energy output is significantly increased, reducing auxiliary heating demand by up to 60%. Additionally, the system-level performance improves with CO₂ emission reductions of up to 33%, highlighting the benefit of utilizing lower-temperature storage regions.

**Relevance**
The study demonstrates that not only high-temperature regions, but also lower-temperature zones within the TES become relevant when coupled with a heat pump. This increases the importance of accurately capturing temperature stratification and mixing processes within the storage. In particular, vertical heat transport mechanisms directly affect the exploitable energy content of the system. Therefore, modeling approaches that account for enhanced mixing or multidimensional effects are essential to assess TES performance in coupled TES–HP systems.
Suggest's study of TES in three segments, system-level, component-level (TES) and hydrogeological level.  


## Dahash et al (2021) - *Techno-economic and exergy analysis of tank and pit thermal energy storage for renewables district heating systems*

**Key Findings**
Tank TES have techno-economic advantages. For a 50m TES depth even with 26cm of insulation the groundwater temperature in 1m radial distance converges to 50°C.

**Relevance**
- MIX Number as stratification measure
- Can be cited in groundwater issues and relevance of insualtion issues.


## Reisenbichler et al (2025) - *Validation of a pit thermal energy storage model: Demonstration of a comprehensive approach*

**Key Findings**
Develops a comprehensive VVUQ framework for PTES models, including measurement, input, and numerical uncertainties. The Modelica-based model is validated against five years of operational data, showing good agreement with errors of ±5% for energy metrics and ±2 K for temperatures. Thermal losses exhibit large uncertainties due to indirect determination. Validation is shown to be inherently context-dependent and limited to specific scenarios. The modeling approach follows the common 1D multi-node fluid representation with simplified physics.

**Relevance**
Provides a key reference for model validation of large-scale TES and is mainly relevant for the general part of this work. It demonstrates that simplified 1D models can achieve good accuracy, while also highlighting limitations regarding mixing and loss prediction. This motivates the inclusion of additional physical effects such as enhanced transport or multidimensional modeling. Furthermore, it underlines that model accuracy must always be assessed with respect to the intended application.

-> cite in general part of literature

## Reisenbichler et al (2023) - *LargeTESModelingToolkit: A Modelica Library for Large-scale Thermal Energy Storage Modeling and Simulation*

**Relevance**
Conference Paper citing the TES package in Modelica.

## Formhals et al (2024) - *Development, validation and demonstration of a new Modelica pit thermal energy storage model for system simulation and optimization*

**Key Findings**
Modelling PTES in Modelica and comparing against COMSOL. 

**Relevance**
Cite in bibliography as newest extension of PTES modelling in nowel frameworks like Modellica.


## Ayele et al. (2021) - *Optimal heat and electric power flows in the presence of intermittent renewable source, heat storage and variable grid electricity tariff*

**Relevance**
Highlights the increasing importance of coupled electricity and heat systems (multi-energy systems) for integrating intermittent renewables and enabling system flexibility. Shows that thermal storage plays a key role in shifting demand and utilizing surplus renewable electricity, achieving over 97% renewable utilization and significant load shifting.

The work underlines that coordinated operation of heat and power networks is essential for decarbonization, as separate optimization leads to suboptimal system performance. Therefore, it supports the relevance of thermo-electric storage coupling as a key concept for future energy systems.

## Johannes et al. (2005) - *Comparison of solar water tank storage modelling solutions*

**Relevance**
Delivers insight into non uniform temperture distribution in a water tank storage along the radial axis, particularly where a steep thermocline in axial direction is present. The tank has an inlet at the top and an outlet approximetely in the middle.


## Untrau et al. (2023) - *A fast and accurate 1-dimensional model for dynamic simulation and optimization of a stratified thermal energy storage*

**Idea**  
Orthogonal collocation replaces the classical spatial discretization of the governing PDE by approximating the temperature field as a linear combination of global basis functions (Lagrange polynomials):
$$
T(z,t) \approx \sum_{i=1}^{N} T_i(t)\, l_i(z)
$$
Spatial derivatives are not computed via finite differences but through precomputed differentiation matrices:
$$
\frac{\partial T}{\partial z} \approx A_{OC}\, T, \quad
\frac{\partial^2 T}{\partial z^2} \approx B_{OC}\, T
$$
This transforms the PDE into a reduced system of ordinary differential equations:
$$
\frac{dT}{dt} = f\bigl(T, A_{OC}T, B_{OC}T\bigr)
$$
Thus, the spatial dimension is discretized globally, leaving only time-dependent states. This significantly reduces the number of degrees of freedom and enables fast simulation and optimization. However, local effects such as sharp thermoclines, buoyancy-driven mixing, or multidimensional transport are only approximated.


**Key Findings**
The orthogonal collocation approach allows accurate representation of stratified temperature profiles with only a small number of collocation points, leading to a strong reduction in computational cost compared to classical multi-node models. The method reduces numerical diffusion and can preserve thermoclines more sharply under suitable conditions.  

At the same time, the accuracy depends strongly on the number and placement of collocation points, and steep gradients may lead to oscillations or reduced fidelity. Physical mixing processes such as free convection are not resolved explicitly and must be incorporated via additional modeling assumptions.  

Overall, the method is particularly suitable for large-scale TES simulations and optimization problems, where low computational cost and reduced model order are essential, but less appropriate for detailed physical analysis of transport phenomena.


**Summary of Thermocline Handling**
Paper delivers a summary of methods to deal with thermoclines.
- First group are "algorithmic" ex-post procedures, where after a time-step the system is screened for inverse thermoclines and effectifely reorganized. Goes back to Franke (1997). Disadvantage is the neglection of mixing effects, and thus an overestimation of exergy.
- The second group is homogenize the temperatures around the inversion, effectively forcing a regression to the mean beahviour.
- Third variant via a modification of the difussion coefficient, where I can add extensive reference to Dahash's papers in the case of large-scale TES. Advantes is the additional coefficient is a continious formulation and can does directly act within the PDE. No ex-post modification is necessary. However, the parametrizaiton of the coefficient demands some sort of experimental/reference data.

**Relevance**
- Reasons that thermoclines can appear especially in the contact to cool surfaces (like the piston or the wall). Also a reference to the paper of Cesaro (2003) is present. The consequence is destratification.
The OC approach delivers interesting alternative for large scale TES calculation, where simple geometry is present. However, in the context of the TPPS with ringgap, moving piston local phenomena cannot be captured accurately by global polynomials. As Reduced-Order Model particulalry interesting for repeated run simulations. Speed is attached to DOF. The degrees of freedom (DOF) correspond to the number of independent temperature states. While classical discretization assigns one DOF per spatial cell, orthogonal collocation reduces this to a small number of global coefficients representing the entire temperature profile.

## 1993 Kleinbach - *Performance Study of One-Dimensional Models for Stratified Thermal Storage Tank*

**Relevance**
Provides a systematic comparison of classical 1D storage tank models (multi-node and plug-flow) as implemented in TRNSYS, with a focus on computational efficiency and prediction accuracy. Particularly relevant as a representative of low-order modeling approaches where free convection is not resolved explicitly but implicitly treated through mixing assumptions. Serves as a baseline for understanding how stratification and destratification are approximated in simplified TES models.

**Key Findings**
The study shows that stratified storage can be reasonably captured using discretized 1D models, where the tank is divided into fully mixed segments and energy balances are solved for each node. However, temperature inversions are not physically modeled but numerically removed by instant mixing of adjacent segments, effectively corresponding to a mass-weighted averaging of temperatures. Plug-flow models improve computational efficiency and better preserve stratification compared to classical multi-node approaches, though they can overpredict energy quantities under certain conditions. The inclusion of plume entrainment aims to mimic mixing processes but does not significantly improve model accuracy while increasing computational cost. Overall, the results highlight that simplified models rely heavily on empirical or numerical treatments of mixing rather than resolving buoyancy-driven flow dynamics.

## Powell, Edgar (2013) - *An adaptive-grid model for dynamic simulation of thermocline thermal energy storage systems*

**Key Idea**
Introduces an adaptive 1D grid that dynamically tracks the thermocline region (i.e., the moving temperature gradient) instead of using a fixed, equidistant discretization. High spatial resolution is applied only in the region of strong temperature gradients, while the remaining parts of the tank are represented by larger, quasi-isothermal control volumes.

**Key Findings**
Models a thermal energy storage with a adaptive 1d grid with coincides with flow. Even when reducing the number of gird points by factor 10 the RMSE stays in comaprative region to standard model (equal spacing).

**Relevance**
Uses as well an flexible difussion coefficient to account for turbulence, inverse thermoclines  and thermal conductivity. Also highlights that analytical solutions are less useful in face of mass flow. Here one can mention that some analytical solutions have been recently developed for specific cases. Figure 6 is a good illustrational example of how inverse thermoclines cause the exergy level to sink. 

Also the adaptive grid can be an interesting option in the discussion. The increase in runtime due to 2D piston can be counteracted with a less fine gitter in total, and a equally fine gitter in water regions where a strong thermocline persists. (Hierzu eventuell eine Grafik.)

The adaptive-grid concept suggests that computational cost can be reduced without sacrificing accuracy by dynamically allocating resolution to regions of strong gradients. This is especially promising for large-scale TES modeling, where sharp thermoclines dominate system behavior.  

However, the approach relies on prior assumptions about the structure of the temperature field (e.g., existence of a single thermocline) and may therefore be less suitable for systems with complex mixing dynamics or multidimensional effects. As such, it represents a complementary strategy to physically detailed models rather than a direct replacement.

## McPherson (2018) - *Deploying storage assets to facilitate variable renewable energy integration: The impacts of grid flexibility, renewable penetration, and market structure*

**Relevance**
Use as citation in introduction for relevance to takle temporal disparities in energy demand and production/supply. 

## Soares (2022) - *Efficient temperature estimation for thermally stratified storage tanks with buoyancy and mixing effects*

**Relevance**
Very handsome introduction. You can oreintate here.

**Key Idea**
The work by Soares et al. builds directly on the smooth 1D model introduced by Lago et al., which was the first to incorporate buoyancy and mixing effects using continuous formulations. :contentReference[oaicite:2]{index=2}  

While Lago et al. focused on the development of the model itself and its validation for large seasonal storage systems, Soares et al. extend this approach by generalizing it to different tank topologies, integrating the model into a unified formulation, and embedding it into a state estimation framework.  

Thus, the paper represents an important step from model development towards practical application in control and optimization contexts.

## Lago (2019) - *A 1-dimensional continuous and smooth model for thermally stratified storage tanks including mixing and buoyancy*

**Relevance**  
The paper introduces a fundamentally different approach to modeling free convection (buoyancy) in 1D thermal energy storage models. Instead of using non-smooth post-processing mixing algorithms, buoyancy effects are directly embedded into the system dynamics using smooth and continuous functions.

This enables the use of gradient-based optimization methods (e.g. automatic differentiation), which are otherwise not applicable to classical TES models due to their non-smooth structure. As a result, the model is particularly relevant for control and optimization applications, where computational efficiency and differentiability are crucial.  

Furthermore, the distinction between slow and fast buoyancy effects provides a more physically consistent representation of mixing processes, especially under charging and discharging conditions.  

For the present work, the paper is relevant as it highlights both the limitations of classical 1D approaches in representing multidimensional transport phenomena (e.g. convection) and a possible modeling strategy to partially compensate for these limitations without increasing spatial dimensionality.


## Benato,  Stoppato (2018) - *Pumped Thermal Electricity Storage: A technology overview*

**Key Findings**


**Relevance**