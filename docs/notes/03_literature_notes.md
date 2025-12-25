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

---

## Schäfer et al. (2021) – *Development of a zero-energy-sauna: Simulation study of thermal energy storage*
**Key idea:**  
Modeling of a zero-energy sauna via a 1D energy conservation model. Includes natural convection through an effective mixing term. Idea is thoroughly explained in Gerle (2019)

**Relevance:**  
My approach for representing free convection in low-velocity water domains.

**Questions**
Explanation of turbulent regions demands further investigation.

## VDI-Wärmeatlas (2019) - *Free Convection*

**Relevance:**  
Give's the theoretical backround for the use of semi-infinite bodies in the calculation of the contact temperature. To be thoroughly derived in model_overview.

---

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

---

## Franke (1997) – *Object-oriented modeling of solar heating systems*

**Key idea:** 
Franke introduces a modular, **object-oriented framework** for simulating solar thermal heating systems. The system is decomposed into physical components (e.g. collectors, pipes, storage), each described by energy and mass balances. Stratified hot water storages are modeled in 1D, with simplified representations of heat losses and mixing.

**Relevance**
Franke provides the methodological foundation later adopted by Cesaro (2003):
- Justification of 1D stratified storage models for system-level simulations
- Clear separation between physical modeling and numerical implementation
- Conceptual basis for source/sink terms representing mixing and losses

→ Relevant mainly as a modeling philosophy reference, not for specific equations.

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


