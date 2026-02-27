# 02 - To Do's, Meetings, Fragen BA TPPS
Hier werden To Do's, Erkentnisse, Fragen sowie Protokolle zu den Treffen mit den Betreuern gesammelt. Oben stehen immer die anstehenden To Do's.

# Open Tasks
- [x] Derive algebraic expressions for the inclusion of the convection term.
- [x] Go through all scripts, update indexing and comments.
- [x] Extensive commenting in code_info of szenario1.m 
- [x] Extensive commenting in code_info of HeatFluidSolid.m 
- [x] To-Do - Aufräumen des dt terms in heattransferSzen und HeatFluidSolid
- [x] Fix difference between own and Dominic Code (until 27.12.2025)
- [x] Understand the made-up of the coupling terms (until 30.12.2025)
- [x] Integrate free convection via mixing term into the existing Matlab code (closed form) (until 03.01.2026)
- [ ] Propadeaticum related literature on free convection big systems, pace of one paper per day (until 06.01.2026)
    - Especially from the Xiang (2022) Paper the Comsol examples of Dahash and several 1D Water + 2D Soil examples should be discussed.
- [x] Understand relevance of Froude, Reynolds and Richardson number as well as their calculation. Calculate for the present case (08.01.2026)
- [ ] Integrate the shift term after free convection (until 23.01.26)
- [x] Singnaling in the time-series plot, where free convection is happening (until 22.01.26)
- [x] Uderstand where the temperature increase in the lower replacement volume comes from in the base model (until 22.01.26)
- [x] set inlet to a fix position (until 25.01.26)
- [x] lambda hybrid in beide richtungen kompatibel (bis 25.01.26)
- [x] Wasservolumen berechnen was pro Zeitschritt zugeführt wird. (bis 26.01.26)

## 02.03.26 Meeting Dominic
- Frage, warum wird beim Init.m file nicht die exakte Fläche der Kreissegmente berechnet?
$$A_{\text{exact}} = \pi \left( r_2^2 - r_1^2 \right)$$
$$A_{\text{FD}} \approx \pi \,\Delta r \left( r_1 + \frac{\Delta r}{4} \right), \quad \text{mit } \Delta r = r_2 - r_1$$
- Loop über Nt2 ist mir nicht gaz klar.

## 23.01.26 Meeting Micha/Dominic
- Freie Konvektion in geschlossener Form als ex-post Update implementiert. Problem: Implementierung benötigt zugeströmte Energeimenge $Q_{add} \propto \Delta t ( T_{\text{upper bypass inlet}} - T_{\text{lower bypass inlet}})$. 
    →  Diese muss vom oberen Volumen entnommen werden. Allerdings wird die zugeführte Energie nicht über einen Massenfluss realisiert sondern über eine Dirichlet-Randbedingung a.k.a "Top-Knoten wird auf 80 °C festgenagelt". Dirichlet kann beliebig viel Energie einspeisen/abziehen – abhängig davon, was der Solver gerade braucht, um den Knoten zu treffen.
- Idee: Füge am Interface Knoten eine Energiequelle hinzu: $$\frac{dT_{\text{interface}}}{dt} += \frac{\text[flow]}{\Delta z}(T_{in}-T_{\text{interface}})$$


## 23.01.26 Meeting Dominic
- Frage zur Kopplung des Ersatzvolumens (Ringspalt) mit dem Kolben. Warum Kolbenoberseite mit Ersatzvolumenoberseite?
- HeattransferSzen 1D->2D mapping Zeile 155
- HeatFluidSolid z.B. Zeile 126-132. Warum unterscheidung zwischen lade/entlade Fall. Beim Laden verlässt Wasser oberes Volumen. Warum wird nicht immer der unterste Wasserknoten genommen?
- Frage auch an mich, wie kommt aktuell warmes Wasser nach unten? Was passiert genau zwischen den Zeitschritten?
- Ask why BCs between piston and water are differently treated than water-insulation and water-soil.
    - Temperature gradients at water-soil and water-insulation are zero, they become set every time the ODE is solved for one iteration step.
- Ersatzvolumen des Ringspalts wird aktuell nicht geplottet
- Warum hat die obere Totzone ebenfalls ein Kuppelvolumen [init.m Zeile 193 -> 2+H(7), aber H(7) ist 1.665m]?

## 25.12.25 To Be Discussed in Literature
Dahash (2017) gives idea of temperature gradient specific thermal conductivity.  
Franke (1997) has to be better understand and derived.

## 12.12.25 - Fragen Dominic
- Warum in HeattransferSzen.m Zeile 110 die Temepratur von replacement_bottom_idx-1, analog Zeile 117 im Ladefall?
- Fragen zu der ganz linken und zweiten von links grafik in TPPS.svg. Links einfach nochmal die Indizes ohne Erdreich und zur Visualisierung des realen Kolbens im Wasser?
- Frage zu timestep.m - Bedeutet die 1 in Zeile 2 von t_hour die Kolbenposition, also ganz oben? 
- In scenario1 liegen unterschiedliche Lade- und Entlade Geschwindigkeiten fürs thermische Laden vor? Zeile 277 in Scenario1_ref. Liegt das daran dass angenommen wird dass in 20 Tagen thermisch komplett be- und entladen wird? Mit anderen Worten die vertikale Strömungsgeschwindigkeit wird auf die Hubhöhe bezogen?

## 02.12.25 - Meeting Micha
- Matlab Lizenz PC
- Frage zur variablen Mischhöhe
- Literatur (Probleme) zur Konvektionsmodellierung via energetischem Ansatz. Wärmeatlas VDI etc. immer fluiddynamischer Ansatz.

## 24.11.25 - Meeting Dominic

**Rückwirkend auf letzten Termin**
- Konkrete Aufgabe bei ortsfester Entnahme/Abgabe bzw. Unterschied zu jetzt?
    - Zurückstellen, Nz(2) Nz(4) sind variabel, Summe konstant. 

**Verständnis der Kopplungsterme**
$$
[\sigma] = \frac{[\lambda] [A_{pist}]*[(T_{pist,0}-T_{pist,1})]}{[\rho_w] [c_w] A_{TPPS} [\Delta z^2]}
= \frac{[m^2 \frac{J}{s m K} K]}{\frac{kg}{m^3}\frac{J}{kg K}m^2 m^2}
= \frac{J m}{\frac{J m s}{K}} = \frac{K}{s}
$$
- Gleichung 6 und 7 im Paper. Der Temperaturstrom über den Rand der nur eine Funktion von z ist, sprich die Temperaturverteilung um einen konstanten Wert pro Gitterpunkt ändert.
- Gleichung 4 und 5 im Paper. Der Temperaturstrom über den Rand der nur **keine** Funktion von z ist, sprich nur an einem Gitterpunkt angreift?
- Gleichung 8, eine gewichteter Mittelwert zwischen T_1 und T_2, wo wird er verwendet?
    - Instante Kontakttemperatur zweier Halb-unendlicher Körper. 

**Exergie HWS**
- Fragen zur Exergie Berechnung
- Frage zur Grafik 10, Grün ist die potentielle Lage Energie?
    - Lage Energie des Kolbens.

**Frage zur Breite des Gitters im Erdreich**
- Im Paper steht "11 nodes starting with an inital spacing of 0.005m and successively doubling the node spacing with each increase"
$$
 \sum_{i=1}^{10} 0.005m* 2^i = 0.005m * 2046 = 10.25m
$$
- Aber das Gitter im Bild zeigt nur eine Breite von 5 Metern im rechten Erdreich?

**Mischzone**
- Im Schäfer (2021) Paper ist die Mischzone fix, zwischen den Einlass und Auslassstutzen. Im TPPS prinzipiell über komplettes Volumen /Wassersäule?

**Code**
- timesteps.m: Zeile 40 Warum Viertelstunden und Stunde?
    - Zwei Zeitrechnungen für Vergleichbarkeit zwischen Comsol und Matlab
- Init.m: Zeile 8 - 4 jahre mit 120 Tagen Initialisierungs lag; Allgemein Init.m Skript nur zur Initialsimulation/Aufladen des Systems - Richtig?
- Szenarion.m: Zeile 51 - Anzahl an Viertelstunden in denen Geladen bzw. Entladen wird ; Zeile 64 - Berechnung der Strömungsgeschwindigkeit für komplette Ladung/Entladung; Zeile 118 Existiert eine Isolierung in radiale Richtung? - Ja!

**Arbeitsweise bis Freitag**
Integration der freien Konvektion ins Modell (rechnerisch) und in den bestehnden Code?

**Fragen an Micha**
- Zugriff auf PC mit der Nummer 17? Zweck Comsol Betrachtung der Initialisierung.

## 22.11.25 – Literature Review - Convective Term in PDE

### Tasks
- Get familiar with natural convection term in Gerle's paper.
- Link it to the work of Häuslein (2024).
- Inspect the source code for scenario 1, understand the .mat files.

---

## 21.11.25 - Project setup & literature management 

### Tasks
- Set up Git infrastructure.
- Created OneDrive structure for BA_bigData.
- Configured Zotero with linked attachments.
- Installed BetterBibTeX.
- Set up automatic export for `bibliography.bib`.

### Insights
- `.bib` is automatically updated only on Uni PC.
- All PDFs must be stored in the OneDrive literature folder, same as big data files.
- Zotero sync replicates collections and metadata across both PCs.

### Next steps
- Start writing structured literature notes.
- Begin analysing the first TPPS equations.
- Create the initial model structure draft.

---
---
## 16.11.25 - Meeting Dominic / Micha
### Festgelegte Ziele und Meilensteine
**Mindestanforderungen**
- Implementierung der freien Konvektion in das bestehende 1D Modell.
    - Ansatz über energetische statt CFD Betrachtung.
- Berücksichtigung des Ringspalts in das TPPS Modell.
    - Ansatz über 1 HS, keine räumliche Auflösung!

**Optional wenn zeitlich Erreichbar**
- Ortsfeste Ein- und Auslassstutzen. Bisher verschieben sie sich mit dem Kolben?
    - Bisher 4 Ein/Auslass, ganz oben und ganz unten. 

**Grober Zeitplan**
- Einarbeitung bis Ende November (Literatur, freie Konvektion + Matlab Code) *[bis Ende November]*
- Formulierung eines freien Konvektionsmodells + Implementierung, sowie Integration fester Ein und Auslassventile *[bis Anfang/Mitte Januar]*
    - Eventuelle Neuformulierung notwendig.
- Implementierung Ringspalt *[Ende Februar]*
- Schreiben *[März]*

**Meilensteine**
 - Vor Weihnachten für ortsfesten Kolben soll die freie Konvektion implementiert werden.
 - Bis Mitte Januar mit beweglichem Kolben.
 - Bis Mitte Februar ist das thermische Modell für den Ringspalt implementiert.

### Weiteres
**Penstock - Druckrohrleitung**
- Bisher Druckrohrleitung zur Hälfte dem oberen und zur Hälfte dem unteren Wasservolumen zugeordnet.

**Literatur**

→ Recherchieren welche Paper Abfallströmungen bei Behältern berücksichtigen bzw. ab wann diese relevant werden. In Literatur /Propadeaticum aufnehmen.

---





