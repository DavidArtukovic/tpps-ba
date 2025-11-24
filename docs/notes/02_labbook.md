# 02 - To Do's, Meetings, Fragen BA TPPS
Hier werden To Do's, Erkentnisse, Fragen sowie Protokolle zu den Treffen mit den Betreuern gesammelt. Oben stehen immer die anstehenden To Do's.

# Open Tasks
- [ ] Derive algebraic expressions for the inclusion of the convection term.
- [ ] Go through all scripts, update indexing and comments.
- [ ] Integrate the free convection into the existing Matlab code.

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

## 22.11.25 – Literature Review - Convective Term in PDE

### Tasks
- Get familiar with natural convection term in Gerle's paper.
- Link it to the work of Häuslein (2024).
- Inspect the source code for scenario 1, understand the .mat files.

---
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




