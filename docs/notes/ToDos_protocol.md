# To Do's und Fragen und offene Punkte BA TPPS
Hier werden Erkentnisse, Fragen sowie Protokolle zu den Treffen mit den Betreuern gesammelt.

## Meeting Dominic / Micha - 16.11.25
### Festgelegte Ziele und Meilensteine
**Mindestanforderungen**
- Implementierung der freien Konvektion in das bestehende 1D Modell.
    - Ansatz über energetische statt CFD Betrachtung.
- Berücksichtigung des Ringspalts in das TPPS Modell.
    - Ansatz über 1 HS, keine räumliche Auflösung!

**Optional wenn zeitlich Erreichbar**
- Ortsfeste Ein- und Auslassstutzen. Bisher verschieben sie sich mit dem Kolben?
    - TODISCUSS

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



## Fragen Dominic - 24.11.25

**Rückwirkend auf letzten Termin**
- Konkrete Aufgabe bei ortsfester Entnahme/Abgabe bzw. Unterschied zu jetzt?

**Verständnis der Kopplungsterme**
$$
[\sigma] = \frac{[\lambda] [A_{pist}]*[(T_{pist,0}-T_{pist,1})]}{[\rho_w] [c_w] A_{TPPS} [\Delta z^2]}
= \frac{[m^2 \frac{J}{s m K} K]}{\frac{kg}{m^3}\frac{J}{kg K}m^2 m^2}
= \frac{J m}{\frac{J m s}{K}} = \frac{K}{s}
$$
- Gleichung 6 und 7 im Paper. Der Temperaturstrom über den Rand der nur eine Funktion von z ist, sprich die Temperaturverteilung um einen konstanten Wert pro Gitterpunkt ändert.
- Gleichung 4 und 5 im Paper. Der Temperaturstrom über den Rand der nur **keine** Funktion von z ist, sprich nur an einem Gitterpunkt angreift?
- Gleichung 8, eine gewichteter Mittelwert zwischen T_1 und T_2, wo wird er verwendet?

**Exergie HWS**
- Fragen zur Exergie Berechnung
- Frage zur Grafik 10, Grün ist die potentielle Lage Energie?

**Frage zur Breite des Gitters im Erdreich**
- Im Paper steht "11 nodes starting with an inital spacing of 0.005m and successively doubling the node spacing with each increase"
$$
 \sum_{i=1}^{10} 0.005m* 2^i = 0.005m * 2046 = 10.25m
$$
- Aber das Gitter im Bild zeigt nur eine Breite von 5 Metern im rechten Erdreich?

**Mischzone**
- Im Schäfer (2021) Paper ist die Mischzone fix, zwischen den Einlass und Auslassstutzen. Im TPPS prinzipiell über komplettes Volumen /Wassersäule?

**Code**





