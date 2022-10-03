# Graph einer Funktion zeichnen
# abgewandelt nach 
# example from https://matplotlib.org/gallery/
# Simple Function Plot 


import matplotlib                # modul um Daten darzustellen
import matplotlib.pyplot as plt
import numpy as np               # modul für Matrixrechnung und math. Funktionen 

# array mit Abszissenwerten (array = vector, matrix, ..) 
# np.arange(Anfangswert, Endwert, Schrittweite) 
x = np.arange(1.4e9, 1.6e9, 1e7)

R = 1.4960e11           # Abstand Sonne Erde
Gm = 3.986004418e14     # Gravitationsparameter Erde
GM = 1.32712442099e20   # Gravitationsparameter Sonne

# array y mit Ordinaten-/Funktionswerten berechnet aus dem array x 
y = 1/(R + x)**2 + Gm/GM * 1/x**2 - (R + x)/R**3

# Zeichenobjekte 
fig = plt.figure()
ax = plt.subplot()

# Funktionsgraph erstellen 
ax.plot(x, y)

# Koordinatenachsen erstellen 
ax.set(xlabel='x  (m)', ylabel='y(x)  (willkürliche Einheiten)',
       title='Funktionsgraph')

# Koordinatenraster 
ax.grid()

# Graph sichern
fig.savefig("Funktionsgraph.png", dpi=300)

# auf dem Bildschirm darstellen
plt.show()
