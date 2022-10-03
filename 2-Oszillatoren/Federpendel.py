# numerische Integration der Bahn eines Federpendels
# mit dem Euler-Cromer Verfahren (symplektisches Eulerverfahren)

import numpy as np
import matplotlib.pyplot as plt


dt = 0.03     # Zeitschritt in Sekunden
tmax = 10     # Zeitraum von 0 bis tmax (s)
n = int(np.round(tmax/dt)+1)  # Anzahl zu berechnende Punkte

t = np.zeros(n)   # array mit Zeiten, nummeriert von 0 bis n-1 
y = np.zeros(n)   # array mit Position
v = np.zeros(n)   # array mit Geschwindigkeit
a = np.zeros(n)   # array mit Beschleunigung

D = 1       # Federkonstante in Newton pro Meter
m = 0.1     # Pendelmasse in kg
y[0] = 0    # Startposition in Meter zur Zeit t=0
v[0] = 1    # Startgeschwindigkeit in m/s


# numerische Integration mit Euler-Cromer 
for i in range(0, n-1):  # für alle array-elemente von Nr. 0 bis n-2 
    t[i+1] = t[i] + dt   # ein Zeitschritt dt vorwärts
    a[i] = -D*y[i]/m     # Beschleunigung zur Zeit t[i]
    v[i+1] = v[i] + a[i]*dt     # Geschwindigkeit zur Zeit t[i+1]
    # y[i+1] = y[i] + v[i]*dt     # Position zur Zeit t[i+1] nach Euler
    y[i+1] = y[i] + v[i+1]*dt  # Position zur Zeit t[i+1] nach Euler-Cromer


# Zeichenobjekte 
fig, ax = plt.subplots()

# scatterplot spezifizieren 
ma = "."          # Punkt (Kreis) als Marke
colors = "green"  # Farbe der Marken
area = 5**2       # Grösse der Marken, Radius in points hoch zwei

ax.scatter(t, y, marker=ma, s=area, c=colors)

# Koordinatenachsen erstellen 
ax.set(xlabel='t  (s)', ylabel='y(t)  (m)')
# Koordinatenraster 
ax.grid()

# Theoriekurve dazu zeichnen
omega = (D/m)**0.5     # Kreisfrequenz der harmonischen Schwingung in 1/s
A = (m*v[0]**2/D)**0.5        # Amplitude der harmonischen Schwingung in Metern 
phi0 = np.arcsin(y[0]/A)     # Startphase in Radiant
z = A * np.sin(omega * t + phi0)   # rechnen mit dem ganzen t-array !

ax.plot(t, z, color="red", marker=None)

# Graph sichern, falls gewünscht
# fig.savefig("Fadenpendel.png")

# auf dem Bildschirm darstellen
plt.show()
