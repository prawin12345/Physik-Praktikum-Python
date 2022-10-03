# numerische Integration einer Keplerbahn
# mit dem Euler-Cromer Verfahren (symplektisches Eulerverfahren)

import numpy as np
import matplotlib.pyplot as plt


dt = 86400       # Zeitschritt in Sekunden
tmax = 8.0957e7  # Zeitraum von 0 bis tmax (s)
n = int(np.round(tmax/dt)+1)  # Anzahl zu berechnende Punkte

t = np.zeros(n)   # array mit Zeiten, nummeriert von 0 bis n-1
x = np.zeros(n)   # arrays mit Positionskomponenten
y = np.zeros(n)
vx = np.zeros(n)  # arrays mit Geschwindigkeitskomponenten
vy = np.zeros(n)
ax = np.zeros(n)  # arrays mit Beschleunigungskomponenten
ay = np.zeros(n)

GM = 1.32712442099e20       # Gravitationsparameter der Sonne in m^3/s^2
x[0] = 4.9468e11    # Anfangsposition eines Planetoiden in Meter
y[0] = 0
vx[0] = 0           # Start im Aphel/Perihel, wenn vx[0]=0 und y[0]=0
vy[0] = 7941.9      # Startgeschwindigkeit des Planetoiden in m/s


# numerische Integration mit Euler-Cromer
for i in range(0, n-1):  # für alle array-elemente von Nr. 0 bis n-2
    t[i+1] = t[i] + dt   # ein Zeitschritt dt vorwärts
    r = np.sqrt(x[i]*x[i]+y[i]*y[i])  # Abstand Sonne-Planet
    ax[i] = -GM*x[i]/r**3     # Beschleunigung zur Zeit t[i]
    ay[i] = -GM*y[i]/r**3
    vx[i+1] = vx[i] + ax[i]*dt     # Geschwindigkeit zur Zeit t[i+1]
    vy[i+1] = vy[i] + ay[i]*dt
    x[i+1] = x[i] + vx[i+1]*dt   # Position zur Zeit t[i+1]
    y[i+1] = y[i] + vy[i+1]*dt

# Zeichenobjekte
fig, ax = plt.subplots()

# scatterplot spezifizieren
ma = "."        # Punkt (Kreis) als Marke
colors = "green"  # Farbe der Marken
area = 5**2     # Grösse der Marken, Radius in points hoch zwei

ax.scatter(x, y, marker=ma, s=area, c=colors)

# Sonne einzeichnen
ax.scatter(0, 0, marker='.', s=8**2, c="royalblue")

# Koordinatenachsen erstellen
ax.set(xlabel='x  (m)', ylabel='y  (m)')
# Koordinatenraster
# ax.grid()

# gleicher Massstab auf beiden Achsen
ax.axis('equal')

a = (np.max(x) - np.min(x))/2
b = (np.max(y) - np.min(y))/2

# eine Ellipse dazu zeichnen
p = b**2/a   # halbes Quermass der Ellipse in m
eps = (1-b**2/a**2)**0.5    # numerische Exzentrizität der Ellipse
phi = np.arange(0, 2*np.pi, 0.01)
xe = p*np.cos(phi)/(1-eps*np.cos(phi))
ye = p*np.sin(phi)/(1-eps*np.cos(phi))

ax.plot(xe, ye, color="red", marker=None)

# Graph sichern, falls gewünscht
# fig.savefig("Kepler-1.png")

# auf dem Bildschirm darstellen
plt.show()
