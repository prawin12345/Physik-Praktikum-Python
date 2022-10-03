# numerische Integration einer Keplerbahn
# mit dem Euler-Cromer Verfahren (symplektisches Eulerverfahren)

# Test des zweiten, keplerschen Gesetzes (Flächensatz)


from numpy import *
import matplotlib.pyplot as plt

dt = 86400       # Zeitschritt in Sekunden
tmax = 8.0957e7  # Zeitraum von 0 bis tmax (s)
n = int(round(tmax/dt)+1)  # Anzahl zu berechnende Punkte

t = zeros(n)   # array mit Zeiten, nummeriert von 0 bis n-1
x = zeros(n)   # arrays mit Positionskomponenten
y = zeros(n)
vx = zeros(n)  # arrays mit Geschwindigkeitskomponenten
vy = zeros(n)
ax = zeros(n)  # arrays mit Beschleunigungskomponenten
ay = zeros(n)
dA = ones(n)   # Flächenstücke

GM = 1.32712442099e20       # Gravitationsparameter der Sonne in m^3/s^2
x[0] = 4.9468e11    # Anfangsposition eines Planetoiden in Meter
y[0] = 0
vx[0] = 0           # Start im Aphel/Perihel, wenn vx[0]=0 und y[0]=0
vy[0] = 7941.9      # Startgeschwindigkeit des Planetoiden in m/s


# numerische Integration mit Euler-Cromer
for i in range(0, n-1):  # für alle array-elemente von Nr. 0 bis n-2
    t[i+1] = t[i] + dt   # ein Zeitschritt dt vorwärts
    r = sqrt(x[i]*x[i]+y[i]*y[i])  # Abstand Sonne-Planet
    ax[i] = -GM*x[i]/r**3     # Beschleunigung zur Zeit t[i]
    ay[i] = -GM*y[i]/r**3
    vx[i+1] = vx[i] + ax[i]*dt     # Geschwindigkeit zur Zeit t[i+1]
    vy[i+1] = vy[i] + ay[i]*dt
    x[i+1] = x[i] + vx[i+1]*dt   # Position zur Zeit t[i+1]
    y[i+1] = y[i] + vy[i+1]*dt
    vec_a = [x[i+1], y[i+1]]
    vec_b = [x[i], y[i]]
    dA[i] = linalg.norm(cross(vec_a, vec_b))/2
    # Flächenstücke berechnen

dA[n-1] = dA[n-2]  # damit das letzte array-element ähnlich gross wird

dA = dA/dA[0] - 1   # skalieren, damit Zahlen nahe Null herauskommen

# Zeichenobjekte
fig, ax = plt.subplots()

# scatterplot spezifizieren
ma = "."        # Punkt (Kreis) als Marke
colors = "green"  # Farbe der Marken
area = 4**2     # Grösse der Marken, Radius in points hoch zwei

ax.scatter(t, dA, marker=ma, s=area, c=colors)  # dA vs t darstellen

# Koordinatenachsen erstellen
ax.set(xlabel='t  (s)', ylabel='dA/dA[0]-1')
# Koordinatenraster
# ax.grid()

# Graph sichern, falls gewünscht
# fig.savefig("Kepler-2.png")

# auf dem Bildschirm darstellen
plt.show()
