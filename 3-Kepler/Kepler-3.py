# numerische Integration einer Keplerbahn
# mit dem Euler-Cromer Verfahren (symplektisches Eulerverfahren)

# Bestimmung der Umlaufzeit und der Halbachsen aus den simulierten Daten
# damit das dritte keplersche Gesetz geprüft werden kann


import numpy as np
import matplotlib.pyplot as plt


def norm(lol):
    return np.linalg.norm(lol)


# Definition einer Funktion
def OrbitalData(X0, VY0):
    global dt, GM  # Konstanten
    global a, b, T  # output
    x1 = X0
    y1 = 0
    vx = 0
    vy = VY0
    T = 0   # kumulierte Zeit
    phi = 0  # kumulierter Polarwinkel
    xmax = 0
    xmin = 0
    ymax = 0
    ymin = 0
    while phi < 2*np.pi:
        x0 = x1
        y0 = y1

        if xmax < x0:
            xmax = x0
        if xmin > x0:
            xmin = x0
        if ymax < y0:
            ymax = y0
        if ymin > y0:
            ymin = y0

        r32 = (x0*x0+y0*y0)**(3/2)
        ax = -GM*x0/r32
        ay = -GM*y0/r32
        vx = vx+ax*dt
        vy = vy+ay*dt
        x1 = x0 + vx*dt
        y1 = y0 + vy*dt
        # Polarwinkel-Inkrement
        vec_a = [x0, y0]
        vec_b = [x1, y1]
        prod = norm(np.cross(vec_a, vec_b))
        # Winkel zw. Ortvektoren (x0,y0) und (x1,y1) AUSRECHNEN !!!!!!
        dphi = np.arcsin(prod/(norm(vec_a)*norm(vec_b)))
        phi = phi + dphi
        T = T + dt

    T = T-dt + (phi-2*np.pi)*dt/dphi  # linear interpolieren
    a = (xmax-xmin)/2    # grosse Halbachse
    b = (ymax-ymin)/2    # kleine Halbachse
    # end of function


dt = 86400     # Zeitschritt in Sekunden

GM = 1.32712442099e20       # Gravitationsparameter der Sonne in m^3/s^2

X0 = 4.9468e11     # Anfangsposition in Meter

VY0 = 7941.9       # Startgeschwindigkeit in m/s

OrbitalData(X0, VY0)  # Funktionsaufruf

print("a = ", a/1.495978707e11, " AE")
print("b = ", b/1.495978707e11, " AE")
print("T^2/a^3 = ", T**2/a**3)

print("lol", 4*np.pi**2/GM)
