# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.bisect.html
# find the root of a function with bisection method 

import numpy as np
from scipy import optimize

R = 1.4960e11           # Abstand Sonne Erde
Gm = 3.986004418e14     # Gravitationsparameter Erde
GM = 1.32712442099e20   # Gravitationsparameter Sonne

# Funktion, von der eine Nullstelle gefunden werden soll 
def f(x):
    global R, Gm, GM
    y = 1/(R + x)**2 - Gm/GM * 1/x**2 - (R + x)/R**3
    return y



# array y mit Ordinaten-/Funktionswerten berechnet aus dem array x 

# Startwerte so, dass f(a)*f(b) < 0 ist
a = 1.47e9
b = 1.52e9

root = optimize.bisect(f, a, b) 

print("Nullstelle: ", root)
# 1501553645.3848324 ~ 1.5e9
print("Funktionswert dort: ", f(root)) 
