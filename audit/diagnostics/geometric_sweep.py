#!/usr/bin/env python3
"""Pure-python check (no CLASS needed): reference-parabola vs geometric-sweep exposure on the
constant-Lambda injection background, and the stored-vs-actual H' mismatch."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
beta = 1/12
Sg = lambda u: 4.5*u*(1-u)*(1-u*(1+3*u)/32)   # closed form at beta=1/12: S_ref*[1-(3/2) beta I_ref u]
print("geometric crossings:", round(brentq(lambda u: Sg(u)-1, .3, .5), 12), round(brentq(lambda u: Sg(u)-1, .5, .7), 12), " (reviewer 0.349349763233 / 0.621228476566)")
h, wb, wc = 0.676255, 0.022539, 0.117409; Ob, Oc = wb/h**2, wc/h**2; Or = (2.469e-5/h**2)*(1+0.2271*3.046); OL = 1-Ob-Oc-Or
def rhs(N, y):
    a = np.exp(N); rc = y[0]; rm = Ob*a**-3+rc; u = OL/(OL+rm); I = .25+.75*u
    return [-3*rc + 4.5*beta*I*u*rm]
N = np.linspace(np.log(1/201), 0, 40000); s = solve_ivp(rhs, (N[0], 0), [Oc*201**3], t_eval=N, rtol=1e-11, atol=1e-14)
a = np.exp(N); rc = s.y[0]; rm = Ob*a**-3+rc; rr = Or*a**-4; rho = rm+OL+rr; H = np.sqrt(rho); u = OL/(OL+rm)
dlnH = np.gradient(np.log(H), N); q_geom = -1-dlnH
S_ref = 4.5*u*(1-u); S_geo = 3*(OL/rho)*(1+q_geom)
D = lambda S: np.clip(1-1/np.maximum(S, 1e-12), 0, None)
i = np.argmin(abs(u-2/3))
print(f"X_H reference parabola : {np.trapezoid(D(S_ref), N):.5f}   (production 0.03588)")
print(f"X_H geometric sweep    : {np.trapezoid(D(S_geo), N):.5f}   (CLASS with interaction_gate_geometric=1: 0.01968)")
print(f"dlnH/dlna at u=2/3     : actual {dlnH[i]:.5f}   stored-formula {(-1.5*(rm+4/3*rr)/rho)[i]:.5f}")
