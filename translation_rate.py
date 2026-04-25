'''
Okay so! A few to dos here:
1. Validate this with actual known translation rates. With dummy variables it seems
    somewhat reasonable (i.e. mRNA decreases with time, protein levels spike and then fall)
    but I would like to validate further.
2. Find actual values that apply to our model. Some are already done, I will keep looking! :D
3. Re-factor some of this code with main.py if desired (i.e. once we get the sequence designed, we can
    auto-plug it in here)
4. Make the model more accurate/complex. ODE logic is based on: https://vixra.org/pdf/2412.0171v1.pdf
    But more recent literature shows that amino acid sequence also impacts translation rate.
    Depends how complicated we want to get, could introduce ML if we want to get fancy.
    We're really assuming first-order kinetics for now and keeping it simple but it's better than nothing.
'''

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def translation_model(t, y, kinit, kexit, kdecay, kbind, kelong, kdegrade):
    M, R, P = y
    
    dR_dt = kinit * M - kexit * R
    dM_dt = -kdecay * M - kbind * R * M
    dP_dt = kelong * R - kdegrade * P
    
    return [dM_dt, dR_dt, dP_dt]

#TO DO: Get actual values if not noted already, otherwise I just made stuff up lol
kinit = 0.2
kexit = 0.5 #TBD - protein dependent, can depend on length?
kdecay = 0.0006 
kbind = 0.1 #TBD - S.cerevisiae has an 80S ribosome, trying to find non-mammalian stats :/
kelong = 9.0
kdegrade = 0.3 #TBD - also protein dependent

#init: approx 80% have translation initiation rates below 0.2 / s, could be more conservative if we want
    #source: https://www.biorxiv.org/content/10.1101/490730v1.full.pdf
#decay: yeast mRNA half life is 20 mins -> convert to rate constant -> 0.0006 
    #source: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100205&ver=23&trm=100205&org=
#elong: elongation rate for yeast is ~6-9AAs per second, more recent value is 9AA
    #source: https://arxiv.org/pdf/2508.14997


#Assume 10 "units" of existing mRNA, no ribosome attachment or protein creation yet.
#[M(0), R(0), P(0)]
y0 = [10, 0, 0]

t_span = (0, 60) #TO DO: Figure out how many mins we want to model
t_eval = np.linspace(0, 60, 200)

solution = solve_ivp(
    translation_model,
    t_span,
    y0,
    args=(kinit, kexit, kdecay, kbind, kelong, kdegrade),
    t_eval=t_eval,
    method='RK45'
)

if not solution.success:
    raise RuntimeError(f"Solver failed: {solution.message}")

M, R, P = solution.y

plt.plot(solution.t, M, label="mRNA (M)")
plt.plot(solution.t, R, label="Ribosomes (R)")
plt.plot(solution.t, P, label="Protein (P)")
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()