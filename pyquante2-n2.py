#%%
from pyquante2 import rhf,uhf,basisset
from pyquante2.basis.data import basis
from pyquante2.geo.molecule import molecule
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%
def n2_molecule(r):
    """creates an N2 molecule using the pyquante.molecule class

    Args:
        r (float): internuclear distance between the molecules

    Returns:
        object: n2 molecule with the given internuclear distance
    """
    n2 = molecule([(7, 0, 0, r/2),
                (7, 0, 0, -r/2)],
                units='Angstrom',
                multiplicity=1,
                name='Nitrogen')
    return n2

# %%
data = pd.DataFrame()
r_angstroms = np.linspace(0.6, 2.5, 25)
data['r'] = r_angstroms
basis_sets = ['sto-3g','6-31g']#, '6-31g**']

for bs in basis_sets:
    """looping through the different basis sets"""
    print(bs) # for progress
    rhf_solver_energies = [] # rhf calculations
    uhf_solver_energies = [] # uhf calculations
    for r in r_angstroms:
        """looping through the internuclear distances"""
        n2 = n2_molecule(r) # build molecule at specific internuclear distance
        bfs = basisset(n2, bs) 

        rhf_solver = rhf(n2, bfs)
        rhf_solver.converge()

        if rhf_solver.converged is True:
            rhf_solver_energies.append(rhf_solver.energy)
        else:
            rhf_solver_energies.append(np.nan)

        uhf_solver = uhf(n2, bfs)
        uhf_solver.converge()

        if uhf_solver.converged is True:
            uhf_solver_energies.append(uhf_solver.energy)
        else:
            uhf_solver_energies.append(np.nan)

    data[f'rhf-{bs}'] = rhf_solver_energies
    data[f'uhf-{bs}'] = uhf_solver_energies
    

#%%
for bs in data.keys()[1:]:
    plt.plot(data['r'], data[bs], label=bs)
plt.legend()
# %%
