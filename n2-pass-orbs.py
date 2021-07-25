#%%
"""checking if I can pass previous orbitals to the next iteration"""
from pyquante2 import rhf, basisset
from pyquante2.geo.molecule import molecule
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyquante2.scf.iterators import SCFIterator as SCFIterator
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

#%%
data = pd.DataFrame(np.linspace(0.6,2.5,25), columns=['r'])

basis_sets = ['sto-3g','6-31g']#, '6-31g**']

# creating dictionaries to go with the basis sets
orbitals = {}
orbital_energy = {}
potential_energy = {}

for bs in basis_sets:
    """looping through the basis sets"""
    orbitals[bs] = [np.zeros((1,1))] # initial orbital 
    # after trying to pass in a number of list and array dimensions
    # I found it need to have a have at least a np.shape(1,1)
    # and be iterable
    orbital_energy[bs] = [np.zeros((1,1))]
    # I don't think I can pass in the orbital energy 
    potential_energy[bs] = []
    # stores the potential after solving 
    for r in data['r']:
        """looping through the internuclear distance"""
        n2 = n2_molecule(r)
        bfs = basisset(n2,bs)
        solver = rhf(n2, bfs)
        #solver.converge(iterator=SCFIterator, c=orbitals[bs][-1])
        solver.converge()

        orbitals[bs].append(solver.orbs)
        orbital_energy[bs].append(solver.orbe)

        if solver.converged is True:
            potential_energy[bs].append(solver.energy)
        else:
            potential_energy[bs].append(np.nan)

#%%
for bs in basis_sets:
    plt.plot(data['r'], potential_energy[bs], '-o', label=bs)
    plt.legend()
#%%
potential_energy['sto-3g']
# %%
r_2 = data['r'][5:]

bs = basis_sets[0]

orbitals_2 = [orbitals[bs][5]]
orbital_energy_2 = [orbital_energy[bs][5]]
potential_energy_2 = [potential_energy[bs][5]]

for r in r_2[1:]:
    """looping through the internuclear distance"""
    print(r)
    n2 = n2_molecule(r)
    bfs = basisset(n2,bs)
    solver = rhf(n2, bfs)
    solver.converge(iterator=SCFIterator, c=orbitals_2[-1])
    #solver.converge()

    orbitals_2.append(solver.orbs)
    orbital_energy_2.append(solver.orbe)

    if solver.converged is True:
        potential_energy_2.append(solver.energy)
    else:
        potential_energy_2.append(np.nan)
# %%
for bs in basis_sets:
    plt.plot(data['r'], potential_energy[bs], '-o', label=bs)
plt.plot(r_2, potential_energy_2, '-o', label=bs+'_2')
plt.legend()
# %%
"""I think what I really need to do is modify hamiltonian to take in energy and c.
That way it checks the E - Eold < tol for Eold !=0 like it's hardcoded to start with."""
