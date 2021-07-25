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
r_angstrom = np.linspace(0.6,2.5,25)

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
    for r in r_angstrom:
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
    plt.plot(r_angstrom, potential_energy[bs], '-o', label=bs)
    plt.legend()

# %%
r_2 = {basis_sets[0]: r_angstrom[5:], basis_sets[1]: r_angstrom[9:]}

orbitals_2 = {basis_sets[0]: [orbitals[basis_sets[0]][5]], basis_sets[1]: [orbitals[basis_sets[1]][9]]}
orbital_energy_2 = {basis_sets[0]: [orbital_energy[basis_sets[0]][5]], basis_sets[1]: [orbital_energy[basis_sets[1]][9]]}
potential_energy_2 = {basis_sets[0]: [potential_energy[basis_sets[0]][5]], basis_sets[1]: [potential_energy[basis_sets[1]][9]]}

for bs in basis_sets:
    print(bs)
    for r in r_2[bs][1:]:
        """looping through the internuclear distance"""
        n2 = n2_molecule(r)
        bfs = basisset(n2,bs)
        solver = rhf(n2, bfs)
        solver.converge(iterator=SCFIterator, c=orbitals_2[bs][-1])
        #solver.converge()

        orbitals_2[bs].append(solver.orbs)
        orbital_energy_2[bs].append(solver.orbe)

        if solver.converged is True:
            potential_energy_2[bs].append(solver.energy)
        else:
            potential_energy_2[bs].append(np.nan)
# %%
for bs in basis_sets:
    plt.plot(r_angstrom, potential_energy[bs], '-o', label=bs)
    plt.plot(r_2[bs], potential_energy_2[bs], '-o', label=bs+' w/ orb')
plt.legend()
# %%
"""I think what I really need to do is modify hamiltonian to take in energy and c.
That way it checks the E - Eold < tol for Eold !=0 like it's hardcoded to start with."""


# %%
