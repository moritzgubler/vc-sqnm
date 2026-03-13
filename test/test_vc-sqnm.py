import numpy as np
import bazant

def _energyandforces(nat, pos, alat):
    epot, forces, stress, deralat = bazant.energyandforces_bazant(alat, pos)
    return epot, forces, deralat

def _rand_vec(nat, sigma):
    import random
    x = np.zeros((3, nat))
    for i in range(nat):
        for j in range(3):
            x[j, i] = random.gauss(0.0, sigma)
    return x

def _tests():
    from ase import io
    import sys
    import sqnm.periodic_sqnm

    filename = sys.argv[1]

    at = io.read(filename)
    pos = at.get_positions()   # (nat, 3)
    lat = at.get_cell().array  # (3, 3), rows are lattice vectors
    nat = at.get_global_number_of_atoms()
    alpha = -.01
    lattice_weight = 2.0

    opt_clean = sqnm.periodic_sqnm.periodic_sqnm(nat, lat, alpha, 10, lattice_weight, 1e-2, 1e-3, use_cupy=False)

    for i in range(50):
        epot, forces, deralat = _energyandforces(nat, pos, lat)
        pos, lat = opt_clean.optimizer_step(pos, lat, epot, forces, deralat)
        opt_clean.lower_bound()
        print(epot, max(np.linalg.norm(forces, axis=1).max(), np.abs(deralat).max()))


if __name__ == "__main__":
    _tests()
