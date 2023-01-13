import numpy as np
import bazant

def _energyandforces(nat, pos, alat):
    import random
    epot, forces, deralat = bazant.energyandforces_bazant(alat, pos, nat)
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
    import time
    import sqnm.periodic_sqnm

    b2a = 0.52917721067

    filename = sys.argv[1]

    at = io.read(filename)
    pos = at.get_positions().T / b2a
    lat = at.get_cell().T / b2a
    nat = at.get_global_number_of_atoms()
    alpha = 2
    lattice_weight = 2.0

    e0 = -1.3670604955980028

    opt_clean = sqnm.periodic_sqnm.periodic_sqnm(nat, lat, alpha, 10, lattice_weight, 1e-2, 1e-3, use_cupy=False)

    for i in range(30):
        epot, forces, deralat = _energyandforces(nat, pos, lat)
        t1 = time.time()
        pos, lat = opt_clean.optimizer_step(pos, lat, epot, forces, deralat)
        print(time.time() - t1)
        est1 = opt_clean.lower_bound()
        print(epot, max(np.linalg.norm(forces, axis=0).max(), np.abs(deralat).max()))


if __name__ == "__main__":
    _tests()