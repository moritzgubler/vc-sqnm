import numpy as np
import bazant

def _energyandforces(nat, pos, alat):
    import random
    epot, forces, stress, deralat = bazant.energyandforces_bazant(alat, pos)
    # bazant returns forces (nat,3) and deralat in row convention;
    # optimizer_step expects (3,nat) and column convention
    return epot, forces.T, deralat.T

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
    pos = at.get_positions()
    lat = at.get_cell().array
    nat = at.get_global_number_of_atoms()
    alpha = -.01
    lattice_weight = 2.0

    e0 = -1.3670604955980028

    opt_clean = sqnm.periodic_sqnm.periodic_sqnm(nat, lat.T, alpha, 10, lattice_weight, 1e-2, 1e-3, use_cupy=False)

    for i in range(30):
        epot, forces, deralat = _energyandforces(nat, pos, lat)
        t1 = time.time()
        pos, lat = opt_clean.optimizer_step(pos.T, lat.T, epot, forces, deralat)
        pos = pos.T
        lat = lat.T




        est1 = opt_clean.lower_bound()
        print(epot, max(np.linalg.norm(forces, axis=0).max(), np.abs(deralat).max()))


if __name__ == "__main__":
    _tests()