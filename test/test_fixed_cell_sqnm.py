import numpy as np
import bazant
import sqnm.free_or_fixed_cell_sqnm


def _energyandforces(pos, alat):
    epot, forces, stress, deralat = bazant.energyandforces_bazant(alat, pos)
    return epot, forces


def _tests():
    from ase import io
    import sys

    filename = sys.argv[1]

    at = io.read(filename)
    pos = at.get_positions()   # (nat, 3)
    lat = at.get_cell().array  # (3, 3), rows are lattice vectors
    nat = at.get_global_number_of_atoms()

    alpha = -.01
    opt = sqnm.free_or_fixed_cell_sqnm.free_sqnm(nat, alpha, 10, 1e-2, 1e-3)

    for i in range(50):
        epot, forces = _energyandforces(pos, lat)
        pos = opt.optimizer_step(pos, epot, forces)
        print(epot, np.linalg.norm(forces, axis=1).max())


if __name__ == "__main__":
    _tests()
