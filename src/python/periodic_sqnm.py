#!/usr/bin/env python3
import numpy as np
import sqnm
import bazant
from ase import io
import sys
import time

class periodic_sqnm:

    def __init__(self, nat, init_lat, initial_step_size, nhist_max, lattice_weigth, alpha_min, eps_subsp):
        self.nat = nat
        self.ndim = 3 * nat + 9
        self.lattice_weight = lattice_weigth
        self.initial_lat = init_lat
        self.initial_lat_inverse = np.linalg.inv(init_lat)
        self.lattice_transformer = np.diag(1 / np.linalg.norm(self.initial_lat, axis=0)) * self.lattice_weight * np.sqrt(nat)
        self.lattice_transformer_inv = np.linalg.inv(self.lattice_transformer)
        self.optimizer = sqnm.SQNM(self.ndim, nhist_max, initial_step_size, eps_subsp, alpha_min)

    def optimizer_step(self, pos, alat, epot, forces, deralat):
        a_inv = np.linalg.inv(alat)

        q = ((self.initial_lat @ a_inv) @ pos).reshape(3 * self.nat)
        df_dq = (- (alat @ self.initial_lat_inverse) @ forces).reshape(3 * self.nat)

        a_tilde = (alat @ self.lattice_transformer).reshape(9)
        df_da_tilde = (- deralat @ self.lattice_transformer).reshape(9)

        q_and_lat = np.concatenate((q, a_tilde))
        dq_and_dlat = np.concatenate((df_dq, df_da_tilde))

        dd = self.optimizer.sqnm_step(q_and_lat, epot, dq_and_dlat)

        q_and_lat = q_and_lat + dd

        q = q_and_lat[:(3 * self.nat)].reshape(3, self.nat)
        a_tilde = q_and_lat[(3*self.nat):].reshape(3, 3)

        alat = a_tilde @ self.lattice_transformer_inv
        pos = (alat @ self.initial_lat_inverse) @ q

        return pos, alat

    def lower_limit(self):
        return self.optimizer.lower_limit()


def _energyandforces(nat, pos, alat):
    epot, forces, deralat = bazant.energyandforces_bazant(alat, pos, nat)
    return epot + 1.3670604955980028, forces, deralat

def _tests():
    b2a = Bohr_Ang = 0.52917721067

    filename = sys.argv[1]

    at = io.read(filename)
    pos = at.get_positions().T / b2a
    lat = at.get_cell().T / b2a
    nat = at.get_global_number_of_atoms()
    alpha = 2


    opt = periodic_sqnm(nat, lat, alpha, 10, 2.0, 1e-2, 1e-4)

    for i in range(100):
        epot, forces, deralat = _energyandforces(nat, pos, lat)
        #print(epot, np.linalg.norm(forces), np.linalg.norm(deralat))
        pos, lat = opt.optimizer_step(pos, lat, epot, forces, deralat)
        print(epot, opt.lower_limit(), np.linalg.norm(forces), np.linalg.norm(deralat))

if __name__ == "__main__":
    _tests()
