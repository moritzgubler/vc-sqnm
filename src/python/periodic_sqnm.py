#!/usr/bin/env python3
import numpy as np
import sqnm


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
        self.fluct = 0.0
        self.l = []
        self.sigsum = 0.0
        self.nsig = 0

    def optimizer_step(self, pos, alat, epot, forces, deralat):


        #fnoise = np.linalg.norm( np.sum(forces, axis=1) ) / np.sqrt(self.nat * 3)
        fnoise = np.abs(( np.sum(forces, axis=1)[0] )) / np.sqrt(self.nat)
        self.l.append(fnoise)
        self.sigsum += np.sum((np.sum(forces, axis=1)**2))
        self.nsig += 3
        print('estim', np.sqrt(self.sigsum / (self.nsig * self.nat)))
        #fnoise = np.linalg.norm( np.sum(forces, axis=1) )**2 / 3
        #print('asdf', np.linalg.norm(np.sum(forces, axis=1)) / np.sqrt(3))
        if self.fluct == 0.0:
            self.fluct = fnoise
        else:
            self.fluct = .8 * self.fluct + .2 * fnoise

        print('fluct, ', self.fluct, fnoise)

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
    import bazant
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

def _cleanup_forces(forces):
    sums = np.sum(forces, axis=1)
    forces[0, :] = forces[0, :] - sums[0] / np.shape(forces)[1]
    forces[1, :] = forces[1, :] - sums[1] / np.shape(forces)[1]
    forces[2, :] = forces[2, :] - sums[2] / np.shape(forces)[1]
    return forces

def _tests():
    from ase import io
    import sys
    import time
    b2a = Bohr_Ang = 0.52917721067

    filename = sys.argv[1]

    at = io.read(filename)
    pos = at.get_positions().T / b2a
    lat = at.get_cell().T / b2a
    nat = at.get_global_number_of_atoms()
    alpha = 2
    lattice_weight = 2.0
    nhist_max = 10


    opt_clean = periodic_sqnm(nat, lat, alpha, nhist_max, lattice_weight, 1e-2, 1e-4)
    opt_noise = periodic_sqnm(nat, lat, alpha, nhist_max, lattice_weight, 1e-2, 1e-4)
    sigma = 0.0001

    posnoise = pos.copy()
    latnoise = lat.copy()

    f = open('estimates.txt', 'w')
    # ground state energy of 64 Si bazant
    e0 = -10.936483964784031

    for i in range(6000):
        epot, forces, deralat = _energyandforces(nat, pos, lat)
        #pos, lat = opt_clean.optimizer_step(pos, lat, epot, forces, deralat)

        # f.write(str(epot - e0) + ' ' + str(abs(opt_clean.lower_limit())) + '\n')

        epotn, forcesn, deralatn = _energyandforces(nat, posnoise, latnoise)
        
        forcesn = forcesn + _rand_vec(nat, sigma)
        deralatn = deralatn + _rand_vec(3, sigma)
        #forcesn = _cleanup_forces(forcesn)
        posnoise, latnoise = opt_noise.optimizer_step(posnoise, latnoise, epotn, forcesn, deralatn)
        print(epot, max(np.max(np.linalg.norm(forces, axis=0)), np.linalg.norm(deralat)),
            epotn, max( np.max(np.linalg.norm(forcesn, axis=0)), np.linalg.norm(deralatn) ))
        
    print(np.sum(opt_noise.l) / 700)

if __name__ == "__main__":
    _tests()
