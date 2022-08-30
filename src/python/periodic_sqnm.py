from mimetypes import init
import numpy as np
import sqnm

class periodic_sqnm:

    def __init__(self, nat, init_lat, initial_step_size, nhist_max, lattice_weigth, alpha_min, eps_subsp):
        self.nat = nat
        self.ndim = 3 * nat + 9
        self.lattice_weight = lattice_weigth
        self.initial_lat = init_lat
        self.initial_lat_inverse = np.linalg.inv(init_lat)
        self.lattice_transformer = np.zeros((3,3))
        for i in range(3):
            self.lattice_transformer[i, i] = 1 / np.linalg.norm(self.initial_lat[:, i])
        self.lattice_transformer = self.lattice_transformer * self.lattice_weight * np.sqrt(nat)
        self.lattice_transformer_inv = np.linalg.inv(self.lattice_transformer)
        self.optimizer = sqnm.SQNM(self.ndim, nhist_max, initial_step_size, eps_subsp, alpha_min)

    def optimizer_step(self, pos, alat, epot, forces, deralat):
        a_inv = np.linalg.inv(alat)

        q = ((self.initial_lat @ a_inv) @ pos).reshape(3 * self.nat)
        df_dq = (- (alat @ self.initial_lat) @ forces).reshape(3 * self.nat)

        a_tilde = (alat @ self.lattice_transformer).reshape(9)
        df_da_tilde = (- deralat @ self.lattice_transformer).reshape(9)

        q_and_lat = np.concatenate((q, a_tilde))
        dq_and_dlat = np.concatenate((df_dq, df_da_tilde))

        dd = self.optimizer.sqnm_step(q_and_lat, epot, dq_and_dlat)

        q_and_lat = q_and_lat + dd

        q = q_and_lat[:(3 * self.nat)]
        a_tilde = q_and_lat[(3*self.nat):]

        alat = a_tilde @ self.lattice_transformer_inv
        pos = (alat @ self.initial_lat_inverse) @ q

        return pos, alat


