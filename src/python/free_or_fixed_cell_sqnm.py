#!/usr/bin/env python3
import numpy as np
import sqnm
import sys


class free_sqnm:
    """
    Implementation of the sqnm method. Adapted to minimization of the PES.
    More informations about the algorithm can be found here: https://arxiv.org/abs/2206.07339
    """

    def __init__(self, nat, initial_step_size, nhist_max, alpha_min, eps_subsp):
        """
        Construct a optimizer object that can be used for free or fixed cell optimization.
         Parameters
        ----------
        nat: int
            Number of atoms
        init_lat: 3*3 numpy matrix 
            Matrix containing initial lattice vectors stored columnwise.
        initial_step_size: double
            initial step size. default is 1.0. For systems with hard bonds (e.g. C-C) use a value between and 1.0 and
            * 2.5. If a system only contains weaker bonds a value up to 5.0 may speed up the convergence.
        nhist_max: int
            Maximal number of steps that will be stored in the history list. Use a value between 3 and 20. Must be <= than 3*nat + 9.
        alpha_min: double
            Lower limit on the step size. 1.e-2 is the default.
        eps_subsp: double 
            Lower limit on linear dependencies of basis vectors in history list. Default 1.e-4.
        """

        self.nat = nat
        self.ndim = 3 * nat
        if self.ndim < nhist_max:
            print("Number of subspace dimensions bigger than number of dimensions.")
            print("Number of subspace dimensions will be reduced")
            nhist_max = self.ndim
        self.optimizer = sqnm.SQNM(self.ndim, nhist_max, initial_step_size, eps_subsp, alpha_min)
        self.fluct = 0.0

    def optimizer_step(self, pos, epot, forces):
        """
        Calculates new atomic coordinates that are closer to the local minimum. Free of fixed cell shape optimization.
        This function should be used the following way:
        1. calculate energies, forces and stress tensor at positions r.
        2. call the step function to update positions r.
        3. repeat.
        Parameters
        ----------
        pos: numpy matrix, dimension(3, nat)
            Input: atomic coordinates, dimension(3, nat). 
            Output: improved coordinates that are calculated based on forces from this and previous iterations.
        epot: double
            Potential energy of current geometry.
        forces: numpy matrix, dimension(3, nat)
            Forces of current geometry
        """

        # check for noise in forces using eq. 23 of vc-sqnm paper
        fnoise = np.linalg.norm(np.sum(forces, axis=1)) / np.sqrt(3 * self.nat)
        if self.fluct == 0.0:
            self.fluct = fnoise
        else:
            self.fluct = .8 * self.fluct + .2 * fnoise
        if self.fluct > 0.2 * np.max( np.abs(forces) ):
            print("""Warning: noise in forces is larger than 0.2 times the largest force component.
            Convergence is not guaranteed.""", file=sys.stderr)
        pos = pos.reshape(3 * self.nat)
        pos = pos + self.optimizer.sqnm_step(pos, epot, -forces.reshape(3 * self.nat))
        pos = pos.reshape((3, self.nat))
        return pos

    def lower_bound(self):
        """ Returns an estimate of a lower bound for the local minumum.
        The estimate is only accurate when the optimization is converged.
        """

        return self.optimizer.lower_bound()


# the rest of this file can be used for testing only

def _energyandforces(nat, pos, alat):
    import bazant
    epot, forces, deralat = bazant.energyandforces_bazant(alat, pos, nat)
    return epot, forces, deralat

def _tests():
    from ase import io
    b2a = Bohr_Ang = 0.52917721067

    filename = sys.argv[1]

    at = io.read(filename)
    pos = at.get_positions().T / b2a
    lat = at.get_cell().T / b2a
    nat = at.get_global_number_of_atoms()
    alpha = 2


    opt = free_sqnm(nat, alpha, 10, 1e-2, 1e-4)

    for i in range(30):
        epot, forces, deralat = _energyandforces(nat, pos, lat)
        print(epot, np.linalg.norm(forces), np.linalg.norm(deralat))
        pos = opt.optimizer_step(pos, epot, forces)

    print('The current energy is: ', epot)
    print('The estimated lower bound of the ground state is:', opt.lower_bound())
    print('The estimated energy error is:', epot - opt.lower_bound())

if __name__ == "__main__":
    _tests()
