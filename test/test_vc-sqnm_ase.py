import sys
from ase import io
import bazant
import sqnm.vcsqnm_for_ase


def _tests():
    filename = sys.argv[1]

    at = io.read(filename)
    
    atlist = []
    atlist.append(at.copy())

    at.calc = bazant.BazantCalculator()

    opt = sqnm.vcsqnm_for_ase.aseOptimizer(at, vc_relax=True, initial_step_size=-0.01,
                                             nhist_max=10, lattice_weigth=2.0,
                                             alpha_min=1e-2, eps_subsp=1e-3)

    for i in range(30):
        opt.step(at)
        print(at.get_potential_energy(), opt._getDerivativeNorm())
        atlist.append(at.copy())

    io.write("trajectory_" + str(len(at)) + ".extxyz", atlist, format="extxyz")


if __name__ == "__main__":
    _tests()
