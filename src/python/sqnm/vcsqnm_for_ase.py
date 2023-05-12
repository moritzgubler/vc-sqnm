#!/usr/bin/env python3

# The variable cell shape optimization method is based on the following 
# paper: https://arxiv.org/abs/2206.07339
# More details about the SQNM optimization method are available here:
# https://comphys.unibas.ch/publications/Schaefer2015.pdf
# Author of this document: Moritz Gubler 

# Copyright (C) 2022 Moritz Gubler
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sqnm.periodic_sqnm
import sqnm.free_or_fixed_cell_sqnm
import sys
import numpy as np
import logging


class aseOptimizer():

    initial_structure = None
    vc_relax = None
    optimizer = None
    nat = 0

    def __init__(self, initial_structure, vc_relax = False, force_tol=1e-2, maximalSteps = 500, initial_step_size = -0.01
            , nhist_max =10, lattice_weigth = 2.0, alpha_min = 1e-3, eps_subsp = 1e-3):
        self.initial_structure = initial_structure
        self.vc_relax = vc_relax
        self.force_tol=force_tol
        self.maximalSteps = maximalSteps

        self.nat = initial_structure.get_global_number_of_atoms()

        self._extractInfo(initial_structure)

        if self.vc_relax:
            self.optimizer = sqnm.periodic_sqnm.periodic_sqnm(self.nat, self.cell, initial_step_size, nhist_max, lattice_weigth, alpha_min, eps_subsp)
        else:
            self.optimizer = sqnm.free_or_fixed_cell_sqnm.free_sqnm(self.nat, initial_step_size, nhist_max, alpha_min, eps_subsp)


    def _extractInfo(self, atoms):
        self.positions = atoms.get_positions().T

        self.forces = atoms.get_forces().T
        self.energy = atoms.get_potential_energy()
        if self.vc_relax:
            self.cell = atoms.get_cell(True).T
            self.stress = atoms.get_stress(voigt=False)
            self.deralat = self._getLatticeDerivative()
    

    def _getLatticeDerivative(self):
        return - np.linalg.det(self.cell) * self.stress @ np.linalg.inv(self.cell).T


    def step(self, atoms):
        self._extractInfo(atoms)

        if self.vc_relax:
            self.positions, self.cell = self.optimizer.optimizer_step(self.positions, self.cell, self.energy, self.forces, self.deralat)
            atoms.set_cell(self.cell.T)
        else:
            self.positions = self.optimizer.optimizer_step(self.positions, self.energy, self.forces)

        atoms.set_positions(self.positions.T)

    
    def _getDerivativeNorm(self):
        forceNorm = np.max(np.abs(self.forces))
        if self.vc_relax:
            forceNorm = max(forceNorm, np.max(np.abs(self._getLatticeDerivative())))
        return forceNorm

    
    def optimize(self):
        # not tested
        i = 0
        while( i < self.maximalSteps and self._getDerivativeNorm() > self.force_tol ):
            logging.INFO("Relaxation step: %d energy: %f norm of forces: %f"% (i, self.initial_structure.get_potential_energy(), np.max(np.abs(self.initial_structure.get_forces()))) )
            self.step(self.initial_structure)
        return self.initial_structure

        