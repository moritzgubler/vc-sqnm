!! The variable cell shape optimization method is based on the following 
!! paper: https://arxiv.org/abs/2206.07339
!! More details about the SQNM optimization method are available here:
!! https://comphys.unibas.ch/publications/Schaefer2015.pdf
!! Author of this document: Moritz Gubler 

! Copyright (C) 2022 Moritz Gubler
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module periodic_optimizer
  use sqnm
  use iso_c_binding
  implicit none

private :: invertalat_lattice_per_opt

  type optimizer_periodic
    real(c_double) :: initial_lattice(3, 3)
    real(c_double) :: initial_lattice_inv(3, 3)
    real(c_double) :: lattice_transformer(3, 3)
    real(c_double) :: lattice_transformer_inv(3, 3)
    type(sqnm_optimizer) :: sqnm_opt
    integer(c_int) :: nat
    integer(c_int) :: ndim
    real(c_double) :: initial_step_size
    real(c_double) :: w
    real(c_double) :: f_sdt_deviation
    contains
    procedure :: initialize_optimizer
    procedure :: optimizer_step
    procedure :: get_lower_energy_bound
  end type optimizer_periodic

contains
  
  subroutine initialize_optimizer(t, nat, init_lat, initial_step_size, nhist_max &
      , lattice_weigth, alpha0, eps_subsp)
    !! This subroutine is used to set up the optimizer obect.
    class(optimizer_periodic) :: t
    integer(c_int), intent(in) :: nat
    !! Number of atoms
    real(c_double), intent(in) :: init_lat(3, 3)
    !! initial lattice vectors a = init_lat(:, 1)
    real(c_double), intent(in) :: initial_step_size
    !! initial step size. default is 1.0. For systems with hard bonds (e.g. C-C) use a value between and 1.0 and
    !! 2.5. If a system only contains weaker bonds a value up to 5.0 may speed up the convergence.
    integer(c_int), intent(in) :: nhist_max
    !! Maximal number of steps that will be stored in the history list. 
    !! Use a value between 3 and 20. Must be <= than 3*nat.
    real(c_double), intent(in) :: lattice_weigth
    !! weight or size of the supercell that is used to transform lattice derivatives. Use a value between 1 and 2. 
    !! Default is 2.
    real(c_double), intent(in) :: alpha0
    !! Lower limit on the step size. 1.e-2 is the default.
    real(c_double), intent(in) :: eps_subsp
    !! Lower limit on linear dependencies of basis vectors in history list. Default 1.e-4.
    !! Increase this parameter if energy or forces contain noise.

    integer(c_int) :: i

    t%nat = nat
    t%ndim = 3*nat + 9
    t%initial_lattice = init_lat
    t%initial_step_size = initial_step_size
    t%w = lattice_weigth
    t%f_sdt_deviation = 0.d0

    call invertalat_lattice_per_opt(t%initial_lattice, t%initial_lattice_inv)
    t%lattice_transformer = 0.d0
    do i = 1, 3
      t%lattice_transformer(i, i) = 1.d0 / norm2(t%initial_lattice(:, i))
    end do
    t%lattice_transformer = t%lattice_transformer * t%w * sqrt(dble(t%nat))
    call invertalat_lattice_per_opt(t%lattice_transformer, t%lattice_transformer_inv)

    call t%sqnm_opt%initialize_sqnm(t%ndim, nhist_max, t%initial_step_size, alpha0, eps_subsp)

  end subroutine initialize_optimizer  

  subroutine optimizer_step(t, r, alat, epot, f, deralat)
    !! This subroutine estimates the coordinates of the closest local minumum
    !! based on the previous points that were visited in the geometry optimization.
    class(optimizer_periodic) :: t
    real(c_double), intent(inout) :: r(3, t%nat)
    !! Positions of the atoms
    real(c_double), intent(inout) :: alat(3, 3)
    !! lattice vectors a = alat(:, 1)
    real(c_double), intent(in) :: epot
    !! potential energy
    real(c_double), intent(in) :: f(3, t%nat)
    !! forces 
    real(c_double), intent(in) :: deralat(3, 3)
    !! negative derivative of the pot, energy w.r. to the lattice vectors

    real(c_double) :: q(3, t%nat)
    real(c_double) :: df_dq(3, t%nat)
    real(c_double) :: a_tilde(3, 3)
    real(c_double) :: df_da_tilde(3, 3)
    real(c_double) :: a_inv(3, 3)
    real(c_double) :: q_and_lat(3, t%nat + 3)
    real(c_double) :: dq_and_dlat(3, t%nat + 3)
    real(c_double) :: dd(3, t%nat + 3)
    real(c_double) :: fnoise

    fnoise = norm2(sum(f, dim=2)) / sqrt(3.d0 * t%nat)
    if ( t%f_sdt_deviation == 0.d0 ) then
      t%f_sdt_deviation = fnoise
    else
      t%f_sdt_deviation = .8d0 * t%f_sdt_deviation + .2d0 * fnoise
    end if

    call invertalat_lattice_per_opt(alat, a_inv)

    ! transform atom coordinates and derivatives
    q = matmul(matmul(t%initial_lattice, a_inv), r)
    df_dq = - matmul(matmul(alat, t%initial_lattice_inv), f)

    ! transform lattice and derivatives
    a_tilde = matmul(alat, t%lattice_transformer)
    df_da_tilde = - matmul(deralat, t%lattice_transformer_inv)

    q_and_lat(:, :t%nat) = q
    q_and_lat(:, t%nat+1:t%nat+3) = a_tilde
    dq_and_dlat(:, :t%nat) = df_dq
    dq_and_dlat(:, t%nat+1:t%nat+3) = df_da_tilde

    call t%sqnm_opt%sqnm_step(q_and_lat, epot, dq_and_dlat, dd)
    !print*, 'dd', norm2(dd)

    q_and_lat = q_and_lat + dd

    q = q_and_lat(:, :t%nat)
    a_tilde = q_and_lat(:, t%nat+1:t%nat+3)

    alat = matmul(a_tilde, t%lattice_transformer_inv)
    r = matmul(matmul(alat, t%initial_lattice_inv), q)
  
  end subroutine optimizer_step

  subroutine invertalat_lattice_per_opt(alat, alatinv)
    !Invert alat matrix
    implicit real*8(a - h, o - z)
    dimension alat(3, 3), alatinv(3, 3)

    div = (alat(1, 1)*alat(2, 2)*alat(3, 3) - alat(1, 1)*alat(2, 3)*alat(3, 2) - &
           alat(1, 2)*alat(2, 1)*alat(3, 3) + alat(1, 2)*alat(2, 3)*alat(3, 1) + &
           alat(1, 3)*alat(2, 1)*alat(3, 2) - alat(1, 3)*alat(2, 2)*alat(3, 1))
    div = 1.d0/div
    alatinv(1, 1) = (alat(2, 2)*alat(3, 3) - alat(2, 3)*alat(3, 2))*div
    alatinv(1, 2) = -(alat(1, 2)*alat(3, 3) - alat(1, 3)*alat(3, 2))*div
    alatinv(1, 3) = (alat(1, 2)*alat(2, 3) - alat(1, 3)*alat(2, 2))*div
    alatinv(2, 1) = -(alat(2, 1)*alat(3, 3) - alat(2, 3)*alat(3, 1))*div
    alatinv(2, 2) = (alat(1, 1)*alat(3, 3) - alat(1, 3)*alat(3, 1))*div
    alatinv(2, 3) = -(alat(1, 1)*alat(2, 3) - alat(1, 3)*alat(2, 1))*div
    alatinv(3, 1) = (alat(2, 1)*alat(3, 2) - alat(2, 2)*alat(3, 1))*div
    alatinv(3, 2) = -(alat(1, 1)*alat(3, 2) - alat(1, 2)*alat(3, 1))*div
    alatinv(3, 3) = (alat(1, 1)*alat(2, 2) - alat(1, 2)*alat(2, 1))*div
  end subroutine invertalat_lattice_per_opt
  
  function get_lower_energy_bound(t) result(lower_bound)
    !! calculates an energy uncertainty (see eq. 20 of vc-sqnm paper)
    class(optimizer_periodic) :: t
    real(c_double) :: lower_bound
  
    lower_bound = t%sqnm_opt%get_lower_bound()
  
  end function get_lower_energy_bound

end module periodic_optimizer