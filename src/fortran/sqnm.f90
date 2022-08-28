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


module sqnm
  use iso_c_binding
  use historylist
  implicit none

  type sqnm_optimizer
    integer(c_int) :: nhistx
    integer(c_int) :: ndim
    real(c_double) :: eps_subsp = 1.d-4
    real(c_double) :: alpha0 = 1.d-2
    type(hist_list) :: x_list
    type(hist_list) :: flist
    real(c_double) :: alpha
    real(c_double), allocatable, dimension(:) :: dir_of_descent
    real(c_double) :: prev_f
    !! previous value of target function
    real(c_double), allocatable, dimension(:) :: prev_df_dx
    !! previous derivative of target function
    real(c_double), allocatable, dimension(:, :) :: s_evec
    real(c_double), allocatable, dimension(:) :: s_eval
    real(c_double), allocatable, dimension(:, :) :: dr_subsp
    real(c_double), allocatable, dimension(:, :) :: df_subsp
    real(c_double), allocatable, dimension(:, :) :: h_evec_subsp
    real(c_double), allocatable, dimension(:, :) :: h_evec
    real(c_double), allocatable, dimension(:) :: h_eval
    real(c_double), allocatable, dimension(:) :: res
    real(c_double), allocatable, dimension(:) :: res_temp
    real(c_double) :: gainratio
    integer :: nhist

    INTEGER :: lwork
    REAL(8), DIMENSION(:), ALLOCATABLE:: work

    contains
    procedure :: initialize_sqnm
    procedure :: sqnm_step
  end type sqnm_optimizer
contains

subroutine initialize_sqnm(t, ndim, nhistx, alpha, alpha0, eps_subsp)
  class(sqnm_optimizer) :: t
  integer(c_int) :: ndim
  integer(c_int) :: nhistx
  real(c_double) :: alpha
  real(c_double) :: alpha0
  real(c_double) :: eps_subsp
  
  t%ndim = ndim
  t%nhistx = nhistx
  t%alpha = alpha
  t%alpha0 = alpha0
  t%eps_subsp = eps_subsp
  call t%x_list%init(ndim, nhistx)
  call t%flist%init(ndim, nhistx)

  allocate(t%s_evec(t%nhistx, t%nhistx), t%s_eval(t%nhistx))
  allocate(t%prev_df_dx(ndim))
  allocate(t%dr_subsp(t%ndim, nhistx))
  allocate(t%df_subsp(t%ndim, nhistx))
  allocate(t%h_evec_subsp(t%nhistx, t%nhistx))
  allocate(t%h_eval(nhistx))
  allocate(t%h_evec(t%ndim, t%nhistx))
  allocate(t%res(nhistx))
  allocate(t%res_temp(ndim))
  allocate(t%dir_of_descent(ndim))

  t%lwork = 100 * t%nhistx
  allocate(t%work(t%lwork))
  
end subroutine initialize_sqnm

subroutine sqnm_step(t, x, f_of_x, df_dx, dir_of_descent)
  class(sqnm_optimizer) :: t
  real(c_double), intent(in) :: x(t%ndim)
  real(c_double), intent(in) :: f_of_x
  real(c_double), intent(in) :: df_dx(t%ndim)
  real(c_double), intent(out) :: dir_of_descent(t%ndim)


  integer :: dim_subsp
  integer :: i, ihist, k, j

  ! lapack variables
  INTEGER :: info

  call t%x_list%add(x)
  call t%flist%add(df_dx)
  t%nhist = t%x_list%get_length()

  if ( t%nhist == 0 ) then !! first step
    t%dir_of_descent = - t%alpha * df_dx
  else
    ! calculate gainratio
    t%gainratio = (f_of_x - t%prev_f) / (.5d0 * dot_product(t%dir_of_descent, t%prev_df_dx))
    if (t%gainratio < 0.5d0 ) t%alpha = max(t%alpha0, t%alpha * 0.65d0)
    if (t%gainratio > 1.05d0) t%alpha = t%alpha * 1.05d0


    ! calculate overlab matrix of basis
    t%s_evec(:t%nhist, :t%nhist) =  matmul(transpose(t%x_list%norm_diff_list(:, :t%nhist)) &
    , t%x_list%norm_diff_list(:, :t%nhist))
    call dsyev('v', 'u', t%nhist, t%s_evec(:t%nhist, :t%nhist), t%nhist, t%s_eval(:t%nhist) &
      , t%work, t%lwork, info)
    if (info /= 0) stop 's dsyev'

    dim_subsp = 0
    do i = 1, t%nhist
      if (t%s_eval(i) / t%s_eval(1) > t%eps_subsp) then
        dim_subsp = dim_subsp + 1
      else
        print*, 'remove dimension'
      end if
    end do
    t%s_eval(1:dim_subsp) = t%s_eval((t%nhist - dim_subsp + 1):)
    t%s_evec(:, 1:dim_subsp) = t%s_evec(:, (t%nhist - dim_subsp + 1):)

    ! compute eq. 11
    t%dr_subsp(:,:dim_subsp) = 0.d0
    t%df_subsp(:,:dim_subsp) = 0.d0
    do i = 1, dim_subsp
      do ihist = 1, t%nhist
        t%dr_subsp(:, i) = t%dr_subsp(:, i) + t%s_evec(ihist, i) * t%x_list%norm_diff_list(:, ihist)
        t%df_subsp(:, i) = t%df_subsp(:, i) + t%s_evec(ihist, i) * t%flist%diff_list(:, ihist) &
          / norm2(t%x_list%diff_list(:, ihist)) 
      end do
      t%dr_subsp(:, i) = t%dr_subsp(:, i) / sqrt(t%s_eval(i))
      t%df_subsp(:, i) = t%df_subsp(:, i) / sqrt(t%s_eval(i))
    end do

    !! compute eq. 13
    t%h_evec_subsp(:dim_subsp, :dim_subsp) = .5d0 * (matmul(transpose(t%df_subsp(:,:dim_subsp)), t%dr_subsp(:,:dim_subsp)) &
        + matmul(transpose(t%dr_subsp(:,:dim_subsp)), t%df_subsp(:,:dim_subsp)))
    call dsyev('v', 'l', dim_subsp, t%h_evec_subsp(:dim_subsp, :dim_subsp)&
      , dim_subsp, t%h_eval(:dim_subsp), t%work, t%lwork, info)
    if (info  /= 0 ) stop 'h_eval dsyev'

    ! compute eq. 15
    t%h_evec = 0.d0
    do i = 1, dim_subsp
      do k = 1, dim_subsp
        t%h_evec(:, i) = t%h_evec(:, i) + t%h_evec_subsp(k, i) * t%dr_subsp(:, k)
      end do
    end do

    ! compute eq. 20
    do j = 1, dim_subsp
      t%res_temp = - t%h_eval(j) * t%h_evec(:, j)
      do k = 1, dim_subsp
        t%res_temp = t%res_temp + t%h_evec_subsp(k, j) * t%df_subsp(:, k)
      end do
      t%res(j) = norm2(t%res_temp)
    end do

    ! modify eigenvalues (eq. 18)
    do i = 1, dim_subsp
      t%h_eval(i) = sqrt(t%h_eval(i)**2 + t%res(i)**2)
    end do

    !print*, 'meigenvaluse', t%h_eval(:dim_subsp)

    ! decompose gradient (eq. 16)
    t%dir_of_descent = df_dx
    do i = 1, dim_subsp
      t%dir_of_descent = t%dir_of_descent &
        - sum(t%h_evec(:, i)*df_dx) * t%h_evec(:, i)
    end do
    t%dir_of_descent = t%dir_of_descent * t%alpha

    ! apply preconditioning to remaining gradien (eq. 21)
    do i = 1, dim_subsp
      t%dir_of_descent = t%dir_of_descent & 
        + sum(df_dx * t%h_evec(:, i)) * t%h_evec(:, i) / t%h_eval(i)
    end do

    t%dir_of_descent = - t%dir_of_descent

  end if

  dir_of_descent = t%dir_of_descent
  t%prev_f = f_of_x
  t%prev_df_dx = df_dx
  
end subroutine sqnm_step

end module sqnm