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
    real(c_double) :: f_prev
    real(c_double), allocatable, dimension(:) :: prev_df_dx
    real(c_double), allocatable, dimension(:, :) :: h_subsp
    real(c_double), allocatable, dimension(:, :) :: h_evec_subsp
    real(c_double), allocatable, dimension(:, :) :: h_evec
    real(c_double), allocatable, dimension(:) :: h_eval
    real(c_double), allocatable, dimension(:) :: res
    real(c_double), allocatable, dimension(:) :: res_temp
    contains
    procedure :: step
  end type sqnm_optimizer
contains

subroutine step(t, x, f_of_x, df_dx, dir_of_descent)
  class(sqnm_optimizer) :: t
  real(c_double), intent(in) :: x(t%ndim)
  real(c_double), intent(in) :: f_of_x
  real(c_double), intent(in) :: df_dx(t%ndim)
  real(c_double), intent(out) :: dir_of_descent(t%ndim)

  
end subroutine step

end module sqnm