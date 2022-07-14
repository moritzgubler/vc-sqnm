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
    real(c_double), allocatable, dimension(:, :) :: s
    !! overlap matrix of historylist
    real(c_double), allocatable, dimension(:, :) :: h_subsp
    real(c_double), allocatable, dimension(:, :) :: h_evec_subsp
    real(c_double), allocatable, dimension(:, :) :: h_evec
    real(c_double), allocatable, dimension(:) :: h_eval
    real(c_double), allocatable, dimension(:) :: res
    real(c_double), allocatable, dimension(:) :: res_temp
    real(c_double) :: gainratio
    integer :: nhist
    contains
    procedure :: initialize
    procedure :: step
  end type sqnm_optimizer
contains

subroutine initialize(t, ndim, nhistx)
  class(sqnm_optimizer) :: t
  integer(c_int) :: ndim
  integer(c_int) :: nhistx
  
  allocate(t%s(nhistx, nhistx))
  allocate(t%prev_df_dx(ndim))

  
end subroutine initialize

subroutine step(t, x, f_of_x, df_dx, dir_of_descent)
  class(sqnm_optimizer) :: t
  real(c_double), intent(in) :: x(t%ndim)
  real(c_double), intent(in) :: f_of_x
  real(c_double), intent(in) :: df_dx(t%ndim)
  real(c_double), intent(out) :: dir_of_descent(t%ndim)

  real(c_double), allocatable, dimension(:, :) :: s_evec
  real(c_double), allocatable, dimension(:) :: s_eval
  real(c_double), allocatable, dimension(:, :) :: dr_subsp
  integer :: dim_subsp
  integer :: i

  ! lapack variables
  INTEGER :: info
  INTEGER :: lwork
  REAL(8), DIMENSION(:), ALLOCATABLE:: work

  call t%x_list%add(x)
  call t%flist%add(df_dx)
  t%nhist = t%x_list%get_length()

  if ( t%nhist == 0 ) then !! first step
    dir_of_descent = - t%alpha * df_dx
    t%dir_of_descent = dir_of_descent
  else
    ! calculate gainratio
    t%gainratio = (f_of_x - t%prev_f) / (.5d0 * dot_product(t%dir_of_descent, t%prev_df_dx))
    if (t%gainratio < 0.5d0 ) t%alpha = max(t%alpha0, t%alpha * 0.65d0)
    if (t%gainratio > 1.05d0) t%alpha = t%alpha * 1.05d0


    allocate(s_evec(t%nhist, t%nhist), s_eval(t%nhist))
    ! calculate overlab matrix of basis
    s_evec =  matmul(transpose(t%x_list%norm_diff_list(:, :t%nhist))&
    , t%x_list%norm_diff_list(:, :t%nhist))
    lwork = 100 * t%nhist
    allocate(work(lwork))
    call dsyev('v', 'l', t%nhist, s_evec, t%nhist, s_eval, work, lwork, info)

    dim_subsp = 0
    do i = 1, t%nhist
      if (s_eval(i) / s_eval(1) > t%eps_subsp) dim_subsp = dim_subsp + 1
    end do

    allocate(dr_subsp(t%ndim, dim_subsp))

  end if
  t%prev_f = f_of_x
  t%prev_df_dx = df_dx
  
end subroutine step

end module sqnm