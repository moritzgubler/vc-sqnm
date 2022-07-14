module historylist
  use iso_c_binding
  implicit none
  type hist_list
  integer(c_int) :: nhistx
  integer(c_int) :: ndim
  integer(c_int), private :: icount
  logical, private :: is_initialized = .false.
  real(c_double), allocatable, dimension(:, :) :: list
  real(c_double), allocatable, dimension(:, :) :: diff_list
  real(c_double), allocatable, dimension(:, :) :: norm_diff_list
  real(c_double), allocatable, dimension(:) :: old_x
  contains
  procedure :: init
  procedure :: add
  procedure :: get_length

  end type hist_list
contains

  subroutine init(t, ndim, nhistx)
    class(hist_list) :: t
    integer(c_int), intent(in) :: ndim
    integer(c_int), intent(in) :: nhistx

    t%icount = 1
    t%ndim = ndim
    t%nhistx = nhistx
    allocate(t%list(ndim, nhistx), t%diff_list(ndim, nhistx))
    allocate(t%old_x(ndim), t%norm_diff_list(ndim, nhistx))
  end subroutine init

  subroutine add(t, x)
    class(hist_list) :: t
    real(c_double), intent(in) :: x(t%ndim)
    integer :: i
   
    if ( t%icount <= t%nhistx) then ! list not yet full
      t%list(:, t%icount) = x
      do i = 2, t%icount
        t%diff_list(:, i - 1) = t%list(:, i) - t%list(:, i - 1)
        t%norm_diff_list(:, i-1) = t%norm_diff_list(:, i-1) / norm2(t%norm_diff_list(:, i-1))
      end do
      t%icount = t%icount + 1
    else ! list is full
      t%old_x = t%list(:, 1)
      do i = 1, t%nhistx - 1
        t%list(:, i) = t%list(:, i + 1)                
      end do
      t%list(:, t%nhistx) = x
      t%diff_list(:, 1) =  t%list(:, 1) - t%old_x
      t%norm_diff_list(:, 1) = t%norm_diff_list(:, 1) / norm2(t%norm_diff_list(:, 1))
      do i = 2, t%nhistx
        t%diff_list(:, i) = t%diff_list(:, i) - t%diff_list(:, i - 1)
        t%norm_diff_list(:, i) = t%norm_diff_list(:, i) / norm2(t%norm_diff_list(:, i))
      end do
    end if
  end subroutine add

  integer function get_length(t)
    class(hist_list) :: t
    get_length = t%icount - 1
  end function get_length

end module historylist