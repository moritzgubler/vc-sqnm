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

module historylist
  use iso_c_binding
  implicit none
  type hist_list
  !! Historylist that is used by the sqnm class.
  !! More informations about the SQNM algorithm can be found here: https://aip.scitation.org/doi/10.1063/1.4905665
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
  procedure :: close_history_list

  end type hist_list
contains

  subroutine init(t, ndim, nhistx)
    !! initializes the historylist object.
    class(hist_list) :: t
    integer(c_int), intent(in) :: ndim
    !! dimension of the optimization problem
    integer(c_int), intent(in) :: nhistx
    !! maximal length of history list

    t%icount = 1
    t%ndim = ndim
    t%nhistx = nhistx
    allocate(t%list(ndim, nhistx), t%diff_list(ndim, nhistx))
    allocate(t%old_x(ndim), t%norm_diff_list(ndim, nhistx))
    t%is_initialized = .true.
  end subroutine init

  subroutine add(t, x)
    !! Add a vector to the history list.
    class(hist_list) :: t
    real(c_double), intent(in) :: x(t%ndim)
    !! Vector to add.
    integer :: i
   
    if ( t%icount <= t%nhistx) then ! list not yet full
      t%list(:, t%icount) = x
      do i = 2, t%icount
        t%diff_list(:, i - 1) = t%list(:, i) - t%list(:, i - 1)
        t%norm_diff_list(:, i-1) = t%diff_list(:, i-1) / norm2(t%diff_list(:, i-1))
      end do
      t%icount = t%icount + 1
    else ! list is full
      t%icount = t%nhistx + 2
      t%old_x = t%list(:, 1)
      do i = 1, t%nhistx - 1
        t%list(:, i) = t%list(:, i + 1)                
      end do
      t%list(:, t%nhistx) = x
      t%diff_list(:, 1) =  t%list(:, 1) - t%old_x
      t%norm_diff_list(:, 1) = t%diff_list(:, 1) / norm2(t%diff_list(:, 1))
      do i = 2, t%nhistx
        t%diff_list(:, i) = t%list(:, i) - t%list(:, i - 1)
        t%norm_diff_list(:, i) = t%diff_list(:, i) / norm2(t%diff_list(:, i))
      end do
    end if
  end subroutine add

  integer function get_length(t)
  !! returns the length of the historylist.
    class(hist_list) :: t
    get_length = t%icount - 2
  end function get_length

  subroutine close_history_list(t)
    !! closes the history list and deallocates the memory.
    class(hist_list) :: t
    if (t%is_initialized) then
      deallocate(t%list, t%diff_list, t%norm_diff_list, t%old_x)
      t%is_initialized = .false.
    end if
  end subroutine close_history_list

end module historylist