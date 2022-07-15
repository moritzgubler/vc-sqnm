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
    contains
    procedure :: initialize_optimizer
    procedure :: optimizer_step
  end type optimizer_periodic

contains
  
  subroutine initialize_optimizer(t, nat, init_lat, initial_step_size, nhist_max &
      , lattice_weigth, alpha0, eps_subsp)
    class(optimizer_periodic) :: t
    integer(c_int), intent(in) :: nat
    real(c_double), intent(in) :: init_lat
    real(c_double), intent(in) :: initial_step_size
    integer(c_int), intent(in) :: nhist_max
    real(c_double), intent(in) :: lattice_weigth
    real(c_double), intent(in) :: alpha0
    real(c_double), intent(in) :: eps_subsp

    integer(c_int) :: i

    t%nat = nat
    t%ndim = 3*nat + 9
    t%initial_lattice = init_lat
    t%initial_step_size = initial_step_size
    t%w = lattice_weigth

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
    class(optimizer_periodic) :: t
    real(c_double), intent(inout) :: r(3, t%nat)
    real(c_double), intent(inout) :: alat(3, 3)
    real(c_double), intent(in) :: epot
    real(c_double), intent(in) :: f(3, t%nat)
    real(c_double), intent(in) :: deralat(3, 3)

    real(c_double) :: q(3, t%nat)
    real(c_double) :: df_dq(3, t%nat)
    real(c_double) :: a_tilde(3, 3)
    real(c_double) :: df_da_tilde(3, 3)
    real(c_double) :: a_inv(3, 3)
    real(c_double) :: q_and_lat(3, t%nat + 3)
    real(c_double) :: dq_and_dlat(3, t%nat + 3)
    real(c_double) :: dd(3, t%nat + 3)

    call invertalat_lattice_per_opt(alat, a_inv)

    ! transform atom coordinates and derivatives
    q = matmul(matmul(t%initial_lattice, a_inv), r)
    df_dq = - matmul(matmul(alat, t%initial_lattice_inv), f)

    ! transform lattice and derivatives
    a_tilde = matmul(alat, t%lattice_transformer)
    df_da_tilde = - deralat * t%lattice_transformer_inv

    q_and_lat(:, :t%nat) = q
    q_and_lat(:, t%nat+1:t%nat+3) = a_tilde
    dq_and_dlat(:, :t%nat) = df_dq
    dq_and_dlat(:, t%nat+1:t%nat+3) = df_da_tilde

    call t%sqnm_opt%sqnm_step(q_and_lat, epot, dq_and_dlat, dd)

    dq_and_dlat = dq_and_dlat + dd

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
  
end module periodic_optimizer