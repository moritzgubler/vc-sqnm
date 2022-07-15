program main
    !! this program reads a periodic structure, puts all atoms back into the periodic box
    !! and writes it to the filename given in the second filename
    use periodic_optimizer
    implicit none
    integer :: nat
    real*8, allocatable, dimension(:,:) :: rxyz
    !! atomic positions always in bohr.
    real*8, dimension(3,3) :: alat
    !! lattice vectors (bohr)
    character(len=250) :: file_in
    !! input/ output file
    character(len=2), allocatable, dimension(:) :: atomnames
    !! chemical name of all the atoms
    character(len=80) :: comment
    !! contents of first comment line
  
    real*8, allocatable, dimension(:,:) :: fxyz
    real(8) :: deralat(3,3), epot, stress(3,3)
    type(optimizer_periodic) optimizer
    integer :: i
    real(8) :: alpha = 2.d0, lattice_weigth = 2.d0, alpha0 = 1.d-2, eps_subsp = 1.d-4
    integer :: nhistx = 10

    call get_command_argument(1,file_in)
    if ( len_trim(file_in) == 0 ) then
      stop "first argument must contain input filename"
    end if
    call get_nat_periodic(file_in, nat)
    allocate(rxyz(3,nat), atomnames(nat), fxyz(3, nat))
    call read_periodic(trim(file_in), nat, rxyz, alat, atomnames, comment)
  
    call optimizer%initialize_optimizer(nat, alat, alpha, nhistx, lattice_weigth, alpha0, eps_subsp)
    
    print*, 'rtest', rxyz(:, 1)
    do i = 1, 5
        call energyandforces_bazant(nat, alat, rxyz, epot, fxyz, deralat, stress)
        print*, 'epot', epot, maxval(abs(fxyz)), maxval(abs(deralat))
        print*, ' '
        call optimizer%optimizer_step(rxyz, alat, epot, fxyz, deralat)
        print*, 'rtest', rxyz(:, 1)
        print*, ' '
    end do


  end program main