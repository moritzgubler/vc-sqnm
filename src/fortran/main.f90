!! This is an example in how the period_optimizer module can be used to optimize
!! both atomic postitions and lattice vectors efficiently.

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

program main
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
    !! forces acting ot the atoms
    real(8) :: deralat(3,3)
    !! negative derivative of the pot, energy w.r. to the lattice vectors
    real(8) :: epot
    !! potential energy
    real(8) :: stress(3,3)
    !! stress tensor
    type(optimizer_periodic) optimizer
    !! optimiyer object
    integer :: i
    !! iteration variable
    real(8) :: alpha = 2.d0
    !! initial step size. default is 1.0. For systems with hard bonds (e.g. C-C) use a value between and 1.0 and
    !! 2.5. If a system only contains weaker bonds a value up to 5.0 may speed up the convergence.
    real(8) :: lattice_weigth = 2.d0
    !! weight or size of the supercell that is used to transform lattice derivatives. Use a value between 1 and 2. 
    !! Default is 2.
    real(8) :: alpha0 = 1.d-2
    !! Lower limit on the step size. 1.e-2 is the default.
    real(8) :: eps_subsp = 1.d-4
    !! Lower limit on linear dependencies of basis vectors in history list. Default 1.e-4.
    !! Increase this parameter if energy or forces contain noise.
    integer :: nhistx = 10
    !! Maximal number of steps that will be stored in the history list. 
    !! Use a value between 3 and 20. Must be <= than 3*nat.

    call get_command_argument(1,file_in)
    if ( len_trim(file_in) == 0 ) then
      print*, "first argument must contain input filename"
      stop
    end if
    call get_nat_periodic(file_in, nat)
    allocate(rxyz(3,nat), atomnames(nat), fxyz(3, nat))
    call read_periodic(trim(file_in), nat, rxyz, alat, atomnames, comment)


    !! initialize the periodic optimizer object.
    call optimizer%initialize_optimizer(nat, alat, alpha, nhistx, lattice_weigth, alpha0, eps_subsp)
    
    print*, "iteration, potential energy, norm of forces, norm of lattice derivatives"

    !! Loop until convergence criteria are reached. A good choice would be the norm of the 
    !! forces and lattice derivatives
    do i = 1, 30
        !! calculate energy, forces and lattice derivatives.
        call energyandforces_bazant(nat, alat, rxyz, epot, fxyz, deralat, stress)

        print"(i3, 2x, g0.7, 2x, g0.2, 2x, g0.2)", i-1, epot, maxval(abs(fxyz)), maxval(abs(deralat))

        !! call this method to get atomic coordinates and lattice vectors that are closer to the 
        !! local mimimun than the previous coordinates. Afterwards evaluate energy and forces again.
        call optimizer%optimizer_step(rxyz, alat, epot, fxyz, deralat)
    end do


  end program main

  !! end of example program. The following subroutines can be used to read and write 
  !! atomic structure files and are taken from: https://github.com/moritzgubler/periodic-io

  subroutine read_ascii(filename, nat, rxyz, alat, atomnames, comment)
    !! reads an ascii file with the specified filename units are assumed to be angstroem.
    !! units are converted to hartree units before they are returned
    implicit none
    character(len=*), intent(in) :: filename
    !! filename of the file which is read
    integer, intent(in) :: nat
    !! number of atoms
    real*8, dimension(3, nat), intent(out) :: rxyz
    !! atom positions (in bohr)
    real*8, dimension(3, 3), intent(out) :: alat
    !! lattice vectors in bohr
    character(len=2), intent(out) :: atomnames(nat)
    !! String containing the chemical symbol of each atom.
    character(len=80), intent(out) :: comment
    !! string that was read from comment line (first line in ascii format)
    character(len=250) :: all_line
    !! sting containing entire line
    integer :: i, io, ios
    real*8 :: alat_temp(2, 3)
    real*8 :: Bohr_Ang = 0.52917721067
    !! Bohr to angstrom conversion factor.
  
    open (newunit=io, file=trim(adjustl(filename)), iostat=ios, status="old")
    if (ios /= 0) then
      print *, "error opening ascii file: "//filename
      stop
    end if
    read (io, "(a80)", iostat=ios) comment
    if (ios /= 0) then
      print *, trim(adjustl(filename)), ios
      print*, "error reading file ascii comment"
      comment = ''
    end if
    i = 1
    alat = 0.d0
    do while ( i < 3 )
      read(io,"(a250)", iostat=ios) all_line
      if (ios /= 0) stop "error reading lattice in moleculario"
      all_line = adjustl(all_line)
      if ( len_trim(all_line) == 0 ) cycle
      if ( index(all_line, "#") == 1 ) cycle
      read(all_line, *, iostat=ios) alat_temp(i,:)
      if (ios /= 0) stop "error parsing lattice in moleculario"
      i = i + 1
    end do
    alat = 0.1
    alat(1, 1) = alat_temp(1, 1)
    alat(1, 2) = alat_temp(1, 2)
    alat(2, 2) = alat_temp(1, 3)
    alat(1, 3) = alat_temp(2, 1)
    alat(2, 3) = alat_temp(2, 2)
    alat(3, 3) = alat_temp(2, 3)
    i = 1
    do while (i <= nat)
      read(io,"(a250)", iostat=ios) all_line
      all_line = adjustl(all_line)
      if ( ios/= 0 ) then
        print*, all_line
        stop "error reading line in read_ascii"
      end if
      if ( len_trim(all_line) == 0 ) cycle
      if ( index(all_line, "#") == 1 ) cycle
      read (all_line, *, iostat=ios) rxyz(1, i), rxyz(2, i), rxyz(3, i), atomnames(i)
      if ( ios/= 0 ) then
        print*, all_line
        stop "error parsing line in read_ascii"
      end if
      i = i + 1
    end do
    close (io)
    alat = alat/Bohr_Ang
    rxyz = rxyz/Bohr_Ang
  end subroutine read_ascii
  
  subroutine get_nat_ascii(filename, nat)
    !! counts the number of atoms in an *.ascii file.
    !! nat doesn't need to be printed on first line.
    implicit none
    character(len=*), intent(in) :: filename
    !! name of the ascii file wich will be read.
    integer, intent(out) :: nat
    !! number of atoms
    ! private variables
    integer :: ios, io
    real*8 :: place(3)
    character(len=250) :: all_line
    character(len=2) :: atname
  
    open (newunit=io, file=trim(adjustl(filename)), iostat=ios, status="old")
    if (ios /= 0) then
      print *, "Error opening ascii file to get number of atoms from. Filename: "//filename
      stop
    end if
    read (io, '(a250)', iostat=ios) all_line
    if (ios /= 0) then
      print *, "ios", ios
      print *, "filename:_", filename
      stop "Error reading file in getnat "!//filename
    end if
    nat = 1
    do while (nat < 3)
      read (io, *, iostat=ios) place(1), place(2), place(3)
      if (ios > 0) then
        cycle
      end if
      if (ios < 0) stop "end of file in get nat ascii"!//filename
      nat = nat + 1
    end do
    nat = 0
    do
      read (io, *, iostat=ios) place(1), place(2), place(3), atname
      if (ios > 0) cycle
      if (ios < 0) exit
      nat = nat + 1
    end do
    close (io)
    !nat = 44
  end subroutine get_nat_ascii
  subroutine read_cif(filename,nat , rxyz, alat, atomnames, comment)
    !! This reader will not be able to read everything. It does not understand the full
    !! cif syntax. But it should understand most of them.
    implicit none
    integer :: nat
    !! number of atoms
    character(len=*) :: filename
    !! filename of the file which is created
    real*8, intent(out), dimension(3, nat) :: rxyz
    !! atom positions (in bohr)
    real*8, dimension(3,3), intent(out) :: alat
    !! lattice vectors in bohr
    character(len=2), dimension(nat), intent(out) :: atomnames
    !! String containing the chemical symbol of each atom.
    character(len=80), intent(out) :: comment
    !! string that was read from first comment line
    integer :: io, ios, lc, k
    character(len=250) :: all_line, num_string
    real*8 :: length_a, length_b, length_c
    real*8 :: alpha, beta, gamma
    logical :: aset = .FALSE., bset = .FALSE., cset = .FALSE.
    logical :: aaset = .FALSE., bbset = .FALSE., ccset = .FALSE.
    integer :: atom_site_count, atom_label_pos, x_pos, y_pos, z_pos
    integer :: nat_read
    integer, parameter :: len_words = 30
    character(len=len_words), dimension(:), allocatable :: words
    integer :: nwords, num_words
    real*8 :: xyzred(3, nat)
    real*8 :: Bohr_Ang = 0.52917721067
    integer :: comment_count
    nat_read = 0
    atom_site_count = 0
    x_pos = -1
    y_pos = -1
    z_pos = -1
    comment_count = 0
  
    open(newunit=io, file=filename, iostat=ios, status="old")
    if ( ios /= 0 ) stop "Error opening file in read_cif"
    lc = 0
    ios = 0
    do
      read(io, fmt="(a250)", iostat=ios) all_line
      if ( ios < 0 ) then
        exit
      end if
      lc = lc + 1
      if ( ios > 0 ) then
        print*, "line", lc, "was formatted wrong and is ignored."
        cycle
      end if
      all_line = adjustl(all_line)
      !check if comment line or ompty
      if ( (index(all_line, "#") == 1) .or. len_trim(all_line) == 0 ) then
        if ( (index(all_line, "#") == 1) ) then ! comment line. read it.
          comment_count = comment_count + 1
          if ( comment_count == 1 ) then
            comment = all_line(2:81)
            comment = adjustl(comment)
          end if
        end if
        cycle
      end if
      k = index(all_line, "_cell_length_a")
      if ( k == 1) then !! line contains cell length a
        num_string = all_line((k+15):len_trim(all_line))
        call remove_uncert_from_numstring(num_string)
        read(num_string, *, iostat=ios) length_a
        if ( ios /= 0 ) then
          print*, 'error reading cell length a'
          print*, trim(all_line)
          stop 'error reading _cell_length_a in read_cif periodicIO'
        end if
        aset = .TRUE.
        cycle
      end if
      k = index(all_line, "_cell_length_b")
      if ( k == 1) then !! line contains cell length b
        num_string = all_line((k+15):len_trim(all_line))
        call remove_uncert_from_numstring(num_string)
        read(num_string, *, iostat=ios) length_b
        if(ios /=0) then
          print*, 'error reading cell length b'
          print*, trim(all_line)
          stop 'error reading _cell_length_b in read_cif periodicIO'
        end if
        bset = .TRUE.
        cycle
      end if
      k = index(all_line, "_cell_length_c")
      if ( k == 1) then !! line contains cell length c
        num_string = all_line((k+15):len_trim(all_line))
        call remove_uncert_from_numstring(num_string)
        read(num_string, *, iostat=ios) length_c
        if(ios /=0) then
          print*, 'error reading cell length c'
          print*, trim(all_line)
          stop 'error reading _cell_length_c in read_cif periodicIO'
        end if
        cset = .TRUE.
        cycle
      end if
      k = index(all_line, "_cell_angle_alpha")
      if ( k == 1) then !! line contains angle alpha
        num_string = all_line((k+18):len_trim(all_line))
        call remove_uncert_from_numstring(num_string)
        read(num_string, *, iostat=ios) alpha
        if ( ios /= 0 ) then
          stop 'error reading cell_angle_alpha in readcif'
        end if
        aaset = .TRUE.
        cycle
      end if
      k = index(all_line, "_cell_angle_beta")
      if ( k == 1) then !! line contains angle beta
        num_string = all_line((k+17):len_trim(all_line))
        call remove_uncert_from_numstring(num_string)
        read(num_string, *, iostat=ios) beta
        if ( ios /= 0 ) then
          stop 'error reading cell_angle_beta in readcif'
        end if
        bbset = .TRUE.
        cycle
      end if
      k = index(all_line, "_cell_angle_gamma")
      if ( k == 1) then !! line contains angle gamma
        num_string = all_line((k+18):len_trim(all_line))
        call remove_uncert_from_numstring(num_string)
        read(num_string, *, iostat=ios) gamma
        if ( ios /= 0 ) then
          stop 'error reading cell_angle_gamma in readcif'
        end if
        ccset = .TRUE.
        cycle
      end if
  
      if ( index(all_line, "_atom_site_") == 1 ) then
        atom_site_count = atom_site_count + 1
        if(index(all_line, "_atom_site_type_symbol") == 1) atom_label_pos = atom_site_count
        if(index(all_line, "_atom_site_fract_x") == 1) x_pos = atom_site_count
        if(index(all_line, "_atom_site_fract_y") == 1) y_pos = atom_site_count
        if(index(all_line, "_atom_site_fract_z") == 1) z_pos = atom_site_count
        cycle
      end if
  
      if ( (x_pos > 0) .and. (y_pos > 0) .and. (z_pos > 0) ) then ! atom line expected after _atom_* lines have been given.
        nat_read = nat_read + 1
        if ( nat_read == 1 ) then
          nwords = atom_site_count
          allocate(words(nwords))
        end if
        if ( nwords /= atom_site_count ) then !! unexpected _atom_site_ line
          stop "nwords /= atom_site_count"
        end if
        call get_num_words(all_line, num_words, 250)
        if ( num_words /= nwords ) then
          print*, trim(all_line)
          stop "numwords /= nwords"
        end if
        call get_words(all_line, nwords, 250, words, len_words)
        if ( len_trim(words(atom_label_pos)) > 2 ) then
          print*, trim(words(atom_label_pos))
          stop "expected atom symbol but got line above"
        end if
        atomnames(nat_read) = trim(words(atom_label_pos))
        read(words(x_pos), *, iostat=ios) xyzred(1, nat_read)
        if ( ios/= 0 ) stop "error reading atom position"
        read(words(y_pos), *, iostat=ios) xyzred(2, nat_read)
        if ( ios/= 0 ) stop "error reading atom position"
        read(words(z_pos), *, iostat=ios) xyzred(3, nat_read)
        if ( ios/= 0 ) stop "error reading atom position"
        if ( nat_read == nat ) then
          exit
        end if
      end if
    end do
    if ( nat /= nat_read ) then
      print*, nat, nat_read
      stop "expected and actual number of atoms do not match in read_cif"
    end if
    length_a = length_a / Bohr_Ang
    length_b = length_b / Bohr_Ang
    length_c = length_c / Bohr_Ang
    call calcAlat
    call frac2cart(nat, alat, xyzred, rxyz)
    if ( comment_count == 0 ) then
      comment = ""
    end if
  contains
    subroutine calcAlat
      implicit none
      real*8 :: conv = 0.01745329251994329576 ! converts degree to radians
      alat = 0
      alat(1,1) = length_a
      alat(1,2) = length_b*cos(gamma*conv)
      alat(2,2) = sqrt(length_b*length_b - alat(1,2)*alat(1,2))
      alat(1,3) = length_c*cos(beta*conv)
      alat(2,3) = (length_b * length_c * cos(alpha*conv) &
                    - alat(1,2) * alat(1,3)) / alat(2,2)
      if(length_c**2 - alat(1,3)**2 - alat(2,3)**2 < 0) then
        stop "square root of negativa number has to be calcuted in read cif."
      end if
      alat(3,3) = sqrt(length_c**2 - alat(1,3)**2 - alat(2,3)**2)
    end subroutine calcAlat
  
    subroutine remove_uncert_from_numstring(numstr)
      !! in some .cif file numbers are written with brackets in the end: 3.141(5)
      !! This subroutine removes these brackets and the number in between them.
      implicit none
      character(len=*), intent(inout) :: numstr
      !! string that may contain brackets that will be removed if present.
      integer :: len_word, numwords, brack_index
  
      len_word = len(numstr)
      call get_num_words(numstr, numwords, len_word)
      if ( numwords /= 1 ) then
        print*, trim(numstr)
        print*, numwords
        stop 'multiple words supplied in remove_encert_from_numstr'
      end if
  
      brack_index = index(numstr, '(')
      if (brack_index <= 0) return
      numstr(brack_index:len_word) = ''
    end subroutine remove_uncert_from_numstring
  
  end subroutine read_cif
  
  subroutine get_nat_cif(filename, nat)
    !! counts the number of atoms in an *.cif file.
    implicit none
    integer, intent(out) :: nat
    !! number of atoms
    character(len=*) :: filename
    !! name of the cif file which will be read
    integer :: io, ios
    character(len=250) :: all_line
    logical :: look_for_atoms, isEle
    integer :: atom_site_pos, nwords
    integer, parameter :: len_words = 30
    character(len=len_words), dimension(:), allocatable :: words
    atom_site_pos = 1
    look_for_atoms = .FALSE.
    open(newunit=io, file=filename, iostat=ios, status="old")
    if ( ios /= 0 ) then
      print*, filename
      stop "error opening cif file to read number of atoms."
    end if
    nat = 0
    do
      read(io, "(a)", iostat=ios) all_line
      if(ios < 0) exit
      if ( ios > 0 ) then
        print*, "reading line in cif file", filename
        stop "error reading line"
      end if
      all_line = adjustl(all_line)
      !call downcase(all_line)
      if (iscommentorempty(all_line)) cycle
      if(index(all_line, "_atom_site_type_symbol") == 1) then
        look_for_atoms = .TRUE.
        allocate(words(atom_site_pos))
        cycle
      end if
  
      if ( .not. look_for_atoms ) then !! look for _atom_site labels
        if ( index(all_line, "_atom_site_") == 1) then
          atom_site_pos = atom_site_pos + 1
          cycle
        end if
      else ! look for atoms.
        call get_num_words(all_line, nwords, 250)
        if ( nwords >= atom_site_pos ) then !! posible atom
          call get_words(all_line, atom_site_pos, 250, words, len_words)
          call isElement(words(atom_site_pos)(1:2), isEle)
          if(isEle) then !! atom found
            nat = nat + 1
          end if
        else !! no atom
          cycle
        end if
      end if
    end do
  
    close(io)
  contains
    function iscommentorempty(line_in) result(bool)
      implicit none
      character(len=*), intent(in) :: line_in
      logical :: bool
      bool = .FALSE.
      if ( len_trim(line_in) == 0 ) then !line is empty
        bool = .TRUE.
      end if
      if ( adjustl(line_in(1:1)) == "#" ) then
        bool = .TRUE. ! line is comment
      end if
      return
    end function iscommentorempty
  
  end subroutine get_nat_cif
  
  subroutine read_dftb(filename, nat, rxyz, alat, atomnames, comment)
    !! reads the .gen format of the dftb+ package. Cell type S (supercell)
    !! and cell type F (supercel with fractional coordintes) are implemented.
    !! Cell type C (cluster) causes the subroutine to stop.
    !! Documentation of file format in Appendix D:
    !! https://github.com/dftbplus/dftbplus/releases/download/21.1/manual.pdf
    use iso_fortran_env
    implicit none
    character(len=*), intent(in) :: filename
    !! filename of the file which is created
    integer, intent(in) :: nat
    !! number of atoms
    real(real64), dimension(3,nat), intent(out) :: rxyz
    !! atom positions (in bohr)
    real(real64), dimension(3,3), intent(out) :: alat
    !! lattice vectors in bohr
    character(len=2), dimension(nat), intent(out) :: atomnames
    !! String containing the chemical symbol of each atom.
    character(len=80), intent(out) :: comment
    !! string that was read from comment line
    integer, parameter :: max_line_length = 250
    character(len=max_line_length) :: all_line
    integer :: io, ios, i, iat
    real(real64) :: Bohr_Ang = 0.52917721067
    character(len=1) :: celltype
    logical :: isfirstcomment
    !! number of different atoms
    integer :: ntypes
    !! Array containing all the different elements. Each element is only listed once.
    !! Only the first ntype entries contain data, the rest is garbage.
    character(len=2), dimension(nat) :: atomtypes
    !! Encodes each element with a number. The element name of iat is accessed the following way:
    !! atomtypes(atomnumbernames(iat))
    integer :: atomnumbernames(nat)
    integer :: natread
    real(real64) :: shift(3)
    real(real64), dimension(:,:), allocatable :: xyzred
  
    isfirstcomment = .TRUE.
    comment = ""
  
    open(newunit=io, file=filename, iostat=ios, status="old")
    if ( ios /= 0 ) then
      print*, "error opening file: ",trim(filename)
      stop "ioerror"
    end if
    ! read first line (nat and cell type)
    do
      read(io, "(a250)", iostat=ios) all_line
      if (ios /= 0) stop "error reading line in read_dftb moleculario"
      if (index(adjustl(all_line), "#") == 1) then !! comment line
        if ( isfirstcomment ) then
          isfirstcomment = .FALSE.
          comment = all_line(3:82)
        end if
        cycle
      end if
      read(all_line, *, iostat=ios) natread, celltype
      if (ios /= 0) stop "error parsing line containg dftb nat celltype moleculario"
      exit
    end do
  
    if ( nat /= natread ) then
      print*, "nat_input, natread :", nat, natread
      stop "nat from argument in read_dftb is not equal to nat read in dftb file. Aborting..."
    end if
  
    if ( celltype /= "S" .and. celltype /= "F" ) then
      print*, "not accepted cell type: ", celltype
      stop "not recognized cell type in read_dftb moleculario"
    end if
  
    ! read atomtypes
    do
      read(io, "(a250)", iostat=ios) all_line
      if (ios /= 0) stop "error reading line in read_dftb moleculario"
      if (index(adjustl(all_line), "#") == 1) cycle !! comment line
      all_line = adjustl(all_line)
      call get_num_words(all_line, ntypes, max_line_length)
      call get_words(all_line, ntypes, max_line_length, atomtypes(1:ntypes), 2)
      exit
    end do
  
    ! read all atoms
    i = 0
    do
      if (nat == i) exit
      read(io, "(a250)", iostat=ios) all_line
      if (ios /= 0) stop "error reading pos in read_dftb moleculario"
      if (index(adjustl(all_line), "#") == 1) cycle !! comment line
      i = i + 1
      read(all_line, *, iostat = ios) iat, atomnumbernames(i), rxyz(1, i), rxyz(2,i), rxyz(3,i)
      if ( ios /= 0 ) then
        stop "error parsing position in read_dftb (moleculario)"
      end if
      if (iat /= i) stop "inconsistent internal state in read_dftb (moleculario)"
      atomnames(i) = atomtypes(atomnumbernames(i))
    end do
  
    ! read shift of origin from lattice cell
    do
      read(io, "(a250)", iostat=ios) all_line
      if (ios /= 0) stop "error reading shift in read_dftb moleculario"
      if (index(adjustl(all_line), "#") == 1) cycle !! comment line
      read(all_line, *, iostat=ios) shift(1), shift(2), shift(3)
      if (ios /=0 ) stop "error parsing shift in read_dftb (moleculario)"
      exit
    end do
  
    ! read lattice cell
    i = 0
    do
      if (i == 3) exit
      read(io, "(a250)", iostat=ios) all_line
      if (ios /= 0) stop "error reading lattice in read_dftb moleculario"
      if (index(adjustl(all_line), "#") == 1) cycle !! comment line
      i = i + 1
      read(all_line, *, iostat=ios) alat(1,i), alat(2,i), alat(3,i)
      if (ios /= 0) stop "error parsing lattice in read_dftb moleculario"
    end do
    if (celltype == "F" ) then
      allocate(xyzred(3,nat))
      xyzred = rxyz
      call frac2cart(nat, alat, xyzred, rxyz)
      deallocate(xyzred)
    end if
  
    !! Positions and lattice vectors are in Angstroem in .gen file format. Convert it
    !! to hartree units.
    rxyz = rxyz / Bohr_Ang
    alat = alat / Bohr_Ang
  
  end subroutine read_dftb
  
  subroutine get_nat_dftb(filename, nat)
    !! counts the number of atoms in an *.gen file.
    implicit none
    integer, intent(out) :: nat
    !! number of atoms
    character(len=*), intent(in) :: filename
      !! name of the gen file wich will be read.
    integer :: io, ios
    character(len=250) :: all_line
  
    open(newunit=io, iostat=ios, file=filename, status="old")
    if (ios /= 0) stop "error opening dftb file to get number of atoms from (moleculario)"
    ! try and read first line and ignore comments.
    do
      read(io, "(a250)", iostat=ios) all_line
      if (ios /=0) stop "error reading line in get_nat_dftb moleculario"
      if (index(adjustl(all_line), "#") == 1) cycle !! comment line
      read(all_line, *, iostat=ios) nat
      if (ios /= 0) stop "error parsing line containg dftb nat moleculario"
      exit
    end do
    close(io)
  end subroutine get_nat_dftb
  subroutine read_in(filename, nat, rxyz, alat, atomnames, comment)
    !! reads a *.in file with the specified filename.
    implicit none
    character(len=*), intent(in) :: filename
    !! filename of the file which is read
    integer, intent(in) :: nat
    !! number of atoms
    real*8, dimension(3,nat), intent(out) :: rxyz
    !! atom positions (in bohr)
    real*8, dimension(3,3), intent(out) :: alat
    !! lattice vectors in bohr
    character(len=2), dimension(nat), intent(out) :: atomnames
    !! String containing the chemical symbol of each atom.
    character(len=80), intent(out) :: comment
    !! string that was read from comment line (first line in ascii format)
    integer, parameter :: max_line_length = 250
    character(len=max_line_length) :: all_line
    character(len=max_line_length) :: dat_line
    integer :: io, ios, i
    real*8 :: Bohr_Ang = 0.52917721067
    integer :: comment_count
  
    comment_count = 0
  
    open(newunit=io, file=filename, iostat=ios, status="old")
    if ( ios /= 0 ) then
      print*, "error opening file: ",trim(filename)
      stop "ioerror"
    end if
    ! read lattice vectors
    i = 1
    lattice_loop: do while ( i <= 3 )
      read(io, "(a250)") all_line
      if ( ios /= 0 ) then
        print*, "error reading lattice vectors in file: ", trim(filename)
        stop "ioerror"
      end if
      all_line = adjustl(all_line)
      if(iscommentorempty(all_line)) then
        if ( index(all_line, "#") == 1 ) then !! comment, read it
          comment_count = comment_count + 1
          if ( comment_count == 1 ) then
            comment = all_line(2:81)
            comment = adjustl(comment)
          end if
        end if
        cycle
      end if
  
      if ( all_line(1:14) == "lattice_vector" ) then ! lattice vector found
        dat_line = adjustl(all_line(15:max_line_length))
        read(dat_line,*, iostat = ios) alat(1,i), alat(2,i), alat(3,i)
        if ( ios /= 0 ) then
          print*, "error in line containing:"
          print*, trim(all_line)
          stop "error reading lattice vectors"
        end if
        i = i + 1
        cycle
      end if
  
      !check if an atom tag is found
      if ( all_line(1:4) == "atom" ) then
        print*, "atom found before all lattice vectors were given"
        stop "format error in .in file"
      end if
      print*, "unable to format line: "//trim(all_line)
      stop "format error"
    end do lattice_loop
    ! read atomic positions
    i = 1
    do while ( i <= nat )
      read(io, "(a250)", iostat=ios) all_line
      if ( ios /= 0 ) then
        print*, "error reading atom in file: ", trim(filename)
        stop "ioerror"
      end if
      all_line = adjustl(all_line)
      if(iscommentorempty(all_line)) cycle
  
      if ( all_line(1:4) == "atom" ) then ! atomic position found
        dat_line = adjustl(all_line(5:max_line_length))
        read(dat_line,*, iostat = ios) rxyz(1,i), rxyz(2,i), rxyz(3,i), atomnames(i)
        if ( ios /= 0 ) then
          print*, "error in line containing:"
          print*, trim(all_line)
          stop "error reading atom"
        end if
        i = i + 1
        cycle
      end if
      if ( all_line(1:9) == "atom_frac" ) then ! reduced atomic position found
        dat_line = adjustl(all_line(5:max_line_length))
        read(dat_line,*, iostat = ios) rxyz(1,i), rxyz(2,i), rxyz(3,i), atomnames(i)
        if ( ios /= 0 ) then
          print*, "error in line containing:"
          print*, trim(all_line)
          stop "error reading atom"
        end if
        rxyz(:,i) = matmul(alat, rxyz(:,i))
        i = i + 1
        cycle
      end if
  
    end do
    close(io)
    rxyz = rxyz / Bohr_Ang
    alat = alat / Bohr_Ang
    if ( comment_count == 0 ) then
      comment = ""
    end if
  contains
    function iscommentorempty(line_in) result(bool)
      implicit none
      character(len=250) :: line_in
      logical :: bool
      bool = .FALSE.
      if ( len_trim(line_in) == 0 ) then !line is empty
        bool = .TRUE.
      end if
      if ( line_in(1:1) == "#" ) then
        bool = .TRUE. ! line is comment
      end if
      return
    end function iscommentorempty
  end subroutine read_in
  
  subroutine get_nat_in(filename, nat)
    !! counts the number of atoms in an *.in file.
    implicit none
    integer, intent(out) :: nat
    !! number of atoms
    character(len=*), intent(in) :: filename
    !! name of the in file wich will be read.
    integer :: io, ios
    character(len=250) :: all_line
    nat = 0
    open(newunit=io, iostat=ios, file=filename, status="old")
    if ( ios /= 0 ) then
      nat = 0
      return
    end if
    ios=0
    do while ( ios == 0 )
      read(io,"(a250)",iostat=ios) all_line
      if ( ios /= 0 ) exit
      all_line = adjustl(all_line)
      if ( all_line(1:4) == "atom" ) then
        nat = nat + 1
      end if
    end do
    if ( ios > 0 ) then !! io error occured
      print*, filename
      stop "io error in getting number of atoms"
    end if
    close(io)
  end subroutine get_nat_in
  
  subroutine read_periodic(filename, nat, rxyz, alat, atomnames, comment)
    !! reads the position, lattice vecors and atomnames from a file. File format is
    !! detected based on the file ending. Supported endings:
    !! .ascii, .cif, .in, .gen, .qe, .vasp and POSCAR (vasp)
    implicit none
    character(len=*), intent(in) :: filename
    !! filename of the file which is created
    integer, intent(in) :: nat
    !! number of atoms
    real*8, dimension(3, nat), intent(out) :: rxyz
    !! atom positions (bohr)
    real*8, dimension(3, 3), intent(out) :: alat
    !! lattice vectors (bohr)
    character(len=2), intent(out) :: atomnames(nat)
    !! String containing the chemical symbol of each atom.
    character(len=80), intent(out) :: comment
    !! content of comment line
    character(len=10) :: filetype
    integer :: l, vaspind
  
    vaspind = index(filename, 'POSCAR')
    l = len_trim(filename)
    !if (vaspind > 0) print*, 'vaspind', l -vaspind
    if ( l - vaspind == 5 .and. vaspind > 0) then !! vasp file
      call read_vasp(filename, nat, rxyz, alat, atomnames, comment)
      return
    end if
  
  
    call get_file_type(filename, filetype)
  
    select case (trim(filetype))
    case ("ascii")
      call read_ascii(filename, nat, rxyz, alat, atomnames, comment)
    case("cif")
      call read_cif(filename, nat, rxyz, alat, atomnames, comment)
    case("in")
      call read_in(filename, nat, rxyz, alat, atomnames, comment)
    case("gen")
      call read_dftb(filename, nat, rxyz, alat, atomnames, comment)
    case("qe")
      call read_quantum_espresso(filename, nat, rxyz, alat, atomnames)
      comment = ''
    case('vasp')
      call read_vasp(filename, nat, rxyz, alat, atomnames, comment)
    case default
      print*, filename
      print*, "filetype, ", filetype
      stop "unknown filetype read_periodic"
    end select
  
  end subroutine read_periodic
  
  subroutine get_nat_periodic(filename, nat)
    !! counts the number of atoms of a vasp file
    implicit none
    character(len=*), intent(in) :: filename
    !! filename of the file which is read
    integer, intent(out) :: nat
    !! number of atoms
    character(len=10) :: filetype
    integer :: l, vaspind
    vaspind = index(filename, 'POSCAR')
    l = len_trim(filename)
    !if (vaspind > 0) print*, 'vaspind', l -vaspind
    if ( l - vaspind == 5 .and. vaspind > 0) then !! vasp file
      call get_nat_vasp(filename, nat)
      return
    end if
  
    call get_file_type(filename, filetype)
  
    select case (trim(filetype))
    case ("ascii")
      call get_nat_ascii(filename, nat)
    case("cif")
      call get_nat_cif(filename, nat)
    case("in")
      call get_nat_in(filename, nat)
    case('qe')
      call get_nat_quantum_espresso(filename, nat)
    case("gen")
      call get_nat_dftb(filename, nat)
    case("vasp")
      call get_nat_vasp(filename, nat)
    case default
      print*, "filetype, ", filetype
      stop "unknown filetype get_nat_periodic"
    end select
  end subroutine get_nat_periodic
  subroutine read_quantum_espresso(filename, nat, rxyz, alat, atomnames)
    !! reads the periodic structure of a quantum espresso file. units are returned in bohr
    !! only ibrav = 0 is currently supported.
    use, intrinsic :: iso_fortran_env
    implicit none
    character(len=*), intent(in) :: filename
    !! filename of the file which is created
    integer, intent(in) :: nat
    !! number of atoms
    real*8, dimension(3, nat), intent(out) :: rxyz
    !! atom positions (bohr)
    real*8, dimension(3, 3), intent(out) :: alat
    !! lattice vectors (bohr)
    character(len=2), intent(out) :: atomnames(nat)
    !! String containing the chemical symbol of each atom.
  
    real(real64) :: Bohr_Ang = 0.52917721067
    integer :: io, ios, i, iend
    character(len=200) :: aline
    integer :: ibrav
    logical :: ibravset
    real(real64) :: convert
    logical :: latread, posread
  
    ibravset = .FALSE.
    latread = .FALSE.
    posread = .FALSE.
  
    open(newunit=io, file=filename, iostat=ios, status="old")
    if ( ios /= 0 ) then
      print*, "error opening file: ",trim(filename)
      stop "ioerror in read_quantum_espresso (molecularIO)"
    end if
  
    floop: do
        read(io, "(a200)", iostat = ios) aline
        if (ios /= 0) exit floop
        if ( index(aline, "ibrav") > 0 .or. index(aline, "IBRAV") > 0) then
          i = index(aline, "=") + 1
          iend = index(aline, ",")
          if (iend <=0 ) iend = 200
          ibravset = .TRUE.
          read(aline(i:iend), *, iostat = ios) ibrav
          if ( ibrav /= 0 ) then
            stop "only ibrav = 0 supported at the moment"
          end if
        end if
  
        if ( index(aline, "CELL_PARAMETERS") > 0 ) then ! cell parameters found
          if(index(aline, "bohr") > 0) then
             convert = 1.0
          else if( index(aline, "angstrom") > 0) then
            convert = 1.0 / Bohr_Ang
          else
            stop "unknown units"
          end if
          do i = 1, 3, 1
            read(io, *, iostat = ios) alat(1, i), alat(2, i), alat(3, i)
            if ( ios /= 0 ) then
              stop "error reading lattice vectors"
            end if
          end do
          alat = alat * convert
          latread = .TRUE.
        end if
  
        ! atoms found
        if ( index(aline, "ATOMIC_POSITIONS") > 0 .or. index(aline, "atomic_positions") > 0 ) then
          ! get units:
          if ( index(aline, "bohr") > 0 ) then
            convert = 1.0
          else if( index(aline, "angstrom") > 0) then
            convert = 1.0 / Bohr_Ang
          else
            print*, "unkown unit on this line:", trim(aline)
            stop "unknown unkwonsn unit"
          end if
          do i = 1, nat, 1
            read(io, *, iostat=ios) atomnames(i), rxyz(1,i), rxyz(2,i), rxyz(3,i)
            if ( ios /= 0 ) then
              stop "error reading positions in qe read moleculario"
            end if
          end do
          rxyz = rxyz * convert
          posread = .TRUE.
        end if
    end do floop
    if ( .not. ibravset ) then
      print*, 'ibrav not set in quantum espresso file. Proceed with care.'
    end if
    if ( .not. latread ) stop "lattice vectors not set"
    if ( .not. posread ) stop "positions not red"
    close(io)
  end subroutine read_quantum_espresso
  
  subroutine get_nat_quantum_espresso(filename, nat)
    !! counts the number of atoms in a quantum espresso file.
    use, intrinsic :: iso_fortran_env
    implicit none
    integer, intent(out) :: nat
    !! number of atoms
    character(len=*) :: filename
    !! filename that will be read.
  
    integer :: io, ios, i, iend
    character(len=200) :: aline
  
    open(newunit=io, file=filename, iostat=ios, status="old")
    if ( ios /= 0 ) then
      print*, "error opening file: ",trim(filename)
      stop "ioerror in read_quantum_espresso (molecularIO)"
    end if
    nat = -1
  
    floop: do
      read(io, "(a200)", iostat = ios) aline
      if (ios /= 0) exit floop
      if ( index(aline, "nat") > 0 .or. index(aline, "NAT") > 0) then
        i = index(aline, "=") + 1
        if (i < 1) stop "no equal sign in nat line qe in file (molecularIO)"
        iend = index(aline, ",")
        if (iend <= 0) iend = 200
          read(aline(i:iend), *, iostat = ios) nat
          if (ios /= 0) stop "error getting nat from quantum espresso file"
        end if
    end do floop
    close(io)
  end subroutine get_nat_quantum_espresso
  subroutine read_vasp(filename, nat, rxyz, alat, atomnames, comment)
    implicit none
    !! reads a vasp file with the specified filename
    !! units are converted to hartree units before they are returned
    character(len=*), intent(in) :: filename
    !! filename of the file which is created
    integer, intent(in) :: nat
    !! number of atoms
    real(8), dimension(3, nat) :: rxyz
    !! atom positions
    real(8), dimension(3, 3) :: alat
    !! lattice vectors
    character(len=2) :: atomnames(nat)
    !! String containing the chemical symbol of each atom.
    character(len=80), intent(out) :: comment
    !! content of the comment line
    character(len=300) :: all_line
    integer :: io, ios, i
    real(8) :: Bohr_Ang = 0.52917721067
    real(8) :: conver
    integer :: ntypes
    character(len=2), allocatable, dimension(:) :: typenames
    integer, allocatable, dimension(:) :: typecount, typesum
    logical :: is_cartesian
    real(8), allocatable, dimension(:, :) :: xyzred
  
    ! open file
    open (newunit=io, file=trim(adjustl(filename)), iostat=ios, status='old')
    if (ios /= 0) then
      print*, 'io', io
      print*, 'ios', ios
      print *, "Error opening file "//trim(filename)
      stop "Error opening file"
    end if
  
    read(io, '(a80)', iostat=ios) comment
    if (ios /= 0 ) stop 'error reading comment in readvasp moleculario'
  
    read(io, *, iostat = ios) conver
    if (ios /= 0 ) stop ' error reading scaling factor readvasp moleculario'
  
    do i = 1, 3, 1
      read(io, *, iostat=ios) alat(1,i), alat(2,i), alat(3,i)
      if(ios /=0 ) stop 'error reading lattice readvasp moleculario'
    end do
  
    read(io, '(a300)', iostat =ios) all_line
    if (ios /=  0 ) stop ' error reading all_line (for atomtypes) in readvasp moleculario'
  
    call get_num_words(all_line, ntypes, 300)
    allocate(typenames(ntypes), typecount(ntypes), typesum(ntypes))
    call get_words(all_line, ntypes, 300, typenames, 2)
  
    read(io, *, iostat=ios) typecount
    if(ios /= 0 ) stop 'error reading typcount readvasp moleculario'
  
    !! check coordinate type (line could contain selective dynamics check that)
    read(io, '(a300)', iostat=ios) all_line
    if (ios /=0 ) stop 'error reading coordinate type readvasp moleculario'
    if ( index(all_line, 'Selective dynamics') > 0 &
        .or. index(all_line, 'selective dynamics') > 0) then
      read(io, '(a300)', iostat=ios) all_line
      if (ios /=0 ) stop 'error reading coordinate type readvasp moleculario'
    end if
  
    is_cartesian = calc_is_cartesian()
  
    do i = 1, nat, 1
      read(io, *, iostat = ios) rxyz(1,i), rxyz(2,i), rxyz(3,i)
      if (ios /= 0) stop 'error reading coordinate line in readvasp moleculario'
    end do
  
    alat = conver*alat / Bohr_Ang
  
    if ( .not. is_cartesian) then
      allocate(xyzred(3,nat))
      xyzred = rxyz
      call frac2cart(nat, alat, xyzred, rxyz)
      deallocate(xyzred)
    else
      rxyz = conver * rxyz / Bohr_Ang
    end if
  
    typesum(1) = 1
    do i = 2 ,ntypes
      typesum(i) = typecount(i-1) + typesum(i - 1)
    end do
  
    do i = 1, ntypes -1, 1
      atomnames(typesum(i): typesum(i+1) -1) = typenames(i)
    end do
    atomnames(typesum(ntypes):nat) = typenames(ntypes)
  
  
  contains
    logical function calc_is_cartesian()
      all_line = adjustl(all_line)
      if (all_line(1:1) == 'c' .or. all_line(1:1) =='C' &
        .or. all_line(1:1) == 'k' .or. all_line(1:1) =='K') then
        calc_is_cartesian = .TRUE.
        return
      end if
      if ( all_line(1:1) =='d' .or. all_line(1:1) == 'D' ) then
        calc_is_cartesian = .FALSE.
        return
      end if
      stop 'error reading coordinate type line in readvasp moleculario'
    end function calc_is_cartesian
  
  end subroutine read_vasp
  
  
  subroutine get_nat_vasp(filename, nat)
    !! counts the number of atoms in a vasp1 file.
    implicit none
    character(len=*), intent(in) :: filename
    !! name of the ascii file wich will be read.
    integer, intent(out) :: nat
    ! private variables
    integer :: ios, io, i, ntypes
    real*8 :: place(3)
    character(len=250) :: all_line
    integer, allocatable, dimension(:) :: typecount
  
    ! open file
    open (newunit=io, file=trim(adjustl(filename)), iostat=ios, status='old')
    if (ios /= 0) then
      print *, "Error opening file"//trim(filename)
      stop "Error opening file"
    end if
  
    read(io, '(a250)', iostat=ios) all_line
    if (ios /= 0 ) stop 'error reading comment in getnatvasp moleculario'
  
    read(io, *, iostat = ios) place(1)
    if (ios /= 0 ) stop ' error reading scaling factor getnatvasp moleculario'
  
    do i = 1, 3, 1
      read(io, *, iostat=ios) place(1), place(2), place(3)
      if(ios /=0 ) stop 'error reading lattice getnatvasp moleculario'
    end do
  
    read(io, '(a250)', iostat=ios) all_line
    if(ios /=0) stop 'error reading typpenames in getnatvasp moleculario'
  
    call get_num_words(all_line, ntypes, 250)
    allocate(typecount(ntypes))
  
    read(io, *, iostat=ios) typecount
    if(ios /=0) stop 'error reading typecount getnatvasp moleculario'
  
    nat = sum(typecount)
    close(io)
  end subroutine get_nat_vasp
  subroutine read_xyz(filename, nat, rxyz, atomnames, comment)
    !! reads xyz file. assumes units of file are in angstrom and converts them to a.u.
    implicit none
    character(len=*), intent(in) :: filename
    !! name of the file that will be read
    integer, intent(inout) :: nat
    !! number of atoms
    real*8, intent(out), dimension(3,nat) :: rxyz
    !! positions of the atoms (bohr)
    real*8 :: Bohr_Ang = 0.52917721067
    !! atomic positions in a.u
    character(len=2), intent(out), dimension(nat) :: atomnames
    !! chemical name of the atoms
    character(len=80), intent(out) :: comment
    !! contents of the first comment line
    integer :: io, ios, i
    open(newunit=io, file=filename, iostat=ios, status="old")
    if(ios/=0) stop "error opening input file"
    read(io,*) nat
    read(io,"(a80)") comment
    do i = 1, nat, 1
      read(io,*, iostat = ios) atomnames(i), rxyz(1,i), rxyz(2,i), rxyz(3,i)
      if ( ios/=0 ) then
        print*, "error reading line: ", i
        stop
      end if
    end do
    close(io)
    rxyz = rxyz / Bohr_Ang
  end subroutine read_xyz
  
  subroutine get_nat_xyz(filename, nat)
    !! counts the number of atoms in a .xyz file
    implicit none
    character(len=*), intent(in) :: filename
    !! filename that will be parsed
    integer, intent(out) :: nat
    !! number of atoms
    integer :: io, ios
    open(newunit=io,file=filename,iostat=ios, status="old")
    if(ios/=0) stop "error opening input file"
    read(io,*) nat
    close(io)
  end subroutine get_nat_xyz
  
  subroutine write_ascii(filename, nat, rxyz, alat, atomnames, comment)
    !! writes the position and lattice vectors to a file in ascii format.
    implicit none
    character(len=*), intent(in) :: filename
    !! filename of the file which is created
    integer, intent(in) :: nat
    !! number of atoms
    real*8, dimension(3, nat), intent(in) :: rxyz
    !! atom positions in bohr
    real*8, dimension(3, 3), intent(in) :: alat
    !! lattice vectors
    character(len=2), intent(in) :: atomnames(nat)
    !! String containing the chemical symbol of each atom.
    character(len=80), intent(in) :: comment
    !! content that will be written to the comment line
  
    ! private variables:
    integer :: io, ios, i
    real*8, dimension(3, nat) :: rxyz_copy
    real*8, dimension(3, 3) :: alat_copy
    real*8 :: Bohr_Ang = 0.52917721067
  
    alat_copy = alat
    rxyz_copy = rxyz
    call alat2ascii(nat, rxyz_copy, alat_copy)
  
    open (newunit=io, file=trim(adjustl(filename)), iostat=ios)
    if (ios /= 0) then
      print *, "Error opening file"//trim(filename)
      stop "Error opening file"
    end if
  
    !Transofrm output to angstroem
    alat_copy = alat_copy*Bohr_Ang
    rxyz_copy = rxyz_copy*Bohr_Ang
  
    write (io, *) trim(comment)
    write (unit=io, fmt=*, iostat=ios) alat_copy(1, 1), alat_copy(1, 2), alat_copy(2, 2)
    if(ios /= 0) stop 'error writing lattice vectors'
    write (unit=io, fmt=*, iostat=ios) alat_copy(1, 3), alat_copy(2, 3), alat_copy(3, 3)
    if(ios /= 0) stop 'error writing lattice vectors'
    do i = 1, nat, 1
      write (unit=io, fmt=*, iostat=ios) rxyz_copy(:, i), atomnames(i)
      if(ios /= 0) stop 'error writing ascii positions'
    end do
    close (io)
  
  contains
    subroutine alat2ascii(nat1, rxyz1, alat_io)
      !! Converts the lattice vectors in ascii format.
     implicit none
     integer, intent(in) :: nat1
     !! Number of atoms
     real*8, dimension(3, nat1), intent(inout) :: rxyz1
     !! Position of atoms.
     real*8, intent(inout), dimension(3, 3) :: alat_io
     !! Lattice vectors that should to be transformed.
  
     ! private variables:
     real*8, dimension(3, 3) :: alat1, t, alat_in_inv
     real*8 :: r1, r2, r3
     integer :: it
     alat1 = 0
     r1 = norm2(alat_io(:, 1))
     r2 = norm2(alat_io(:, 2))
     r3 = norm2(alat_io(:, 3))
     alat1(1, 1) = r1
     alat1(1, 2) = dot_product(alat_io(:, 1), alat_io(:, 2))/r1
     alat1(1, 3) = dot_product(alat_io(:, 1), alat_io(:, 3))/r1
     alat1(2, 2) = sqrt(r2**2 - alat1(1, 2)**2)
     alat1(2, 3) = (dot_product(alat_io(:, 2), alat_io(:, 3)) - alat1(1, 2)*alat1(1, 3)) &
                  /alat1(2, 2)
     alat1(3, 3) = sqrt(r3**2 - alat1(1, 3)**2 - alat1(2, 3)**2)
     call matinv3(alat_io, alat_in_inv)
     t = matmul(alat1, alat_in_inv)
     do it = 1, nat1, 1
       rxyz1(:, it) = matmul(t, rxyz1(:, it))
     end do
     alat_io = alat1
     !call back2cell(nat1, rxyz1, alat1)
    end subroutine alat2ascii
  
    subroutine matinv3(A, B)
      !! Performs a direct calculation of the inverse of a 3×3 matrix.
      implicit none
      real*8, intent(in), dimension(3, 3) :: A
      !! Input Matrix
      real*8, intent(out), dimension(3, 3) :: B
      !! Inverse matrix
      real*8 :: detinv
  
      ! Calculate the inverse determinant of the matrix
      detinv = 1/(A(1, 1)*A(2, 2)*A(3, 3) - A(1, 1)*A(2, 3)*A(3, 2) &
                  - A(1, 2)*A(2, 1)*A(3, 3) + A(1, 2)*A(2, 3)*A(3, 1) &
                  + A(1, 3)*A(2, 1)*A(3, 2) - A(1, 3)*A(2, 2)*A(3, 1))
  
      ! Calculate the inverse of the matrix
      B(1, 1) = +detinv*(A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2))
      B(2, 1) = -detinv*(A(2, 1)*A(3, 3) - A(2, 3)*A(3, 1))
      B(3, 1) = +detinv*(A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))
      B(1, 2) = -detinv*(A(1, 2)*A(3, 3) - A(1, 3)*A(3, 2))
      B(2, 2) = +detinv*(A(1, 1)*A(3, 3) - A(1, 3)*A(3, 1))
      B(3, 2) = -detinv*(A(1, 1)*A(3, 2) - A(1, 2)*A(3, 1))
      B(1, 3) = +detinv*(A(1, 2)*A(2, 3) - A(1, 3)*A(2, 2))
      B(2, 3) = -detinv*(A(1, 1)*A(2, 3) - A(1, 3)*A(2, 1))
      B(3, 3) = +detinv*(A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))
    end subroutine matinv3
  end subroutine write_ascii
  subroutine write_cif(filename, nat, rxyz, alat, atomnames, comment)
    !! writes a peridic structure to a file with the cif standard.
    implicit none
    integer :: nat
    !! number of atoms
    character(len=*) :: filename
    !! name of the file
    real*8, intent(in), dimension(3, nat) :: rxyz
    !! positions of the atoms (units in bohr)
    real*8, dimension(3,3), intent(in) :: alat
    !! lattice vectors (units in bohr)
    character(len=2), dimension(nat), intent(in) :: atomnames
    !! chemical names of the atoms
    character(len=80), intent(in) :: comment
    !! content that will be written to the comment line
    integer :: io, ios, i
    real*8 :: length_a, length_b, length_c
    real*8 :: alpha, beta, gamma
    real*8 :: xyzred(3,nat)
    real*8 :: Bohr_Ang = 0.52917721067
    real*8 :: alat_convert(3,3)
  
    alat_convert = Bohr_Ang * alat
    call cif_angles()
    call cart2frac(nat, alat, rxyz, xyzred)
  
    open(newunit=io, file=filename, iostat=ios)
    if ( ios /= 0 ) then
      print*, "cif filename: ", trim(filename)
      stop "error opening file in write_cif"
    end if
    if(len_trim(comment) > 0) then
      write(io, *) "# "//trim(comment)
    end if
    write(io,*) "data_cif"
    write(io, "(a)") "_space_group_name_H-M_alt 'P 1'"
    write(io,"(a, 1x, f20.15)") "_cell_length_a", length_a
    write(io,"(a, 1x, f20.15)") "_cell_length_b", length_b
    write(io,"(a, 1x, f20.15)") "_cell_length_c", length_c
    write(io,"(a, 1x, f20.15)") "_cell_angle_alpha", alpha
    write(io,"(a, 1x, f20.15)") "_cell_angle_beta", beta
    write(io,"(a, 1x, f20.15)") "_cell_angle_gamma", gamma
  
    write(io, "(a)") "loop_"
    write(io, "(a)") " _atom_site_type_symbol"
    write(io, "(a)") " _atom_site_label"
    write(io, "(a)") " _atom_site_fract_x"
    write(io, "(a)") " _atom_site_fract_y"
    write(io, "(a)") " _atom_site_fract_z"
  
    do i = 1, nat, 1
      write(io, "(2x, a2, 1x, a, i3.3, 1x, 3f18.14)") atomnames(i)&
            , trim(atomnames(i)), i, xyzred(1,i), xyzred(2,i), xyzred(3,i)
    end do
    close(io)
  contains
    subroutine cif_angles
      implicit none
      real*8, parameter :: converts = 57.295779513
      !! converts radian to degree
      length_a = norm2(alat_convert(:,1))
      length_b = norm2(alat_convert(:,2))
      length_c = norm2(alat_convert(:,3))
      alpha = acos(dot_product(alat_convert(:,2), alat_convert(:,3)) / (length_b * length_c))
      beta = acos(dot_product(alat_convert(:,1), alat_convert(:,3)) / (length_a * length_c))
      gamma = acos(dot_product(alat_convert(:,1), alat_convert(:,2)) / (length_a * length_b))
      alpha = converts * alpha
      beta = converts * beta
      gamma = converts * gamma
    end subroutine cif_angles
  end subroutine write_cif
  subroutine write_dftb(filename, nat, rxyz, alat, atomnames, comment)
    !! writes periodic structure to a file in the dftb+ (.gen) format
    use iso_fortran_env
    implicit none
    !! filename of the file which is created
    character(len=*), intent(in) :: filename
    !! number of atoms
    integer, intent(in) :: nat
    !! atom positions in bohr
    real(real64), dimension(3, nat), intent(in) :: rxyz
    !! lattice vectors
    real(real64), dimension(3, 3), intent(in) :: alat
    !! String containing the chemical symbol of each atom.
    character(len=2), intent(in) :: atomnames(nat)
    !! comment string that will be written to file
    character(len=80), intent(in) :: comment
    !! number of different atoms
    integer :: ntypes
    !! Array containing all the different elements. Each element is only listed once.
    !! Only the first ntype entries contain data, the rest is garbage.
    character(len=2), dimension(nat) :: atomtypes
    !! Encodes each element with a number. The element name of iat is accessed the following way:
    !! atomtypes(atomnumbernames(iat))
    integer :: atomnumbernames(nat)
  
    integer :: u, ityp, iat, ios, i, j
  
    real(real64), parameter :: Bohr_Ang = 0.52917721067
    real(real64) :: pos_ang(3,nat), lat_ang(3,3)
  
    pos_ang = rxyz * Bohr_Ang
    lat_ang = alat * Bohr_Ang
  
    call back2cell(nat, pos_ang, lat_ang)
  
    ! calculate ntypes and atomtypes
    ntypes = 1
    atomtypes(1) = atomnames(1)
    do i = 2, nat, 1
      do j = 1, ntypes, 1
        if ( atomnames(i) == atomtypes(j) ) then
          goto 123
        end if
      end do
      ntypes = ntypes + 1
      atomtypes(ntypes) = atomnames(i)
      123 continue
    end do
  
    ! calculate atomnumbernames
    do i = 1, nat, 1
      do j = 1, ntypes, 1
        if ( atomnames(i) == atomtypes(j) ) then
          atomnumbernames(i) = j
        end if
      end do
    end do
  
    open(newunit=u, file=filename, iostat=ios)
    if ( ios /= 0 ) stop "error opening dftb output file in moleculario"
  
    if ( len_trim(comment) > 0 ) then
        write(u, "(a1, 1x, a)", iostat=ios) "#", trim(comment)
        if (ios /= 0) stop "error writing comment in write dftb in moleculario"
    end if
  
    ! first non comment line contains number of atoms and format type. Here 'S' is
    ! used for a periodic supercell
    write(u,'(i5.5,a)', iostat=ios) nat, " S"
    if (ios /= 0) stop "error writing to dftb  output file in moleculario"
  
    write(u,*, iostat=ios) (trim(adjustl(atomtypes(ityp)))//" ", ityp=1,ntypes)
    if (ios /= 0) stop "error writing atomtypes to dftb  output file in moleculario"
  
    do iat = 1, nat
        write(u,'(i5,1x,i5,3(1x,es25.15))', iostat=ios) iat, atomnumbernames(iat), pos_ang(:,iat)
        if (ios /= 0) stop "error writing positions to dftb  output file in moleculario"
    end do
  
    ! write origin of coordinate system
    write(u,'(3(1x,es25.15))', iostat=ios) 0.d0,0.d0,0.d0
    if (ios /= 0) stop "error writing zeros to dftb output file in moleculario"
  
    write(u,'(3(1x,es25.15))', iostat=ios) lat_ang(:, 1)
    if (ios /= 0) stop "error writing lattice vector 1 to dftb output file in moleculario"
  
    write(u,'(3(1x,es25.15))', iostat=ios) lat_ang(:, 2)
    if (ios /= 0) stop "error writing lattice vector 2 to dftb output file in moleculario"
  
    write(u,'(3(1x,es25.15))', iostat=ios) lat_ang(:, 3)
    if (ios /= 0) stop "error writing lattice vector 3 to dftb output file in moleculario"
  
    close(u)
  end subroutine write_dftb
  subroutine write_in(filename, nat, rxyz, alat, atomnames, comment)
    !! writes a periodic structure to a . in file
    implicit none
    character(len=*) :: filename
    !! output file
    integer, intent(in) :: nat
    !! number of atoms
    real*8, intent(in), dimension(3,nat) :: rxyz
    !! atomic positions in bohr
    character(len=2), intent(in), dimension(nat) :: atomnames
    !! atomnames
    real*8, dimension(3,3), intent(in) :: alat
    !! lattice vectors in bohr
    character(len=80), intent(in) :: comment
    !! content that will be written to comment line
    integer :: io, ios, i
    real*8 :: Bohr_Ang = 0.52917721067
    real*8 :: alat_convert(3,3), xyz_convert(3,nat)
  
    alat_convert = alat * Bohr_Ang
    xyz_convert = rxyz * Bohr_Ang
  
    open(newunit=io,file=filename,iostat=ios)
    if(ios /= 0) stop "error opening output file"
    if ( len_trim(comment) > 0 ) then
      write(io, *) "# "//trim(comment)
    end if
    write(io,*) "lattice_vector", alat_convert(1,1), alat_convert(2,1), alat_convert(3,1)
    write(io,*) "lattice_vector", alat_convert(1,2), alat_convert(2,2), alat_convert(3,2)
    write(io,*) "lattice_vector", alat_convert(1,3), alat_convert(2,3), alat_convert(3,3)
    do i = 1, nat, 1
      write(io,*) "atom", xyz_convert(1,i), xyz_convert(2,i), xyz_convert(3,i), atomnames(i)
    end do
    close(io)
  end subroutine write_in
  subroutine write_periodic(filename, nat, rxyz, alat, atomnames, comment)
    !! writes a periodic structure to a file specified in filename variable.
    !! format is detected in ending of filename.
    !! allowed values:
    !! .asci, .cif, .gen, .in, .qe, .vasp and POSCAR (vasp)
    character(len=*), intent(in) :: filename
    !! filename of the file which is created
    integer, intent(in) :: nat
    !! number of atoms
    real*8, dimension(3, nat), intent(in) :: rxyz
    !! atom positions
    real*8, dimension(3, 3), intent(in) :: alat
    !! lattice vectors
    character(len=2), intent(in) :: atomnames(nat)
    !! String containing the chemical symbol of each atom.
    character(len=80), intent(in) :: comment
    !! content that will be written to comment line.
    integer :: l, vaspind
    character(len=10) :: filetype
  
    vaspind = index(filename, 'POSCAR')
    l = len_trim(filename)
  
    if ( l - vaspind == 5 .and. vaspind > 0) then !! vasp file
      call write_vasp(filename, nat, rxyz, alat, atomnames, comment)
      return
    end if
  
    call get_file_type(filename, filetype)
  
    select case (trim(filetype))
    case ("ascii")
      call write_ascii(filename, nat, rxyz, alat, atomnames, comment)
    case("cif")
      call write_cif(filename, nat, rxyz, alat, atomnames, comment)
    case("in")
      call write_in(filename, nat, rxyz, alat, atomnames, comment)
    case("gen")
      call write_dftb(filename, nat, rxyz, alat, atomnames, comment)
    case("qe")
      call write_quantum_espresso(filename, nat, rxyz, alat, atomnames)
    case('vasp')
      call write_vasp(filename, nat, rxyz, alat, atomnames, comment)
    case default
      print*, "filetype,", filetype
      print*, "filename, ", trim(filename)
      stop "unknown filetype write_periodic"
    end select
  
  end subroutine write_periodic
  
  subroutine write_quantum_espresso(filename, nat, rxyz, alat, atomnames)
    !! This subroutine appends cell parameters and atomic positions to the file specified in
    !! the filename variable. It is useful to create input files for quantum espresso
    use, intrinsic :: iso_fortran_env
    implicit none
    character(len=*), intent(in) :: filename
    !! filename of the file which is created
    integer, intent(in) :: nat
    !! number of atoms
    real(8), dimension(3, nat), intent(in) :: rxyz
    !! atom positions in bohr
    real(8), dimension(3, 3), intent(in) :: alat
    !! lattice vectors
    character(len=2), intent(in) :: atomnames(nat)
    !! String containing the chemical symbol of each atom.
  
    integer :: ntypes
    character(len=2), dimension(:), allocatable :: atom_types
    integer :: io, ios
    integer :: i
  
    call get_atom_types
  
    open(file=filename, newunit=io, iostat=ios, position="append")
    if (ios /= 0) stop "error reading file in write_quantum_espresso (libmoleulario)"
  
    write(io, *) "CELL_PARAMETERS bohr"
    do i = 1, 3, 1
      write(io,*) alat(1, i), alat(2, i), alat(3, i)
    end do
  
    write(io,*) "ATOMIC_POSITIONS bohr"
    do i = 1, nat, 1
      write(io,*) atomnames(i), rxyz(1,i), rxyz(2,i), rxyz(3,i)
    end do
  
    close(io)
  contains
    !! deprecated but should work. The subroutine write_quantum_espresso does and
    !! should not call this subroutine.
    subroutine get_atom_types
      integer :: ii, jj
      character(len=2), dimension(nat) :: an_copy
  
      ntypes = 0
      an_copy = atomnames
      do ii = 1, nat, 1
        if (len_trim(an_copy(ii)) /= 0) then ! new atom species found. delete all the following from list
          ntypes = ntypes + 1
          do jj = ii + 1, nat, 1
            if (an_copy(jj) == an_copy(ii)) then
              an_copy(jj) = ""
            end if
          end do
        end if
      end do
      if (allocated(atom_types)) deallocate(atom_types)
      allocate (atom_types(ntypes))
      jj = 0
      do ii = 1, nat, 1
        if (an_copy(ii) /= "") then
          jj = jj + 1
          atom_types(jj) = an_copy(ii)
        end if
      end do
    end subroutine get_atom_types
  end subroutine write_quantum_espresso
  subroutine write_vasp(filename, nat, rxyz, alat, atomnames, comment)
    !! writes periodic geometry  in  vasp format
    !! Multiple instances of an atom type must be contigouos. This subroutine
    !! takes care of this. ( C C C C H H Si Si ...)
    !! https://www.vasp.at/wiki/index.php/POSCAR
    implicit none
    character(len=*), intent(in) :: filename
    !! filename of the file which is created
    integer, intent(in) :: nat
    !! number of atoms
    real(8), dimension(3, nat), intent(in) :: rxyz
    !! atom positions in bohr
    real(8), dimension(3, 3), intent(in) :: alat
    !! lattice vectors (bohr)
    character(len=2), intent(in) :: atomnames(nat)
    !! String containing the chemical symbol of each atom.
    character(len=80), intent(in) :: comment
    !! content that will be written to the comment line
  
    ! private variables:
    integer :: io, ios, i
    real(8) :: Bohr_Ang = 0.52917721067
    real(8) :: p_out(3,nat), a_out(3,3)
    !! number of different atoms
    integer :: ntypes
    !! Array containing all the different elements. Each element is only listed once.
    !! Only the first ntype entries contain data, the rest is garbage.
    character(len=2), dimension(nat) :: atomtypes
    !! Encodes each element with a number. The element name of iat is accessed the following way:
    !! atomtypes(atomnumbernames(iat))
    integer :: atomnumbernames(nat), itype
    integer, allocatable, dimension(:) :: typepos, typecount
  
    if ( nat <= 0 ) then
      stop 'number of atoms is non positive in write_vasp'
    end if
  
    ! calculate ntypes and atomtypes
    ntypes = 1
    atomtypes(1) = atomnames(1)
    do i = 2, nat, 1
      do itype = 1, ntypes, 1
        if ( atomnames(i) == atomtypes(itype) ) then
          goto 123
        end if
      end do
      ntypes = ntypes + 1
      atomtypes(ntypes) = atomnames(i)
      123 continue
    end do
  
    ! calculate atomnumbernames
    allocate(typepos(ntypes), typecount(ntypes))
    typecount = 0
    do i = 1, nat, 1
      do itype = 1, ntypes, 1
        if ( atomnames(i) == atomtypes(itype) ) then
          atomnumbernames(i) = itype
          typecount(itype) = typecount(itype) + 1
        end if
      end do
    end do
    typepos(1) = 1
    do i = 2 ,ntypes
      typepos(i) = typecount(i-1) + typepos(i - 1)
    end do
  
    !! sort array
    do i = 1, nat, 1
      p_out(:, typepos(atomnumbernames(i)) ) = rxyz(:,i)
      typepos(atomnumbernames(i)) = typepos(atomnumbernames(i)) + 1
    end do
  
    deallocate(typepos)
  
  
    ! open file
    open (newunit=io, file=trim(adjustl(filename)), iostat=ios)
    if (ios /= 0) then
      print *, "Error opening file"//trim(filename)
      stop "Error opening file"
    end if
  
    !! write comment to file
    write(io, '(a)', iostat=ios) trim(comment)
    if(ios /= 0) then
      print*, 'error writing to file in write_vasp moleculario'
      print*, trim(filename)
      stop 'error writing to file'
    end if
  
    ! convert to angstroem
    p_out = p_out * Bohr_Ang
    a_out = alat * Bohr_Ang
  
    write(io, '(a6)', iostat=ios) '  1.00'
    if(ios/=0) stop 'error in write_vasp'
  
    !! write lattice
    do i = 1, 3, 1
      write(io, *, iostat=ios) a_out(1, i), a_out(2, i), a_out(3, i)
      if(ios /= 0) stop 'error writing lattice in write_vasp'
    end do
  
    ! write atomtypes
    write(io, *, iostat=ios) atomtypes(1:ntypes)
    if (ios /= 0 ) stop ' error writing atomtypes in moleculario'
  
    ! write atomcounts
    write(io, *, iostat=ios) typecount
    if (ios /= 0) stop 'error writing typecount in moleculario'
  
    ! write cartesian flag
    write(io, '(a)', iostat= ios) 'Cartesian'
    if (ios /= 0) stop 'error writing cartesian flag in moleculario'
  
    ! write positions
    do i = 1, nat, 1
      write(io, *, iostat  = ios) p_out(1,i), p_out(2,i), p_out(3,i)
      if ( ios /= 0) stop 'error writing position vasp moleculario'
    end do
  
    close(io, iostat=ios)
    if (ios /= 0 ) stop 'error closing vasp file moleculario'
  
  end subroutine write_vasp
  subroutine write_xyz(filename, nat, rxyz, atomnames, comment)
    !! writes xyz file. assumes units of rxyz are in bohr and writes them to file
    !! in a.u.
    implicit none
    character(len=*), intent(in) :: filename
    !! name of the output file
    integer, intent(in) :: nat
    !! number of atoms
    real*8, intent(in), dimension(3,nat) :: rxyz
    !! atomic positions in a.u
    character(len=2), intent(in), dimension(nat) :: atomnames
    !! chemical symbol of the atoms
    character(len=80), intent(in) :: comment
    !! content that will be written to the comment line
    integer :: io, ios, i
    real*8 :: Bohr_Ang = 0.52917721067
    real*8 :: xyz_convert(3,nat)
  
    xyz_convert = Bohr_Ang * rxyz
    open(newunit=io,file=filename,iostat=ios)
    if(ios/=0) stop "error opening output file"
    write(io,*) nat
    write(io,*) trim(comment)
    do i = 1, nat, 1
      write(io,*, iostat = ios) atomnames(i), xyz_convert(1,i), xyz_convert(2,i), xyz_convert(3,i)
      if ( ios/=0 ) then
        print*, "error writing line: ", i
        stop
      end if
    end do
    close(io)
  end subroutine write_xyz
  subroutine isElement(atom_in, is_element)
    !! returns wether chemical symbol is an element or not
    implicit none
    character(len=2), intent(in) :: atom_in
    !! String containing chemical symbol
    logical, intent(out) :: is_element
    !! logical, true if string is element, false if not
    select case (trim(atom_in))
    case('H')
       is_element = .true.
    case('He')
       is_element = .true.
    case('Li')
       is_element = .true.
    case('Be')
       is_element = .true.
    case('B')
       is_element = .true.
    case('C')
       is_element = .true.
    case('N')
       is_element = .true.
    case('O')
       is_element = .true.
    case('F')
       is_element = .true.
    case('Ne')
       is_element = .true.
    case('Na')
       is_element = .true.
    case('Mg')
       is_element = .true.
    case('Al')
       is_element = .true.
    case('Si')
       is_element = .true.
    case('P')
       is_element = .true.
    case('S')
       is_element = .true.
    case('Cl')
       is_element = .true.
    case('Ar')
       is_element = .true.
    case('K')
       is_element = .true.
    case('Ca')
       is_element = .true.
    case('Sc')
       is_element = .true.
    case('Ti')
       is_element = .true.
    case('V')
       is_element = .true.
    case('Cr')
       is_element = .true.
    case('Mn')
       is_element = .true.
    case('Fe')
       is_element = .true.
    case('Co')
       is_element = .true.
    case('Ni')
       is_element = .true.
    case('Cu')
       is_element = .true.
    case('Zn')
       is_element = .true.
    case('Ga')
       is_element = .true.
    case('Ge')
       is_element = .true.
    case('As')
       is_element = .true.
    case('Se')
       is_element = .true.
    case('Br')
       is_element = .true.
    case('Kr')
       is_element = .true.
    case('Rb')
       is_element = .true.
    case('Sr')
       is_element = .true.
    case('Y')
       is_element = .true.
    case('Zr')
       is_element = .true.
    case('Nb')
       is_element = .true.
    case('Mo')
       is_element = .true.
    case('Tc')
       is_element = .true.
    case('Ru')
       is_element = .true.
    case('Rh')
       is_element = .true.
    case('Pd')
       is_element = .true.
    case('Ag')
       is_element = .true.
    case('Cd')
       is_element = .true.
    case('In')
       is_element = .true.
    case('Sn')
       is_element = .true.
    case('Sb')
       is_element = .true.
    case('Te')
       is_element = .true.
    case('I')
       is_element = .true.
    case('Xe')
       is_element = .true.
    case('Cs')
       is_element = .true.
    case('Ba')
       is_element = .true.
    case('La')
       is_element = .true.
       !     case('Ce')
       !     case('Pr')
       !     case('Nd')
       !     case('Pm')
       !     case('Sm')
       !     case('Eu')
       !     case('Gd')
       !     case('Tb')
       !     case('Dy')
       !     case('Ho')
       !     case('Er')
       !     case('Tm')
       !     case('Yb')
    case('Lu')
       is_element = .true.
    case('Hf')
       is_element = .true.
    case('Ta')
       is_element = .true.
    case('W')
       is_element = .true.
    case('Re')
       is_element = .true.
    case('Os')
       is_element = .true.
    case('Ir')
       is_element = .true.
    case('Pt')
       is_element = .true.
    case('Au')
       is_element = .true.
    case('Hg')
       is_element = .true.
    case('Tl')
       is_element = .true.
    case('Pb')
       is_element = .true.
    case('Bi')
       is_element = .true.
       !     case('Po')
       !     case('At')
    case('Rn')
       is_element = .true.
      case default
        is_element = .FALSE.
    end select
  end subroutine isElement
  
  subroutine cart2frac(nat, alat, rxyz, xyzred)
    !! converts cartesian coordinates rxyz to reduced coordinates xyzred
    implicit none
    integer, intent(in) :: nat
    !! Number of Atoms
    real*8, intent(in), dimension(3, 3) :: alat
    !! Lattice Vectors.
    real*8, intent(in), dimension(3, nat) :: rxyz
    !! Position of the Atoms in cartesian coorinates.
    real*8, intent(out), dimension(3, nat) :: xyzred
    !! Position of the Atoms in reduced coordinates.
  
    !private variables
    real*8, dimension(3, 3) :: alatinv
    !! inverse of the lattice matrix
    real*8 :: div
    !! inverse volume of the lattice matrix
  
    div = alat(1, 1)*alat(2, 2)*alat(3, 3) - alat(1, 1)*alat(2, 3)*alat(3, 2) - &
          alat(1, 2)*alat(2, 1)*alat(3, 3) + alat(1, 2)*alat(2, 3)*alat(3, 1) + &
          alat(1, 3)*alat(2, 1)*alat(3, 2) - alat(1, 3)*alat(2, 2)*alat(3, 1)
    if ( abs(div) < 1.d-7 ) then
      stop "cell is singular in cart2frac"
    end if
    div = 1.d0/div
    alatinv(1, 1) = (alat(2, 2)*alat(3, 3) - alat(2, 3)*alat(3, 2))
    alatinv(1, 2) = -(alat(1, 2)*alat(3, 3) - alat(1, 3)*alat(3, 2))
    alatinv(1, 3) = (alat(1, 2)*alat(2, 3) - alat(1, 3)*alat(2, 2))
    alatinv(2, 1) = -(alat(2, 1)*alat(3, 3) - alat(2, 3)*alat(3, 1))
    alatinv(2, 2) = (alat(1, 1)*alat(3, 3) - alat(1, 3)*alat(3, 1))
    alatinv(2, 3) = -(alat(1, 1)*alat(2, 3) - alat(1, 3)*alat(2, 1))
    alatinv(3, 1) = (alat(2, 1)*alat(3, 2) - alat(2, 2)*alat(3, 1))
    alatinv(3, 2) = -(alat(1, 1)*alat(3, 2) - alat(1, 2)*alat(3, 1))
    alatinv(3, 3) = (alat(1, 1)*alat(2, 2) - alat(1, 2)*alat(2, 1))
    alatinv = alatinv * div
    xyzred = matmul(alatinv, rxyz)
  end subroutine cart2frac
  
  subroutine frac2cart(nat, alat, xyzred, rxyz)
    !! Converts reduced coordinates xyzred to cartesian coordinates rxyz
    implicit none
    integer, intent(in) :: nat
    !! Number of atoms.
    real*8, intent(in), dimension(3, 3) :: alat
    !! Lattice Vecors
    real*8, dimension(3, nat), intent(in) :: xyzred
    !! Position of the atoms in reduced coordinates.
    real*8, dimension(3, nat), intent(out) :: rxyz
    !! Position of the atoms in cartesian coordinates.
    rxyz = matmul(alat, xyzred)
  end subroutine frac2cart
  
  subroutine back2cell(nat, rxyz, alat)
    !! Translates atoms outside the cell back into the cell.
    implicit none
    integer, intent(in) :: nat
    !! Number of atoms
    real*8, dimension(3, nat), intent(inout) :: rxyz
    !! Positions of the atoms.
    real*8, dimension(3, 3), intent(in) :: alat
    !! Lattice vectors of the atoms.
  
    ! private variables
    real*8, dimension(3, nat) :: xyzred
    !! Position of the atoms in cartesian coordinates.
  
    call cart2frac(nat, alat, rxyz, xyzred)
    xyzred = modulo(xyzred, 1.d0)
    call frac2cart(nat, alat, xyzred, rxyz)
  end subroutine back2cell
  subroutine get_num_words(string, nwords, len_string)
    !! gets the number of words (characters separated by one or multiple spaces) in a string
    implicit none
    integer, intent(out) :: nwords
    !! number of words
    integer, intent(in) :: len_string
    !! length of input string
    character(len=len_string), intent(in) :: string
    !! input string were number of words should be counted.
    integer :: ichar
    !! iteration variable of characters of string
    logical :: in_word
  
    if ( len_trim(string) == 0 ) then
      nwords = 0
      return
    end if
    ichar = 1
    ! advance until first word starts
    do while ( string(ichar:ichar) == " " .and. ichar <= len_string)
      ichar = ichar + 1
    end do
    nwords = 0
    ichar = 1
    in_word = .FALSE.
    do while ( ichar <= len_string ) !! find the atom site column
      if ( (string(ichar:ichar) == " ") .and. in_word ) then
        in_word = .FALSE.
        nwords = nwords + 1
      else
        if(string(ichar:ichar) /= " ") in_word = .TRUE.
      end if
      ichar = ichar + 1
    end do
  end subroutine get_num_words
  
  
  subroutine get_words(string, nwords, len_string, words, len_words)
    !! splits string seperated by spaces into words.
    implicit none
    !! number of words that are in string (use getnumword to calculate this)
    integer, intent(in) :: nwords
    !! length of the string
    integer, intent(in) :: len_string
    !! string that shoul be split
    character(len=len_string), intent(in) :: string
    !! max length of each word
    integer, intent(in) :: len_words
    !! string array containing all the words
    character(len=len_words), dimension(nwords) :: words
    integer :: ichar, iword, wordend
    ichar = 1
    words = ""
    do iword = 1, nwords, 1
      do while ( string(ichar:ichar) == " " .and. ichar <= len_string)
        ichar = ichar + 1
      end do
      wordend = index(string(ichar:len_string), " ")
      if ( wordend <= 0 ) then
        words(iword) = string(ichar:len_string)
        return
      end if
      words(iword) = string(ichar:(ichar + wordend-1))
      ichar = ichar + wordend - 1
    end do
  end subroutine get_words
  
  subroutine get_file_type(filename, filetype)
    !! returns the filetype of file with filename.
    !! example: test.txt -> txt; test.dat -> dat
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(out) :: filetype
    integer :: sep_pos, len_file
    len_file = len_trim(filename)
    do sep_pos = len_file, 1, -1
      if ( filename(sep_pos:sep_pos) == "." ) then
        exit
      end if
    end do
    if ( len_file == 0 ) then
      stop "empy filename in get_file_type"
    end if
    if ( sep_pos <= 0 ) then
      filetype = ""
      return
    end if
    if ( sep_pos == len_file ) then
      filetype = ""
      return
    end if
    filetype = filename((sep_pos + 1):len_file)
  end subroutine get_file_type
  
  subroutine get_basename(filename, bname)
    implicit none
    !! filename where basename should be extracted from
    character(len=*), intent(in) :: filename
    !! basename of filename (stripped of ending and prepending path)
    character(len=*), intent(out) :: bname
    integer :: lenfile, i
    !! position of point in filename
    integer :: point_pos
    !! position of last "/".
    integer :: slash_pos
    lenfile = len_trim(filename)
    if ( lenfile == 0) then
      stop "empty string in get basename"
    end if
    do point_pos = lenfile, 1, -1
      if ( filename(point_pos:point_pos) == "." ) then
        exit
      end if
    end do
    if ( point_pos <= 0 ) then
      stop "filename does not contain a dot in moleculario library."
    end if
    slash_pos = 0
    do i = 1, point_pos, 1
      if ( filename(i:i) == "/" ) then
        slash_pos = i
      end if
    end do
  
    bname = filename(slash_pos+1:point_pos - 1)
  end subroutine get_basename
  
  subroutine sym2amu(sym, amu)
    !! returns the atomic mass of atom with chemical symbol sym
    !! mass in amu taken from https://www.angelo.edu/faculty/kboudrea/periodic/structure_numbers.htm
    implicit none
    real*8  :: amu
    !! atomic mass unit
    character(len=2) :: sym
    !! chemical symbol
    select case (trim(sym))
    case ('H')
      amu = 1.00797
    case ('He')
      amu = 4.00260
    case ('Li')
      amu = 6.941
    case ('Be')
      amu = 9.01218
    case ('B')
      amu = 10.81
    case ('C')
      amu = 12.011
    case ('N')
      amu = 14.0067
    case ('O')
      amu = 15.9994
    case ('F')
      amu = 18.998403
    case ('Ne')
      amu = 20.179
    case ('Na')
      amu = 22.98977
    case ('Mg')
      amu = 24.305
    case ('Al')
      amu = 26.98154
    case ('Si')
      amu = 28.0855
    case ('P')
      amu = 30.97376
    case ('S')
      amu = 32.06
    case ('Cl')
      amu = 35.453
    case ('Ar')
      amu = 39.948
    case ('K')
      amu = 39.0983
    case ('Ca')
      amu = 40.08
    case ('Sc')
      amu = 44.9559
    case ('Ti')
      amu = 47.90
    case ('V')
      amu = 50.9415
    case ('Cr')
      amu = 51.996
    case ('Mn')
      amu = 54.9380
    case ('Fe')
      amu = 55.847
    case ('Co')
      amu = 58.9332
    case ('Ni')
      amu = 58.70
    case ('Cu')
      amu = 63.546
    case ('Zn')
      amu = 65.38
    case ('Ga')
      amu = 69.72
    case ('Ge')
      amu = 72.59
    case ('As')
      amu = 74.9216
    case ('Se')
      amu = 78.96
    case ('Br')
      amu = 79.904
    case ('Kr')
      amu = 83.80
    case ('Rb')
      amu = 85.4678
    case ('Sr')
      amu = 87.62
    case ('Y')
      amu = 88.9059
    case ('Zr')
      amu = 91.22
    case ('Nb')
      amu = 92.9064
    case ('Mo')
      amu = 95.94
    case ('Tc')
      amu = 98
    case ('Ru')
      amu = 101.07
    case ('Rh')
      amu = 102.9055
    case ('Pd')
      amu = 106.4
    case ('Ag')
      amu = 107.868
    case ('Cd')
      amu = 112.41
    case ('In')
      amu = 114.82
    case ('Sn')
      amu = 118.69
    case ('Sb')
      amu = 121.75
    case ('Te')
      amu = 127.60
    case ('I')
      amu = 126.9045
    case ('Xe')
      amu = 131.30
    case ('Cs')
      amu = 132.9054
    case ('Ba')
      amu = 137.33
    case ('La')
      amu = 138.9055
    case('Ce')
      amu = 140.12
    case('Pr')
      amu = 140.9077
    case('Nd')
      amu = 144.24
    case('Pm')
      amu = 145
    case('Sm')
      amu = 150.4
    case('Eu')
      amu = 151.96
    case('Gd')
      amu = 157.25
    case('Tb')
      amu = 158.9254
    case('Dy')
      amu = 162.50
    case('Ho')
      amu = 164.9304
    case('Er')
      amu = 167.26
    case('Tm')
      amu = 168.9342
    case('Yb')
      amu = 173.04
    case ('Lu')
      amu = 174.967
    case ('Hf')
      amu = 178.49
    case ('Ta')
      amu = 180.9479
    case ('W')
      amu = 183.85
    case ('Re')
      amu = 186.207
    case ('Os')
      amu = 190.2
    case ('Ir')
      amu = 192.22
    case ('Pt')
      amu = 195.09
    case ('Au')
      amu = 196.9665
    case ('Hg')
      amu = 200.59
    case ('Tl')
      amu = 204.37
    case ('Pb')
      amu = 207.2
    case ('Bi')
      amu = 208.9804
    case('Po')
      amu = 209
    case('At')
      amu = 210
    case ('Rn')
      amu = 222
    case('Fr')
      amu = 223
    case('Ra')
      amu = 226.0254
    case('Ac')
      amu = 227.0278
    case('Th')
      amu = 232.0381
    case('Pa')
      amu = 231.0359
    case('U')
      amu = 238.029
    case('Np')
      amu = 237.0482
    case('Pu')
      amu = 242
    case('Am')
      amu = 243
    case('Cm')
      amu = 247
    case default
      print *, " Not recognized atomic type "//sym
      stop
    end select
  
  end subroutine sym2amu
  