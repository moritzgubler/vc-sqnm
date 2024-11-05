import numpy as np
import numba

# @numba.jit(nopython=True)
def nnlist(nat, nnbrx, alat, cutoff, rxyz):
#   implicit none
#   integer, parameter :: nwork=1000
#   real(8) :: rxyz(3, nat), rel(5, nat*nnbrx)
#   real(8) :: alat(3, 3)
#   integer :: lsta(2, nat), lstb(nnbrx*nat)
#   real(8) :: alatalat(3, 3), eigalat(3), workalat(nwork)
#   integer :: ixyzmax, info, ind, iat, jat, ix, iy, iz, nat, nnbrx, i, j
#   real(8) :: cutoff2, cutoff, dist2, relx, rely, relz, tt, ttinv, xj, yj, zj

    # if (alat(1, 1)*alat(2, 2)*alat(3, 3) .eq. 0.d0) then ! no periodic boundary condition
    if np.abs(np.linalg.det(alat)) < 1e-10:
        ixyzmax = 0
    else:  #! periodic boundary conditions
        alatalat = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                alatalat[i, j] = alat[0, i] * alat[0, j] + alat[1, i] * alat[1, j] + alat[2, i] * alat[2, j]

    #   call dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
    #   !   write(*,*) !  'alat !  EVals',eigalat !  write(*,*) 'ixyzmax',int(sqrt(1.d0/eigalat(1))*radius_cutoff)
    #   ! ixyzmax determines over how many periodiv images one has to search to fill the sphere with atoms
    #   ixyzmax = int(sqrt(1.d0/eigalat(1))*cutoff) + 1
    # end if
        eigval, eigvec = np.linalg.eigh(alatalat)
        ixyzmax = int(np.sqrt(1.0 / eigval[0]) * cutoff) + 1
    cutoff2 = cutoff * cutoff

#   !  write(*,*) 'ixyzmax ',ixyzmax

    ind = 0
    lsta = np.zeros((2, nat), dtype=np.int32)
    lstb = np.zeros(nnbrx * nat, dtype=np.int32)
    rel = np.zeros((5, nat * nnbrx))
#   do iat = 1, nat
    for iat in range(nat):
        lsta[0, iat] = ind

    # do jat = 1, nat
        for jat in range(nat):

    #   do ix = -ixyzmax, ixyzmax
    #     do iy = -ixyzmax, ixyzmax
    #       do iz = -ixyzmax, ixyzmax
            for ix in range(-ixyzmax, ixyzmax + 1):
               for iy in range(-ixyzmax, ixyzmax + 1):
                   for iz in range(-ixyzmax, ixyzmax + 1):
                        xj = rxyz[0, jat] + ix * alat[0, 0] + iy*alat[0, 1] + iz*alat[0, 2]
                        yj = rxyz[1, jat] + ix * alat[1, 0] + iy*alat[1, 1] + iz*alat[1, 2]
                        zj = rxyz[2, jat] + ix * alat[2, 0] + iy*alat[2, 1] + iz*alat[2, 2]
                        relx = xj - rxyz[0, iat]
                        rely = yj - rxyz[1, iat]
                        relz = zj - rxyz[2, iat]
                        dist2 = relx**2 + rely**2 + relz**2

                        if dist2 > 1.0e-20 and dist2 <= cutoff2:
                            ind = ind + 1
                            if (ind >= nnbrx*nat):
                                print('enlarge nnbrx')
                                # quit()
                            lstb[ind] = jat
                            tt = np.sqrt(dist2)
                            ttinv = 1.0 / tt
                            rel[0, ind] = relx*ttinv
                            rel[1, ind] = rely*ttinv
                            rel[2, ind] = relz*ttinv
                            rel[3, ind] = tt
                            rel[4, ind] = ttinv
    #         end if
    #       end do
    #     end do
    #   end do
    # end do
        lsta[1, iat] = ind
#   end do

#   !  do iat=1,nat
#   !  write(*,'(i3,1x,20(1x,i2))') iat,(lstb(j),j=lsta(1,iat),lsta(2,iat))
#   !  write(*,'(i3,1x,20(1x,e9.2))') iat,(rel(4,j),j=lsta(1,iat),lsta(2,iat))
#   !  enddo
    return lsta, lstb, rel
# end subroutine nnlist

# @numba.jit(nopython=True)
def energyandforces_bazant(nat, alat0, rxyz0):
#
#  !   INTERFACE
#  !   ---------
#  !    Internally things are computed in Angstroem and eV, The input (rxyz0,alat0) is however in bohr and
#  !     the output energy in Hartree, fxyz and deralat  in Hartree/Bohr
#
#  !     nat : number of particles
#  !     rxyz0 : array (3,N) of positions
#  !     alat0 : lattice vecors (each column is one vector)
#  !     Etot : returned energy in eV
#  !     fxyz : returned forces array (3,N) in eV/Angstroms
#  !     alatder  :  derivative of energy with respsct to altice vectors
#
#  !     neighbors(p_nbrs(i)),...,neighbors(p_nbrs(i+1)) are the
#  !     atoms which are "neighbors" of atom i, using a standard Verlet
#  !     neighbor list. These are a global arrays, not passed, that are declared
#  !     in "edip_neighbors_include.h". This way of storing atomi! positions
#  !     is not unique, and will require a customized patch to the main MD program.
#
#  !     The parameters of the potential initialized in input_EDIP_params() are global
#  !     variables declared in "edip_pot_include.h".
#
#  !   PARAMETERS
#  !   ----------
#
#  !    par_cap_A,par_cap_B,par_rh,par_a,par_sig
#  !    par_lam,par_gam,par_b,par_c,par_delta
#  !    par_mu,par_Qo,par_palp,par_bet,par_alp
#
#  !    5.6714030     2.0002804     1.2085196     3.1213820     0.5774108
#  !    1.4533108     1.1247945     3.1213820     2.5609104    78.7590539
#  !    0.6966326   312.1341346     1.4074424     0.0070975     3.1083847
#
#  !    Connection between these parameters and those given in the paper,
#  !    Justo et al., Phys. Rev. B 58, 2539 (1998):
#
#  !    A((B/r)**rh-palp*exp(-bet*Z*Z)) = A'((B'/r)**rh-exp(-bet*Z*Z))
#
#  !    so in the paper (')
#  !    A' = A*palp
#  !    B' = B * palp**(-1/rh)
#  !    eta = detla/Qo
#
#  !    Non-adjustable parameters for tau(Z) from Ismail & Kaxiras, 1993,
#  !    also in Bazant, Kaxiras, Justo, PRB (1997):
#
#  !    u1 = -0.165799;
#  !    u2 = 32.557;
#  !    u3 = 0.286198;
#  !    u4 = 0.66;
#  implicit none
#
#  !  ------------------------- VARIABLE DECLARATIONS -------------------------
#  integer:: nat
#  integer, parameter :: nnbrx = 50 !number of geighbors
#  real*8 :: Ha_eV, Bohr_Ang
#
#  integer i, j, k, l, n, n2, n3, nz, iat
#  real*8 :: dx, dy, dz, r, asqr
#  real*8 :: rinv, rmainv, xinv, xinv3, den, Z, fZ
#  real*8 :: dV2j, dV2ijx, dV2ijy, dV2ijz, pZ, dp
#  real*8 :: temp0, temp1
#  real*8 :: Qort, muhalf, u5
#  real*8 :: rmbinv, winv, dwinv, tau, dtau, lcos, x, H, dHdx, dhdl
#  real*8 :: dV3rij, dV3rijx, dV3rijy, dV3rijz
#  real*8 :: dV3rik, dV3rikx, dV3riky, dV3rikz
#  real*8 :: dV3l, dV3ljx, dV3ljy, dV3ljz, dV3lkx, dV3lky, dV3lkz
#  real*8 :: dV2dZ, dxdZ, dV3dZ
#  real*8 :: dEdrl, dEdrlx, dEdrly, dEdrlz
#  real*8 :: bmc, cmbinv
#  real*8 :: fjx, fjy, fjz, fkx, fky, fkz
#
    nnbrx = 50
#  real*8 s2_t0(nnbrx), s2_t1(nnbrx), s2_t2(nnbrx), s2_t3(nnbrx), s2_dx(nnbrx), s2_dy(nnbrx), s2_dz(nnbrx), s2_r(nnbrx), &
#    s3_g(nnbrx), s3_dg(nnbrx), s3_rinv(nnbrx), s3_dx(nnbrx), s3_dy(nnbrx), s3_dz(nnbrx), s3_r(nnbrx), &
#    sz_df(nnbrx), sz_sum(nnbrx), sz_dx(nnbrx), sz_dy(nnbrx), sz_dz(nnbrx), sz_r(nnbrx)

    s2_t0 = np.zeros(nnbrx)
    s2_t1 = np.zeros(nnbrx)
    s2_t2 = np.zeros(nnbrx)
    s2_t3 = np.zeros(nnbrx)
    s2_dx = np.zeros(nnbrx)
    s2_dy = np.zeros(nnbrx)
    s2_dz = np.zeros(nnbrx)
    s2_r = np.zeros(nnbrx)
    s3_g = np.zeros(nnbrx)
    s3_dg = np.zeros(nnbrx)
    s3_rinv = np.zeros(nnbrx)
    s3_dx = np.zeros(nnbrx)
    s3_dy = np.zeros(nnbrx)
    s3_dz = np.zeros(nnbrx)
    s3_r = np.zeros(nnbrx)
    sz_df = np.zeros(nnbrx)
    sz_sum = np.zeros(nnbrx)
    sz_dx = np.zeros(nnbrx)
    sz_dy = np.zeros(nnbrx)
    sz_dz = np.zeros(nnbrx)
    sz_r = np.zeros(nnbrx)

#  integer num2(nnbrx), num3(nnbrx), numz(nnbrx)
    num2 = np.zeros(nnbrx, dtype=np.int32)
    num3 = np.zeros(nnbrx, dtype=np.int32)
    numz = np.zeros(nnbrx, dtype=np.int32)
#
#  integer nj, nk, nl
#  !   indices for the store arrays
#  real*8 ::  virial, virial_xyz(3)
#
#  real*8 ::    par_cap_A, par_cap_B, par_rh, par_a, par_sig
#  real*8 ::    par_lam, par_gam, par_b, par_c, par_delta
#  real*8 ::    par_mu, par_Qo, par_palp, par_bet, par_alp
#  real*8 ::    u1, u2, u3, u4
#  real*8 :: par_bg
#  real*8 :: par_eta
#  real*8 :: cutoff
#  real*8 :: delta_safe
#
#  ! My variables
#  integer :: lsta(2, nat), lstb(nnbrx*nat)
#  real*8  :: rel(5, nat*nnbrx)
#  real*8 :: etot, ener_iat, ener, ener2
#  real*8 :: alat0(3, 3), alat(3, 3), deralat(3, 3)
#  real*8 :: rxyz0(3, nat), rxyz(3, nat), fxyz(3, nat), alatinv(3, 3)
#  real*8 :: si1_sj1, si2_sj2, si3_sj3
#  real(8), intent(out) :: stress(3, 3)
#  real(8) :: vol

#  !         End my variables
#  ! EDIP parameters
#  !         taken from Justo et al., Phys. Rev. B 58, 2539 (1998).

    par_cap_A = 5.6714030
    par_cap_B = 2.0002804
    par_rh = 1.2085196
    par_a = 3.1213820
    par_sig = 0.5774108
    par_lam = 1.4533108
    par_gam = 1.1247945
    par_b = 3.1213820
    par_c = 2.5609104
    par_delta = 78.7590539
    par_mu = 0.6966326
    par_Qo = 312.1341346
    par_palp = 1.4074424
    par_bet = 0.0070975
    par_alp = 3.1083847
  
    par_bg = par_a
    par_eta = par_delta/par_Qo
    cutoff = par_a  
    delta_safe = 0.2
  
    u1 = -0.165799
    u2 = 32.557
    u3 = 0.286198
    u4 = 0.66
  #   !end parameters
  
    Ha_eV = 27.211399
    Bohr_Ang = 0.529177
  #   !Do some preparation including the construction of the pair list
  #   do j = 1, 3
  #   do i = 1, 3
  #     alat(i, j) = alat0(i, j)*Bohr_Ang
  #   end do
  #   end do
    alat = alat0 * Bohr_Ang

#   do iat = 1, nat
#     for i in range(nat):
#     rxyz(1, iat) = rxyz0(1, iat)*Bohr_Ang
#     rxyz(2, iat) = rxyz0(2, iat)*Bohr_Ang
#     rxyz(3, iat) = rxyz0(3, iat)*Bohr_Ang
#   end do
    rxyz = rxyz0 * Bohr_Ang
    # call back2cell(nat, rxyz, alat)

    alatinv = np.linalg.inv(alat)

    # call nnlist(nat, nnbrx, alat, cutoff, rxyz, lsta, lstb, rel)
    lsta, lstb, rel = nnlist(nat, nnbrx, alat, cutoff, rxyz)
    # call invertalat(alat, alatinv)

    # !Allocation of temporary arrays

    # !          L_x_div_2 = L_x/2.0D0
    # !          L_y_div_2 = L_y/2.0D0
    # !          L_z_div_2 = L_z/2.0D0

    # !          do i=1, N_own
    # fxyz = 0.0
    fxyz = np.zeros((3, nat))
    ener = 0.0
    ener2 = 0.0

    virial = 0.0
    # virial_xyz(:) = 0.0d0
    virial_xyz = np.zeros(3)
    # deralat = 0.d0
    deralat = np.zeros((3, 3))

    # !   COMBINE COEFFICIENTS

    asqr = par_a * par_a
    Qort = np.sqrt(par_Qo)
    muhalf = par_mu*0.50
    u5 = u2*u4
    bmc = par_b - par_c
    cmbinv = 1.0 / (par_c - par_b)

#   !  --- LEVEL 1: OUTER LOOP OVER ATOMS ---

#   do i = 1, nat
    for i in range(nat):

        #!   RESET COORDINATION AND NEIGHBOR NUMBERS
        ener_iat = 0.0
        Z = 0.0
        n2 = 0 # used to be 1
        n3 = 0 # used to be 1
        nz = 0 # used to be 1

        # !  --- LEVEL 2: LOOP PREPASS OVER PAIRS ---

        # !            do n=p_nbrs(i), p_nbrs(i+1)-1
        # !              j = neighbors(n)
        # do n = lsta(1, i), lsta(2, i)
        for n in range(lsta[0, i], lsta[1, i] + 1):
            j = lstb[n]

            #   !   PARTS OF TWO-BODY INTERACTION r<par_a

            num2[n2] = j
            #   !                rinv = 1.0/r
            #   !                dx = dx * rinv
            #   !                dy = dy * rinv
            #   !                dz = dz * rinv
            dx = -rel[0, n]
            dy = -rel[1, n]
            dz = -rel[2, n]
            r = rel[3, n]
            rinv = rel[4, n]

            rmainv = 1.0 / (r - par_a)
            s2_t0[n2] = par_cap_A * np.exp(par_sig*rmainv)
            s2_t1[n2] = (par_cap_B*rinv)**par_rh
            s2_t2[n2] = par_rh*rinv
            s2_t3[n2] = par_sig*rmainv*rmainv
            s2_dx[n2] = dx
            s2_dy[n2] = dy
            s2_dz[n2] = dz
            s2_r[n2] = r
            n2 = n2 + 1

            if (n2 > nnbrx):
                print('WARNING1 enlarge nnbrx')
                # quit()

            # !!Additional part from stefan
            # !! coordination number calculated with soft cutoff between first and
            # !! second nearest neighbor
            # !        if (r.le.2.36d0) then
            # !        coord_iat=coord_iat+1.d0
            # !        else if (r.ge.3.83d0) then
            # !        else
            # !        xarg=(r-2.36d0)*(1.d0/(3.83d0-2.36d0))
            # !        coord_iat=coord_iat+(2*xarg+1.d0)*(xarg-1.d0)**2
            # !        endif
            # !-----------------------------

            # !   RADIAL PARTS OF THREE-BODY INTERACTION r<par_b

            if (r < par_bg):
                num3[n3] = j
                rmbinv = 1.0 / (r - par_bg)
                temp1 = par_gam*rmbinv
                temp0 = np.exp(temp1)
                s3_g[n3] = temp0
                s3_dg[n3] = -rmbinv*temp1*temp0
                s3_dx[n3] = dx
                s3_dy[n3] = dy
                s3_dz[n3] = dz
                s3_rinv[n3] = rinv
                s3_r[n3] = r
                n3 = n3 + 1

                # !       if(fixZ .eq. 0) then
                # !Additional part from Stefan
                if (n3 >= nnbrx):
                    print('WARNING2 enlarge nnbrx')

                # !   COORDINATION AND NEIGHBOR FUNCTION par_c<r<par_b

                if (r < par_b):
                    if (r < par_c):
                        Z = Z + 1.0
                    else:
                        xinv = bmc / (r - par_c)
                        xinv3 = xinv*xinv*xinv
                        den = 1.0 / (1.0 - xinv3)
                        temp1 = par_alp*den
                        fZ = np.exp(temp1)
                        Z = Z + fZ
                        numz[nz] = j
                        sz_df[nz] = fZ * temp1 * den * 3.0 * xinv3 * xinv * cmbinv
                        # !   df/dr
                        sz_dx[nz] = dx
                        sz_dy[nz] = dy
                        sz_dz[nz] = dz
                        sz_r[nz] = r
                        nz = nz + 1
                        # !Additional part from Stefan
                        if (nz > nnbrx):
                            print('WARNING3 enlarge nnbrx')
                            # quit()
                    # end if
                    #   !  r < par_C
                # end if
                # !  r < par_b
            # end if
            #   !  fixZ .eq. 0
            #   !                end if
            #   !  r < par_bg
            #   !              end if
            #   !  rsqr < asqr
            #   !              end if
            #   !  dz < par_a
            #   !              end if
            #   !  dy < par_a
            #   !              end if
            #   !  dVz < par_a
        # end do

        # !            if(fixZ .ne. 0) then
        # !
        # !              Z = tricks_Zfix
        # !              pZ = par_palp*dexp(-par_bet*Z*Z)
        # !              dp = 0.0

        # !            else

        # !   ZERO ACCUMULATION ARRAY FOR ENVIRONMENT FORCES

        # do nl = 1, nz - 1
        #   sz_sum(nl) = 0.0d0
        # end do
        sz_sum[0:nz - 1] = 0.0

        # !   ENVIRONMENT-DEPENDENCE OF PAIR INTERACTION

        temp0 = par_bet * Z
        pZ = par_palp * np.exp(-temp0 * Z)
        # !   bond order
        dp = -2.0*temp0*pZ
        # !   derivative of bond order


        # !  --- LEVEL 2: LOOP FOR PAIR INTERACTIONS ---

        # do nj = 1, n2 - 1
        for nj in range(n2 - 1):

            temp0 = s2_t1[nj] - pZ

            # !   two-body energy V2(rij,Z)

            ener_iat = ener_iat + temp0*s2_t0[nj]

            # two-body forces

            dV2j = -s2_t0[nj]*(s2_t1[nj]*s2_t2[nj] + temp0*s2_t3[nj])
            # !   dV2/dr
            dV2ijx = dV2j*s2_dx[nj]
            dV2ijy = dV2j*s2_dy[nj]
            dV2ijz = dV2j*s2_dz[nj]
            fxyz[0, i] = fxyz[0, i] + dV2ijx
            fxyz[1, i] = fxyz[1, i] + dV2ijy
            fxyz[2, i] = fxyz[2, i] + dV2ijz
            j = num2[nj]
            fxyz[0, j] = fxyz[0, j] - dV2ijx
            fxyz[1, j] = fxyz[1, j] - dV2ijy
            fxyz[2, j] = fxyz[2, j] - dV2ijz

            # ! dV2/dr contribution to virial

            virial_xyz[0] = virial_xyz[0] - s2_r[nj] * (dV2ijx * s2_dx[nj])
            virial_xyz[1] = virial_xyz[1] - s2_r[nj] * (dV2ijy * s2_dy[nj])
            virial_xyz[2] = virial_xyz[2] - s2_r[nj] * (dV2ijz * s2_dz[nj])
            virial = virial - s2_r[nj] * (dV2ijx*s2_dx[nj] + dV2ijy*s2_dy[nj] + dV2ijz*s2_dz[nj])

            #   !Cell gradient part
            #   !My own implementation
            si1_sj1 = alatinv[0, 0]*s2_dx[nj] + alatinv[0, 1]*s2_dy[nj] + alatinv[0, 2]*s2_dz[nj]
            si2_sj2 = alatinv[1, 0]*s2_dx[nj] + alatinv[1, 1]*s2_dy[nj] + alatinv[1, 2]*s2_dz[nj]
            si3_sj3 = alatinv[2, 0]*s2_dx[nj] + alatinv[2, 1]*s2_dy[nj] + alatinv[2, 2]*s2_dz[nj]
            deralat[0, 0] = deralat[0, 0] - s2_r[nj] * dV2ijx * si1_sj1
            deralat[0, 1] = deralat[0, 1] - s2_r[nj] * dV2ijx * si2_sj2
            deralat[0, 2] = deralat[0, 2] - s2_r[nj] * dV2ijx * si3_sj3
            deralat[1, 0] = deralat[1, 0] - s2_r[nj] * dV2ijy * si1_sj1
            deralat[1, 1] = deralat[1, 1] - s2_r[nj] * dV2ijy * si2_sj2
            deralat[1, 2] = deralat[1, 2] - s2_r[nj] * dV2ijy * si3_sj3
            deralat[2, 0] = deralat[2, 0] - s2_r[nj] * dV2ijz * si1_sj1
            deralat[2, 1] = deralat[2, 1] - s2_r[nj] * dV2ijz * si2_sj2
            deralat[2, 2] = deralat[2, 2] - s2_r[nj] * dV2ijz * si3_sj3

            #   !              if(fixZ .eq. 0) then

            #   !  --- LEVEL 3: LOOP FOR PAIR COORDINATION FORCES ---

            dV2dZ = - dp * s2_t0[nj]
            # do nl = 1, nz - 1
            #   sz_sum(nl) = sz_sum(nl) + dV2dZ
            # end do
            sz_sum[0:nz - 1] = sz_sum[0:nz - 1] + dV2dZ

            #   !              end if
            #   !  fixZ
        # end do
        # !Commented out by Stefan
        # !            if(fixZ .ne. 0) then
        # !              winv = Qort*dexp(-muhalf*Z)
        # !              dwinv = 0.0
        # !              temp0 = dexp(-u4*Z)
        # !              tau = u1+u2*temp0*(u3-temp0)
        # !              dtau = 0.0
        # !            else
        # !----------------------------------

        # !   COORDINATION-DEPENDENCE OF THREE-BODY INTERACTION

        winv = Qort * np.exp(-muhalf * Z)
        # !   inverse width of angular function
        dwinv = -muhalf*winv
        # !   its derivative
        temp0 = np.exp(-u4 * Z)
        tau = u1 + u2*temp0*(u3 - temp0)
        # !   -cosine of angular minimum
        dtau = u5*temp0*(2.0*temp0 - u3)
        # !   its derivative
        # !            end if

        # !  --- LEVEL 2: FIRST LOOP FOR THREE-BODY INTERACTIONS ---

        # do nj = 1, n3 - 2
        for nj in range(n3 - 2):

            j = num3[nj]

            # !  --- LEVEL 3: SECOND LOOP FOR THREE-BODY INTERACTIONS ---

            # do nk = nj + 1, n3 - 1
            for nk in range(nj + 1, n3):
                k = num3[nk]
                # !   angular function h(l,Z)
                lcos = s3_dx[nj] * s3_dx[nk] + s3_dy[nj] * s3_dy[nk] + s3_dz[nj] * s3_dz[nk]
                x = (lcos + tau) * winv
                temp0 = np.exp(-x*x)

                H = par_lam * (1.0 - temp0 + par_eta * x * x)
                dHdx = 2.0 * par_lam * x * (temp0 + par_eta)

                dhdl = dHdx*winv

                # !   three-body energy

                temp1 = s3_g[nj] * s3_g[nk]
                ener_iat = ener_iat + temp1*H

                # !   (-) radial force on atom j

                dV3rij = s3_dg[nj] * s3_g[nk] * H
                dV3rijx = dV3rij*s3_dx[nj]
                dV3rijy = dV3rij*s3_dy[nj]
                dV3rijz = dV3rij*s3_dz[nj]
                fjx = dV3rijx
                fjy = dV3rijy
                fjz = dV3rijz

                # !   (-) radial force on atom k

                dV3rik = s3_g[nj]*s3_dg[nk]*H
                dV3rikx = dV3rik*s3_dx[nk]
                dV3riky = dV3rik*s3_dy[nk]
                dV3rikz = dV3rik*s3_dz[nk]
                fkx = dV3rikx
                fky = dV3riky
                fkz = dV3rikz

                # !   (-) angular force on j

                dV3l = temp1*dhdl
                dV3ljx = dV3l * (s3_dx[nk] - lcos*s3_dx[nj])*s3_rinv[nj]
                dV3ljy = dV3l * (s3_dy[nk] - lcos*s3_dy[nj])*s3_rinv[nj]
                dV3ljz = dV3l * (s3_dz[nk] - lcos*s3_dz[nj])*s3_rinv[nj]
                fjx = fjx + dV3ljx
                fjy = fjy + dV3ljy
                fjz = fjz + dV3ljz

                # !   (-) angular force on k

                dV3lkx = dV3l*(s3_dx[nj] - lcos*s3_dx[nk])*s3_rinv[nk]
                dV3lky = dV3l*(s3_dy[nj] - lcos*s3_dy[nk])*s3_rinv[nk]
                dV3lkz = dV3l*(s3_dz[nj] - lcos*s3_dz[nk])*s3_rinv[nk]
                fkx = fkx + dV3lkx
                fky = fky + dV3lky
                fkz = fkz + dV3lkz

                # !   apply radial + angular forces to i, j, k

                fxyz[0, j] = fxyz[0, j] - fjx
                fxyz[1, j] = fxyz[1, j] - fjy
                fxyz[2, j] = fxyz[2, j] - fjz
                fxyz[0, k] = fxyz[0, k] - fkx
                fxyz[1, k] = fxyz[1, k] - fky
                fxyz[2, k] = fxyz[2, k] - fkz
                fxyz[0, i] = fxyz[0, i] + fjx + fkx
                fxyz[1, i] = fxyz[1, i] + fjy + fky
                fxyz[2, i] = fxyz[2, i] + fjz + fkz

                # !   dV3/dR contributions to virial

                virial = virial - s3_r[nj] * (fjx*s3_dx[nj] + fjy*s3_dy[nj] + fjz*s3_dz[nj])
                virial = virial - s3_r[nk] * (fkx*s3_dx[nk] + fky*s3_dy[nk] + fkz*s3_dz[nk])
                virial_xyz[0] = virial_xyz[0] - s3_r[nj] * (fjx * s3_dx[nj])
                virial_xyz[1] = virial_xyz[1] - s3_r[nj] * (fjy * s3_dy[nj])
                virial_xyz[2] = virial_xyz[2] - s3_r[nj] * (fjz * s3_dz[nj])
                virial_xyz[0] = virial_xyz[0] - s3_r[nk] * (fkx * s3_dx[nk])
                virial_xyz[1] = virial_xyz[1] - s3_r[nk] * (fky * s3_dy[nk])
                virial_xyz[2] = virial_xyz[2] - s3_r[nk] * (fkz * s3_dz[nk])

                # !Cell gradient part
                # !My own implementation
                si1_sj1 = alatinv[0, 0] * s3_dx[nj] + alatinv[0, 1] * s3_dy[nj] + alatinv[0, 2] * s3_dz[nj]
                si2_sj2 = alatinv[1, 0] * s3_dx[nj] + alatinv[1, 1] * s3_dy[nj] + alatinv[1, 2] * s3_dz[nj]
                si3_sj3 = alatinv[2, 0] * s3_dx[nj] + alatinv[2, 1] * s3_dy[nj] + alatinv[2, 2] * s3_dz[nj]
                deralat[0, 0] = deralat[0, 0] - s3_r[nj] * fjx * si1_sj1
                deralat[0, 1] = deralat[0, 1] - s3_r[nj] * fjx * si2_sj2
                deralat[0, 2] = deralat[0, 2] - s3_r[nj] * fjx * si3_sj3
                deralat[1, 0] = deralat[1, 0] - s3_r[nj] * fjy * si1_sj1
                deralat[1, 1] = deralat[1, 1] - s3_r[nj] * fjy * si2_sj2
                deralat[1, 2] = deralat[1, 2] - s3_r[nj] * fjy * si3_sj3
                deralat[2, 0] = deralat[2, 0] - s3_r[nj] * fjz * si1_sj1
                deralat[2, 1] = deralat[2, 1] - s3_r[nj] * fjz * si2_sj2
                deralat[2, 2] = deralat[2, 2] - s3_r[nj] * fjz * si3_sj3

                # !Cell gradient part
                # !My own implementation
                si1_sj1 = alatinv[0, 0]*s3_dx[nk] + alatinv[0, 1] * s3_dy[nk] + alatinv[0, 2] * s3_dz[nk]
                si2_sj2 = alatinv[1, 0]*s3_dx[nk] + alatinv[1, 1] * s3_dy[nk] + alatinv[1, 2] * s3_dz[nk]
                si3_sj3 = alatinv[2, 0]*s3_dx[nk] + alatinv[2, 1] * s3_dy[nk] + alatinv[2, 2] * s3_dz[nk]
                deralat[0, 0] = deralat[0, 0] - s3_r[nk] * fkx * si1_sj1
                deralat[0, 1] = deralat[0, 1] - s3_r[nk] * fkx * si2_sj2
                deralat[0, 2] = deralat[0, 2] - s3_r[nk] * fkx * si3_sj3
                deralat[1, 0] = deralat[1, 0] - s3_r[nk] * fky * si1_sj1
                deralat[1, 1] = deralat[1, 1] - s3_r[nk] * fky * si2_sj2
                deralat[1, 2] = deralat[1, 2] - s3_r[nk] * fky * si3_sj3
                deralat[2, 0] = deralat[2, 0] - s3_r[nk] * fkz * si1_sj1
                deralat[2, 1] = deralat[2, 1] - s3_r[nk] * fkz * si2_sj2
                deralat[2, 2] = deralat[2, 2] - s3_r[nk] * fkz * si3_sj3

                # !                if(fixZ .eq. 0) then

                # !   prefactor for 4-body forces from coordination
                dxdZ = dwinv*(lcos + tau) + winv*dtau
                dV3dZ = temp1*dHdx*dxdZ

                # !  --- LEVEL 4: LOOP FOR THREE-BODY COORDINATION FORCES ---

                # do nl = 1, nz - 1
                #   sz_sum(nl) = sz_sum(nl) + dV3dZ
                # end do
                sz_sum[0:nz - 1] = sz_sum[0:nz - 1] + dV3dZ
                #   !                end if
            # end do
        # end do

        # !            if(fixZ .eq. 0) then

        # !  --- LEVEL 2: LOOP TO APPLY COORDINATION FORCES ---

        # do nl = 1, nz - 1
        for nl in range(nz - 1):

            dEdrl = sz_sum[nl] * sz_df[nl]
            dEdrlx = dEdrl*sz_dx[nl]
            dEdrly = dEdrl*sz_dy[nl]
            dEdrlz = dEdrl*sz_dz[nl]
            fxyz[0, i] = fxyz[0, i] + dEdrlx
            fxyz[1, i] = fxyz[1, i] + dEdrly
            fxyz[2, i] = fxyz[2, i] + dEdrlz
            l = numz[nl]
            fxyz[0, l] = fxyz[0, l] - dEdrlx
            fxyz[1, l] = fxyz[1, l] - dEdrly
            fxyz[2, l] = fxyz[2, l] - dEdrlz

            # !   dE/dZ*dZ/dr contribution to virial

            virial = virial - sz_r[nl] * (dEdrlx*sz_dx[nl] + dEdrly*sz_dy[nl] + dEdrlz*sz_dz[nl])
            virial_xyz[0] = virial_xyz[0] - sz_r[nl]*(dEdrlx * sz_dx[nl])
            virial_xyz[1] = virial_xyz[1] - sz_r[nl]*(dEdrly * sz_dy[nl])
            virial_xyz[2] = virial_xyz[2] - sz_r[nl]*(dEdrlz * sz_dz[nl])

            # !Cell gradient part
            # !My own implementation
            si1_sj1 = alatinv[0, 0]*sz_dx[nl] + alatinv[0, 1]*sz_dy[nl] + alatinv[0, 2]*sz_dz[nl]
            si2_sj2 = alatinv[1, 0]*sz_dx[nl] + alatinv[1, 1]*sz_dy[nl] + alatinv[1, 2]*sz_dz[nl]
            si3_sj3 = alatinv[2, 0]*sz_dx[nl] + alatinv[2, 1]*sz_dy[nl] + alatinv[2, 2]*sz_dz[nl]
            deralat[0, 0] = deralat[0, 0] - sz_r[nl] * dEdrlx * si1_sj1
            deralat[0, 1] = deralat[0, 1] - sz_r[nl] * dEdrlx * si2_sj2
            deralat[0, 2] = deralat[0, 2] - sz_r[nl] * dEdrlx * si3_sj3
            deralat[1, 0] = deralat[1, 0] - sz_r[nl] * dEdrly * si1_sj1
            deralat[1, 1] = deralat[1, 1] - sz_r[nl] * dEdrly * si2_sj2
            deralat[1, 2] = deralat[1, 2] - sz_r[nl] * dEdrly * si3_sj3
            deralat[2, 0] = deralat[2, 0] - sz_r[nl] * dEdrlz * si1_sj1
            deralat[2, 1] = deralat[2, 1] - sz_r[nl] * dEdrlz * si2_sj2
            deralat[2, 2] = deralat[2, 2] - sz_r[nl] * dEdrlz * si3_sj3

        # end do

        # !           end if
        ener = ener + ener_iat
        ener2 = ener2 + ener_iat**2
    # end do

    #   !        call getvol(alat,vol)
    #   !          if(vol.lt.0.d0) then
    #   !            write(77,*) nat
    #   !            write(77,*) alat(:,1)
    #   !            write(77,*) alat(:,2)
    #   !            write(77,*) alat(:,3)
    #   !            do i=1,nat
    #   !              write(77,*)  rxyz(:,i),"Si"
    #   !            enddo
    #   !
    #   !          endif

    etot = ener

    etot = etot/Ha_eV
    # do iat = 1, nat
    #   fxyz(1, iat) = -fxyz(1, iat)/Ha_eV*Bohr_Ang
    #   fxyz(2, iat) = -fxyz(2, iat)/Ha_eV*Bohr_Ang
    #   fxyz(3, iat) = -fxyz(3, iat)/Ha_eV*Bohr_Ang
    # end do
    fxyz = -fxyz / Ha_eV * Bohr_Ang

    # do j = 1, 3
    # do i = 1, 3
    #   deralat(i, j) = deralat(i, j)/Ha_eV*Bohr_Ang
    # end do
    # end do
    deralat = deralat / Ha_eV * Bohr_Ang


    # ! Stress
    # ! Transform the forces on the lattice vectors into the stress tensore according to the paper Tomas Bucko and Jurg Hafner
    # !        do i=1,3
    # !           tmplat(:,i)=latvec(i,:)
    # !        enddo
    # vol = (alat0(1, 1)*alat0(2, 2)*alat0(3, 3) - alat0(1, 1)*alat0(2, 3)*alat0(3, 2) - &
    #   alat0(1, 2)*alat0(2, 1)*alat0(3, 3) + alat0(1, 2)*alat0(2, 3)*alat0(3, 1) + &
    #   alat0(1, 3)*alat0(2, 1)*alat0(3, 2) - alat0(1, 3)*alat0(2, 2)*alat0(3, 1))

    vol = np.abs(np.linalg.det(alat0))

    # stress=-matmul(deralat,transpose(alat0))/vol

    stress = - deralat @ alat0.T / vol
    return etot, fxyz, deralat, stress

#   !stress = stress / Ha_eV *Bohr_Ang**3
#   !        strten(1) = stress(1,1)
#   !        strten(2) = stress(2,2)
#   !        strten(3) = stress(3,3)
#   !        strten(6) = stress(2,1)
#   !        strten(5) = stress(3,1)
#   !        strten(4) = stress(3,2)
#   !!This is not very clear yet...
#   !        strten=strten/Ha_eV*Bohr_Ang**3
# end subroutine energyandforces_bazant


# subroutine invertalat(alat, alatinv)

def test():
    from ase.io import read
    import bazant_fortran as bf

    import time

    atoms = read('test.ascii')
    nat = len(atoms)

    import ase.units as units

    alat = atoms.get_cell().T
    positions = atoms.get_positions().T

    # convert to atomic units
    alat /= units.Bohr
    positions /= units.Bohr

    t1 = time.time()
    etot, fxyz, deralat = bf.energyandforces_bazant(alat, positions)
    t2 = time.time()

    print(etot)

    t3 = time.time()
    etot, fxyz, deralat, stress = energyandforces_bazant(nat, alat, positions)
    t4 = time.time()

    print(etot)

    t5 = time.time()
    etot, fxyz, deralat, stress = energyandforces_bazant(nat, alat, positions)
    t6 = time.time()

    print(etot)

    print('Fortran time:', t2 - t1)
    print('Python time:', t4 - t3)
    print('Python time (compiled):', t6 - t5)
    print('Speedup:', (t4 - t3) / (t2 - t1), (t6 - t5) / (t2 - t1))

if __name__ == '__main__':
    test()