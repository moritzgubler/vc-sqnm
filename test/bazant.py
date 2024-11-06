import numpy as np
try:
    from numba import njit
except ImportError:
    def njit(f):
        return f

@njit()
def nnlist(nnbrx, alat, cutoff, rxyz):
    """
    Calculate the neighbor list for a given set of atoms and lattice vectors

    Parameters

    nnbrx : int Maximum number of neighbors
    alat : np.ndarray(3, 3) Lattice vectors, each row is one vector
    cutoff : float Cutoff radius
    rxyz : np.ndarray(nat, 3) Atomic positions
    """

    nat = rxyz.shape[0]

    # if cell volume is zero, we are in a non-periodic system
    if np.abs(np.linalg.det(alat)) < 1e-10:
        ixyzmax = 0
    else:  #! periodic boundary conditions
        eigval, eigvec = np.linalg.eigh(alat.T @ alat)
        ixyzmax = int(np.sqrt(1.0 / eigval[0]) * cutoff) + 1
        print(np.sqrt(1.0 / eigval[0]) * cutoff, np.sqrt(1.0 / eigval[1]) * cutoff, np.sqrt(1.0 / eigval[2]) * cutoff)
    cutoff2 = cutoff * cutoff

    ind = 0
    lsta = np.empty((2, nat), dtype=np.int32)
    lstb = np.empty(nnbrx * nat, dtype=np.int32)
    rel = np.empty((nat * nnbrx, 5))
    for iat in range(nat):
        lsta[0, iat] = ind
        for jat in range(nat):
            for ix in range(-ixyzmax, ixyzmax + 1):
               for iy in range(-ixyzmax, ixyzmax + 1):
                   for iz in range(-ixyzmax, ixyzmax + 1):
                        xj = rxyz[jat, 0] + ix * alat[0, 0] + iy*alat[1, 0] + iz*alat[2, 0]
                        yj = rxyz[jat, 1] + ix * alat[0, 1] + iy*alat[1, 1] + iz*alat[2, 1]
                        zj = rxyz[jat, 2] + ix * alat[0, 2] + iy*alat[1, 2] + iz*alat[2, 2]
                        relx = xj - rxyz[iat, 0]
                        rely = yj - rxyz[iat, 1]
                        relz = zj - rxyz[iat, 2]
                        dist2 = relx * relx + rely * rely + relz * relz

                        if dist2 > 1.0e-20 and dist2 <= cutoff2:
                            if (ind >= nnbrx*nat):
                                raise ValueError('enlarge nnbrx')
                            lstb[ind] = jat
                            tt = np.sqrt(dist2)
                            ttinv = 1.0 / tt
                            rel[ind, 0] = relx*ttinv
                            rel[ind, 1] = rely*ttinv
                            rel[ind, 2] = relz*ttinv
                            rel[ind, 3] = tt
                            rel[ind, 4] = ttinv
                            ind = ind + 1
        lsta[1, iat] = ind - 1
    return lsta, lstb, rel

@njit()
def energyandforces_bazant(alat0, rxyz0):
    """
    Calculate the energy and forces for a given set of atoms and lattice vectors using the Bazant EDIP potential.
    Internally things are computed in Angstroem and eV, The input (rxyz0,alat0) is however in bohr and
    the output energy in Hartree, fxyz and deralat  in Hartree/Bohr

    Parameters

    alat0 : np.ndarray(3, 3) Lattice vectors, each row is one vector, units in bohr
    rxyz0 : np.ndarray(nat, 3) Atomic positions, units in bohr

    Returns
    ener : float Energy in Hartree
    fxyz : np.ndarray(nat, 3) Forces in Hartree/Bohr
    stress: np.ndarray(3, 3) Stress tensor in Hartree/Bohr^3

    """


    # Connection between these parameters and those given in the paper,
    # Justo et al., Phys. Rev. B 58, 2539 (1998):
    #
    # A((B/r)**rh-palp*exp(-bet*Z*Z)) = A'((B'/r)**rh-exp(-bet*Z*Z))
    #
    # so in the paper (')
    # A' = A*palp
    # B' = B * palp**(-1/rh)
    # eta = detla/Qo
    #
    # Non-adjustable parameters for tau(Z) from Ismail & Kaxiras, 1993,
    # also in Bazant, Kaxiras, Justo, PRB (1997):
    #
    # u1 = -0.165799;
    # u2 = 32.557;
    # u3 = 0.286198;
    # u4 = 0.66;


    #  integer i, j, k, l, n, n2, n3, nz, iat
    i = 0
    j = 0
    k = 0
    l = 0
    n = 0
    n2 = 0
    n3 = 0
    nz = 0
    iat = 0
    #  real*8 :: dx, dy, dz, r, asqr
    dx = 0.0
    dy = 0.0
    dz = 0.0
    r = 0.0
    asqr = 0.0

    #  real*8 :: rinv, rmainv, xinv, xinv3, den, Z, fZ
    rinv = 0.0
    rmainv = 0.0
    xinv = 0.0
    xinv3 = 0.0
    den = 0.0
    Z = 0.0
    fZ = 0.0
    #  real*8 :: dV2j, dV2ijx, dV2ijy, dV2ijz, pZ, dp
    dV2j = 0.0
    dV2ijx = 0.0
    dV2ijy = 0.0
    dV2ijz = 0.0
    pZ = 0.0
    dp = 0.0
    #  real*8 :: temp0, temp1
    temp0 = 0.0
    temp1 = 0.0
    #  real*8 :: Qort, muhalf, u5
    Qort = 0.0
    muhalf = 0.0
    u5 = 0.0
    #  real*8 :: rmbinv, winv, dwinv, tau, dtau, lcos, x, H, dHdx, dhdl
    rmbinv = 0.0
    winv = 0.0
    dwinv = 0.0
    tau = 0.0
    dtau = 0.0
    lcos = 0.0
    x = 0.0
    H = 0.0
    dHdx = 0.0
    dhdl = 0.0
    #  real*8 :: dV3rij, dV3rijx, dV3rijy, dV3rijz
    dV3rij = 0.0
    dV3rijx = 0.0
    dV3rijy = 0.0
    dV3rijz = 0.0
    #  real*8 :: dV3rik, dV3rikx, dV3riky, dV3rikz
    dV3rik = 0.0
    dV3rikx = 0.0
    dV3riky = 0.0
    dV3rikz = 0.0
    #  real*8 :: dV3l, dV3ljx, dV3ljy, dV3ljz, dV3lkx, dV3lky, dV3lkz
    dV3l = 0.0
    dV3ljx = 0.0
    dV3ljy = 0.0
    dV3ljz = 0.0
    dV3lkx = 0.0
    dV3lky = 0.0
    dV3lkz = 0.0
    #  real*8 :: dV2dZ, dxdZ, dV3dZ
    dV2dZ = 0.0
    dxdZ = 0.0
    dV3dZ = 0.0
    #  real*8 :: dEdrl, dEdrlx, dEdrly, dEdrlz
    dEdrl = 0.0
    dEdrlx = 0.0
    dEdrly = 0.0
    dEdrlz = 0.0
    #  real*8 :: bmc, cmbinv
    bmc = 0.0
    cmbinv = 0.0
    #  real*8 :: fjx, fjy, fjz, fkx, fky, fkz
    fjx = 0.0
    fjy = 0.0
    fjz = 0.0
    fkx = 0.0
    fky = 0.0
    fkz = 0.0


    nnbrx = 50
    nat = rxyz0.shape[0]

    # work arrays
    s2_t0 =   np.zeros(nnbrx)
    s2_t1 =   np.zeros(nnbrx)
    s2_t2 =   np.zeros(nnbrx)
    s2_t3 =   np.zeros(nnbrx)
    s2_dx =   np.zeros(nnbrx)
    s2_dy =   np.zeros(nnbrx)
    s2_dz =   np.zeros(nnbrx)
    s2_r =    np.zeros(nnbrx)
    s3_g =    np.zeros(nnbrx)
    s3_dg =   np.zeros(nnbrx)
    s3_rinv = np.zeros(nnbrx)
    s3_dx =   np.zeros(nnbrx)
    s3_dy =   np.zeros(nnbrx)
    s3_dz =   np.zeros(nnbrx)
    s3_r =    np.zeros(nnbrx)
    sz_df =   np.zeros(nnbrx)
    sz_sum =  np.zeros(nnbrx)
    sz_dx =   np.zeros(nnbrx)
    sz_dy =   np.zeros(nnbrx)
    sz_dz =   np.zeros(nnbrx)
    sz_r =    np.zeros(nnbrx)

    num2 = np.zeros(nnbrx, dtype=np.int32)
    num3 = np.zeros(nnbrx, dtype=np.int32)
    numz = np.zeros(nnbrx, dtype=np.int32)

    # EDIP parameters
    # taken from Justo et al., Phys. Rev. B 58, 2539 (1998).

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
    alat = alat0 * Bohr_Ang
    rxyz = rxyz0 * Bohr_Ang
    
    alatinv = np.linalg.inv(alat)

    fracpos = rxyz @ alatinv
    fracpos = fracpos - np.floor(fracpos)
    rxyz = fracpos @ alat

    # call nnlist(nat, nnbrx, alat, cutoff, rxyz, lsta, lstb, rel)
    lsta, lstb, rel = nnlist(nnbrx, alat, cutoff, rxyz)

    fxyz = np.zeros((nat, 3))
    ener = 0.0
    ener2 = 0.0

    virial = 0.0
    virial_xyz = np.zeros(3)
    deralat = np.zeros((3, 3))

    # COMBINE COEFFICIENTS
    asqr = par_a * par_a
    Qort = np.sqrt(par_Qo)
    muhalf = par_mu*0.50
    u5 = u2*u4
    bmc = par_b - par_c
    cmbinv = 1.0 / (par_c - par_b)

    #   --- LEVEL 1: OUTER LOOP OVER ATOMS ---

    for i in range(nat):

        #   RESET COORDINATION AND NEIGHBOR NUMBERS
        ener_iat = 0.0
        Z = 0.0
        n2 = 0 # used to be 1
        n3 = 0 # used to be 1
        nz = 0 # used to be 1

        #  --- LEVEL 2: LOOP PREPASS OVER PAIRS ---

        # do n = lsta(1, i), lsta(2, i)
        for n in range(lsta[0, i], lsta[1, i] + 1):
            j = lstb[n]

            #   !   PARTS OF TWO-BODY INTERACTION r<par_a

            num2[n2] = j
            dx = - rel[n, 0]
            dy = - rel[n, 1]
            dz = - rel[n, 2]
            r  =   rel[n, 3]
            rinv = rel[n, 4]

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
                raise ValueError('enlarge nnbrx')
            #   RADIAL PARTS OF THREE-BODY INTERACTION r<par_b
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

                if (n3 >= nnbrx):
                    raise ValueError('enlarge nnbrx')

                #   COORDINATION AND NEIGHBOR FUNCTION par_c<r<par_b

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
                        sz_dx[nz] = dx
                        sz_dy[nz] = dy
                        sz_dz[nz] = dz
                        sz_r[nz] = r
                        nz = nz + 1
                        if (nz > nnbrx):
                            raise ValueError('enlarge nnbrx')
        #   ZERO ACCUMULATION ARRAY FOR ENVIRONMENT FORCES

        sz_sum[0:nz - 1] = 0.0

        #   ENVIRONMENT-DEPENDENCE OF PAIR INTERACTION
        temp0 = par_bet * Z
        pZ = par_palp * np.exp(-temp0 * Z)
        #   bond order
        dp = -2.0*temp0*pZ
        #   derivative of bond order


        #  --- LEVEL 2: LOOP FOR PAIR INTERACTIONS ---

        for nj in range(n2):
            temp0 = s2_t1[nj] - pZ

            #   two-body energy V2(rij,Z)
            ener_iat = ener_iat + temp0*s2_t0[nj]

            # two-body forces
            dV2j = -s2_t0[nj]*(s2_t1[nj]*s2_t2[nj] + temp0*s2_t3[nj])
            #   dV2/dr
            dV2ijx = dV2j*s2_dx[nj]
            dV2ijy = dV2j*s2_dy[nj]
            dV2ijz = dV2j*s2_dz[nj]
            fxyz[i, 0] = fxyz[i, 0] + dV2ijx
            fxyz[i, 1] = fxyz[i, 1] + dV2ijy
            fxyz[i, 2] = fxyz[i, 2] + dV2ijz
            j = num2[nj]
            fxyz[j, 0] = fxyz[j, 0] - dV2ijx
            fxyz[j, 1] = fxyz[j, 1] - dV2ijy
            fxyz[j, 2] = fxyz[j, 2] - dV2ijz

            #  dV2/dr contribution to virial
            virial_xyz[0] = virial_xyz[0] - s2_r[nj] * (dV2ijx * s2_dx[nj])
            virial_xyz[1] = virial_xyz[1] - s2_r[nj] * (dV2ijy * s2_dy[nj])
            virial_xyz[2] = virial_xyz[2] - s2_r[nj] * (dV2ijz * s2_dz[nj])
            virial = virial - s2_r[nj] * (dV2ijx*s2_dx[nj] + dV2ijy*s2_dy[nj] + dV2ijz*s2_dz[nj])

            #   Cell gradient part
            #   My own implementation
            si1_sj1 = alatinv[0, 0]*s2_dx[nj] + alatinv[1, 0]*s2_dy[nj] + alatinv[2, 0]*s2_dz[nj]
            si2_sj2 = alatinv[0, 1]*s2_dx[nj] + alatinv[1, 1]*s2_dy[nj] + alatinv[2, 1]*s2_dz[nj]
            si3_sj3 = alatinv[0, 2]*s2_dx[nj] + alatinv[1, 2]*s2_dy[nj] + alatinv[2, 2]*s2_dz[nj]
            deralat[0, 0] = deralat[0, 0] - s2_r[nj] * dV2ijx * si1_sj1
            deralat[1, 0] = deralat[1, 0] - s2_r[nj] * dV2ijx * si2_sj2
            deralat[2, 0] = deralat[2, 0] - s2_r[nj] * dV2ijx * si3_sj3
            deralat[0, 1] = deralat[0, 1] - s2_r[nj] * dV2ijy * si1_sj1
            deralat[1, 1] = deralat[1, 1] - s2_r[nj] * dV2ijy * si2_sj2
            deralat[2, 1] = deralat[2, 1] - s2_r[nj] * dV2ijy * si3_sj3
            deralat[0, 2] = deralat[0, 2] - s2_r[nj] * dV2ijz * si1_sj1
            deralat[1, 2] = deralat[1, 2] - s2_r[nj] * dV2ijz * si2_sj2
            deralat[2, 2] = deralat[2, 2] - s2_r[nj] * dV2ijz * si3_sj3

            #     --- LEVEL 3: LOOP FOR PAIR COORDINATION FORCES ---

            dV2dZ = - dp * s2_t0[nj]
            sz_sum[0:nz - 1] = sz_sum[0:nz - 1] + dV2dZ
        #   COORDINATION-DEPENDENCE OF THREE-BODY INTERACTION

        winv = Qort * np.exp(-muhalf * Z)
        #   inverse width of angular function
        dwinv = -muhalf*winv
        #   its derivative
        temp0 = np.exp(-u4 * Z)
        tau = u1 + u2*temp0*(u3 - temp0)
        #   -cosine of angular minimum
        dtau = u5*temp0*(2.0*temp0 - u3)
        #   its derivative
        #            end if

        #  --- LEVEL 2: FIRST LOOP FOR THREE-BODY INTERACTIONS ---

        for nj in range(n3 - 1):
            j = num3[nj]

            #  --- LEVEL 3: SECOND LOOP FOR THREE-BODY INTERACTIONS ---

            for nk in range(nj + 1, n3):
                k = num3[nk]
                #   angular function h(l,Z)
                lcos = s3_dx[nj] * s3_dx[nk] + s3_dy[nj] * s3_dy[nk] + s3_dz[nj] * s3_dz[nk]
                x = (lcos + tau) * winv
                temp0 = np.exp(-x*x)

                H = par_lam * (1.0 - temp0 + par_eta * x * x)
                dHdx = 2.0 * par_lam * x * (temp0 + par_eta)

                dhdl = dHdx*winv

                #   three-body energy
                temp1 = s3_g[nj] * s3_g[nk]
                ener_iat = ener_iat + temp1*H

                #   (-) radial force on atom j
                dV3rij = s3_dg[nj] * s3_g[nk] * H
                dV3rijx = dV3rij*s3_dx[nj]
                dV3rijy = dV3rij*s3_dy[nj]
                dV3rijz = dV3rij*s3_dz[nj]
                fjx = dV3rijx
                fjy = dV3rijy
                fjz = dV3rijz

                #   (-) radial force on atom k
                dV3rik = s3_g[nj]*s3_dg[nk]*H
                dV3rikx = dV3rik*s3_dx[nk]
                dV3riky = dV3rik*s3_dy[nk]
                dV3rikz = dV3rik*s3_dz[nk]
                fkx = dV3rikx
                fky = dV3riky
                fkz = dV3rikz

                #    (-) angular force on j
                dV3l = temp1*dhdl
                dV3ljx = dV3l * (s3_dx[nk] - lcos*s3_dx[nj])*s3_rinv[nj]
                dV3ljy = dV3l * (s3_dy[nk] - lcos*s3_dy[nj])*s3_rinv[nj]
                dV3ljz = dV3l * (s3_dz[nk] - lcos*s3_dz[nj])*s3_rinv[nj]
                fjx = fjx + dV3ljx
                fjy = fjy + dV3ljy
                fjz = fjz + dV3ljz

                #   (-) angular force on k
                dV3lkx = dV3l*(s3_dx[nj] - lcos*s3_dx[nk])*s3_rinv[nk]
                dV3lky = dV3l*(s3_dy[nj] - lcos*s3_dy[nk])*s3_rinv[nk]
                dV3lkz = dV3l*(s3_dz[nj] - lcos*s3_dz[nk])*s3_rinv[nk]
                fkx = fkx + dV3lkx
                fky = fky + dV3lky
                fkz = fkz + dV3lkz

                #   apply radial + angular forces to i, j, k
                fxyz[j, 0] = fxyz[j, 0] - fjx
                fxyz[j, 1] = fxyz[j, 1] - fjy
                fxyz[j, 2] = fxyz[j, 2] - fjz
                fxyz[k, 0] = fxyz[k, 0] - fkx
                fxyz[k, 1] = fxyz[k, 1] - fky
                fxyz[k, 2] = fxyz[k, 2] - fkz
                fxyz[i, 0] = fxyz[i, 0] + fjx + fkx
                fxyz[i, 1] = fxyz[i, 1] + fjy + fky
                fxyz[i, 2] = fxyz[i, 2] + fjz + fkz

                #   dV3/dR contributions to virial
                virial = virial - s3_r[nj] * (fjx*s3_dx[nj] + fjy*s3_dy[nj] + fjz*s3_dz[nj])
                virial = virial - s3_r[nk] * (fkx*s3_dx[nk] + fky*s3_dy[nk] + fkz*s3_dz[nk])
                virial_xyz[0] = virial_xyz[0] - s3_r[nj] * (fjx * s3_dx[nj])
                virial_xyz[1] = virial_xyz[1] - s3_r[nj] * (fjy * s3_dy[nj])
                virial_xyz[2] = virial_xyz[2] - s3_r[nj] * (fjz * s3_dz[nj])
                virial_xyz[0] = virial_xyz[0] - s3_r[nk] * (fkx * s3_dx[nk])
                virial_xyz[1] = virial_xyz[1] - s3_r[nk] * (fky * s3_dy[nk])
                virial_xyz[2] = virial_xyz[2] - s3_r[nk] * (fkz * s3_dz[nk])

                # Cell gradient part
                # My own implementation
                si1_sj1 = alatinv[0, 0] * s3_dx[nj] + alatinv[1, 0] * s3_dy[nj] + alatinv[2, 0] * s3_dz[nj]
                si2_sj2 = alatinv[0, 1] * s3_dx[nj] + alatinv[1, 1] * s3_dy[nj] + alatinv[2, 1] * s3_dz[nj]
                si3_sj3 = alatinv[0, 2] * s3_dx[nj] + alatinv[1, 2] * s3_dy[nj] + alatinv[2, 2] * s3_dz[nj]
                deralat[0, 0] = deralat[0, 0] - s3_r[nj] * fjx * si1_sj1
                deralat[1, 0] = deralat[1, 0] - s3_r[nj] * fjx * si2_sj2
                deralat[2, 0] = deralat[2, 0] - s3_r[nj] * fjx * si3_sj3
                deralat[0, 1] = deralat[0, 1] - s3_r[nj] * fjy * si1_sj1
                deralat[1, 1] = deralat[1, 1] - s3_r[nj] * fjy * si2_sj2
                deralat[2, 1] = deralat[2, 1] - s3_r[nj] * fjy * si3_sj3
                deralat[0, 2] = deralat[0, 2] - s3_r[nj] * fjz * si1_sj1
                deralat[1, 2] = deralat[1, 2] - s3_r[nj] * fjz * si2_sj2
                deralat[2, 2] = deralat[2, 2] - s3_r[nj] * fjz * si3_sj3

                # Cell gradient part
                # My own implementation
                si1_sj1 = alatinv[0, 0]*s3_dx[nk] + alatinv[1, 0] * s3_dy[nk] + alatinv[2, 0] * s3_dz[nk]
                si2_sj2 = alatinv[0, 1]*s3_dx[nk] + alatinv[1, 1] * s3_dy[nk] + alatinv[2, 1] * s3_dz[nk]
                si3_sj3 = alatinv[0, 2]*s3_dx[nk] + alatinv[1, 2] * s3_dy[nk] + alatinv[2, 2] * s3_dz[nk]
                deralat[0, 0] = deralat[0, 0] - s3_r[nk] * fkx * si1_sj1
                deralat[1, 0] = deralat[1, 0] - s3_r[nk] * fkx * si2_sj2
                deralat[2, 0] = deralat[2, 0] - s3_r[nk] * fkx * si3_sj3
                deralat[0, 1] = deralat[0, 1] - s3_r[nk] * fky * si1_sj1
                deralat[1, 1] = deralat[1, 1] - s3_r[nk] * fky * si2_sj2
                deralat[2, 1] = deralat[2, 1] - s3_r[nk] * fky * si3_sj3
                deralat[0, 2] = deralat[0, 2] - s3_r[nk] * fkz * si1_sj1
                deralat[1, 2] = deralat[1, 2] - s3_r[nk] * fkz * si2_sj2
                deralat[2, 2] = deralat[2, 2] - s3_r[nk] * fkz * si3_sj3

                #   prefactor for 4-body forces from coordination
                dxdZ = dwinv*(lcos + tau) + winv*dtau
                dV3dZ = temp1*dHdx*dxdZ

                #   --- LEVEL 4: LOOP FOR THREE-BODY COORDINATION FORCES ---
                sz_sum[0:nz - 1] = sz_sum[0:nz - 1] + dV3dZ
                #                 end if

        #           if(fixZ .eq. 0) then

        #  --- LEVEL 2: LOOP TO APPLY COORDINATION FORCES ---

        for nl in range(nz - 1):

            dEdrl = sz_sum[nl] * sz_df[nl]
            dEdrlx = dEdrl*sz_dx[nl]
            dEdrly = dEdrl*sz_dy[nl]
            dEdrlz = dEdrl*sz_dz[nl]
            fxyz[i, 0] = fxyz[i, 0] + dEdrlx
            fxyz[i, 1] = fxyz[i, 1] + dEdrly
            fxyz[i, 2] = fxyz[i, 2] + dEdrlz
            l = numz[nl]
            fxyz[l, 0] = fxyz[l, 0] - dEdrlx
            fxyz[l, 1] = fxyz[l, 1] - dEdrly
            fxyz[l, 2] = fxyz[l, 2] - dEdrlz

            # dE/dZ*dZ/dr contribution to virial
            virial = virial - sz_r[nl] * (dEdrlx*sz_dx[nl] + dEdrly*sz_dy[nl] + dEdrlz*sz_dz[nl])
            virial_xyz[0] = virial_xyz[0] - sz_r[nl]*(dEdrlx * sz_dx[nl])
            virial_xyz[1] = virial_xyz[1] - sz_r[nl]*(dEdrly * sz_dy[nl])
            virial_xyz[2] = virial_xyz[2] - sz_r[nl]*(dEdrlz * sz_dz[nl])

            si1_sj1 = alatinv[0, 0]*sz_dx[nl] + alatinv[1, 0]*sz_dy[nl] + alatinv[2, 0]*sz_dz[nl]
            si2_sj2 = alatinv[0, 1]*sz_dx[nl] + alatinv[1, 1]*sz_dy[nl] + alatinv[2, 1]*sz_dz[nl]
            si3_sj3 = alatinv[0, 2]*sz_dx[nl] + alatinv[1, 2]*sz_dy[nl] + alatinv[2, 2]*sz_dz[nl]
            deralat[0, 0] = deralat[0, 0] - sz_r[nl] * dEdrlx * si1_sj1
            deralat[1, 0] = deralat[1, 0] - sz_r[nl] * dEdrlx * si2_sj2
            deralat[2, 0] = deralat[2, 0] - sz_r[nl] * dEdrlx * si3_sj3
            deralat[0, 1] = deralat[0, 1] - sz_r[nl] * dEdrly * si1_sj1
            deralat[1, 1] = deralat[1, 1] - sz_r[nl] * dEdrly * si2_sj2
            deralat[2, 1] = deralat[2, 1] - sz_r[nl] * dEdrly * si3_sj3
            deralat[0, 2] = deralat[0, 2] - sz_r[nl] * dEdrlz * si1_sj1
            deralat[1, 2] = deralat[1, 2] - sz_r[nl] * dEdrlz * si2_sj2
            deralat[2, 2] = deralat[2, 2] - sz_r[nl] * dEdrlz * si3_sj3

        ener = ener + ener_iat
        ener2 = ener2 + ener_iat**2

    etot = ener
    etot = etot/Ha_eV
    fxyz = -fxyz / Ha_eV * Bohr_Ang
    deralat = deralat / Ha_eV * Bohr_Ang
    vol = np.abs(np.linalg.det(alat0))

    # formula in fortran form of memory layout
    # stress = - deralat @ alat0.T / vol
    stress = - deralat.T @ alat0 / vol
    return etot, fxyz, stress

def test():
    from ase.io import read
    import bazant_fortran as bf

    import time

    atoms = read('test64.ascii')
    nat = len(atoms)

    atoms.rattle(0.01)

    import ase.units as units

    alat = atoms.get_cell().T
    positions = atoms.get_positions().T

    # convert to atomic units
    alat /= units.Bohr
    positions /= units.Bohr

    t1 = time.time()
    etot1, fxyz1, stress1 = bf.energyandforces_bazant(alat, positions)
    t2 = time.time()

    print(etot1)

    alat = alat.T
    positions = positions.T


    t3 = time.time()
    etot, fxyz, stress = energyandforces_bazant(alat, positions)
    t4 = time.time()

    print(etot)

    print("difference in energy:", etot - etot1)
    print("difference in forces:", np.linalg.norm(fxyz - fxyz1.T))
    print('difference in stress:', np.linalg.norm(stress - stress1))

    t5 = time.time()
    etot, fxyz, stress = energyandforces_bazant(alat, positions)
    t6 = time.time()

    print(etot)

    print('Fortran time:', t2 - t1)
    print('Python time:', t4 - t3)
    print('Python time (compiled):', t6 - t5)
    print('Speedup:', (t4 - t3) / (t2 - t1), (t6 - t5) / (t2 - t1))

if __name__ == '__main__':
    test()