
  module MDPARAM
     implicit none
     integer, parameter :: N = 64, Nsite = 3
     real*8 :: REQ = 0.94, & !A
               THETAEQ = 104.0, & !deg
               K1 = 547.5, & !(kcal/mol)ang^-2
               KTHETA = 49.9, & !(kcal/mol)rad^-2
               A = 650000, & !(kcal/mol)ang^12
               B = 625.47, & !(kcal/mol)ang^6
               EL = 1.6021766208D-19, & !C
               DT = 1.0D-15, & !sec
               TEMP0 = 373, & !K
               EPS = 8.85418782D-12, & !permittivity of free space,(C^2)(J^-1)(m^-1)
               DENS = 1.0 !g/cm^3
     real*8, parameter ::  BOLZ = 1.380658D-23, &
                          AVO = 6.0221367D+23, &
                          BOHR = 0.529177249, &
                          HARTREE = 4.3597482D-18, &
                          AUTIME = 4.134136844D+16, &
                          EMASS = 9.1093897D-28, &
                          CAL = 4.814D3, &
                          PI = 3.14159265358979
     real*8, dimension(Nsite) :: MASS = (/16.0, 1.0, 1.0/), & !g/mol
                                 Q = (/-0.82, 0.41, 0.41/), & !e
                                 C_cmass
     real*8, dimension(3, Nsite) :: C = reshape( &
                                    (/0.000000, 0.000000, 0.000000, &
                                      0.940000, 0.000000, 0.000000, &
                                      -0.226567, 0.912287, 0.000000/), (/3, Nsite/))
     real*8 :: RBOX, RCUT, WMASS, KC
     integer :: MAXSTP = 10000, &
                NPRINT = 10, &
                CPRINT = 500, &
                NEQUIL = 10000, &
                NINTERVAL = 100
     integer :: IOUT = 1, IIOUT = 2
     character(LEN=20) :: FOUT = 'water_init.xyz', FFOUT = 'water.xyz'
  end module MDPARAM

!*******************************************************************

  module VARIABLES
     use MDPARAM
     implicit none
     real*8 :: E, T, U,Uintera,Uinter, TEMP
     real*8, dimension(3, Nsite, N) :: V, R, F,Fintera,Finter
  end module VARIABLES

!*******************************************************************

  program MDMAIN
     use MDPARAM
     use VARIABLES
     implicit none
     integer :: ISTEP, I, J
     call INIT
     call F
end program MDMAIN
!*******************************************************************

  subroutine INIT
     use MDPARAM
     use VARIABLES
     implicit none
     integer :: IX, IY, IZ, I, J, IMAX
     real*8 :: DR, SMALL, RAND, SMALLNEW

!atomic units
     WMASS = 0.0D0
     do I = 1, Nsite
        WMASS = WMASS + MASS(I)
     end do
     do I = 1, 3
         do J = 1, 3
            C_cmass(I)=C_cmass(I) + C(I,J) * MASS(J) 
         enddo
     enddo
     C_cmass = C_cmass / WMASS
     DENS = (DENS/WMASS)*AVO*(BOHR*1.0D-8)**3
     MASS = (MASS/AVO)/EMASS
     WMASS = (WMASS/AVO)/EMASS
     RBOX = (N/DENS)**(1.0D0/3.0D0)
     RCUT = RBOX/2.0D0
     C = C/BOHR
     DT = DT*AUTIME
     REQ = REQ/BOHR
     THETAEQ = THETAEQ/(180.0D0/PI)
     K1 = K1*(BOHR**2)*CAL/HARTREE/AVO
     KTHETA = KTHETA*CAL/HARTREE/AVO
     A = A/(BOHR**12)*CAL/HARTREE/AVO
     B = B/(BOHR**6)*CAL/HARTREE/AVO

     IMAX = INT((N - 1)**(1.0D0/3.0D0)) + 1
     DR = RBOX/IMAX
     SMALL = DR*0.2D0
     J = 0
     do IX = 1, IMAX
        do IY = 1, IMAX
           do IZ = 1, IMAX
              J = J + 1
              if (J > N) exit
              call ROTATION
              call RANDOM_NUMBER(RAND)
              R(1, 1, J) = C(1, 1) + IX*DR + SMALL*(2.0D0*RAND - 1.0D0)
              R(2, 1, J) = C(2, 1) + IY*DR + SMALL*(2.0D0*RAND - 1.0D0)
              R(3, 1, J) = C(3, 1) + IZ*DR + SMALL*(2.0D0*RAND - 1.0D0)
              R(1, 2, J) = C(1, 2) + IX*DR + SMALL*(2.0D0*RAND - 1.0D0)
              R(2, 2, J) = C(2, 2) + IY*DR + SMALL*(2.0D0*RAND - 1.0D0)
              R(3, 2, J) = C(3, 2) + IZ*DR + SMALL*(2.0D0*RAND - 1.0D0)
              R(1, 3, J) = C(1, 3) + IX*DR + SMALL*(2.0D0*RAND - 1.0D0)
              R(2, 3, J) = C(2, 3) + IY*DR + SMALL*(2.0D0*RAND - 1.0D0)
              R(3, 3, J) = C(3, 3) + IZ*DR + SMALL*(2.0D0*RAND - 1.0D0)
!        print *,DR,IX*DR,IY*DR,IZ*DR
           end do
        end do
     end do
     V = 0.0D0
  end subroutine INIT
!*******************************************************************
  subroutine ROTATION
     use MDPARAM
     implicit none
     real*8 :: PHI, PSI, THETA,C_rotated(3,Nsite)
     integer :: I
     call RANDOM_NUMBER(PHI)
     call RANDOM_NUMBER(PSI)
     call RANDOM_NUMBER(THETA)
     PHI = PHI*2.0D0*PI
     PSI = PSI*2.0D0*PI
     THETA = THETA*PI !0<THETA<PI
     do I = 1, Nsite
        C_rotated(1, I) = (cos(PSI)*cos(PHI) - cos(THETA)*sin(PHI)*sin(PSI))  *C(1, I) &
                     + (cos(PSI)*sin(PHI) + cos(THETA)*cos(PHI)*sin(PSI))*C(2, I) &
                     + (sin(PSI)*sin(THETA))                             *C(3, I)
        C_rotated(2, I) = (-sin(PSI)*cos(PHI) - cos(THETA)*sin(PHI)*cos(PSI))  *C(1, I) &
                     + (-sin(PSI)*sin(PHI) + cos(THETA)*cos(PHI)*cos(PSI))*C(2, I) &
                     + (cos(PSI)*sin(THETA))                              *C(3, I)
        C_rotated(3, I) = (sin(THETA)*sin(PHI))               *C(1, I) &
                     + (-sin(THETA)*cos(PHI))            *C(2, I) &
                     + cos(THETA)                        *C(3, I)
     end do
    C = C_rotated
  end subroutine ROTATION
!*******************************************************************
