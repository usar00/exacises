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
     real*8, dimension(3, Nsite) :: C_original = reshape( &
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
     print *, R(1,1,1)
end program MDMAIN
!*******************************************************************

