module MDPARAM

     implicit none
     integer, parameter :: N = 2, Nsite = 3
     double precision :: REQ = 0.94d0, & !A
               THETAEQ = 104.0d0, & !deg
               KBOND = 547.5d0, & !(kcal/mol)ang^-2
               KTHETA = 49.9d0, & !(kcal/mol)rad^-2
               A = 650000d0, & !(kcal/mol)ang^12
               B = 625.47d0, & !(kcal/mol)ang^6
               EL = 1.6021766208D-19, & !C
               DT = 1.0D-15, & !sec
               TEMP0 = 373d0, & !K
               EPS = 8.85418782D-12, & !permittivity of free space,(C^2)(J^-1)(m^-1)
               DENS = 1.0 !g/cm^3
     double precision, parameter ::  BOLZ = 1.380658D-23, &
                          AVO = 6.0221367D+23, &
                          BOHR = 0.529177249d0, &
                          HARTREE = 4.3597482D-18, &
                          AUTIME = 4.134136844D+16, &
                          EMASS = 9.1093897D-28, &
                          CAL = 4.814D3, &
                          PI = 3.14159265358979d0
     double precision, dimension(Nsite) :: MASS = (/16.0d0, 1.0d0, 1.0d0/), & !g/mol
                                 Q = (/-0.82d0, 0.41d0, 0.41d0/), & !e
                                 C_cmass
     double precision, dimension(3, Nsite) :: C = reshape( &
                                    (/0.000000d0, 0.000000d0, 0.000000d0, &
                                      0.940000d0, 0.000000d0, 0.000000d0, &
                                      -0.226567d0, 0.912287d0, 0.000000d0/), (/3, Nsite/)),&
                                      C_rotated
     double precision :: RBOX, RCUT, WMASS, KC
     integer :: MAXSTP = 10000, &
                NPRINT = 10, &
                CPRINT = 500, &
                NEQUIL = 10000, &
                NINTERVAL = 100
     integer :: IOUT = 1, IIOUT = 2
     end module MDPARAM
  module VARIABLES
     use MDPARAM
     implicit none
     double precision :: E, T, U, Uintera, Uinter, TEMP
     double precision, dimension(3, Nsite, N) :: V, R, F, Fintera, Finter
  end module VARIABLES

!*******************************************************************
!*******************************************************************
!*******************************************************************
    module common_toolbox
    contains
!*******************************************************************
   function DISTANCE(C1, C2)
     implicit none
     integer :: i
     double precision :: C1(3), C2(3),DISTANCE
     DISTANCE = 0.D00
     do i = 1, 3
        DISTANCE = DISTANCE + (C1(i) - C2(i))**2.0D00
     enddo
     DISTANCE = sqrt(DISTANCE)
  end function DISTANCE
!*******************************************************************
  subroutine CENTER_OF_MASS(C, M, Nsite, C_cmass, wmass)
     implicit none
     integer :: Nsite, i, j
     double precision :: C(3, Nsite), M(Nsite), C_cmass(3, Nsite)
     double precision, optional :: Wmass
     C_cmass = 0.D0
     wmass = 0.D0
     do i = 1, Nsite
        do j = 1, 3
           C_cmass(j, i) = C_cmass(j, i) + C(j, i)*M(i)
        enddo
        wmass = wmass + M(i)
     enddo
     C_cmass = C_cmass/wmass
  end subroutine CENTER_OF_MASS
!*******************************************************************
  end module common_toolbox
!*******************************************************************
program main
    use MDPARAM , only :C,MASS,Nsite,C_cmass,wmass
    use VARIABLES, only :C_rotated
    use common_toolbox,only :CENTER_OF_MASS
    implicit none
    integer :: i,j
    call CENTER_OF_MASS(C,MASS,Nsite,C_cmass,wmass)
    do i=1,Nsite
        do j=1,3
            C(j,i) = C(j,i) -C_cmass(j) 
        enddo
    enddo
    C_rotated = C
    print *,"C_rotated angle", ANGLE(C_rotated(1,1),C_rotated(1,2),C_rotated(1,3))
    do i=1,3
        print *,C_rotated(:,i)
    enddo
    call ROTATION
    print *,"C_rotated1 angle", ANGLE(C_rotated(1,1),C_rotated(1,2),C_rotated(1,3))
    do i=1,3
        print *,C_rotated(:,i)
    enddo
    call ROTATION
    print *,"C_rotated2=",ANGLE(C_rotated(1,1),C_rotated(1,2),C_rotated(1,3))
    do i=1,3
        print *,C_rotated(:,i)
    enddo
    call ROTATION
    print *,"C_rotated3=",ANGLE(C_rotated(1,1),C_rotated(1,2),C_rotated(1,3))
    do i=1,3
        print *,C_rotated(:,i)
    enddo
    contains

     subroutine ROTATION
        use MDPARAM, only : pi,C,Nsite
        use VARIABLES, only : C_rotated
        implicit none
        double precision :: PHI, PSI, THETA
        integer :: I
!        call RANDOM_NUMBER(PHI)
!        call RANDOM_NUMBER(PSI)
!        call RANDOM_NUMBER(THETA)
        PHI = 0.5d0*2.D0*PI
        PSI =0.d0*2.D0*PI
        THETA = 0.d0*PI !0<THETA<PI
        do I = 1, Nsite
           C_rotated(1, I) = (cos(PSI)*cos(PHI) - cos(THETA)*sin(PHI)*sin(PSI))*C(1, I) &
                             + (cos(PSI)*sin(PHI) + cos(THETA)*cos(PHI)*sin(PSI))*C(2, I) &
                             + (sin(PSI)*sin(THETA))*C(3, I)
           C_rotated(2, I) = (-sin(PSI)*cos(PHI) - cos(THETA)*sin(PHI)*cos(PSI))*C(1, I) &
                             + (-sin(PSI)*sin(PHI) + cos(THETA)*cos(PHI)*cos(PSI))*C(2, I) &
                             + (cos(PSI)*sin(THETA))*C(3, I)
           C_rotated(3, I) = (sin(THETA)*sin(PHI))*C(1, I) &
                             + (-sin(THETA)*cos(PHI))*C(2, I) &
                             + cos(THETA)*C(3, I)
        end do
     end subroutine ROTATION
!*******************************************************************
     function ANGLE(C1, C2, C3)
        use common_toolbox, only : DISTANCE 
        implicit none
        double precision,intent(in) ::C1(3), C2(3), C3(3)
        double precision :: ANGLE 
        double precision :: cosTHETA,&
                   R12(3), R13(3), &
                  R12_distance, R13_distance
        integer :: i
        cosTHETA = 0.0D00
        R12(3) = 0.0D00
        R13(3) = 0.0D00
        do i = 1, 3 !i represent xyz
           R12(i) = C2(i) - C1(i)
           R13(i) = C3(i) - C1(i)
        enddo
        do i = 1, 3
        enddo
        R12_distance = DISTANCE(C1(1),C2(1))
        R13_distance = DISTANCE(C1(1),C3(1))
!******* caluculate triangle function of theta
        do i = 1, 3
           cosTHETA = cosTHETA + R12(i)*R13(i)
        enddo
        cosTHETA = cosTHETA/(R12_distance*R13_distance)
        !if(cosTHETA > 1.0D00) print *,"cosTHETA has an erroer",cosTHETA
        ANGLE = acos(cosTHETA)
        end function ANGLE
!*******************************************************************
end program main
