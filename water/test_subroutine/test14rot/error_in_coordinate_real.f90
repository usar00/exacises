module MDPARAM
    implicit none
     integer, parameter :: N = 2, Nsite = 3
     double precision,parameter :: REQ = 0.94, & !A
               THETAEQ = 104.0, & !deg
               KBOND = 547.5, & !(kcal/mol)ang^-2
               KTHETA = 49.9, & !(kcal/mol)rad^-2
               A = 650000, & !(kcal/mol)ang^12
               B = 625.47, & !(kcal/mol)ang^6
               EL = 1.6021766208D-19, & !C
               DT = 1.0D-15, & !sec
               TEMP0 = 373, & !K
               EPS = 8.85418782D-12, & !permittivity of free space,(C^2)(J^-1)(m^-1)
               DENS = 1.0 !g/cm^3
     double precision, parameter ::  BOLZ = 1.380658D-23, &
                          AVO = 6.0221367D+23, &
                          BOHR = 0.529177249, &
                          HARTREE = 4.3597482D-18, &
                          AUTIME = 4.134136844D+16, &
                          EMASS = 9.1093897D-28, &
                          CAL = 4.814D3, &
                          PI = 3.14159265358979
     double precision, dimension(Nsite),parameter :: MASS = (/16.0, 1.0, 1.0/), & !g/mol
                                 Q = (/-0.82, 0.41, 0.41/)
     double precision, dimension(3, Nsite),parameter :: C = reshape( &
                                    (/0.000000, 0.000000, 0.000000, &
                                      0.940000, 0.000000, 0.000000, &
                                      -0.226567, 0.912287, 0.000000/), (/3, Nsite/))
end module MDPARAM
!*******************************************************************

  module VARIABLES
     use MDPARAM
     implicit none
     double precision,dimension(3,Nsite) :: C_rotated
     double precision :: E, T, U, Uintera, Uinter, TEMP
     double precision, dimension(3, Nsite, N) :: V, R, F, Fintera, Finter
  end module VARIABLES

!*******************************************************************
!*******************************************************************
    module common_toolbox
    contains
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
  end module common_toolbox
!*******************************************************************
program main
    use MDPARAM
    use VARIABLES
    implicit none
    print *,"C angle", ANGLE(C(1,1),C(1,2),C(1,3))
    call ROTATION
    print *,"C1 angle", ANGLE(C(1,1),C(1,2),C(1,3))
    call ROTATION
    print *,"C2=",ANGLE(C(1,1),C(1,2),C(1,3))
    call ROTATION
    print *,"Cx=",ANGLE(C(1,1),C(1,2),C(1,3))
    contains

     subroutine ROTATION
        use MDPARAM, only : pi
        use VARIABLES, only : C_rotated
        implicit none
        double precision :: PHI, PSI, THETA
        integer :: I
        call RANDOM_NUMBER(PHI)
        call RANDOM_NUMBER(PSI)
        call RANDOM_NUMBER(THETA)
        PHI = PHI*2.D0*PI
        PSI = PSI*2.D0*PI
        THETA = THETA*PI !0<THETA<PI
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
        double precision :: cosTHETA = 0.0D00,&
                   R12(3), R13(3), &
                  R12_distance, R13_distance
        integer :: i
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
