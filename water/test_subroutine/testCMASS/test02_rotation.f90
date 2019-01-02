  module MDPARAM
     implicit none
     integer, parameter :: N = 64, Nsite = 3
     real*8 :: REQ = 0.94, & !A
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

program main
     implicit none
     integer, parameter :: N = 64, Nsite = 3
     integer :: i,j
     real*8 :: WMASS
     real*8, dimension(Nsite) :: MASS = (/16.0, 1.0, 1.0/), & !g/mol
                                 Q = (/-0.82, 0.41, 0.41/), & !e
                                 C_cmass
     real*8, dimension(3, Nsite) :: C = reshape( &
                                    (/0.000000, 0.000000, 0.000000, &
                                      0.940000, 0.000000, 0.000000, &
                                      -0.226567, 0.912287, 0.000000/), (/3, Nsite/))
     WMASS = 0.0D0
     do I = 1, Nsite
        WMASS = WMASS + MASS(I)
     end do
     do I = 1, Nsite
         do J = 1, 3
            C_cmass(I)=C_cmass(I) + C(I,J) * MASS(J) 
         enddo
     enddo
     C_cmass = C_cmass / WMASS
     do I = 1, Nsite
     C(:,I) = C(:,I) - C_cmass
     print *,C_cmass
     print *,C
end program main
!*******************************************************************
  subroutine ROTATION
     use MDPARAM
     implicit none
     real*8 :: PHI, PSI, THETA,PI=acos(-1.0D0),C_rotated(3,Nsite)
     integer :: I
     call RANDOM_NUMBER(PHI)
     call RANDOM_NUMBER(PSI)
     call RANDOM_NUMBER(THETA)
     PHI = PHI*2.0D0*PI
     PSI = PSI*2.0D0*PI
     THETA = THETA*PI !0<THETA<PI
     C_rotated = C
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
!*******************************************************************
  subroutine FORCR_intra(Fintera,Uintera)
     use MDPARAM
     use VARIABLES
     implicit none
     integer ::i,j,k
     real*8, dimension(3, Nsite) ::DTDX123
     real*8 :: R1(3), R2(3),R3(3),A213, R12,R13, DRDX12(3, 2),DRDX13(3, 2),Ua
        do i=1,N
            do j=1,3
                R1(j)= R(j,1,i)
                R2(j)= R(j,2,i)
                R3(j)= R(j,3,i)
            enddo
            call BOND(R1(1),R2(1),R12,DRDX12(1,1)) 
            call BOND(R1(1),R3(1),R13,DRDX13(1,1)) 
            call ANGLE(R1(1),R2(1),R3(1),A213,DTDX(1,1)) 
            do j=1,3
                    Fintera(j,1,i)= 2.0D0*KTHETA*(A213 - THETAEQ)*DTDX123(j,1)& 
                                    +2.0D0*KBOND*(R12 - REQ)*DRDX12(j,1)&
                                    +2.0D0*KBOND*(R13 - REQ)*DRDX13(j,1)
                    Fintera(j,2,i)= 2.0D0*KTHETA*(A213 - THETAEQ)*DTDX123(j,2)& 
                                    +2.0D0*KBOND*(R12 - REQ)*DRDX12(j,2)&
                    Fintera(j,3,i)= 2.0D0*KTHETA*(A213 - THETAEQ)*DTDX123(j,2)& 
                                    +2.0D0*KBOND*(R13 - REQ)*DRDX13(j,3)
            end do
            Uintera = Uintera + KTHETA*(THETA - THETAEQ)**2 +KBOND*(R12 - REQ)**2 + KBOND*(R13 - REQ)**2 
        end do
  end subroutine FORCR_intra


!*******************************************************************
subroutine FORCE_inter(Finter,Uinter)
    use MDPARAM
    use VARIABLES
    implicit none
    integer :: i,j,k,Isite,Jsite,xyz
    real*8 :: R_Isite(3,3),R_Jsite(3,3),R_IJ(3,3,3),R_IJsq,Rij_distance,Fi(3,3)
    do i=1,N
        do j=1, i-1
                do xyz=1,3
                    R_Isite(xyz,1) = R(xyz,1,i)
                    R_Jsite(xyz,1) = R(xyz,1,j)
                    R_IJ(xyz)=R_Jsite(xyz,1)-R_Isite(xyz,1)
                    R_IJ(xyz)=R_Jsite(xyz,1)-R_Isite(xyz,1) - ANINT(R_Jsite(xyz,1)-R_Isite(xyz,1),RBOX)*RBOX
                    R_IJsq =R_IJsq +  R_IJ(xyz) **2
                enddo
                if(R_IJsq < RCUT **2) then
                    do xyz=1,3
                        do k=1,3
                        R_Isite(xyz,k)=R(xyz,k,i)
                        R_Jsite(xyz,k)=R(xyz,k,j) - ANINT(R_Jsite(xyz,1)-R_Isite(xyz,1),RBOX)*RBOX
                        enddo
                    enddo
                    do Isite=1,Nsite
                        do Isite=1,Nsite
                            if(Isite==1 .and. Jsite==1) then
                                do xyz =1,3
                                Fi(xyz,1) =Fi(xy,z)  (12.0D0*A*RIJ(XYZ)/RSQ**14 &
                                                                 - 6.0D0*B*RIJ(XYZ)/RSQ**8 + Q(ISITE)*Q(JSITE)*RIJ(XYZ)/RSQ**3)
    


                enddo
                    R_IJ= 
!*******************************************************************
subroutine BOND(C1, C2, R, DRDX,Ua)
   implicit none
   real*8 :: C1(3), C2(3), R, DRDX(3, 2)
   integer :: i
   DRDX(3, 2) = 0.0D00
   R = 0.0D00
   do i = 1, 3
      R = R + (C1(i) - C2(i))**2.0D00
   enddo
   R = sqrt(R)
!*******
   do i = 1, 3
      DRDX(i, 1) = -C1(i)/R
      DRDX(i, 2) = C2(i)/R
   enddo
end subroutine BOND
!*******************************************************************
subroutine DISTANCE(C1,C2,R)
    implicit none
    integer :: i
    real*8 :: C1(3),C2(3),R
    R = 0.D00
    do i=1,3
        R=R+C1(i)**2+C2(i)**2
     enddo
    R = sqrt(R)
end subroutine DISTANCE
!*******************************************************************
subroutine CENTER_OF_MASS(C,M,N,Cmass,wmass)
    implicit none
    integer :: N,i,j
    real*8 :: C(3,N),M(N),Cmass(3,N)
    real*8,optional :: Wmass
    Cmass = 0.D0
    wmass = 0.D0
    do i =1,N
        do j=1,3
            Cmass(j,i) = C(j,i) * M(i)
        enddo
        wmass = wmass + M(i) 
        enddo
    Cmass = Cmass / wmass
end subroutine CENTER_OF_MASS

!*******************************************************************
subroutine ANGLE(C1, C2, C3, THETA, DTDX)
   implicit none
   real*8 ::C1(3), C2(3), C3(3), THETA, cosTHETA = 0.0D00, sinTHETA, &
             DTDX(3, 3), R12(3), R13(3), &
             R12_distance, R13_distance
   integer :: i
   DTDX(3, 3) = 0.0D00
   R12(3) = 0.0D00
   R13(3) = 0.0D00
   R12_distance = 0.0D00
   R13_distance = 0.0D00
!*******
!******* caluculate vector
   do i = 1, 3 !i represent xyz
      R12(i) = C2(i) - C1(i)
      R13(i) = C3(i) - C1(i)
   enddo
   do i = 1, 3
      R12_distance = R12_distance + R12(i)**2
      R13_distance = R13_distance + R13(i)**2
   enddo
   R12_distance = sqrt(R12_distance)
   R13_distance = sqrt(R13_distance)
!*******
!******* caluculate triangle function of theta
   do i = 1, 3
      cosTHETA = cosTHETA + R12(i)*R13(i)
   enddo
   cosTHETA = cosTHETA/(R12_distance*R13_distance)
   !if(cosTHETA > 1.0D00) print *,"cosTHETA has an erroer",cosTHETA
   sinTHETA = 1.0D00 - cosTHETA**2.0D00
   THETA = acos(cosTHETA)
   !if(abs(sinTHETA - sin(THETA)) > 0.000001D00) print*,"THETA error",THETA
!******
!****** caluculate DTDX
   do i = 1, 3
      DTDX(i, 2) = -(R13(i) - cosTHETA*(R12_distance/R13_distance)*R12(i)) &
                   /(sinTHETA*R12_distance*R13_distance)
      DTDX(i, 3) = -(R12(i) - cosTHETA*(R13_distance/R12_distance)*R13(i)) &
                   /(sinTHETA*R12_distance*R13_distance)
      DTDX(i, 1) = -DTDX(i, 2) - DTDX(i, 3)
   enddo
end subroutine ANGLE
!******
