
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
     F = 0.0D0
	 Fintera=0.0D0
	 Finter=0.0D0
     U = 0.0D0
	 Uintera=0.0D0
	 Uinter=0.0D0
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

  subroutine FORCE
     use MDPARAM
     use VARIABLES
     implicit none
     F = 0.0D0
	 Fintera=0.0D0
	 Finter=0.0D0
     U = 0.0D0
	 Uintera=0.0D0
	 Uinter=0.0D0
     call FORCR_intra !intramolcule
!     call FORCE_inter(Finter,Uinter) !intermolcule
	 F=Fintera+Finter
	 U=Uintera+Uinter
  end subroutine FORCE

!*******************************************************************
!*******************************************************************
!*******************************************************************
  subroutine FORCR_intra
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
!*******************************************************************
  subroutine FORCETEST
     use MDPARAM
     use VARIABLES
     implicit none
     integer :: I, J, K
     real*8 :: D = 1.0D-5, U0, Up2, Up, Un, Un2, R0
     real*8, dimension(3, Nsite, N) :: F0, F1
     real*8 :: DAVG = 0.0D0

     call FORCE
     U0 = U
     F0 = F
     do I = 1, N
        do J = 1, Nsite
           do K = 1, 3
              R0 = R(K, J, I)
              R(K, J, I) = R0 + D + D
              call FORCE
              Up2 = U
              R(K, J, I) = R0 + D
              call FORCE
              Up = U
              R(K, J, I) = R0 - D
              call FORCE
              Un = U
              R(K, J, I) = R0 - D - D
              call FORCE
              Un2 = U
              F1(K, J, I) = -(8.0D0*(Up - Un) - (Up2 - Un2))/(D*1.2D1)
              R(K, J, I) = R0
              DAVG = DAVG + (F1(K, J, I) - F0(K, J, I))**2
!       print *,i,j,k,Up,Un,Up2,Un2,F0(K,J,I),F1(K,J,I)
           end do
           print '(I5,6E14.6)', I, (F0(K, J, I), K=1, 3), (F1(K, J, I) - F0(K, J, I), K=1, 3)
        end do
     end do
     DAVG = SQRT(DAVG/(3*N))
     print *, 'DAVG=', DAVG
     return
  end subroutine FORCETEST
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
        R=R+(C1(i) - C2(i))**2.0D00
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
            Cmass(j,i) = Cmass(j,i)+ C(j,i) * M(i)
        enddo
        wmass = wmass + M(i) 
        enddo
    Cmass = Cmass / wmass
end subroutine CENTER_OF_MASS
