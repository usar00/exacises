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
     call INIT
!     open (IOUT, file=FOUT, status='replace', access='sequential', form='formatted')
!     call INITOUTPUT
 !    open (IIOUT, file=FFOUT, status='replace', access='sequential', form='formatted')
     call FORCE
!  call FORCETEST
     do ISTEP = 1, NEQUIL
        call MOVEA
        call FORCE
        call MOVEB
        if (mod(ISTEP, NINTERVAL) == 0) call VSCALE
     end do

     do ISTEP = 0, MAXSTP
        call MOVEA
        call FORCE
        call MOVEB
        if (mod(ISTEP, CPRINT) == 0) call OUTPUT(ISTEP)
        if (mod(ISTEP, NPRINT) == 0) call COORDOUT
     end do
!     close (IOUT, status='keep')
!     close (IIOUT, status='keep')
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
            C_cmass(I)=C_cmass(I) + C_original(I,J) * MASS(J) 
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
     call FORCR_intra(Fintera,Uintera) !intramolcule
     call FORCE_inter(Finter,Uinter) !intermolcule
	 F=Fintera+Finter
	 U=Uintera+Uinter
  end subroutine FORCE

!*******************************************************************

  subroutine FORCR_intra
     use MDPARAM
     use VARIABLES
     implicit none
     integer :: I, J, K
     real*8, dimension(3, Nsite) :: RIJ
                  + 2.0D0*KTHETA*(THETA - THETAEQ) &
                  *(RIJ(K, 2)*cos(THETA)/(R1**2) - RIJ(K, 3)/(R1*R2))/sin(THETA) ! d(FORCR_intra)/d(H1)
           FH2O = 2.0D0*K1*(R2 - REQ)*RIJ(K, 3)/R2 &
                  + 2.0D0*KTHETA*(THETA - THETAEQ) &
                  *(RIJ(K, 3)*cos(THETA)/(R2**2) - RIJ(K, 2)/(R1*R2))/sin(THETA) ! d(FORCR_intra)/d(H2)
           F(K, 2, I) = F(K, 2, I) - FH1O
           F(K, 3, I) = F(K, 3, I) - FH2O
           F(K, 1, I) = F(K, 1, I) + FH1O + FH2O
        end do
!    print *,i,THETA,KTHETA,THETAEQ,KTHETA*(THETA-THETAEQ)**2
        U = U + K1*(R1 - REQ)**2 + K1*(R2 - REQ)**2 + KTHETA*(THETA - THETAEQ)**2
     end do
  end subroutine FORCR_intra

!*******************************************************************

  subroutine FORCE_inter
     use MDPARAM
     use VARIABLES
     implicit none
     integer :: I, J, ISITE, JSITE, XYZ
     real*8, dimension(3) :: RIJ
     real*8 :: RSQ
     do I = 2, N
        do J = 1, I - 1
           do ISITE = 1, Nsite
              do JSITE = 1, Nsite
                 do XYZ = 1, 3
                    RIJ(XYZ) = R(XYZ, ISITE, I) - R(XYZ, JSITE, J)
                    RIJ(XYZ) = RIJ(XYZ) - RBOX*ANINT(RIJ(XYZ)/RBOX)
                 end do
                 RSQ = RIJ(1)**2 + RIJ(2)**2 + RIJ(3)**2
                 if (RSQ < (RCUT**2)) then
                    RSQ = SQRT(RSQ)
!            if(I==2) print *,RSQ
                    U = U + Q(ISITE)*Q(JSITE)/RSQ
                    if (ISITE == 1 .and. JSITE == 1) then
                       U = U + A/RSQ**12 - B/RSQ**6
                       do XYZ = 1, 3
                          F(XYZ, ISITE, I) = F(XYZ, ISITE, I) + (12.0D0*A*RIJ(XYZ)/RSQ**14 &
                                                                 - 6.0D0*B*RIJ(XYZ)/RSQ**8 + Q(ISITE)*Q(JSITE)*RIJ(XYZ)/RSQ**3)
                          F(XYZ, JSITE, J) = F(XYZ, JSITE, J) - (12.0D0*A*RIJ(XYZ)/RSQ**14 &
                                                                 - 6.0D0*B*RIJ(XYZ)/RSQ**8 + Q(ISITE)*Q(JSITE)*RIJ(XYZ)/RSQ**3)
                       end do
                    else
                       do XYZ = 1, 3
                          F(XYZ, ISITE, I) = F(XYZ, ISITE, I) + Q(ISITE)*Q(JSITE)*RIJ(XYZ)/RSQ**3
                          F(XYZ, JSITE, J) = F(XYZ, JSITE, J) - Q(ISITE)*Q(JSITE)*RIJ(XYZ)/RSQ**3
                       end do
                    end if
                 end if
              end do
           end do
        end do
     end do
  end subroutine FORCE_inter

!*******************************************************************

  subroutine MOVEA
     use MDPARAM
     use VARIABLES
     implicit none
     integer :: I, J, K
     do I = 1, N
        do J = 1, Nsite
           do K = 1, 3
              R(K, J, I) = R(K, J, I) + V(K, J, I)*DT + (F(K, J, I)*DT**2)/(2.0D0*MASS(J))
              V(K, J, I) = V(K, J, I) + F(K, J, I)*DT/(2.0D0*MASS(J))
           end do
        end do
     end do
  end subroutine MOVEA

!*******************************************************************

  subroutine MOVEB
     use MDPARAM
     use VARIABLES
     implicit none
     integer :: I, J, K
     T = 0.0D0
     do I = 1, N
        do J = 1, Nsite
           do K = 1, 3
              V(K, J, I) = V(K, J, I) + F(K, J, I)*DT/(2.0D0*MASS(J))
              T = T + (V(K, J, I)**2*MASS(J))/2.0D0
           end do
        end do
     end do
  end subroutine MOVEB

!*******************************************************************

  subroutine VSCALE
     use MDPARAM
     use VARIABLES
     implicit none
     real*8 :: FAC
     TEMP = T*2.0*HARTREE/BOLZ/(3*N*Nsite - 3)
     FAC = SQRT(TEMP0/TEMP)
     V = V*FAC
     T = T*(FAC**2)
     E = T + U
!  print *,TEMP0,TEMP
  end subroutine VSCALE

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

  subroutine ROTATION
     use MDPARAM
     implicit none
     real*8 :: PHI, PSI, THETA
     integer :: I
     call RANDOM_NUMBER(PHI)
     call RANDOM_NUMBER(PSI)
     call RANDOM_NUMBER(THETA)
     PHI = PHI*2*PI
     PSI = PSI*2*PI
     THETA = THETA*PI !0<THETA<PI
     do I = 1, Nsite
        C(1, I) = (cos(PSI)*cos(PHI) - cos(THETA)*sin(PHI)*sin(PSI))  *C_original(1, I) &
                     + (cos(PSI)*sin(PHI) + cos(THETA)*cos(PHI)*sin(PSI))*C_original(2, I) &
                     + (sin(PSI)*sin(THETA))                             *C_original(3, I)
        C(2, I) = (-sin(PSI)*cos(PHI) - cos(THETA)*sin(PHI)*cos(PSI))  *C(1, I) &
                     + (-sin(PSI)*sin(PHI) + cos(THETA)*cos(PHI)*cos(PSI))*C(2, I) &
                     + (cos(PSI)*sin(THETA))                              *C(3, I)
        C(3, I) = (sin(THETA)*sin(PHI))               *C(1, I) &
                     + (-sin(THETA)*cos(PHI))            *C(2, I) &
                     + cos(THETA)                        *C(3, I)
     end do
  end subroutine ROTATION

!*******************************************************************
  subroutine BOND(C1, C2, R, DRDX)
     implicit none
     real*8 :: C1(3), C2(3), R, DRDX(3, 2)

  end subroutine BOND
!*******************************************************************
  subroutine ANGLE(C1, C2, C3, DTDX)
     implicit none
     real*8 ::C1(3), C2(3), C3(3), DTDX(3, 3)

  end subroutine ANGLE
!*******************************************************************

  subroutine INITOUTPUT
     use MDPARAM
     use VARIABLES
     implicit none
     integer :: I
     write (IOUT, *) N*Nsite
     write (IOUT, *) 'Water molecule'
     do I = 1, N
        write (IOUT, *) 'O', R(1, 1, I), R(2, 1, I), R(3, 1, I)
        write (IOUT, *) 'H', R(1, 2, I), R(2, 3, I), R(3, 3, I)
        write (IOUT, *) 'H', R(1, 3, I), R(2, 3, I), R(3, 3, I)
     end do
  end subroutine INITOUTPUT

!*******************************************************************

  subroutine OUTPUT(ISTEP)
     use VARIABLES
     use MDPARAM
     implicit none
     integer :: ISTEP
     E = T + U
     if (ISTEP == 0) write (*, 9000)
     write (*, 9010) ISTEP, DT*ISTEP, E, T, U
9000 format('#STEP#', '    TIME', '     TOTAL ENERGY', '     KINETIC', &
            '       POTENTIAL''     unit: a.u.')
9010 format(I7, E10.3, 3E15.6)
  end subroutine OUTPUT

!*******************************************************************

  subroutine COORDOUT
     use MDPARAM
     use VARIABLES
     implicit none
     integer :: I
     write (IIOUT, *) N*Nsite
     write (IIOUT, *) 'Water molecule'
     do I = 1, N
        write (IIOUT, *) 'O', R(1, 1, I), R(2, 1, I), R(3, 1, I)
        write (IIOUT, *) 'H', R(1, 2, I), R(2, 3, I), R(3, 3, I)
        write (IIOUT, *) 'H', R(1, 3, I), R(2, 3, I), R(3, 3, I)
     end do
  end subroutine COORDOUT

!*******************************************************************

!*******************************************************************

!*******************************************************************

!*******************************************************************

!*******************************************************************

!*******************************************************************
