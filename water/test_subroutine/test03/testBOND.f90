program main
   implicit none
   real*8 :: C(3, 3), R, DRDX(3, 2), DTDX(3, 3), Kr, Kt, THETA,Ua,Ub
   integer :: i, j
   C(:, 1) = (/0.0D00, 0.0D00, 0.0D00/)
   C(:, 2) = (/1.0D00, 1.0D00, 0.0D00/)
   C(:, 3) = (/-1.0D00, 1.0D00, 0.0D00/)
   print *, "C=", C
   call BOND(C(1, 1), C(2, 1), R, DRDX,Ua)
   print *, "R=", R
   print *, "DRDX=", DRDX
   call ANGLE(C(1, 1), C(1, 2), C(1, 3), THETA, DTDX,Ub)
   print *, "theta=", theta
   print *, "DTDX=", DTDX
!*******
!******
!*******
   print *, "************"
   print *, "***ROUND2***"
   print *, "************"
   do i = 1, 3
      do j = 1, 3
         C(i, j) = C(i, j) + 1.0D00
      enddo
   enddo
   call BOND(C(1, 1), C(2, 1), R, DRDX)
   print *, "R=", R
   print *, "DRDX=", DRDX
   call ANGLE(C(1, 1), C(1, 2), C(1, 3), THETA, DTDX)
   print *, "theta=", theta
   print *, "DTDX=", DTDX

!*******
!*******
!*******
   print *, "************"
   print *, "***ROUND3***"
   print *, "************"
   call ROTATION(C(1, 1))
   call ROTATION(C(1, 2))
   call ROTATION(C(1, 3))
   call ANGLE(C(1, 1), C(1, 2), C(1, 3), THETA, DTDX)
   print *, "theta=", theta
   print *, "DTDX=", DTDX
end program main!

!*******
!*******
!*******
!*******
!*******
!*******
subroutine BOND(C1, C2, R, DRDX,Ua)
   implicit none
   real*8 :: C1(3), C2(3), R, DRDX(3, 2),Ua
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
!******
!******
subroutine FORCETEST




subroutine ROTATION(C)
   implicit none
   real*8 :: PHI, PSI, THETA, C(3), C_new(3), PI
   PI = acos(-1.0D00)
   PHI = 0.223D00*2*PI
   PSI = 0.378D00*2*PI
   THETA = 0.181D00*PI !0<THETA<PI
   C_new(1) = (cos(PSI)*cos(PHI) - cos(THETA)*sin(PHI)*sin(PSI))*C(1) &
              + (cos(PSI)*sin(PHI) + cos(THETA)*cos(PHI)*sin(PSI))*C(2) &
              + (sin(PSI)*sin(THETA))*C(3)
   C_new(2) = (-sin(PSI)*cos(PHI) - cos(THETA)*sin(PHI)*cos(PSI))*C(1) &
              + (-sin(PSI)*sin(PHI) + cos(THETA)*cos(PHI)*cos(PSI))*C(2) &
              + (cos(PSI)*sin(THETA))*C(3)
   C_new(3) = (sin(THETA)*sin(PHI))*C(1) &
              + (-sin(THETA)*cos(PHI))*C(2) &
              + cos(THETA)*C(3)
   C = C_new
!    print *, "C=" , C
end subroutine ROTATION

!******
!******
!******
