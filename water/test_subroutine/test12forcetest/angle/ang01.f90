program main
 implicit none
 double precision,dimension(3) :: c1,c2,c3
 double precision ::  theta,r12,t13
 integer :: i,k

    c1= (/0.d0, 0.d0, 0.d0/)
    c2= (/0.d0, 1.d0, 0.d0/)
    c3= (/1.d0, 0.d0, 0.d0/)
    print *,"distance d12, d13", distance(c1(1),c2(1)),distance(c1(1),c3(1))
    print *, "angle ",angle(c1,c2,c3)
 contains

!**************************************************************************
function distance(c1,c2)
    implicit none
    double precision,dimension(3) :: c1,c2
    double precision :: distance
    integer :: i
    distance = 0.d0
    do i = 1, 3
        distance =distance + (c1(i)-c2(i))**2
    enddo
    distance = SQRT(distance)
end function distance

!**************************************************************************
function angle(c1,c2,c3)
    implicit none
    double precision,dimension(3) :: c1,c2,c3,v12,v13
    double precision :: angle ,distance12,distance13
    integer :: i
    angle =0.d0
    do i = 1,3
        v12(i)= c2(i)- c1(i)
        v13(i)=c3(i) - c1(i)
    enddo
    
    distance12=distance(c1,c2)
    distance13=distance(c1,c3)
    do i =1,3
        angle =angle + v12(i) * v13(i)
    enddo
    angle = angle /(distance12*distance13)
    angle = acos(angle)
 end function angle
!**************************************************************************
!*******************************************************************
  subroutine ROTATION
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
 end program main
 
