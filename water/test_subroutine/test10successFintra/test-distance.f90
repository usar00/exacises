program main
    implicit none
    double precision :: C1(3),C2(3),d
    C1= 0.D0
    C2 =1.D0
    d = DISTANCE(C1,C2)
    print *, d

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
end program main
