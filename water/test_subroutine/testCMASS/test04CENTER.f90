program main
    implicit none
    real*8,dimension(3,5) :: C,M
    real*8 :: R, Wmass,Cmass(3)
    integer :: i,j,Nsite
    do i = 1, 5 
        do j =1,3
            C(j,i) = dble(i) - 3.D0
            enddo
            enddo
    M=1.D0
    Nsite = 5
    call DISTANCE(C(1,1),C(1,2),R)
    print *, "R=",R
    call CENTER_OF_MASS(C,M,Nsite,Cmass,wmass)
    print *,"Cmass=",Cmass
    print *,"wmass=",wmass
end program main
!*******************************************************************
subroutine DISTANCE(C1,C2,R)
    implicit none
    integer :: i
    real*8 :: C1(3),C2(3),R
    R = 0.D00
    do i=1,3
        R=R+(C1(i) - C2(i))**2
     enddo
    R = sqrt(R)
end subroutine DISTANCE
!*******************************************************************
subroutine CENTER_OF_MASS(C,M,N,Cmass,wmass)
    implicit none
    integer :: N,i,j
    real*8 :: C(3,N),M(N),Cmass(3)
    real*8 :: Wmass
    Cmass = 0.D0
    wmass = 0.D0
    do i =1,N
        do j=1,3
            Cmass(j) = Cmass(j) + C(j,i) * M(i)
        enddo
        wmass = wmass + M(i) 
        enddo
    Cmass = Cmass / wmass
end subroutine CENTER_OF_MASS

!*******************************************************************
