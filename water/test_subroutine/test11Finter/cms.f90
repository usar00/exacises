program main
    integer :: Nsite = 3, i, k
    real :: C(3,3),M(3),Cmass(3),wmass
    do i = 1,Nsite
        do k= 1,3
        C(k,i)= k+i
        enddo
        M(i)=1
    enddo
    call CENTER_OF_MASS(C(1,1),M(1),Nsite,Cmass(1),wmass)
end program main
!*******************************************************************
subroutine CENTER_OF_MASS(C,M,N,Cmass,wmass)
    implicit none
    integer :: N,i,j
    real*8 :: C(3,N),M(N),Cmass(3)
    real*8,optional :: Wmass
    Cmass = 0.D0
    wmass = 0.D0
    do i =1,N
        do j=1,3
            Cmass(j) = Cmass(j)+ C(j,i) * M(i)
        enddo
        wmass = wmass + M(i) 
        enddo
    Cmass = Cmass / wmass
end subroutine CENTER_OF_MASS
!*******************************************************************
