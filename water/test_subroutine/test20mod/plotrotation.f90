module mdparam
        double precision :: pi = acos(-1.0),&
                            c(3) = (/1.d0, 0.d0, 0.d0/)
        integer :: nsite = 1
end module
module variables
        double precision :: c_rotated(3)
end module variables

program main
        use variables
        implicit none
        integer :: n_times
        do n_times = 1,1000
            call rotation
            print * ,c_rotated
        enddo

contains

     subroutine rotation
        use mdparam, only : pi,c,nsite
        use variables, only : c_rotated
        implicit none
        double precision :: phi, psi, theta
        integer :: i_atoms
        call random_number(phi)
        call random_number(psi)
        call random_number(theta)
        phi =2.d0 *phi*pi
        psi =2.d0 *psi*pi
        theta = theta*pi !0<theta<pi
        do i_atoms = 1, nsite
           c_rotated(1) = (cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi))*c(1) &
                             + (cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi))*c(2) &
                             + (sin(psi)*sin(theta))*c(3)
           c_rotated(2) = (-sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi))*c(1) &
                             + (-sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi))*c(2) &
                             + (cos(psi)*sin(theta))*c(3)
           c_rotated(3) = (sin(theta)*sin(phi))*c(1) &
                             + (-sin(theta)*cos(phi))*c(2) &
                             + cos(theta)*c(3)
        end do
     end subroutine rotation

end program main
