program main
     implicit none
     integer, parameter :: N = 64, Nsite = 3
     integer :: i,j
     real*8 :: WMASS
     real*8, dimension(Nsite) :: MASS = (/16.0, 1.0, 1.0/), & !g/mol
                                 Q = (/-0.82, 0.41, 0.41/), & !e
                                 C_cmass
     real*8, dimension(3, Nsite) :: C_original = reshape( &
                                    (/0.000000, 0.000000, 0.000000, &
                                      0.940000, 0.000000, 0.000000, &
                                      -0.226567, 0.912287, 0.000000/), (/3, Nsite/))
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
     print *,C_cmass
end program main
