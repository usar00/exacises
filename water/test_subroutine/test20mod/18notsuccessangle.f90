module mdparam

     implicit none
     integer, parameter :: n =64, nsite = 3
     double precision :: req = 0.94d0, & !a
               thetaeq = 104.0d0, & !deg
               kbond = 547.5d0, & !(kcal/mol)ang^-2
               ktheta = 49.9d0, & !(kcal/mol)rad^-2
               a = 650000d0, & !(kcal/mol)ang^12
               b = 625.47d0, & !(kcal/mol)ang^6
               el = 1.6021766208d-19, & !c
               dt = 1.0d-15, & !sec
               temp0 = 373d0, & !k
               eps = 8.85418782d-12, & !permittivity of free space,(c^2)(j^-1)(m^-1)
               dens = 1.0 !g/cm^3
     double precision, parameter ::  bolz = 1.380658d-23, &
                          avo = 6.0221367d+23, &
                          bohr = 0.529177249d0, &
                          hartree = 4.3597482d-18, &
                          autime = 4.134136844d+16, &
                          emass = 9.1093897d-28, &
                          cal = 4.814d3, &
                          pi = 3.14159265358979d0
     double precision, dimension(nsite) :: mass = (/16.0d0, 1.0d0, 1.0d0/), & !g/mol
                                 q = (/-0.82d0, 0.41d0, 0.41d0/), & !e
                                 c_cmass
     double precision, dimension(3, nsite) :: c = reshape( &
                                    (/0.000000d0, 0.000000d0, 0.000000d0, &
                                      0.940000d0, 0.000000d0, 0.000000d0, &
                                      -0.226567d0, 0.912287d0, 0.000000d0/), (/3, nsite/)),&
                                      c_rotated
     double precision :: rbox, rcut, wmass, kc
     integer :: maxstp = 10000, &
                nprint = 10, &
                cprint = 500, &
                nequil = 10000, &
                ninterval = 100
     integer :: iout = 1, iiout = 2
     end module mdparam
  module variables
     use mdparam
     implicit none
     double precision :: e, t, u, uintera, uinter, temp
     double precision, dimension(3, nsite, n) :: v, r, f, fintera, finter
  end module variables

!*******************************************************************
!*******************************************************************
!*******************************************************************
    module module_dist_angle_cms
    contains
!*******************************************************************
   function distance(c1, c2)
     implicit none
     integer :: i
     double precision :: c1(3), c2(3),distance
     distance = 0.d00
     do i = 1, 3
        distance = distance + (c1(i) - c2(i))**2.0d00
     enddo
     distance = sqrt(distance)
  end function distance
!*******************************************************************
  subroutine center_of_mass(c, m, nsite, c_cmass, wmass)
     implicit none
     integer :: nsite, i_xyz, j_nsite
     double precision :: c(3, nsite), m(nsite), c_cmass(3)
     double precision, optional :: wmass
     c_cmass = 0.d0
     wmass = 0.d0
     do i_xyz = 1,3 
        do j_nsite = 1,nsite
           c_cmass(i_xyz) = c_cmass(i_xyz) + c(i_xyz, j_nsite)*m(i_xyz)
        enddo
        wmass = wmass + m(i_xyz)
     enddo
     c_cmass = c_cmass/wmass
     print *,"c_cmass",c_cmass
  end subroutine center_of_mass
!*******************************************************************
  end module module_dist_angle_cms
!*******************************************************************
!*******************************************************************
!*******************************************************************
program main
    call init
    call forcetest
end program main 
!*******************************************************************
!*******************************************************************
!*******************************************************************
!    use mdparam , only :c,mass,nsite,c_cmass,wmass
!    use variables, only :c_rotated
!    use module_dist_angle_cms,only :center_of_mass
!    implicit none
!    integer :: i,j
!    call center_of_mass(c,mass,nsite,c_cmass,wmass)
!    do i=1,nsite
!        do j=1,3
!            c(j,i) = c(j,i) -c_cmass(j) 
!        enddo
!    enddo
!    c_rotated = c
!    print *,"c_rotated angle", angle(c_rotated(1,1),c_rotated(1,2),c_rotated(1,3))
!    do i=1,3
!        print *,c_rotated(:,i)
!    enddo
!    call rotation
!    print *,"c_rotated1 angle", angle(c_rotated(1,1),c_rotated(1,2),c_rotated(1,3))
!    do i=1,3
!        print *,c_rotated(:,i)
!    enddo
!    call rotation
!    print *,"c_rotated2=",angle(c_rotated(1,1),c_rotated(1,2),c_rotated(1,3))
!    do i=1,3
!        print *,c_rotated(:,i)
!    enddo
!    call rotation
!    print *,"c_rotated3=",angle(c_rotated(1,1),c_rotated(1,2),c_rotated(1,3))
!    do i=1,3
!        print *,c_rotated(:,i)
!    enddo

subroutine init
     use mdparam
     use variables
     use module_dist_angle_cms, only :center_of_mass 
     implicit none
     integer :: number_of_dr_in_x, number_of_dr_in_y, number_of_dr_in_z, i, j, max_number_of_dr
     integer :: j_nsite,i_xyz
     double precision :: dr, small_dr, ramdom_from_0_to_1
!atomic units
     dens = (dens/wmass)*avo*(bohr*1.0d-8)**3
     mass = (mass/avo)/emass
     wmass = (wmass/avo)/emass
     rbox = (n/dens)**(1.0d0/3.0d0)
     rcut = rbox/2.0d0
     c = c/bohr
     dt = dt*autime
     req = req/bohr
     thetaeq = thetaeq/(180.0d0/pi)
     kbond = kbond*(bohr**2)*cal/hartree/avo
     ktheta = ktheta*cal/hartree/avo
     a = a/(bohr**12)*cal/hartree/avo
     b = b/(bohr**6)*cal/hartree/avo

     max_number_of_dr = int((n - 1)**(1.0d0/3.0d0)) + 1
     dr = rbox/max_number_of_dr
     small_dr = dr*0.2d0
     j = 0
     call center_of_mass(c,mass,nsite,c_cmass,wmass)
     print *,"c_cmass",c_cmass
     do i_xyz= 1,3
         do j_nsite =1,nsite
            c(i_xyz,j_nsite) = c(i_xyz,j_nsite) - c_cmass(i_xyz)
         enddo
     enddo
     print *,"c",c
     do number_of_dr_in_x = 1, max_number_of_dr
        do number_of_dr_in_y = 1, max_number_of_dr
           do number_of_dr_in_z = 1, max_number_of_dr
              j = j + 1
              if (j > n) exit
              call rotation
              call random_number(ramdom_from_0_to_1)
              r(1, 1, j) = c_rotated(1, 1) + number_of_dr_in_x*dr + small_dr*(2.0d0*ramdom_from_0_to_1 - 1.0d0)
              r(2, 1, j) = c_rotated(2, 1) + number_of_dr_in_y*dr + small_dr*(2.0d0*ramdom_from_0_to_1 - 1.0d0)
              r(3, 1, j) = c_rotated(3, 1) + number_of_dr_in_z*dr + small_dr*(2.0d0*ramdom_from_0_to_1 - 1.0d0)
              r(1, 2, j) = c_rotated(1, 2) + number_of_dr_in_x*dr + small_dr*(2.0d0*ramdom_from_0_to_1 - 1.0d0)
              r(2, 2, j) = c_rotated(2, 2) + number_of_dr_in_y*dr + small_dr*(2.0d0*ramdom_from_0_to_1 - 1.0d0)
              r(3, 2, j) = c_rotated(3, 2) + number_of_dr_in_z*dr + small_dr*(2.0d0*ramdom_from_0_to_1 - 1.0d0)
              r(1, 3, j) = c_rotated(1, 3) + number_of_dr_in_x*dr + small_dr*(2.0d0*ramdom_from_0_to_1 - 1.0d0)
              r(2, 3, j) = c_rotated(2, 3) + number_of_dr_in_y*dr + small_dr*(2.0d0*ramdom_from_0_to_1 - 1.0d0)
              r(3, 3, j) = c_rotated(3, 3) + number_of_dr_in_z*dr + small_dr*(2.0d0*ramdom_from_0_to_1 - 1.0d0)
!        print *,dr,number_of_dr_in_x*dr,number_of_dr_in_y*dr,number_of_dr_in_z*dr
           end do
        end do
     end do
     v = 0.0d0
     f = 0.0d0
     fintera = 0.0d0
     finter = 0.0d0
     u = 0.0d0
     uintera = 0.0d0
     uinter = 0.0d0

    contains

!*******************************************************************
     subroutine rotation
        use mdparam, only : pi,c,nsite
        use variables, only : c_rotated
        implicit none
        double precision :: phi, psi, theta
        integer :: i
        call random_number(phi)
        call random_number(psi)
        call random_number(theta)
        phi = phi*pi
        psi =psi*pi
        theta = theta*pi !0<theta<pi
        do i = 1, nsite
           c_rotated(1, i) = (cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi))*c(1, i) &
                             + (cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi))*c(2, i) &
                             + (sin(psi)*sin(theta))*c(3, i)
           c_rotated(2, i) = (-sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi))*c(1, i) &
                             + (-sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi))*c(2, i) &
                             + (cos(psi)*sin(theta))*c(3, i)
           c_rotated(3, i) = (sin(theta)*sin(phi))*c(1, i) &
                             + (-sin(theta)*cos(phi))*c(2, i) &
                             + cos(theta)*c(3, i)
        end do
     end subroutine rotation
!*******************************************************************
     function angle(c1, c2, c3)
        use module_dist_angle_cms, only : distance 
        implicit none
        double precision,intent(in) ::c1(3), c2(3), c3(3)
        double precision :: angle 
        double precision :: costheta,&
                   r12(3), r13(3), &
                  r12_distance, r13_distance
        integer :: i_xyz
        costheta = 0.0d00
        r12(3) = 0.0d00
        r13(3) = 0.0d00
        do i_xyz = 1, 3 !i represent xyz
           r12(i_xyz) = c2(i_xyz) - c1(i_xyz)
           r13(i_xyz) = c3(i_xyz) - c1(i_xyz)
        enddo
        r12_distance = distance(c1(1),c2(1))
        r13_distance = distance(c1(1),c3(1))
!******* caluculate triangle function of theta
        do i_xyz = 1, 3
           costheta = costheta + r12(i_xyz)*r13(i_xyz)
        enddo
        costheta = costheta/(r12_distance*r13_distance)
        !if(costheta > 1.0d00) print *,"costheta has an erroer",costheta
        angle = acos(costheta)
        end function angle
!*******************************************************************
end subroutine init
!*******************************************************************
!*******************************************************************
!*******************************************************************
  subroutine force
     use mdparam
     use variables
     implicit none
     f = 0.0d0
     fintera = 0.0d0
     finter = 0.0d0
     u = 0.0d0
     uintera = 0.0d0
     uinter = 0.0d0
     call forcr_intra 
!     call force_inter(finter,uinter) !intermolcule
     f = fintera + finter
     u = uintera + uinter
  end subroutine force

!*******************************************************************
!*******************************************************************
!*******************************************************************
  subroutine forcr_intra
     use mdparam
     use variables
     implicit none
     integer ::number_of_molecule,xyz,k
     real*8, dimension(3, nsite) ::dtdx123
     real*8 :: r1(3), r2(3),r3(3),a213, r12,r13, drdx12(3, 2),drdx13(3, 2)
        do number_of_molecule=1,n
            do xyz=1,3
                r1(xyz)= r(xyz,1,number_of_molecule)
                r2(xyz)= r(xyz,2,number_of_molecule)
                r3(xyz)= r(xyz,3,number_of_molecule)
            enddo
            call bond_distance_differencial(r1(1),r2(1),r12,drdx12(1,1)) 
            call bond_distance_differencial(r1(1),r3(1),r13,drdx13(1,1)) 
            call angle_differencial(r1(1),r2(1),r3(1),a213,dtdx123(1,1)) 
            do xyz=1,3
                    fintera(xyz,1,number_of_molecule)=-(2.0d0*ktheta*(a213 - thetaeq)*dtdx123(xyz,1)& 
                                    +2.0d0*kbond*(r12 - req)*drdx12(xyz,1)&
                                    +2.0d0*kbond*(r13 - req)*drdx13(xyz,1))
                    fintera(xyz,2,number_of_molecule)=-(2.0d0*ktheta*(a213 - thetaeq)*dtdx123(xyz,2)& 
                                    +2.0d0*kbond*(r12 - req)*drdx12(xyz,2))
                    fintera(xyz,3,number_of_molecule)=-(2.0d0*ktheta*(a213 - thetaeq)*dtdx123(xyz,3)& 
                                    +2.0d0*kbond*(r13 - req)*drdx13(xyz,2))
            end do
            uintera = uintera + ktheta*(a213- thetaeq)**2 +kbond*(r12 - req)**2 + kbond*(r13 - req)**2 
        end do
  contains
!*******************************************************************
     subroutine bond_distance_differencial(c1, c2, r, drdx)
        use module_dist_angle_cms, only : distance
        implicit double precision(a-h,o-z),integer(i-n)
        double precision,intent(in) :: c1(3), c2(3) 
        double precision,intent(out) :: r, drdx(3, 2)
        drdx(3, 2) = 0.0d00
        r = distance(c1(1),c2(1)) 
        do i_xyz = 1, 3
           drdx(i_xyz, 1) = -c1(i_xyz)/r
           drdx(i_xyz, 2) = c2(i_xyz)/r
        enddo
     end subroutine bond_distance_differencial
!*******************************************************************
     subroutine angle_differencial(c1, c2, c3, theta, dtdx)
        use module_dist_angle_cms, only : distance
        implicit double precision(a-h,o-z),integer(i-n)
        double precision,intent(in) ::c1(3), c2(3), c3(3)
        double precision,intent(out) :: theta, dtdx(3, 3)
        double precision :: costheta, sintheta, &
                   r12(3), r13(3), &
                  r12_distance, r13_distance
        costheta = 0.0d00
        dtdx(3, 3) = 0.0d00
        r12(3) = 0.0d00
        r13(3) = 0.0d00
!******* caluculate vector
        do i_xyz = 1, 3 
           r12(i_xyz) = c2(i_xyz) - c1(i_xyz)
           r13(i_xyz) = c3(i_xyz) - c1(i_xyz)
        enddo
        r12_distance = distance(c1(1),c2(1))
        r13_distance = distance(c1(1),c3(1))
!******* caluculate triangle function of theta
        do i_xyz = 1, 3
           costheta = costheta + r12(i_xyz)*r13(i_xyz)
        enddo
        costheta = costheta/(r12_distance*r13_distance)
        sintheta = sqrt(1.0d00 - costheta**2)
        theta = acos(costheta)
!****** caluculate dtdx
        do i_xyz = 1, 3
           dtdx(i_xyz, 2) = -(r13(i_xyz) - costheta*(r12_distance/r13_distance)*r12(i_xyz)) &
                        /(sintheta*r12_distance*r13_distance)
           dtdx(i_xyz, 3) = -(r12(i_xyz) - costheta*(r13_distance/r12_distance)*r13(i_xyz)) &
                        /(sintheta*r12_distance*r13_distance)
           dtdx(i_xyz, 1) = -dtdx(i_xyz, 2) - dtdx(i_xyz, 3)
        enddo
     end subroutine angle_differencial
!*******************************************************************
  end subroutine forcr_intra

!*******************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************
  subroutine forcetest
     use mdparam
     use variables
     implicit none
     integer :: i, j, k
     double precision :: d = 1.0d-8, u0, up2, up, un, un2, r0
     double precision, dimension(3, nsite, n) :: f0, f1
     double precision :: davg = 0.0d0

     call force
     u0 = u
     f0 = f
     do i = 1, n
        do j = 1, nsite
           do k = 1, 3
              r0 = r(k, j, i)
              r(k, j, i) = r0 + d + d
              call force
              up2 = u
              r(k, j, i) = r0 + d
              call force
              up = u
              r(k, j, i) = r0 - d
              call force
              un = u
              r(k, j, i) = r0 - d - d
              call force
              un2 = u
              f1(k, j, i) = -(8.0d0*(up - un) - (up2 - un2))/(d*1.2d1)
              r(k, j, i) = r0
              davg = davg + (f1(k, j, i) - f0(k, j, i))**2
              print *,"i,j,k=",i,j,k,"f0,f1=",f0(k,j,i),f1(k,j,i)
           end do
!            print '(i5,6e14.6)', i, (f0(k, j, i), k=1, 3), (f1(k, j, i) - f0(k, j, i), k=1, 3)
        end do
     end do
     davg = sqrt(davg/(3*n))
     print *, 'davg=', davg
        do j=1,nsite
            print *, "r",r(:,j,1)
        enddo
     return
  end subroutine forcetest
!*******************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************
