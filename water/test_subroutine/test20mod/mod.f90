program main
    implicit none
    double precision :: a,b,c
    a = 1.1d0
    b = 0.9d0
    c = mod(a,b)
    print *,"a,b,mod(a,b)",a,b,c
    a = 2.1d0
    b = 0.9d0
    c = mod(a,b)
    print *,"a,b,mod(a,b)",a,b,c
    a = 3.1d0
    b = 0.9d0
    c = mod(a,b)
    print *,"a,b,mod(a,b)",a,b,c
    a = 1.1d0
    b = 0.8d0
    c = mod(a,b)
    print *,"a,b,mod(a,b)",a,b,c
    a = 2.1d0
    b = 0.8d0
    c = mod(a,b)
    print *,"a,b,mod(a,b)",a,b,c
    a = 3.1d0
    b = 0.8d0
    c = mod(a,b)
    print *,"a,b,mod(a,b)",a,b,c
end program main
