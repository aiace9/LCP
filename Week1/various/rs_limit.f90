      program rs_limit
!
! determines the machine precision (eps) in SINGLE PRECISION
!
! the MACHINE PRECISION is the MAX POSITIVE REAL NUMBER eps such that,
! added to the "unity" (one) stored in the memory, does not affect such number,
! i.e. one = one + eps
!
! In order to determine "eps" we try adding to "one" 
! a number smaller and smaller up to reach that condition
!
      real :: eps, factor, one     
      integer :: i, N              
!
!
      eps = 1.0                   
      factor = 2.0                
      N = 100                     
      do i=1,N                       !\   DO LOOP
         eps = eps/factor            ! |  add to "one" a quantity "eps"
         one = 1.0 + eps             !  > smaller and smaller and compare
         print*,  i, one, eps        ! |  "one" with "one + eps"
      end do                         !/

      print*,' using the intrinsic function epsilon(x):' , epsilon(one)

      stop
      end program rs_limit
      
