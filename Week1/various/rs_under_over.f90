      program rs_under_over
!
! test UNDERFLOW and OVERFLOW limits
! using standard precision for REAL numbers 
! which is typically 32 bits (= 4 bytes), i.e. SINGLE PRECISIONE 
!

      real :: under, over, factor   ! declaration of variables
      integer :: i, N               ! and constants
!
!
      under = 1.0                   ! initialize variable and constants
      over  = 1.0                   !
      factor = 2.0                  ! 
      N      = 500
      do i=1,N                       !\   
         under = under/factor        ! |  DO LOOP
         over  =  over*factor        !  > multiply (divide) iteratively
         print* ,   i, under, over   ! |  in order to obtain a progressively
      end do                         !/   larger (smaller) quantity
      stop
      end program rs_under_over
