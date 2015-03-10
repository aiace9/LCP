      program rd_under_over
!
! test UNDERFLOW and OVERFLOW limits
! using improved precision for REAL numbers 
! i.e. DOUBLE PRECISION for those architectures where
! the standard is 32 bits (= 4 bytes)
!
      integer, parameter :: dp=selected_real_kind(13)
      real(kind=dp) :: under, over, factor     ! declaration of variables
      integer :: i, N                   ! and constants
!
!
      under = 1.0_dp                   ! inizialize variables and constants
      over  = 1.0_dp                   !
      factor = 2.0_dp                  ! 
      N      = 1100
      do i=1,N                       !\   
         under = under/factor        ! |  DO LOOP
         over  =  over*factor        !  > multiply (divide) iteratively
         print* ,   i, under, over   ! |  in order to obtain
      end do                         !/   a larger (smaller) quantity
      print* ,   dp 
      stop  
      end program rd_under_over
