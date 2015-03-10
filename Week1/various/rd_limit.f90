      program rd_limit
!                                                                                
! determines the machine precision (eps) in DOUBLE PRECISION 
!                                                                               
! the MACHINE PRECISION is the MAX POSITIVE REAL NUMBER eps such that,          
! added to the "unity" (one) stored in the memory, does not affect such number, 
! i.e. one = one + eps                                                          
!                                                                               
! In order to determine "eps" we try adding to "one"                            
! a number smaller and smaller up to reach that condition                        
!
      integer, parameter :: dp=selected_real_kind(13)
      real(kind=dp) :: eps, factor, one     
      integer :: i, N                
!
      eps = 1.e-25_dp
	print*,eps             
      factor = 2.0_dp                
      N = 100                    
      do i=1,N                       !\   DO LOOP
         eps = eps/factor            ! |  add to  "one" a quantity "eps"
         one = 1.0_dp + eps             !  > smaller and smaller and compare
         print*,  i, one, eps        ! |  "one" to "one + eps"
      end do                         !/

      print*,' using the intrinsic function epsilon(x):' , epsilon(eps)

      stop
      end program rd_limit
      
