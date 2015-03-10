      program i_min_max
!
! calculate the minimum and maximum possible integer number 
! for a 32-bit architecture using 4 bytes for integer representation
!
      implicit none
      integer :: imin, imax, iexp
!
      print*," max positive integer:"
      iexp = 1
      imax = 1
      do 
         imax = 2**iexp - 1
         print*, " iexp=",iexp,"     2**iexp - 1=",imax
	 if (imax <= 0 ) go to 10
         imax = 2**iexp
         print*, " iexp=",iexp,"     2**iexp    =",imax
	 if (imax <= 0 ) go to 10
         iexp = iexp + 1
         end do
10	continue
!
      print*, " min negative integer:"
      iexp = 1
      imin = -1
      do 
         imin = - 2**iexp + 1
         print*, " iexp=",iexp,"    -2**iexp + 1=",imin
         if (imin >= 0)exit
         imin = - 2**iexp
         print*, " iexp=",iexp,"    -2**iexp    =",imin
         iexp = iexp + 1
         if (imin >= 0)exit
         end do
!C
      stop
      end program i_min_max
