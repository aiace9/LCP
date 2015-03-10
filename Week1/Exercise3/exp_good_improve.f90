program expgood
! variable declaration:
  !     x
  !     accuracy limit: min
  !
  implicit none
  integer, parameter :: dp=selected_real_kind(13)
  real (kind = dp) ::  element, sum, x, min, inc, max
  integer :: n
  open(unit=7,file="exp-good-imp.dat") !,position="append",action="write")
  write(unit=7,fmt=*) "x, n, sum, exp(-x), abs(sum-exp(-x))/sum" ! FMT MEANS FORMAT, * = free
  
  open(unit=8,file="exp-good-util-imp.dat") !,position="append",action="write")
  write(unit=8,fmt=*) "n, sum, element, exp(-x), " ! FMT MEANS FORMAT, * = free
  !
  ! execute
  !
  min=1.e-10
  
  x=200.0_dp
  inc=0.1_dp
  max = 200.0_dp
  !write(6,*)' enter x:'  ! or write (*,*) or write(6,*) or write(unit=6,fmt=*)
  !read(5,*) x  ! or read(*,*) since unit=5 is standard input , i.e. keyboard
 
  
  do 
	  element = 1
	  sum = 1
	  if (x>max) exit
	  
	  do  n=1, 1000000
		 element = element*(x)/n
     print*, sum,element
		 sum = sum + element
     write(unit=8,fmt=*) n, sum, element, element/sum
		 
     if((element <  epsilon(sum)*sum) .and. (sum /= 0)) then
			write(*,*) x, n, 1./sum, exp(-x), abs(1./sum-exp(-x))*sum 
			write(unit=7,fmt=*) x, n, 1./sum, exp(-x), abs(1./sum-exp(-x))*sum 
			exit !or go to 10 . ATTENTION: exit is DIFFERENT FROM cycle
		 endif
	  enddo
	  x = x + inc
  enddo
! 10 continue
  close(7)
  close(8)
  !      stop "data saved in exp-good.dat"      
 stop  
end program expgood





