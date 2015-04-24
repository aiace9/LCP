program integration_test
	use omp_lib
	use integration
	implicit none
	integer, parameter :: dp = selected_real_kind(14,200)
	real(kind=dp) :: a,b,c,d,int_true
	real(kind=dp) ::  int_stimated, int_stimated2, sigma1, sigma2
	integer :: n

    a = 0.d0
    b = 1.0d0
    c = 0.d0
    d = 4.d0
    int_true = acos(-1.0_dp)
    n = 1d6
    !$ print*, ('compiled with openmp')
    print 10, int_true
 10 format("TI", es22.14)

    int_stimated = MC(f, a, b, c, d, n, sigma1)

    print 11, int_stimated, sigma1
 11 format("MC", es22.14, " sigma", es22.14)



	
contains

    real(kind=dp) function f(x)
        implicit none
        real(kind=dp), intent(in) :: x 
        
        f = 4 * sqrt(1 - x**2)
    end function f

end program integration_test
