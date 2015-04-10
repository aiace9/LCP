program integration_test
	use omp_lib
	use integration
	implicit none
	integer, parameter :: dp = selected_real_kind(14,200)
	real(kind=dp) :: a,b,int_true
	real(kind=dp) ::  int_stimated, int_stimated2, sigma1, sigma2
	integer :: n

    a = 0.d0
    b = 1.d0
    int_true = 0.746824
    n = 100000
    !$ print*, ('compiledi with openmp')
    print 10, int_true
 10 format("TI", es22.14)

    int_stimated = MC(f, a, b, n, sigma1)

    print 11, int_stimated, sigma1
 11 format("MC", es22.14, " sigma", es22.14)

    int_stimated2 = MC_imp(f, a, b, n, sigma2)

    print 12, int_stimated2, sigma2
 12 format("IS", es22.14, " sigma", es22.14)


	
contains

    real(kind=dp) function f(x)
        implicit none
        real(kind=dp), intent(in) :: x 
        
        f = exp(- x**2)
    end function f

end program integration_test
