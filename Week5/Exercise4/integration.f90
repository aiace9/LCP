module integration

	use omp_lib
	implicit none
	public :: MC, MC_imp
	private :: dp
	integer, parameter :: dp = selected_real_kind(14,200)
	
	contains

	function MC (f, a, b, c, d, n, sigma)
		!f is the function to integrate
		!a-b c-d is the area base x high where points are trow
		!n is the number of points
		real(kind = dp) :: MC
		real(kind = dp), external :: f
		real(kind = dp), intent(in) :: a,b,c,d
		integer, intent(in) :: n
		real(kind = dp), intent(out) :: sigma

		logical :: debug
		real(kind = dp) :: rnd(2)
		real(kind = dp) :: med2, med
		integer :: int_value
		integer :: i
		
		!$ call omp_set_num_threads(4)
		debug = .true.

		int_value = 0
		med = 0
		med2= 0
		!$omp parallel do reduction ( + : int_value) privete(rnd(2)) 

		if ( debug )  open(unit=1, file='name.dat')

		do i = 1, n, 1
			call random_number(rnd)
			!range fix
			rnd(1) = rnd(1) * ( b - a ) + a
			rnd(2) = rnd(2) * ( d - c ) + c
			rnd(1) = f(rnd(1))

			if ( rnd(1) >= rnd (2) ) then
				int_value = int_value + 1
				med =  med + rnd(1)
				med2 = med2 + rnd(1)**2
			end if
			if ( debug )write(unit=1, fmt=*) i, (int_value / real(i)) - acos(-1.0_dp)!subtract the real value here 

		end do

		if ( debug ) close(unit=1)
		
		sigma = sqrt(med2/ real(n) + (med / real(n))**2)
		MC = int_value / real(n) *(b-a)*(d-c)

	end function MC

	function MC_imp (f, a, b, n,sigma)
		!important sample with an exponential

		real(kind = dp) :: MC_imp
		real(kind = dp), external :: f
		real(kind = dp), intent(in) :: a,b
		real(kind = dp), intent(out) ::sigma
		integer, intent(in) :: n

		logical :: debug
		real(kind = dp) :: rnd
		real(kind = dp) :: int_value
		real(kind = dp) :: med2, med
		real(kind = dp) :: c
		integer :: i
		debug = .false.
		!$ call omp_set_num_threads(4)

		int_value = 0.0
		med = 0
		med2= 0
		c = 1.0

		if ( debug )  open(unit=1, file='name.dat')

		!$omp parallel do privete(rnd reduction ( + : int_value)
		do i = 1, n, 1
			do
				call random_number(rnd)
				!distribution convertion
				if ( rnd > 0 ) then
					rnd = -(1.0/c) * log(rnd)
					if ( rnd >= a .and. rnd <= b ) exit					
				end if
			end do
			!summ					
			int_value = int_value +f(rnd)/ (c * exp(-rnd))
			med =  med + rnd
			med2 = med2 + rnd**2
			if ( debug )write(unit=1, fmt=*) i, (int_value / real(i)) * (c*(exp(-a) - exp(-b))) !subtract the real value here 
		end do

		if ( debug ) close(unit=1)

		sigma = sqrt(med2/ real(n) + (med / real(n))**2)

		MC_imp = int_value / real(n) * (c*(exp(-a) - exp(-b)))

	end function MC_imp

end module integration
