program istogram_module_test
	 !where is indicated "valid only in this case" is because
	 !there is hidden theorical calculation behind a value that is done for THIS distribution
	use	istogram_module
	
	implicit none
	integer, parameter :: dp = selected_real_kind(14,200)

	integer :: N, points
	integer :: i
	real(kind = dp), dimension(:), allocatable :: rnd
	real(kind = dp) :: x, med, med2, sigma2, z4, z2, z

	!error flags
	integer :: err

	!histogram
	call init_ist(0d0,5d2,1d-4)
	
	!input
	N=500
	points =100000

	!internal initialization
	med = 0
	med2 = 0
	sigma2 = 0

	allocate(rnd(N), stat=err)
	if (err /= 0) print *, "rnd: Allocation request denied"

	!points generation
	do i = 1, points, 1
		call random_number(rnd)
		x = sum(rnd, dim=1)

		call push_ist(x)

	    !numerical check
	    med = med + x
	    med2 = med2 + x**2
	    z = (x - N *0.5) / (N * 1/12.0) !valid only in this case
	    z4 = z4 + z**4
	    z2 = z2 + z**2
	    end do

	!histogram writinig
	call save_data_ist()

  	!more numerical stuff

  	med = med/points
  	med2 = med2/points
  	z2= ( z2 / N )**2
  	z4= ( z4 / N )
  	sigma2 = med2 - med**2

  	print*, 'final check'
  	print*, '---'
  	print*, 'numerical mean value = ', med
  	print*, 'analitical mean value= ', N * 0.5 !valid only in this case
  	print*, '---'
  	print*, 'numerical sigma2 = ', sigma2
  	print*, 'analitical sigma2= ', N * 1/12.0!valid only in this case
  	print*, '---'
  	print*, 'those number must be similar for some theorical reasons'
  	print*, z4
  	print*, 3 * z2
  	
	
	if (allocated(rnd)) deallocate(rnd, stat=err)
	if (err /= 0) print *, "rnd: Deallocation request denied"
	

end program istogram_module_test
