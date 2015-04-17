program gauss_dist
	 !where is indicated "valid only in this case" is because
	 !there is hidden theorical calculation behind a value that is done for THIS distribution
	implicit none
	integer, parameter :: dp = selected_real_kind(14,200)

	integer :: N, points
	integer :: i, j
	real(kind = dp), dimension(:), allocatable :: rnd
	real(kind = dp) :: x, med, med2, sigma2, z4, z2, z, c

	!histogram
	integer, dimension(:), allocatable :: hist
	real(kind = dp) :: min, max, step
	integer :: nbin, ibin
	integer :: ios

	!error flags
	integer :: err

	!histogram
	max = 1000
	min = 0
	step = 1
	nbin= int((max-min)/step)
	allocate(hist(nbin))
	hist = 0
	open(unit=1, file='hist.dat', iostat=ios, status="unknown", action="write")
	if ( ios /= 0 ) stop "Error opening file name"
	
	!input
	N=500
	points =1000
	c = 1.0
	!internal initialization
	med = 0
	med2 = 0
	sigma2 = 0

	allocate(rnd(N), stat=err)
	if (err /= 0) print *, "rnd: Allocation request denied"

	!points generation
	do i = 1, points, 1
		call random_number(rnd)

			!distribution convertion
			do j = 1, N, 1
				do
					if (rnd(j)==0) then
						call random_number(rnd(j))
					else
						exit
					end if
				end do
			end do

		rnd = -(1.0/c) * log(rnd)				

		x = sum(rnd, dim=1)

		!histogram population
		ibin = floor ( x / step) + 1
		if (ibin <= nbin )then
      		hist(ibin) = hist(ibin) + 1
	    end if

	    !numerical check
	    med = med + x
	    med2 = med2 + x**2
	    z = (x - N * 1/c) / (N * 1/c**2) !valid only in this case
	    z4 = z4 + z**4
	    z2 = z2 + z**2
	    end do

	!histogram writinig
	do ibin= 1 ,nbin
    	write(unit=1,fmt=*)(ibin-0.5)*step - min ,hist(ibin)/float(points)/step
  	end do
  	close(unit=1, iostat=ios)
  	if ( ios /= 0 ) stop "Error closing file unit 1"

  	!more numerical stuff

  	med = med/points
  	med2 = med2/points
  	z2= ( z2 / N )**2
  	z4= ( z4 / N )
  	sigma2 = med2 - med**2

  	print*, 'final check'
  	print*, '---'
  	print*, 'numerical mean value = ', med
  	print*, 'analitical mean value= ', N * 1/c !valid only in this case
  	print*, '---'
  	print*, 'numerical sigma2 = ', sigma2
  	print*, 'analitical sigma2= ', N * 1/c**2!valid only in this case
  	print*, '---'
  	print*, 'those number must be similar for some theorical reasons'
  	print*, z4
  	print*, 3 * z2
  	
	
	if (allocated(rnd)) deallocate(rnd, stat=err)
	if (err /= 0) print *, "rnd: Deallocation request denied"
	

end program gauss_dist
