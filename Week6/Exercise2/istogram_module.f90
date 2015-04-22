module istogram_module

	implicit none
	
	private :: dp, hist, a, b, step, nbin


	public :: push_ist, init_ist, save_data_ist

	!privat declaretions
	integer, parameter :: dp = selected_real_kind(14,200)
	
	integer, save, dimension(:), allocatable :: hist
	real(kind = dp), save :: a, b, step
	integer,save :: nbin


contains
	subroutine init_ist(min,max,delta)
		real(kind = dp), intent(in) :: min,max,delta
		integer :: err

		b = max
	 	a = min
		step = delta
		nbin= int((b-a)/step)
		allocate(hist(nbin), stat=err)
		if (err /= 0) print *, "hist: Allocation request denied"
		hist = 0		
	end subroutine init_ist

	subroutine push_ist(x)
		real(kind = dp), intent(in) :: x
		integer :: ibin
		logical :: debug = .false.

		ibin = floor ( x / step) + 1

		if (ibin <= nbin )then
      		hist(ibin) = hist(ibin) + 1
	    end if

	    if ( debug ) print*, 'D: hist(ibin)', hist(ibin)
	end subroutine push_ist

	subroutine save_data_ist()
		integer :: myunit, ibin, points
		integer :: ios, err
		logical :: debug = .false.

		open(newunit=myunit, file='ist.dat', iostat=ios, status="unknown", action="write")
		if ( ios /= 0 ) stop "Error opening file ist.dat"

		points =  sum(hist, dim=1)!calculete the total points
		if ( debug ) print*, 'D: number of points', points

		
		do ibin= 1 ,nbin
    		write(unit=myunit,fmt=*)(ibin-0.5)*step - a ,hist(ibin)/float(points)/step
  		end do

		
		if (allocated(hist)) deallocate(hist, stat=err)
		if (err /= 0) print *, "hist: Deallocation request denied"

		close(unit=myunit, iostat=ios)
		if ( ios /= 0 ) stop "Error closing file ist.dat"
	end subroutine save_data_ist

end module istogram_module
