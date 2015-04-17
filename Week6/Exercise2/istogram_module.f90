module istogram_module

	implicit none
	
	private :: dp, hist, a, b, step, nbin


	public :: push_ist, init_ist, save_data_ist

	!privat declaretions
	integer, parameter :: dp = selected_real_kind(14,200)
	integer, dimension(:), allocatable :: hist
	save
	real(kind = dp) :: a, b, step
	save
	integer :: nbin
	save


contains
	subroutine init(min,max,delta)
		real(kind = dp), intent(in) :: min,max,delta
		integer :: err

		b = max
	 	a = min
		step = delta
		nbin= int((b-a)/step)
		allocate(hist(bin), stat=err)
		if (err /= 0) print *, "hist: Allocation request denied"
		hist = 0		
	end subroutine init

	subroutine push_ist(x)
		real(kind = dp), intent(in) :: x
		integer :: ibin

		ibin = floor ( x / step) + 1

		if (ibin <= nbin )then
      		hist(ibin) = hist(ibin) + 1
	    end if
	end subroutine push_ist

	subroutine save_data()
		integer :: myunit, ibin, points
		integer :: ios, err

		open(newunit=myunit, file='ist.dat', iostat=ios, status="unknown", action="write")
		if ( ios /= 0 ) stop "Error opening file ist.dat"

		points = !calculete the total points
		
		do ibin= 1 ,nbin
    		write(unit=1,fmt=*)(ibin-0.5)*step - min ,hist(ibin)/float(points)/step
  		end do

		
		if (allocated(hist)) deallocate(hist, stat=err)
		if (err /= 0) print *, "hist: Deallocation request denied"

		close(unit=myunit, iostat=ios)
		if ( ios /= 0 ) stop "Error closing file unit", myunit
	end subroutine save_data

end module istogram_module
