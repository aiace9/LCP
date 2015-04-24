module autocorrelation

	implicit none
	
	public :: autocorrelation_push, autocorrelation_init, autocorrelation_end
	
	private :: autovector

	integer parameter :: dp = selected_real_kind(14,200)
	real(dp), dimension(:), allocatable,save :: autovector_old, autovector
	real(dp), save :: x,x2,xij
	integer, save:: counter,lenght
	
	subroutine autocorrelation_init(l)
		integer, intent(in) :: l
		integer :: err
		lenght = l
		allocate(autovector, stat=err)
		if (err /= 0) print *, "autovector: Allocation request denied"

		allocate(autovector_old, stat=err)
		if (err /= 0) print *, "autovector_old: Allocation request denied"
		
		autovector = 0
		autovector_old =0
	end subroutine autocorrelation_init

	subroutine autocorrelation_push(val)
		real(dp), intent(in) :: val

		autovector_old=autovector
		do i = 1, lenght-1,1
			autovector(i)= autovector_old(i+1)
		end do
		autovector(i+1) = val

		counter = counter + 1
		if (counter >= lenght) then
			x = x + val
			x2 = x2 + val
			xij = xij + autovector(1) * autovector(i+1)
		endif

	end subroutine autocorrelation

	subroutine autocorrelation_pull(auto)
		real(dp), intent(out) :: auto
		integer :: supp
		supp = counter - lenght
		if ( supp > 0 ) then
			auto = (xij/supp -(x/supp)**2)/ (x2/supp-(x/supp)**2)
		else
			auto = NaN
		endif
	end subroutine autocorrelation_pull

	subroutine autocorrelation_end()

		if (allocated(autovector)) deallocate(autovector, stat=err)
		if (err /= 0) print *, "autovector: Deallocation request denied"

	end subroutine autocorrelation_end

end module autocorrelation
