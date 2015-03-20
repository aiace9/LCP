program strange_disy
  implicit none
  real :: delta,x
  integer :: i,n,nbin,ibin, sizer
  integer :: ios
  integer, dimension(:), allocatable :: histo, seed
  logical :: debug = .false.
  print*, " Generates random numbers distributed as 1/pi - (1 - x^2)^1/2"
  call random_seed(sizer)
  allocate(seed(sizer))
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  read(*,*)seed
  call random_seed(put=seed)
  print *," length of the sequence >"
  read *, n	
  print *," Insert number of bins in the histogram MUST BE EVEN>"
  read *, nbin
  
  delta = 2./nbin
  allocate (histo(nbin))
  histo = 0

  do i = 1,n
    call strangeDist(x)
    ibin = int (x/delta) + 1

    if ( debug ) print*, ibin
    
    if (ibin <= nbin/2 .and. ibin >= -nbin/2)then
      histo( ibin + nbin / 2 ) = histo( ibin + nbin / 2) + 1
      if ( debug ) print*, histo( ibin + nbin / 2)
    end if

  end do
  
  open (unit=7,file="strangedist.dat",status="replace",action="write")
  
  do ibin= 1 ,nbin
     write(unit=7,fmt=*)(ibin-nbin/2-0.5)*delta,histo(ibin)
  end do

  close(unit=7, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 7"
  
contains

  subroutine strangeDist(x)
    real , parameter:: pi = 3.14
    REAL, intent (out) :: x
    REAL :: r
    
    call random_number(r)

    x = sin( pi * (2 * r - 1 ) )
    
  END subroutine strangeDist


end program strange_disy
