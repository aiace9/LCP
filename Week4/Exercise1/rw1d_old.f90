! rw1d.f90
! A simple random walk program  in 1D.

program rw1d
  
  implicit none
  integer ::         N  ! number of steps
  integer, dimension(:), allocatable :: seed
  integer :: icount1, icount2, icount_rate, ix, irun, istep, nruns, sizer
  real, dimension(:), allocatable :: rnd          ! array of random numbers
  real :: pleft,jump,delta, deltath
  integer, dimension(:), allocatable :: x_N, x2_N ! sum of deviations and 
                                                  ! squares over the runs
  integer, dimension(:), allocatable :: P_N     ! final positions, sum over runs
  real, dimension(:), allocatable :: P_N_anal, P_N_anal_sto     ! final positions, sum over runs

  integer :: ios
  integer :: i,j

  print *, "Enter number of steps >"
  read *, N
  print *, "Enter number of runs >"
  read *, nruns
  print*, "Enter P left [0-1] >"
  read *, pleft
  call random_seed(sizer)
  allocate(seed(sizer))
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  read(*,*)seed
  call random_seed(put=seed)


  allocate(rnd(N))
  allocate(x_N(N))
  allocate(x2_N(N))
  allocate(P_N(-N:N))
  allocate(P_N_anal(-N:N))
  allocate(P_N_anal_sto(-N-1:N+1))
  x_N  = 0
  x2_N = 0
  P_N  = 0
  
  jump = 1


  !open files
  open(unit=2, file='walk.dat', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file name"

  open(unit=3, file='delta-deltath.dat', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file name"

  open(unit=4, file='delta.dat', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file name"

  !create the anal(cit.)itical distributions P_N
  !initialization:
  P_N_anal = 0.0
  P_N_anal(0) = 1.0
  P_N_anal_sto = 0.0

  do i = 1, N, 1
    !create the storical arry
    do j = -N, N, 1
      P_N_anal_sto(j) = P_N_anal(j)
    end do
    !update the new line
    do j = -N, N, 1
      P_N_anal(j) = (1 / 2.0) * P_N_anal_sto(j-1) +(1 / 2.0) * P_N_anal_sto(j+1)
    end do
  end do


  do irun = 1, nruns
    ix = 0 ! initial position of each run
    call random_number(rnd) ! get a sequence of random numbers
    do istep = 1, N
      if (rnd(istep) < pleft) then ! random move
        ix = ix - jump ! left
      else
        ix = ix + jump ! right
      end if
      x_N (istep) = x_N (istep) + ix !has the position of all the walkers!!
      x2_N(istep) = x2_N(istep) + ix**2
    end do
    
    P_N(ix) = P_N(ix) + 1 ! accumulate (only for istep = N)
  end do

  !write walk.dat to see only the positions occupied by one walker you must chose nruns = 1
  !and then plot the first and the second column.
  do istep = 1, N, 1
    write(unit=2, fmt=*, iostat=ios) istep, real(x_N(istep))/nruns, real(x2_N(istep))/nruns, &
     (real(x2_N(istep))/nruns - (real(x_N(istep))/nruns)**2) 
    if ( ios /= 0 ) stop "Write error in file unit 2" 
  end do

 
  !write delta-th e delta. dat
  do istep = 1, N, 1
    delta = (real(x2_N(istep))/nruns - (real(x_N(istep))/nruns)**2)
    deltath = (4 * (1-pleft) * pleft * real(istep) * jump**2) - 1.0

    write(unit=3, fmt=*, iostat=ios) istep, abs(delta/deltath - 1.0)
    if ( ios /= 0 ) stop "Write error in file unit 2"

    write(unit=4, fmt=*, iostat=ios) istep, delta
    if ( ios /= 0 ) stop "Write error in file unit 2" 
  end do

  !chek how the simulation wents.....
  print*, " "
  print*, "significant results"
  print*,"# N=",N,"  nruns=",nruns
  print*,"# <x_N>   = ",real(x_N(N))/nruns, "th value =", N*((1.0-pleft)-pleft)*jump 
  print*,"# <x^2_N> = ",real(x2_N(N))/nruns, "th value =", N * jump**2
  print*,"# <x^2_N> -  <x_N>^2 = ",real(x2_N(N))/nruns - (real(x_N(N))/nruns)**2 


  
  open(1,file="P_N.dat",STATUS="REPLACE", ACTION="WRITE")
  write(1,*)"# N=",N,"  nruns=",nruns
  write(1,*)"# <x_N>   = ",real(x_N(N))/nruns 
  write(1,*)"# <x^2_N> = ",real(x2_N(N))/nruns
  write(1,*)"# <x^2_N> - <x_N>^2 = ",real(x2_N(N))/nruns-(real(x_N(N))/nruns)**2 
  write(1,*)" "

  write(1,*)"# N, mean deviations, mean squared deviations, sigma^2"
  do ix = - N, N
  write(1,*)ix,real(P_N(ix))/nruns, P_N_anal(ix)
  end do

  !closing file
  close(unit=1, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 1"

  close(unit=2, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 2"
  
  close(unit=3, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 3"

  close(unit=4, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 3"




  deallocate(rnd, x_N, x2_N, P_N, seed)

  stop
end program rw1d