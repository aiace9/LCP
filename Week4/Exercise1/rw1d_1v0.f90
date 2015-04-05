! rw1d.f90
! A simple random walk program  in 1D.

program rw1d
  
  implicit none
  integer, parameter :: dp=selected_real_kind(12)
  integer ::         N  ! number of steps
  integer, dimension(:), allocatable :: x_runner, x2_runner, seed
  integer :: icount1, icount2, icount_rate, irun, istep, nruns, sizer, jump
  real (kind = dp), dimension(:), allocatable :: rnd          ! array of random numbers
  integer, dimension(:), allocatable :: ix
  real (kind = dp) :: pleft,delta, deltath
  integer, dimension(:), allocatable :: x_N, x2_N ! sum of deviations and 
                                                  ! squares over the runs
  integer, dimension(:), allocatable :: P_N     ! final positions, sum over runs
  real (kind = dp), dimension(:), allocatable :: P_N_anal, P_N_anal_sto     ! final positions, sum over runs

  integer :: ios
  integer :: i,j

  print *, "Enter number of steps >"
  read *, N
  print *, "Enter number of runner >"
  read *, nruns
  print*, "Enter P left [0-1] >"
  read *, pleft
  call random_seed(sizer)
  allocate(seed(sizer))
  print *,'Here the seed has ',sizer,' components; insert them (or print "/") >'  
  read(*,*)seed
  call random_seed(put=seed)

  !allocates various quanitites
  allocate(rnd(nruns))        !random numbers
  allocate(x_runner(nruns))   !positions vector
  allocate(x2_runner(nruns))  !disps vector
  allocate(ix(nruns))         !move for each runner
  !allocates distribution stuff
  allocate(P_N(-N:N))
  allocate(P_N_anal(-N:N))
  allocate(P_N_anal_sto(-N-1:N+1))
  !inizialization
  x_runner = 0
  x2_runner = 0
  P_N  = 0
  !unit jump! the code is not yet developed for this quantity
  jump = 1


  !open files
  open(unit=1, file='P_N.dat', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file name"

  open(unit=2, file='walks.dat', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file name"

  open(unit=3, file='avgwalk.dat', iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file name"
  

  !create the analitical distributions P_N
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

  !this algorithm work but only becouse the central limit is valid in this case

  !do irun = 1, nruns
  !  ix = 0 ! initial position of each run
  !  call random_number(rnd) ! get a sequence of random numbers
  !  do istep = 1, N
  !    if (rnd(istep) < pleft) then ! random move
  !      ix = ix - jump ! left
  !    else
  !      ix = ix + jump ! right
  !    end if
  !    x_N (istep) = x_N (istep) + ix !has the position of all the walkers!!
  !    x2_N(istep) = x2_N(istep) + ix**2
  !  end do
  !  ------how the vector X_N came out------
  !  irun = 1
  !  1 step       - 2 steps      - .... - n  steps
  !  irun = 2
  !  2 steps      - 4 steps      - .... - 2n steps
  !   .....       -  .....       - .... -  .....
  !  irun = nruns
  !  nruns steps - 2*nruns seps - .... -  nruns*n steps
  !   
  !  ----- this is correct, in the way it will be used i think ----
  !
  !  P_N(ix) = P_N(ix) + 1 ! accumulate (only for istep = N)
  !end do

  !correct version:
  do istep = 1, N
    call random_number(rnd) ! get a sequence of random numbers
    do irun = 1, nruns      
      if (rnd(irun) < pleft) then ! random move
        ix(irun) = - jump ! left
      else
        ix(irun) = + jump ! right
      end if
    end do

    x_runner = x_runner  + ix !has the position of all the walkers!!
    !x2_runner = x2_runner + ix**2

    !file: step, pos, pos^2 for each runner
    if (.false.) then
      write(unit=2, fmt=*, iostat=ios) istep, x_runner
      if ( ios /= 0 ) stop "Write error in file unit 2"
    endif

    !file: 
    !step,  position,  sqeare position, variance , theoretical variance, %  of difference
    !avreaged over the runners
    
    delta = ( Sum(x_runner**2)/real(nruns) - (Sum(x_runner)/real(nruns))**2)
    deltath = (4.0 * (1-pleft) * pleft * istep * jump**2) 
    
    write(unit=3, fmt=*, iostat=ios) istep, Sum(x_runner)/real(nruns) , Sum(x_runner**2)/real(nruns), &
      delta, deltath, abs(delta/deltath - 1.0) * 100
    if ( ios /= 0 ) stop "Write error in file unit 3"
  
  end do

  do irun = 1, nruns, 1    
    P_N(x_runner(irun)) = P_N(x_runner(irun)) + 1 
  end do

  write(1,*)"# N, mean deviations, mean squared deviations, sigma^2"
  
  do i = - N, N
    write(1,*) i,real(P_N(i))/nruns, P_N_anal(i)
  end do


  !chek how the simulation wents.....
  print*, " "
  print*, "significant results"
  print*,"# N=",N,"  nruns=",nruns
  print*,"# <x_N>   = ",real(sum(x_runner))/real(nruns), "th value =", N*((1.0-pleft)-pleft)*jump 
  print*,"# <x^2_N> = ",real(sum(x_runner**2))/real(nruns), "th value =", real(N) * jump**2
  print*,"# <x^2_N> -  <x_N>^2 = ", delta, "th value=", deltath





  !closing file
  close(unit=1, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 1"

  close(unit=2, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 2"
  
  close(unit=3, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 3"




  deallocate(rnd)        
  deallocate(x_runner)   
  deallocate(x2_runner)  
  deallocate(ix)         
  deallocate(P_N)
  deallocate(P_N_anal)
  deallocate(P_N_anal_sto)
  deallocate(seed)

  stop
end program rw1d