! Program to simulate the growth of a 2D crystal that forms
! under diffusion-limited aggregation
! 
! G. Hart
! NAU
! March 2002
!
program dla2d
use random_stuff
implicit none

integer, parameter :: Nw = 200  ! Number of walkers (particles in the crystal)
logical, dimension(-Nw:Nw,-Nw:Nw) :: occupied ! Grid where the xtal grows
logical :: stuck ! Did the current walker get stuck yet?
integer :: mass ! number of particles in the cluster inside a given radius
integer, dimension(2) :: newpos, prevpos ! current and previous position of the current walker
integer :: i, j, idist ! general loop counters
real :: radius ! outer radius of the crystal plus a little
real :: distance ! distance of the walker from the origin
real :: theta, rndstep ! random numbers for starting and stepping the walkers, respectively
real :: twopi


call set_random_seed(0)
radius = 5
twopi = 8*atan(1.0)
occupied(:,:) = .false.! Initialize the array
occupied(0,0) = .true. ! Make the origin occupied (this is the seed crystal)

do ! Start a walker at a random postion outside the crystal
   call random_number(theta)
   theta = theta * twopi
   newpos = nint( radius*(/cos(theta),sin(theta)/)) ! Start a new walker
   if(occupied(newpos(1),newpos(2)))cycle ! Already occupied, try again
   prevpos = newpos
   do 
      newpos = prevpos
      call random_number(rndstep)
      select case(int(rndstep*4)+1)
         case(1)
            newpos(1) = newpos(1) -1
         case(2)
            newpos(1) = newpos(1) +1
         case(3)
            newpos(2) = newpos(2) -1
         case(4)
            newpos(2) = newpos(2) +1
         end select
      if(any(abs(newpos) > Nw)) exit ! Walker stepped out of the box. Start a new one
      if(occupied(newpos(1),newpos(2))) then ! walker gets stuck to the crystal
         occupied(prevpos(1),prevpos(2)) = .true. ! Add the walker to the crystal
         distance = sqrt(real(dot_product(prevpos,prevpos)))
         if(distance > (radius-5)) radius = distance + 5 ! Make the starting circle larger if necessary
         exit
      endif
      prevpos = newpos ! Walker made a valid move (didn't get stuck or wander away). Update and keep it moving
   enddo
   if(radius > Nw) exit
enddo

! Write occupied sites to disk
open(10,file="dla2d.data",status="replace")
do i = -Nw, Nw
   do j = -Nw, Nw
      if(occupied(i,j)) write(10,*) i,j
   enddo
enddo
close(10)

! Do the m(r) analysis and write results to disk
open(11,file="dlamass.data")
do idist = 2, int(0.75*distance)
   mass = 0
   do i = -Nw, Nw
      do j = -Nw, Nw
         if(occupied(i,j)) then
            if(idist**2 >= i**2 + j**2) mass = mass + 1
         endif
      enddo
   enddo
   write(11,'(2i10)') idist, mass
enddo
close(11)

end program dla2d
