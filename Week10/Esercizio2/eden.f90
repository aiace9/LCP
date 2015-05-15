!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! eden.f90
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module common
  
  implicit none
  public::load, init, edengen
  integer, parameter, public :: d=2
  integer, public ::  Lx, Ly, nmcs, posx, c, xmax
  real,public :: rnd
  integer, public, dimension(1)::seed

contains

  !grid parameters

  subroutine load()

    print*, "L>"
    read*, Lx

    print*, "nmcs>"
    read*, nmcs
    Ly=nmcs

    print*, "seed>"
    read*, seed

  end subroutine load

  !Initialize the lattice

  subroutine init(grid, Lx,Ly,s, v)
    integer,intent(inout) :: Lx,Ly,v
    integer :: i
    integer, dimension(Lx,Ly),  intent(inout) :: grid
    integer, dimension(2,v),  intent(inout) :: s

    grid = 0
    s = 0

    do i=1,Lx
       grid(i,1) = 1
       s(1,i)=i
       s(2,i)=2
    end do

    !do i= 1, nmcs
    !   w(i) = 0.0_d
    !   hmed(i) = 0.0_d
    !end do

  end subroutine init

  !eden model

  subroutine edengen(grid,Lx,Ly,v,s)

    integer,intent(inout) :: v,Lx,Ly
    integer :: i,ccp,j
    integer, dimension(Lx,Ly),  intent(inout) :: grid
    integer, dimension(2,v),  intent(inout) :: s
    integer, dimension(2) :: loc

    call random_seed (put = seed)

    loc = minloc(s)
    xmax = loc(2) - 1
    print*,"xmax = ",xmax

    call random_number(rnd)
    posx=int(rnd*xmax)+1
    c=0
    print*,"posx=",posx, "s(1:2,pox)=",s(:,posx)


    grid(s(1,posx),s(2,posx))=1

    if (s(1,posx)==1) then

       ccp=Lx
    else
       ccp=s(1,posx)-1
    end if

    if (grid(ccp,s(2,posx))==0) then

       do i=1,xmax
          if ((s(1,i)/=ccp) .and. (s(2,i)/=s(2,posx))) then
             s(1,posx)=ccp
             s(2,posx)=s(2,posx)
             c=1
          end if
       end do
    end if

    if (s(1,posx)==Lx) then
       ccp=1
    else
       ccp=s(1,posx)+1
    end if

    if (grid(ccp,s(2,posx))==0) then

       do i=1,xmax
          if ((s(1,i)/=ccp) .and. (s(2,i)/=s(2,posx))) then
             if (c==0) then
                s(1,posx)=ccp
                s(2,posx)=s(2,posx)
                c=1
             else
                s(1,xmax+1)=ccp
                s(2,xmax+1)=s(2,posx)
                xmax=xmax+1
             end if
          end if
       end do
    end if

    if (s(2,posx)==Ly) then
       ccp=1
    else
       ccp=s(2,posx)+1
    end if

    if (grid(s(1,posx),ccp)==0) then
       print*,"s(1,posx)=",s(1,posx)
       print*,"ccp=",ccp
       do i=1,xmax
          if ((s(1,i)/=s(1,posx)) .and. (s(2,i)/=ccp)) then
             print*,"si"
             if (c==0) then
                s(1,posx)=s(1,posx)
                s(2,posx)=ccp
                c=1
             else
                s(1,xmax+1)=s(1,posx)
                s(2,xmax+1)=ccp
                xmax=xmax+1
             end if
          end if
       end do
    end if

    !   do i=1,2
    !   print*,(s(i,j), j=1,v)
    !   end do
    !   do i=1,Ly
    !   print*,(grid(j,i), j=1,Lx)
    !   end do

    if (grid(s(1,posx),s(2,posx)-1)==0) then
       print*,"s4"
       do i=1,xmax
          if ((s(1,i)/=s(1,posx)) .and. (s(2,i)/=s(2,posx)-1)) then
             if (c==0) then
                s(1,posx)=s(1,posx)
                s(2,posx)=s(2,posx)-1
                c=1
             else
                s(1,xmax+1)=s(1,posx)
                s(2,xmax+1)=s(2,posx)-1
                xmax=xmax+1
             end if
          end if
       end do
    end if

  end subroutine edengen

end module common

program eden

  use common
  implicit none
  integer::i,v,j
  integer, dimension(:,:), allocatable :: grid
  real, dimension (:), allocatable :: w, hmed
  integer, dimension(:,:), allocatable :: s

  open (unit=1, file="eden1.dat", status= "replace", action="write")

  call load ()

  v=lx*ly

  allocate(grid(Lx,Ly))
  allocate(w(nmcs))
  allocate(hmed(nmcs))
  allocate(s(2,v))

  call init(grid,Lx,Ly, s, v)
  print*,"v=lx*ly=",v
  print*,"nmcs=",nmcs
  do i=1, nmcs
 !    print*," imcs=",i,"    grid:"
 !    do j=ly,1,-1
 !       print*,grid(1:lx,j)
 !    end do
     call edengen(grid,Lx,Ly,v,s)
  end do

  write(unit=1,fmt=*) s

  close(unit=1)

  deallocate(grid,w,hmed,s)

end program eden




