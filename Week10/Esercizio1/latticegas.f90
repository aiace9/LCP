program latticegas

  implicit none

  logical, allocatable :: lattice(:,:) 
  integer, allocatable :: x(:),y(:)
  double precision, allocatable :: dx(:),dy(:)

  integer  :: Nsteps,Np,L 
  integer  :: istep,isubstep, q 
  integer  :: dir,i,j,nfail,njumps
  integer :: seed1
  integer, dimension(12) :: seed

  integer,  parameter  ::  MAXINT=1000000000  ! Variables  for  counting
  ! allowed    directions   
  integer :: free(4),nfree   
  integer :: dxtrial(4),dytrial(4) 
  integer :: xnew(4),ynew(4)

  real, dimension(2) :: rnd(2)
  real :: rnd1
  double  precision :: dxsum,dysum,dxsqsum,dysqsum 
  double  precision :: t,deltat,drsqave,D,a,help

  !error label
  integer :: err,ios

  ! Set average time  between jumps and jump length Units is  s and cm 
  ! although actually this is not needed for the simulation 
  !deltat=1d-9; ! 1 ns 
  !a=2e-8; ! 2  A
  deltat=1.0; ! 1 ns 
  a=1.0; ! 2  A

  print*," # of steps>"
  read*, Nsteps  
  print*," # of particles>"
  read*, Np 
  print*," Side of the lattice (cm?)>"
  read*, L
  print*," seed>"
  read*, seed1
  do q = 1, 12, 1
      seed(q)= seed1 + 37 * (q-1)
  end do

  call random_seed(put=seed) 
  print *,'Doing lattice gas walk to',Nsteps,'MC steps, initial seed',seed1
  print *,'using',Np,' particles  on a',L,'^2 square lattice'
  print*, 'Density = ', Np/(L*L*1.0)

  if (Np >= L*L) then
     print *,'Number of particles > number of sites !!!' 
     STOP 'Too small lattice' 
  endif
  
  allocate(lattice(0:L-1,0:L-1),stat=err)
  if (err /= 0) print *, "array: Allocation request denied"
  allocate(x(Np),stat=err)
  if (err /= 0) print *, "array: Allocation request denied"
  allocate(y(Np),stat=err)
  if (err /= 0) print *, "array: Allocation request denied"
  allocate(dx(Np),stat=err)
  if (err /= 0) print *, "array: Allocation request denied"
  allocate(dy(Np),stat=err)
  if (err /= 0) print *, "array: Allocation request denied"

  ! Mark all positions as empty 
  lattice = .false.
  !do i=0,L-1
  !   do j=0,L-1
  !      lattice(i,j) = .false. 
  !   enddo
  !enddo

  ! enumeration of directions: 1 left 2 right 3 up 4 down 
  dxtrial(1)=+1; dytrial(1)= 0;   
  dxtrial(2)=-1; dytrial(2)= 0;   
  dxtrial(3)= 0; dytrial(3)=+1; 
  dxtrial(4)= 0; dytrial(4)=-1;

  nfail=0; njumps=0; 
  ! Generate particles on lattice 
  
  open(unit=7, file='D_t.dat', iostat=ios)
  if ( ios /= 0 ) stop "Error opening file D_t.dat"

  do i=1,Np
     do ! Loop until empty position found
        !  To   be  on  safe  side,   check  that  upper   limit  not  reached
        call random_number(rnd)
        x(i)=int(rnd(1)*L);  if (x(i)>=L) x(i)=L-1;  
        y(i)=int(rnd(2)*L);  if (y(i)>=L) y(i)=L-1; 
        if (lattice(x(i),y(i))) then
           ! Position already filled, loop to find new trial 
           cycle 
        else
           lattice(x(i),y(i))=.true. 
           !  Success, go  to next particle  
           exit 
        endif
     enddo
     dx(i)=0.0d0; dy(i)=0.0d0; 
  enddo

  t=0.0; 
  do istep=0,Nsteps-1 ! Loop over MC steps
     do isubstep=1,Np ! Do all particles on average once every MC step
        ! Pick one particle at random 
        call random_number(rnd1)
        i=int(rnd1*Np)+1; if (i>Np) i=Np;

        ! Find possible directions, store it in free() 
        nfree=0 
        do j=1,4
           xnew(j)=x(i)+dxtrial(j);
           !verifico le PBC
           if  (xnew(j) >= L) xnew(j)=0; if (xnew(j)<0) xnew(j)=L-1; 
           ynew(j)=y(i)+dytrial(j); 
           !verifico le PBC
           if  (ynew(j) >= L) ynew(j)=0; if (ynew(j)<0) ynew(j)=L-1;

           if (.not. lattice(xnew(j),ynew(j))) then
              ! Success: position free 
              nfree=nfree+1 
              free(nfree)=j 
           endif
        enddo

        ! If no possible directions, get new particle 
        if (nfree == 0) then
           nfail=nfail+1 
           cycle 
        endif
        njumps=njumps+1

        !  Pick  one of  the  possible directions  randomly  
        !  Note that  the dir>nfree  check here  really is  needed!  
        call random_number(rnd1)
        dir=int(rnd1*nfree)+1; if (dir>nfree) dir=nfree 
        j=free(dir)

        ! Now x(i),y(j) is old position and xnew(j),ynew(j) new 
        ! Double check that new site really is free 
        if (lattice(xnew(j),ynew(j))) then
           print *,'ERROR:    THIS   SHOULD   BE   IMPOSSIBLE'   
           print *,i,j,dir,nfree  
           print *,free  
           print *,x(i),y(i),xnew(j),ynew(j) 
           STOP 'ERROR  new  site  bug'  
        endif
        !Empty  old  position  and  fill  new
        lattice(x(i),y(i))=.false. 
        lattice(xnew(j),ynew(j))=.true.

        x(i)=xnew(j); y(i)=ynew(j);          
        dx(i)=dx(i)+dxtrial(j); dy(i)=dy(i)+dytrial(j); 
     enddo

     if (mod(istep,1000) == 0) then
        ! Calculate  and print intermediate results
        ! Get total displacement from dx,dy 
        dxsum=0.0d0; dysum=0.0d0; 
        dxsqsum=0.0d0; dysqsum=0.0d0; 
        do i=1,Np
           dxsum=dxsum+dx(i);   dysum=dysum+dy(i);   
           dxsqsum=dxsqsum+dx(i)*dx(i);
           dysqsum=dysqsum+dy(i)*dy(i); 
        enddo
        drsqave=(dxsqsum+dysqsum)/(1.0*Np)

        if (t>0.0) then
           !  Get  diffusion coefficient  by  proper scaling  
           D=drsqave*a*a/(4*t)
           !write(*,fmt='(3(a,1pe10.2))')&
           !     'At ',t,' <dR^2>=',drsqave*a*a,' D=',D,' cm^2/s'
           write(unit=7, fmt=*, iostat=ios) t,D
           if ( ios /= 0 ) stop "Write error in file unit 7"
        endif

     endif

     t=t+deltat  

  enddo

  close(unit=7, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 7"

  !  Get   total  displacement   from  dx,dy
  dxsum=0.0d0; dysum=0.0d0; 
  dxsqsum=0.0d0; dysqsum=0.0d0; 
  do i=1,Np
     dxsum=dxsum+dx(i);   dysum=dysum+dy(i);   
     dxsqsum=dxsqsum+dx(i)*dx(i); dysqsum=dysqsum+dy(i)*dy(i);     
  enddo
  print *,'dxsum',dxsum,' dysum',dysum 
  print *,'dxsqsum',dxsqsum,' dysqsum',dysqsum

  drsqave=(dxsqsum+dysqsum)/(1.0*Np)  
  print *,'drsqave',drsqave 
  print *,'Number of  failed jumps',nfail,' number of  successes',njumps 
  ! Get diffusion  coefficient  by  proper scaling  
  D=drsqave*a*a/(4*t) 
  write(*,fmt='(a,f6.4,a)')'For Np/L^2=',real(Np)/L**2,' :'
  write(*,fmt='(3(a,1pe10.2))')&
       'at ',t,' <dR^2>=',drsqave*a*a,' D=',D,' cm^2/s'

end  program latticegas
