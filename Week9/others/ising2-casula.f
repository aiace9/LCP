      program main
      implicit none
      integer n, l, m, dim, opt, time
      real t, j, k, prob(-4:4), energy, magnetisation
      parameter(k=1.,j=1.,dim=20)
      integer spin(dim,dim)
      print*, 'enter temperature'
      read*, t
      print*, 'enter opt'
      read*, opt
      print*, 'enter time'
      read*, time
      open(5,file='ising.txt',status='unknown')
      open(6,file='isinganalisys.txt',status='unknown')
      call srand(30)
      call initial(spin,opt,dim)
      prob(-4)=exp(-4*j/(k*t))
      prob(-2)=exp(-2*j/(k*t))
      do 10 n=1,time
         call analisys(j,dim,spin,energy,magnetisation)
         call flip(dim,spin,prob)
         write(6,200) n-1,energy,magnetisation
 200     format(i4,3x,f5.2,3x,f5.2)
 10   continue   
      write(5,100) ((spin(l,m),l=1,20),m=1,20)
 100  format(20(1x,i2))      
      close(5)
      close(6)
      stop
      end
c
c
      subroutine initial(spin,opt,dim)
      integer  x, y, opt, dim
      integer spin(dim,dim)
      real*4 rand
      if(opt.eq.1) then
         do 10 x=1,dim
            do 20 y=1,int(0.5*dim+0.5)
               spin(x,y)=1
 20         continue
            do 30 y=int(0.5*dim+0.5),dim
               spin(x,y)=-1
 30         continue
 10      continue
c
      else if(opt.eq.2) then
         do 40 x=1,dim
            do 50 y=1,dim
               if(rand().lt.0.5) then
                  spin(x,y)=1
               else
                  spin(x,y)=-1
               endif
 50         continue
 40      continue
c
      else if(opt.eq.3) then
         do 60 x=1,dim
            do 70 y=1,dim
               if(mod(x,2)+mod(y,2).eq.1) then
                  spin(x,y)=1
               else
                  spin(x,y)=-1
               endif
 70         continue
 60      continue
c
      else if(opt.eq.4) then
         do 80 x=1,dim
            do 90 y=1,dim
               spin(x,y)=1
 90         continue
 80      continue
      endif   
      return
      end
c
c
      subroutine flip(dim,spin,prob)
      integer dim, x, y, n
      integer spin(dim,dim), spinflip, test
      real prob(-4:4)
      real*4 rand
      do 10 n=1,dim*dim
         x=1+int(dim*rand())
         y=1+int(dim*rand())
         spinflip=-1*spin(x,y)
         if(x.gt.1 .and. x.lt.dim) then
            if(y.gt.1 .and. y.lt.dim) then
               test=spinflip*(spin(x-1,y)+spin(x,y+1)+spin(x,y-1)+
     &              spin(x+1,y))
            else if(y.eq.1) then
               test=spinflip*(spin(x-1,y)+spin(x,y+1)+spin(x,dim)+
     &              spin(x+1,y))
            else if(y.eq.dim) then
               test=spinflip*(spin(x-1,y)+spin(x,1)+spin(x,y-1)+
     &              spin(x+1,y))
            endif
c
         else if(x.eq.1) then
            if(y.gt.1 .and. y.lt.dim) then
               test=spinflip*(spin(dim,y)+spin(x,y+1)+spin(x,y-1)+
     &              spin(x+1,y))
            else if(y.eq.1) then
               test=spinflip*(spin(dim,y)+spin(x,y+1)+spin(x,dim)+
     &              spin(x+1,y))
            else if(y.eq.dim) then
               test=spinflip*(spin(dim,y)+spin(x,1)+spin(x,y-1)+
     &              spin(x+1,y))
            endif
c            
         else if(x.eq.dim) then
            if(y.gt.1 .and. y.lt.dim) then
               test=spinflip*(spin(x-1,y)+spin(x,y+1)+spin(x,y-1)+
     &              spin(1,y))
            else if(y.eq.1) then
               test=spinflip*(spin(x-1,y)+spin(x,y+1)+spin(x,dim)+
     &              spin(1,y))
            else if(y.eq.dim) then
               test=spinflip*(spin(x-1,y)+spin(x,1)+spin(x,y-1)+
     &              spin(1,y))
            endif
         endif
c     
         if(test.ge.0) then
            spin(x,y)=spinflip
         else
            if(rand().lt.prob(test)) then
               spin(x,y)=spinflip
            endif
         endif
 10   continue
      return
      end
c
      subroutine analisys(j,dim,spin,energy,magnetisation)
      integer dim, x, y, n, m, spinflip, test
      integer spin(dim,dim)
      real energy, magnetisation, j
      energy=0.
      magnetisation=0.
      do 10 x=1,dim
         do 20 y=1,dim
            magnetisation=spin(x,y)+magnetisation
            spinflip=spin(x,y)
            if(x.gt.1 .and. x.lt.dim) then
               if(y.gt.1 .and. y.lt.dim) then
                  test=spinflip*(spin(x-1,y)+spin(x,y+1)+spin(x,y-1)+
     &                 spin(x+1,y))
               else if(y.eq.1) then
                  test=spinflip*(spin(x-1,y)+spin(x,y+1)+spin(x,dim)+
     &                 spin(x+1,y))
               else if(y.eq.dim) then
                  test=spinflip*(spin(x-1,y)+spin(x,1)+spin(x,y-1)+
     &                 spin(x+1,y))
               endif
c
            else if(x.eq.1) then
               if(y.gt.1 .and. y.lt.dim) then
                  test=spinflip*(spin(dim,y)+spin(x,y+1)+spin(x,y-1)+
     &                 spin(x+1,y))
               else if(y.eq.1) then
                  test=spinflip*(spin(dim,y)+spin(x,y+1)+spin(x,dim)+
     &                 spin(x+1,y))
               else if(y.eq.dim) then
                  test=spinflip*(spin(dim,y)+spin(x,1)+spin(x,y-1)+
     &                 spin(x+1,y))
               endif
c            
            else if(x.eq.dim) then
               if(y.gt.1 .and. y.lt.dim) then
                  test=spinflip*(spin(x-1,y)+spin(x,y+1)+spin(x,y-1)+
     &                 spin(1,y))
               else if(y.eq.1) then
                  test=spinflip*(spin(x-1,y)+spin(x,y+1)+spin(x,dim)+
     &                 spin(1,y))
               else if(y.eq.dim) then
                  test=spinflip*(spin(x-1,y)+spin(x,1)+spin(x,y-1)+
     &                 spin(1,y))
               endif
            endif            
            energy=0.5*test*j+energy
 20      continue
 10   continue
      magnetisation=magnetisation/(dim*dim)
      energy=-energy/(dim*dim)
      return
      end
            
            







      

      

