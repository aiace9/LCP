      program ising2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c algoritmo di Metropolis per calcolare <E>, <M>, nell'insieme canonico
c (fisso T,N,V) con un modello di Ising 2D
c
c Qui e' posto: K_B = 1
c               J   = 1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c il file 'ising2.dat' e' scritto nel formato opprtuno per usare GNUPLOT:
c vedasi i seguenti comandi:
c 
c gnuplot> set view 0,0
c gnuplot> set nosurface
c gnuplot> set contour
c gnuplot> set cntrparam levels discrete -1,1
c gnuplot> splot "ising2.dat"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c La stima di E e M e' implementata in diversi modi:
c a) con Econf (updating di E per ogni configurazione, cioe' dopo ogni
c                 singola mossa MC)
c b) con Estep (updating di E per ogni step MC, cioe' dopo 1 mossa su
c                 tutti gli spin)
c NB: a) e b) NON sono perfettamente equivalenti !
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      parameter (maxL=32,maxN=maxL*maxL)
      dimension  ispin(maxL,maxL),w(-4:4)
      data ispin/maxN*0/
      Estepcum  = 0.
      Estep2cum = 0.
      Econfcum  = 0.
      Econf2cum = 0.
      Mconf     = 0
      M2conf     = 0

      ecum   = 0.
      e2cum  = 0.
      xmcum  = 0.
      xm2cum = 0.
      open(unit=8,file='ising2.dat',status='unknown')
      call start(N,L,T,nmcs,nmcseq,ispin,E,M,w)
      if(L.gt.maxL)call error(L,maxL)
      print*,' |magn.| iniziale netta per spin = ',abs(float(M))/N
      print*,' energia iniziale per spin = ',E/N
      do 2001 j = 1,L
 2001    write(*,'(32i2)')(ispin(i,j),i=1,L)
      do 100 imcs = 1,nmcs
         call Metrop(N,L,ispin,E,M,w,ratio,Econf,Mconf)
         Estep = 0.
         do i=1,L
            in = i+1
            if(in.gt.L)in=1
            do j=1,L
               jn = j+1
               if(jn.gt.L)jn=1
               Estep = Estep - ispin(i,j)*(ispin(i,jn)+ispin(in,j))
            enddo
         enddo
      Mconfcum  = Mconfcum + abs(Mconf)
      Mconf2cum = Mconf2cum + Mconf * Mconf
      Mconfave  = Mconfcum/float(N)/float(imcs)
      Mconf2ave = Mconf2cum/float(N*N)/float(imcs)
      sigmam2conf = Mconf2ave - Mconfave * Mconfave
      
      Estepcum  = Estepcum + Estep
      Estep2cum = Estep2cum + Estep * Estep
      Estepave  = Estepcum/N/float(imcs)
      Estep2ave = Estep2cum/(N*N)/float(imcs)
      sigma2step = Estep2ave - Estepave * Estepave

      Econfcum  = Econfcum + Econf
      Econf2cum = Econf2cum + Econf * Econf
      Econfave  = Econfcum/N/float(imcs)
      Econf2ave = Econf2cum/(N*N)/float(imcs)
      sigma2conf = Econf2ave - Econfave * Econfave
c stampa gli ultimi 20 MCsteps
c      if (imcs.gt.nmcs-20)then
c         print *, 
c     & imcs,E/N,Estep/N,Estepcum/N/float(imcs),Econfcum/N/float(imcs),
c     & Mconfcum/N/float(imcs)
c         print *, sqrt(sigma2step/(imcs-1)), sqrt(sigma2conf/(imcs-1)),
c     & sqrt(sigmam2conf/(imcs-1))
c         end if
         if(imcs.gt.nmcseq)call data(E,M,ecum,e2cum,xmcum,xm2cum,ratio)
 100     continue
         print *, 
     & imcs,E/N,Estep/N,Estepcum/N/float(imcs),Econfcum/N/float(imcs),
     & Mconfcum/N/float(imcs)
         print *, sqrt(sigma2step/(imcs-1)), sqrt(sigma2conf/(imcs-1)),
     & sqrt(sigmam2conf/(imcs-1))
         call output(N,nmcs,nmcseq,ecum,e2cum,xmcum,xm2cum,ratio)
      do 201 j = 1,L
 201     write(*,'(32i2)')(ispin(i,j),i=1,L)
      do 111 j = 1,L
      do 112 i = 1,L
         write(8,*)ispin(i,j)
 112  continue
         write(8,*)
 111  continue
         close(8)
         stop
         end

      subroutine start(N,L,T,nmcs,nmcseq,ispin,E,M,w)
      dimension ispin(32,32),w(-4:4)
      write(6,*)' dimensione lineare del reticolo (un intero, L) = '
      read(5,*) L
      N = L*L ! numero di spin
      write(6,*)' temperatura (T) = '
      read(5,*) T
      write(6,*)' numero totale di steps Monte Carlo per spin (nmcs) = '
      read(5,*) nmcs
      write(6,*)' numero di steps Monte Carlo per spin di equilibratura 
     &           (nmcseq) = '
      read(5,*) nmcseq
      write(6,*)' n. intero > 0 per inizializzare la seq. random = '
      read(5,*) iseed0
      call ranset(iseed0)
c genera la configurazione iniziale
      print*,' inizio: spin down,a caso,up,meta`,antif (-1,0,1,2,3)? >'
      read(*,*)iconfig
      M = 0.
      E = 0.
      if ( iconfig .eq. 0) then
      do 200 j = 1,L
         do 100 i = 1,L
            if( ranf() .lt. 0.5 ) then
               ispin(i,j) = 1
               else
               ispin(i,j) = -1
               end if
               M = M + ispin(i,j)  ! magnetizzazione iniziale netta
 100           continue
 200           continue
      else if (iconfig .eq. -1 .or. iconfig .eq. 1) then   !CONFIG. ORDINATA
      do 201 j = 1,L
         do 101 i = 1,L
            ispin(i,j) = iconfig
 101        continue
 201        continue
               M = N * iconfig
      else if (iconfig .eq. 2) then  !ALTRA CONFIG. ORDINATA
      do 211 j = 1,L
         do 111 i = 1,L/2
            ispin(i,j) = 1
 111        continue
         do 311 i = L/2+1,L
            ispin(i,j) = -1
 311        continue
 211        continue
               M = 0
      else if (iconfig .eq. 3) then  !ALTRA CONFIG. ORDINATA: ANTIFERROM.
      do 212 j = 1,L
         do 112 i = 1,L
            ispin(i,j) = 1
            if( mod(i,2)+mod(j,2).eq.1 ) ispin(i,j) = -1 
 112     continue
 212     continue
               M = 0
       else 
                  print*,' caso non previsto'
                  stop
       end if
c calcola l'energia iniziale sommando su tutte le coppie di NN
c (dato uno spin, solo con il NN sopra e il NN destra, 
c e non sotto e a sin., cosi' contiamo le interazioni una volta sola)
      do 300 j = 1,L
         if (j .eq. L) then
            iup = 1
            else
               iup = j + 1
               endif
               do 400 i = 1,L
                  if (i .eq. L) then
                     iright = 1
                     else
                        iright = i + 1
                        endif
                        sum = ispin(i,iup) + ispin(iright,j)
                        E = E - ispin(i,j)*sum ! ENERGIA TOTALE iniziale
 400                    continue
 300                    continue
c calcola le probabilita' di transizione secondo la distrib. di Boltzmann
c (exp(-deltaE/KT)
c Avendo posto la forza delle interazioni (il param. J)=1, 
c le possibili variazioni di energia per uno spin flip sono -8,-4,0,+4,+8: 
c Argomento dell'exp  = cambiamento in energia
      e4 = exp(-4/T)
      e8 = e4 * e4
c argomento di w = isum = somma degli spin 4 NN (puo' essere -4,-2,0,+2,+4)
c Fissiamo il caso in cui lo spin centrale sia +1: allora
c se i 4 NN sono tutti +1, il suo flip corrisponde a deltaE = +8, e
c quindi associamo la prob. exp(-8/T):
      w(4)  = e8    
c e via dicendo per gli altri casi possibili 
c (se isum=0 inutile sprecar tempo, w(0)=1 q.que sia T)
      w(-4) = e8 !1/e8 ! cervellotico, ma tanto poi la usiamo solo per 
c             variaz. di E POSITIVE (le negative son comunque accettate)
      w(2)  = e4
      w(-2) = e4 !1/e4
      return
      end



      subroutine Metrop(N,L,ispin,E,M,w,ratio,Econf,Mconf)
      dimension ispin(32,32),w(-4:4)
      Mconf = 0
      Econf = 0.
      do 100 is = 1,N
c prende uno ispin a caso (i,j sono le coordinate) 
c int: returns the largest integer whose absolute value
c does not exceed the absolute value of the argument and has the same
c sign as the argument.
 11      r1 = ranf()
         if(r1.eq.1.)go to 11
 12      r2 = ranf()
         if(r2.eq.1.)go to 12
         i = int (L*r1 + 1.)
         j = int (L*r2 + 1.)
         if(i.gt.L.or.j.gt.L)then
            print*,' ATTENZIONE! i,j,L=',i,j,L
         end if
c determina i valori degli spin NN usando le PBC
         call periodic(i,j,L,ispin,isum)
         if (ispin(i,j)*isum .le. 0) then  ! se la somma dei 4 spin NN a
                                          ! (i,j) ha segno opposto a
                                          ! spin(i,j) allora il deltaE
                                          ! per spin flip di (i,j) e' < 0
                                          ! e la mossa e' accettata
            call accept(i,j,M,E,isum,ispin,ratio)
            else if ( ranf() .lt. 
     &                w(-ispin(i,j)*isum) ) then
               call accept(i,j,M,E,isum,ispin,ratio) ! altrimenti, viene
                                                    ! accettata o no
                                                    ! secondo la
                                                    ! probabilita' di
                                                    ! transizione exp(-dE/T)
c cosi' l'argomento di w e' sempre <0 , =-2 o =-4
               end if
               Econf = Econf - ispin(i,j)*isum
               Mconf = Mconf + ispin(i,j)
 100           continue

      Econf = Econf / 2.  ! perche' ho contato le coppie 2 volte !
c      if(Econf.ne.E)print*,' energie:',Econf,E
               
      return
      end
               

      subroutine periodic(i,j,L,ispin,isum)
c determina la somma dei valori sui NN usando le PBC
      dimension ispin(32,32)
      ir = i+1
      il = i-1
      ju = j+1
      jd = j-1
      if(ir.gt.L)ir=1
      if(il.le.0)il=L
      if(ju.gt.L)ju=1
      if(jd.le.0)jd=L
      isum = ispin(ir,j) + ispin(il,j) + ispin(i,ju) + ispin(i,jd)
      return
      end

      subroutine periodic_old(i,j,L,ispin,isum)
c determina la somma dei valori sui NN usando le PBC
      dimension ispin(32,32)
      if (i .eq. 1) then
         left = ispin(L,j)
         else
            left = ispin(i-1,j)
            endif
      if (i .eq. L) then
         iright = ispin(1,j)
         else
            iright = ispin(i+1,j)
            endif
      if (j .eq. 1) then
         idown = ispin(i,L)
         else
            idown = ispin(i,j-1)
            endif
      if (j .eq. L) then
         iup = ispin(i,1)
         else
            iup = ispin(i,j+1)
            endif
      isum = left + iright + iup + idown
      return
      end
      


      subroutine accept(i,j,M,E,isum,ispin,ratio)
c entrano M,E della config. prec. ed escono quelli della nuova confign.
      dimension ispin(32,32)
      ispin(i,j) = -ispin(i,j)
      ratio = ratio + 1        ! poi verra' rinormalizzata per
                               ! det. l'acc. ratio
      M = M + 2*ispin(i,j)      ! qui ci va il 2 perche' e' la variaz.
      E = E - 2*ispin(i,j)*isum ! idem qui
      return
      end


      subroutine data(E,M,ecum,e2cum,xmcum,xm2cum,ratio)
c accumula dati dopo ogni step Monte Carlo per spin
      ecum = ecum + E
      e2cum = e2cum + E*E
      xmcum = xmcum + abs(M)
      xm2cum = xm2cum + M*M
      return
      end

      subroutine ranset(idum)
c initialize variables in /cran/ (from ran3(idum) of numerical recipes)
      implicit none
      integer idum,mbig,mseed,mz,ma,mj,mk,i,ii,k,inext,inextp
      real*8 fac
      dimension ma(55)
      common /cran/fac,ma,mbig,mz,inext,inextp
      mbig=1000000000
      mseed=161803398
      mz=0
      fac=1.d0/mbig
c
      mj=mseed-iabs(idum)
      mj=mod(mj,mbig)
      ma(55)=mj
      mk=1
      do i=1,54
         ii=mod(21*i,55)
         ma(ii)=mk
         mk=mj-mk
         if(mk.lt.mz)mk=mk+mbig
         mj=ma(ii)
      enddo
      do k=1,4
         do i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
         enddo
      enddo
      inext=0
      inextp=31
      return
      end
      function ranf()
c generates a uniform deviate in (0,1) (from ran3 on numerical recipes)
      implicit none
      integer mbig,mz,ma,mj,inext,inextp
      real*8 ranf,fac
      dimension ma(55)
      common /cran/fac,ma,mbig,mz,inext,inextp
c
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ranf=1.d0-mj*fac
c
      return
      end


      subroutine output(N,nmcs,nmcseq,ecum,e2cum,xmcum,xm2cum,ratio)
      znorm = 1./float((nmcs-nmcseq)*N)
c medie per spin
      ratio = ratio*znorm
      eave  = ecum*znorm
      e2ave = e2cum*znorm/N
      xmave  = xmcum*znorm
      xm2ave = xm2cum*znorm/N
      sigmae = sqrt((e2ave - eave*eave)/(nmcs-nmcseq-1))
      sigmam = sqrt((xm2ave - xmave*xmave)/(nmcs-nmcseq-1))
      write(6,*)' acceptance ratio              = ',ratio
      write(6,*)' en. media per spin            = ',eave
      write(6,*)' en. quadratica media per spin = ',e2ave
      write(6,*)' errore sull`en.               = ',sigmae
      write(6,*)' |magnetizzazione| media          = ',xmave
      write(6,*)' magnetizzazione quadratica media = ',xm2ave
      write(6,*)' errore sulla magnetizzazione     = ',sigmam
      return
      end


      subroutine error(n,m)
      print*,' lato troppo grande (troppi spins)! (',n,'.ge.',m,')'
      stop
      end







