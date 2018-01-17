      MODULE my_param
      IMPLICIT NONE
      SAVE

      INTEGER  ::  N,ISTART,MAXBIN,MAXNAY,Ndensity,ncfg,prifreq
      REAL*8  ::  NSTEP, NCOUNT, Nmax, Nequil, Limit, Ndold,ilat
      REAL*8  ::  Kompress,Nrecord,Nresult,Ntrial
      REAL*8  ::  nconfig
      INTEGER, DIMENSION(:), ALLOCATABLE :: hist
      INTEGER, DIMENSION(:), ALLOCATABLE :: jc
      INTEGER, DIMENSION(:), ALLOCATABLE :: nay
      INTEGER, DIMENSION(:), ALLOCATABLE :: kw
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: nab
      REAL*8, PARAMETER :: temp = 1.0e0
      REAL*8, PARAMETER :: pi = 3.141592654e0
      REAL*8 :: dia,d0,dsq,tbig,tim,psum,zed,range,diaset
      REAL*8 :: equilibrated,tcount,zlast,timeSum,fract,gcon
      REAL*8 :: delr,ds,runfor,diastart,zedset,zmax
      REAL*8, DIMENSION(:), ALLOCATABLE :: x
      REAL*8, DIMENSION(:), ALLOCATABLE :: y
      REAL*8, DIMENSION(:), ALLOCATABLE :: z
      REAL*8, DIMENSION(:), ALLOCATABLE :: vx
      REAL*8, DIMENSION(:), ALLOCATABLE :: vy
      REAL*8, DIMENSION(:), ALLOCATABLE :: vz
      REAL*8, DIMENSION(:), ALLOCATABLE :: x0
      REAL*8, DIMENSION(:), ALLOCATABLE :: y0
      REAL*8, DIMENSION(:), ALLOCATABLE :: z0
      REAL*8, DIMENSION(:), ALLOCATABLE :: xd
      REAL*8, DIMENSION(:), ALLOCATABLE :: yd
      REAL*8, DIMENSION(:), ALLOCATABLE :: zd
      REAL*8, DIMENSION(:), ALLOCATABLE :: told
      REAL*8, DIMENSION(:), ALLOCATABLE :: tnew
      REAL*8, DIMENSION(:), ALLOCATABLE :: gra
      REAL*8, DIMENSION(:), ALLOCATABLE :: grat
      REAL*8, DIMENSION(0:12) :: density
      REAL*8, DIMENSION(9)  :: collkind
      REAL*8, DIMENSION(0:21,5) :: ss
      REAL*8:: avgbond
      INTEGER, DIMENSION(100) :: iq6dist
      REAL*8 :: RQ6NEIGHB, gq6,pram,pramsum,gq6sum,transord
      END MODULE my_param


!*********************************************************************
      subroutine initrnd(isum)
!
!------------------------------------------------------------------
!   initialize a random number generator using the hp-ux system time
!
!   compile with the +E1 option
!
!------------------------------------------------------------------
!   isum  is initial random number seed returned on
!         the call to initrnd
!------------------------------------------------------------------
!
      INTEGER dd,mm,yy
      CHARACTER*8 CLOCK
      CHARACTER*1 CI
!
!      call idate(mm,dd,yy)
      n1=dd/10
      n2=dd-n1*10
!
      call ctime(dd,clock)
      i=0
      ii=0
      isum=0
    3 i=i+1
      if(i.eq.9) go to 4
      if(i.eq.3) then
         ici=n1
        else if(i.eq.6) then
         ici=n2
        else
         ci=clock(9-i:9-i)
         ici=ichar(ci)-48
      endif
      isum=isum+ici*(10**(i))/10
      go to 3
    4 continue

      RETURN
      END
!***************

      subroutine SPHERES90(np,phi)

      
      USE my_param
      
      IMPLICIT NONE
      
      REAL*8 :: SIDE, SINV,zero,zero1,zero2,delta,rseeder
      real*8 :: randumb,dummy,rnn1,rxij,ryij,rzij,rNNsum
      REAL*8 :: RIJRT,RIJSQ,rNNmax,rNNmin
      INTEGER :: I,iseeder,J,np
      double precision :: phi

      n=np
      fract=phi

      istart=1
      prifreq=0

      pramsum=0.0e0
      gq6sum=0.0e0

      
!     READ INPUT DATA

!      write(*,*)'enter the number of molecules'
!      read(*,*)n

!      n=500
      
!      write(*,*)'enter the # of configurations to be saved'
!      write(*,*)'**(separated by 25N sphere collisions)**'
!      read(*,*)nconfig
      nconfig=1

      nconfig=nconfig
!      write(*,*)'enter the final volume fraction'
!      read(*,*)fract
      
!      WRITE(*,'('' enter Q6 neighbor cut-off'')')
!      READ (*,*) RQ6NEIGHB

!      write(*,*)'method for generating first config'
!      WRITE(*,'(''RSA and compress(0), READ in (1) or FCC (2)'')')
!      READ (*,*) ilat
      ilat=0

!      WRITE(*,'(''ENTER dimensionless compression rate'')')
!      write(*,*)' between .02 and .2 to suppress ordering'
!      write(*,*)'during compressions (best ds=0.2)'
!      read(*,*) ds
      ds=0.2
!     as defined in Lubachevsky-Stillinger

!      ds=.02

!      write(*,*)'Save all or just last configuration'
!      write(*,*)'for all enter 0, for one enter 1'
!      write(*,*) '**routine writes out coordinates and velocities**'
!      read(*,*)ncfg
      ncfg =0
      
!       write(*,*) 'enter the bin width'
!       read(*,*) delr
!      delr=.01e0
      delr=0.1

!      write(*,*)'enter 4-digit seed'
!      read(*,*)iseeder

      call initrnd(iseeder)
      iseeder=mod(iseeder,10000)
      
!      write(*,*)'enter limit for compression cycle'
!      read(*,*)limit
      limit=50000

      side = (real(n,kind(fract))/(fract*6.0E0/PI))**(1.0e0/3.0e0)
      MAXBIN = INT(SIDE/(2.0E0*DELR))
      
      write(*,*)MAXBIN,'bins for g(r)'

      do i=1,iseeder
         rseeder=randumb(dummy)
      end do

!     Allocate Arrays
      
      if (allocated(x)) deallocate(x,y,z)
      if (allocated(vx)) deallocate(vx,vy,vz)
      if (allocated(x0)) deallocate(x0,y0,z0)
      if (allocated(xd)) deallocate(xd,yd,zd)
      if (allocated(told)) deallocate(told,tnew)
      if (allocated(jc)) deallocate(jc,kw,nay,nab)
      if (allocated(gra)) deallocate(gra,grat,hist)
      ALLOCATE(x(1:N))
      ALLOCATE(y(1:N))
      ALLOCATE(z(1:N))
      ALLOCATE(vx(1:N))
      ALLOCATE(vy(1:N))
      ALLOCATE(vz(1:N))
      ALLOCATE(x0(1:N))
      ALLOCATE(y0(1:N))
      ALLOCATE(z0(1:N))
      ALLOCATE(xd(1:N))
      ALLOCATE(yd(1:N))
      ALLOCATE(zd(1:N))
      ALLOCATE(told(1:N))
      ALLOCATE(tnew(1:N))
      ALLOCATE(jc(1:N))
      ALLOCATE(kw(1:N))
      ALLOCATE(nay(1:N))
      ALLOCATE(nab(1:N,40))
      ALLOCATE(gra(1:MAXBIN))
      ALLOCATE(grat(1:MAXBIN))
      ALLOCATE(hist(1:MAXBIN))
      
!     initialize g(r) tally
      grat=0.0e0
      iq6dist=0

      call lattice
      
      call ResetSums
!      open(unit=23,status='old',file='pressure.dat')
!      open(unit=24,status='old',file='q6vst')
!      open(unit=25,status='old',file='q6dist')
!      open(unit=26,status='old',file='transorder')
!     prepare to compress
      
      kompress=-2.0e0
      equilibrated=-2.0e0
      limit=limit*real(N,kind(fract))
      
!     max number ofmoves to compress for
!      ds=0.01e0
      runfor = 2.0e0*real(N,kind(fract))

!...... compress to zedset
        delta=(zedset-zed)/25.0e0
        zlast=zed
	write(*,*)zed,zedset
        do while (zed.lt.zedset)
                kompress=2.0e0
                limit=limit-Ncount
                call ResetSums
                call advance
                kompress=2.0e0
                limit=limit-Ncount
                call ResetSums
                call advance
                
                if (zed.gt.zlast+delta) then
                   write(*,*)' compressing zed ds:',&
                       zed*sqrt(2.)*PI/6.,ds
                   zlast=zed
                end if
                
                if (zed.gt.zedset) then
                   zed=zedset
                   zmax=zedset
                end if
                if (Ncount.gt.limit) then
                   write(*,*)'jammed'
                   zmax=zed
                   exit
                end if
             end do   
!     ......... decompress to zedset
             if (zed.gt.zedset) zed=zedset

             fract=zed/real(1.35047447,kind(zed))
             ds=0.0e0
             kompress=-2.0e0
             side = (real(n,kind(fract))&
            /(fract*6.0E0/PI))**(1.0e0/3.0e0)
!........... prepare to equilibrate
!...........zero msd for particles
            diaset=(6.0e0*FRACT/PI/&
            real(N,kind(zed)))**(1.0e0/3.0e0) 
            RQ6NEIGHB =  RQ6NEIGHB * diaset
            avgbond=0.0E0

!            open(unit=33,status='unknown',file='meansq')

            
             Nequil=1000.0*real(N,kind(fract))
           
             runfor=Nequil
             call ResetSums
             do i =1,N
                xd(i)=x(i)
                yd(i)=y(i)
                zd(i)=z(i)
             end do

            
!................equilibrate for Nequil events


            

!             write(*,*) 'into equilibrate'



             do while (Ncount.lt.Nequil)
                call advance
             end do
!...............prepare to run 
             equilibrated=2.0e0
             if (nconfig.eq.1) then
                call saveconfig
                write(*,*)'done. go back to wrap program.' ! Look!!!
                return
            end if
             Nmax=25.0e0*nconfig*real(N,kind(fract))
             runfor=Nmax
             call ResetSums 
!..............run for Nmax collisions
             write(*,*) ' into main run'

             do while (Ncount.lt.Nmax) 

                call advance
                
             end do 

!...............write out g(r) info

 	open(unit=17,status='unknown',file='gofr')
	
	DO 400 I=1,MAXBIN
           
           
           zero = 0.0e0
           zero1 = .999e0  
           zero2 = 1.0e0
           
           WRITE(17,395)(I-0.5)*delr,GRAt(I)/real(nconfig)
!           if((((i-0.5)*delr).gt.(1.-delr)).and
!     :          .(((i-0.5)*delr).lt.1.001))then   
!              write(17,395)zero1,zero
!              write(17,395)zero2,gcon
!           endif
           
           
           
 395       format('',F7.3,F15.4,F15.4,F15.4)
 400    CONTINUE
        close(unit=17)
!        close(unit=33)
!        close(unit=23)

        rNNsum=0.0e0
        rNNmax=0.0e0
        rNNmin=100.0e0
        do i=1,N
           rNN1=100.0e0
           do j=1,N
              if (i.ne.j) then
                 
                 RXIJ = X(I)-X(J)
                 RYIJ = Y(I)-Y(J)
                 RZIJ = Z(I)-Z(J)
                 RXIJ = RXIJ - ANINT(RXIJ)
                 RYIJ = RYIJ - ANINT(RYIJ) 
                 RZIJ = RZIJ - ANINT(RZIJ)
                 RIJSQ=RXIJ*RXIJ +RZIJ*RZIJ + RYIJ*RYIJ
                 
!     scale to particle size
                 RIJRT=(RIJSQ)**(0.5e0)*side
                 if (RIJRT.lt.rNN1) rNN1=RIJRT
              end if
           end do
           rNNsum=rNNsum+rNN1
           if (rNN1.gt.rNNmax) rNNmax=rNN1
           if (rNN1.lt.rNNmin) rNNmin=rNN1
        end do
        
        rNNsum=rNNsum/real(n,kind(fract))

!     WRITE OUT Q6 info...

!        open(unit=25,status='unknown',file='q6dist')        
!       do i=1,100
!           write(25,*)(-0.5E0+real(i,KIND(zero)))&
!          /100.0e0,100.0E0*iq6dist(i)&
!          /real(n,KIND(zero))/real(nconfig,KIND(zero))
!        end do
!        close(unit=25)
        
!        write(*,*) 'avg. coordination=',avgbond/
!     :       real(nconfig,KIND(ZERO))
!        write(*,*) 'Global Q6=',gq6sum/
!     :       real(nconfig,KIND(ZERO))
!
!        write(*,*) 'Nb^0.5*Q6=',gq6sum/
!     :       real(nconfig,KIND(ZERO))*(real(N)*avgbond/
!     :       real(nconfig,KIND(ZERO)))**0.5
!
!        write(*,*) 'T=',pramsum/
!     :       real(nconfig,KIND(ZERO))

        
!        write(*,*) 'mean nearest neighbor dist=',rNNsum
!        write(*,*) 'largest nearest neighbor dist=',rNNmax
!        write(*,*) 'minimum nearest neighbor dist=',rNNmin
!        close(24)
!        close(unit=26)

      return 
      end subroutine spheres90

      
      REAL*8 FUNCTION RANDUMB ( DUMMY )
      
!     *********************************************************
!     ** RETURNS A UNIFORM RANDUMB VARIATE IN THE RANGE 0 TO 1. 
!     **                                                       
!     **                 ***************                       
!     **                 **  WARNING  **                       
!     **                 ***************                       
!     **                                                       
!     ** GOOD RANDUMB NUMBER GENERATORS ARE MACHINE SPECIFIC.   
!     ** PLEASE USE THE ONE RECOMMENDED FOR YOUR MACHINE.      
!     *********************************************************
      
      IMPLICIT NONE
      INTEGER  ::   L, C, Mo
      PARAMETER ( L = 1029, C = 221591, Mo = 1048576 )
      
      INTEGER    :: SEED
      REAL*8     ::   DUMMY
      SAVE       :: SEED
      DATA        SEED / 0 /
      
!     *******************************************************
      
      SEED = MOD ( SEED * L + C, Mo )
      RANDUMB = REAL ( SEED ) / Mo
      
      RETURN
      END
      
      
!     *****************************************
      REAL*8 Function Sptime(dx,dy,dz,dvx,dvy,dvz,a2,b2,dsq)

      IMPLICIT NONE
      REAL*8 dx,dy,dz,dvx,dvy,dvz,a2,b2,a,b,q,c,dsq


!     using this function slows execution by.lt.5%
!     using the same function when a2=b2=0 slows.lt.1%
      sptime=1.0E+09
      dx=dx-anint(dx)
      dy=dy-anint(dy)
      dz=dz-anint(dz)
      a=dvx*dvx+dvy*dvy+dvz*dvz+a2
      b=2*(dx*dvx+dy*dvy+dz*dvz)+b2
      if (b.lt.0.0e0 .or. a.lt.0.0e0) then 
!     surfaces moving closer
         c=dx*dx+dy*dy+dz*dz-dsq
         if (c.gt.0.0e0) then
            q=b*b-4.0e0*a*c
            if (q.ge.0.0e0)  sptime=-(b+sqrt(q))/(2.0e0*a)
         else  
! they overlap so force a collision
            sptime=0.0e0
         end if
      end if
      return
      end
      
!     *******************************************
      subroutine update(xi,yi,zi,vxi,vyi,vzi,toldi,tnow)
      IMPLICIT NONE
      REAL*8 tlag,tnow,toldi,xi,yi,zi,vxi,vyi,vzi
!     using this subroutine outside Bangt slows execution by =89 5%
      tlag=tnow-toldi
      xi=xi+tlag*vxi
      yi=yi+tlag*vyi
      zi=zi+tlag*vzi
      
      toldi=tnow
      return
      end

!     *******************************************
      Subroutine Bangt(i)
      USE my_param
      IMPLICIT NONE

      REAL*8 dtim,vxi,vyi,vzi,xi,yi,zi,a2,b2,t,sptime
      REAL*8 delx,dely,delz,delvx,delvy,delvz
      INTEGER I,L,j

      dtim=tim-1E-07
      if (told(i).lt.dtim) &
       call update(x(i),y(i),z(i),vx(i),vy(i),vz(i),told(i),tim)
      xi=x(i)
      yi=y(i)
      zi=z(i)
      vxi=vx(i)
      vyi=vy(i)
      vzi=vz(i)
      tnew(i)=tbig
        jc(i)=0  
!     identifier for the sphere that i collides with
        kw(i)=9  
!     the kind of collision (sphere-sphere=4, none=9)
        a2=-ds*ds
        b2=-2.0e0*ds*dia
        do  l=1, nay(i)
           j=nab(i,l)
           if (told(j).lt.dtim) & 
            call update(x(j),y(j),z(j),vx(j),vy(j),vz(j),told(j),tim)
           
           delx =xi-x(j)
           dely =yi-y(j)
           delz =zi-z(j)
           delvx =vxi-vx(j)
           delvy =vyi-vy(j)
           delvz =vzi-vz(j)


           t=sptime(delx,dely,delz,delvx,delvy,delvz,a2,b2,dsq)
           t=t+tim
!          ***

           if (t.lt.tnew(i)) then
              tnew(i)=t
              jc(i)=j
              kw(i)=4
           end if
        end do
!     write(*,*) i,jc(i),kw(i),tnew(i)
        return
        end


!     *****************************************
      subroutine change(i,j)   
!     to change velocities after a collision

      USE my_param
      IMPLICIT NONE
      INTEGER I,J
      REAL*8 f,rsq,dx,dy,dz,b


      call update(x(i),y(i),z(i),vx(i),vy(i),vz(i),told(i),tim)
      if (kw(i).eq. 4 ) then         
!     a pair collides
         call update(x(j),y(j),z(j),vx(j),vy(j),vz(j),told(j),tim)
!     Change velocities of i and j
         dx=x(i)-x(j)
         dy=y(i)-y(j)
         dz=z(i)-z(j)
         dx=dx-anint(dx)
         dy=dy-anint(dy)
         dz=dz-anint(dz)
         rsq=dx**2.0e0+dy**2.0e0+dz**2.0e0
         
         b=dx*(vx(i)-vx(j))+dy*(vy(i)-vy(j))+dz*(vz(i)-vz(j))
         
         if (Kompress.gt.0.0e0) then  
!     if compressing,increase the size now
            dia=diastart+ds*tim
            dsq=dia*dia
            zed=(dia/d0)**3.0e0
         end if
         f=-b/rsq+ds/dia
         psum=psum-b
         vx(i)=vx(i)+f*dx
         vy(i)=vy(i)+f*dy
         vz(i)=vz(i)+f*dz
         vx(j)=vx(j)-f*dx
         vy(j)=vy(j)-f*dy
         vz(j)=vz(j)-f*dz
      end if
      return
      end
      



!     ************************************************
      subroutine lattice
      
      USE my_param
      IMPLICIT NONE
      
      INTEGER :: i,j,k,l,ix,iy,iz
      REAL*8  :: r2, xij, yij, zij, xi, yi, zi,randumb,dummy
      REAL*8  :: fractest,rx
      
!     IF READING COORDINATES AND VELOCITIES...

      if(ilat.eq.1) then
         open(unit=34,status='unknown',file='config')
         read(34,*)fractest
 
         
         do i=1,N
            read(34,*)x(i),y(i),z(i)
         end do
         do i=1,N
            read(34,*)vx(i),vy(i),vz(i)
         end do
         close(34)
         zed=fractest*real(1.35047447,kind(zed))
         density(0)=fract*real(1.35047447,kind(zed))
         
!     if N=1372 it wont freeze below zed=0.73
         
         Ndensity=0
         Ndold=0
         zedset=density(0)
!     write(*,*)zedset*1.414213562e0*PI/6.0e0
         d0=(sqrt(2.0)/real(N,kind(fract)))**(1./3.)
         
         dia=d0*zed**(1.0e0/3.0e0)
         dsq=dia*dia
         write(*,*) "zed  zedset"
         write(*,*)zed, zedset
         
      else if(ilat.eq.0) then
!     set density for rsa start-up at a vol fract approx 30%
         zed = .405e0
!     zed = fract/fractcp
         
!     assign randumb velocities
         
         do i=1,N
            vx(i)=randumb(dummy)
            vy(i)=randumb(dummy)
            vz(i)=randumb(dummy)
         end do
         
         density(0)=6.0e0*fract/PI/1.414213562e0
         
!     if N=1372 it wont freeze below zed=0.73
         
         Ndensity=0
         Ndold=0
         zedset=density(0)
!     write(*,*)zedset*1.414213562e0*PI/6.0e0
         d0=(sqrt(2.0)/real(N,kind(fract)))**(1./3.)
         
         dia=d0*zed**(1.0e0/3.0e0)
         dsq=dia*dia
         write(*,*) "zed  zedset",zed, zedset
         
         do i = 1,n 
 10         x(i)=(randumb(dummy)-0.50e0)
            y(i)=(randumb(dummy)-0.50e0)
            z(i)=(randumb(dummy)-0.50e0)   
            
            do j=1,i
               if(j.ne.i) then
                  xi=x(i) 
                  yi=y(i)
                  zi=z(i)
                  xij=x(j)-xi
                  yij=y(j)-yi
                  zij=z(j)-zi
                  xij = xij - anint(xij)
                  yij = yij - anint(yij)
                  zij = zij - anint(zij)
                  r2=xij*xij + yij*yij + zij*zij
                  if (r2.lt.dsq) then   
                     go to 10
                  end if 
               end if
            end do
!            write(*,*)i,'got in'
         end do
         
      else

         do i=1,N
!     assign random velocities
            vx(i)=randumb(dummy)
            vy(i)=randumb(dummy)
            vz(i)=randumb(dummy)
         end do   
         
         density(0)=fract*real(1.35047447,kind(zed))
         Ndensity=0
         Ndold=0
         zedset=density(0)
         zed=zedset
         d0=(sqrt(2.0e0)/real(N,kind(zed)))**(1./3.)
         dia=d0*zed**(1./3.)
         dsq=dia*dia
         write(*,*) "zed  zedset",zed, zedset
         
         l=(N*2)**(1./3.)
         k=0
         rx=1./float(l)
         do ix=1,l
            do iy=1,l
               do iz=1,l
                  j=ix+iy+iz
                  if (mod(j,2).eq.0) then
                     k=k+1
                     x(k)=rx*ix
                     y(k)=rx*iy
                     z(k)=rx*iz
                  end if
               end do
            end do
         end do
      end if
      
      
      
      
      
      return
      end
      
!     *************************************
      subroutine ResetSums
      
      USE my_param
      IMPLICIT NONE
      
      tim=0.0e0
      timeSum=0.0e0
      psum=0.0e0
      Ncount=0.0e0
      Nstep=0.0e0
      Nrecord=0.0e0
      tcount=0.0e0
      dia=d0*zed**(1.0e0/3.0e0)
      dsq=dia*dia
      Ntrial=0.0e0
      told=0.0e0
      
      return
      end
      
      
!     ******************************************
      subroutine advance
      
      USE my_param
      IMPLICIT NONE
      
      INTEGER :: k,mt,i,lasti,lastj,i2
      REAL*8 :: Ncoll,NPRINT, NabRemake,t,nabtest,randumb,dummy


      tbig=1.0E+09
      tim=0.0e0
      diastart=dia
      Ncoll=0.0e0
      Nprint=runfor/20.0e0
      NabRemake = -2.0e0
      NabTest=real(N,kind(fract))/10.0e0
      call rescale
      call nablist
      
!     set up impending event times tnew(i)
      do k=1, N       
         call bangt(k)
      end do
      
      mt=N+1
      do while (NabRemake.lt.0.0e0.and.mt.gt.0.0e0&
     .and.Ncoll.lt.runfor)
         t=tbig-10.0e0
         mt=0
         
!     find shortest time t=tnew(i)
         
         do i=1,N                         
            if (tnew(i).lt.t) then
               t=tnew(i)
               mt=i
            end if
         end do
         if (mt.gt.0) then
            
!.....to save a config every 25n collisions(used to apply:z>0.77)
            if(((equilibrated.gt.0.0e0).and.(zed.gt.0.30e0)).and.&
           (mod(Ncount,(25.0e0*real(n,kind(fract)))).lt.1)) &
           then
               
!     update the positions of all particles by a
!     randumb fraction of the time
!     to the next collision
            
               
            tim = tim + (t-tim)*randumb(dummy)      
            
!     0< ran <1
            
            do i=1,N
               call update(x(i),y(i),z(i),vx(i),vy(i),&
              vz(i),told(i),tim)
               told(i)=tim
            end do
            call saveconfig
            call q6r

!            call torder 
            pramsum=pramsum+transord


            gq6sum=gq6sum+gq6
!            write(24,*) INT(Ncount/N), gq6

!            write(26,*) INT(Ncount/N), transord

         end if
!     ...........................................
            
            
!     identify pair and change velocities
            lasti=mt
            lastj=jc(lasti)
            if (t.gt.tim) tim=t  
!     any negative times are set to zero
            Ncount=Ncount+1.0e0
            Ncoll=Ncoll+1.0e0
            call change(lasti,lastj)
!     write(*,*) 'i kwi jci t Ncount'
!     write(*,*)lasti,kw(lasti),jc(lasti),t,Ncount
            
!     revise list of impending event times for lasti
!     lastj and any others which involve them
            
            jc(lastj)=lasti         
!     to ensurethat j is processed
            do i2=1, nay(lasti)
               k=nab(lasti,i2)
               if (jc(k).eq.lasti) call bangt(k)
            end do
            
            if (lastj.gt.0) then
               jc(lasti)=lastj         
!     to ensure that i is processed
               do i2=1, nay(lastj)
                  k=nab(lastj,i2)
                  if (jc(k).eq.lastj) call bangt(k)
               end do
            else
               call bangt(lasti)
            end if
!...................................................................
            
            if (Kompress.lt.0.0e0 .and. mod(Ncount,&
           Nprint).eq.0) call printout
            
            
            if (Ncoll.gt.Nabtest ) call &
           Test(Nabtest,NabRemake,Ncoll)
            
         end if          
         
      end do
      
!.................................................................
      Nstep=Nstep+1.0e0
      timeSum=timeSum+tim
      do i=1,N
         call update(x(i),y(i),z(i),vx(i),vy(i),vz(i),told(i),tim)
         told(i)=0.0e0
      end do
      return
      end
!     ******************************************
      subroutine test(Nabtest,NabRemake,Ncoll) 

!     to test if nablist needs to be remade
      USE my_param
      IMPLICIT NONE
      
      INTEGER :: i
      REAL*8 :: rsq,distance,dmax,nabtest,ncoll,nabremake

      NabRemake=-2.0e0
      distance=(range-dia)/2.0e0
      dmax=0.0e0
      do i=1,N
         rsq=(x(i)-x0(i))**2+(y(i)-y0(i))**2+(z(i)-z0(i))**2
         if (rsq.gt.dmax) dmax=rsq
      end do
      dmax=sqrt(dmax)
      if (dmax.gt.distance*0.9e0) then
           NabRemake=2.0e0             
!     forces remake of nablist
        else if (dmax.gt.0.0e0) then
           Nabtest=(distance*Ncoll/dmax+Ncoll)/2.0e0
           if (Nabtest.gt.(Ncoll+real(N,kind(fract)))) &
          Nabtest=Ncoll+real(N, kind(fract))
        end if
        return
        end
      
!       ***************************************************
        subroutine rescale     
! to rescale velocites and set momenta to zero
        USE my_param
        IMPLICIT NONE
        
        REAL*8 :: vsq,dvx,dvy,dvz,factor
        INTEGER :: i

        vsq=0.0e0
        dvx=0.0e0
        dvy=0.0e0
        dvz=0.0e0
!     sum momenta (mass = 1)
        do i=1,N
           dvx=dvx+vx(i)
           dvy=dvy+vy(i)
           dvz=dvz+vz(i)
        end do
        dvx=dvx/real(N,kind(fract))
        dvy=dvy/real(N,kind(fract))
        dvz=dvz/real(N,kind(fract))
      
!     set momenta to zero
        do i=1,N
           vx(i)=vx(i)-dvx
           vy(i)=vy(i)-dvy
           vz(i)=vz(i)-dvz
           vsq=vsq+vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
        end do

!     rescale to temperature (kT = 1)
        
        factor=sqrt(3.0e0*temp*real(N,kind(fract))/vsq)     
        do i=1,N
           vx(i)=vx(I)*factor
           vy(i)=vy(I)*factor
           vz(i)=vz(I)*factor
        end do
        if (abs(factor-1.0e0).gt.0.1e0) write(*,*)' factor=',factor
        return
        end

!     ******************************************
      subroutine nablist
      USE my_param
      IMPLICIT NONE
      
      INTEGER :: i,j
      REAL*8 :: rsq,rasq,dx,dy,dz,xi,yi,zi

      do i=1,N
         Nay(i)=0
         x0(i)=x(i)
         y0(i)=y(i)
         z0(i)=z(i)
      end do
      
      range=(3.0e0-2.2e0*zed)*dia
      if (zed.gt.0.6e0)  range=1.4e0*dia
      rasq=range**2
      do i=1,N-1
         xi=x(i)
         yi=y(i)
         zi=z(i)
         do j=i+1,N
            dx=xi-x(j)
            dx=dx-anint(dx)
            if (abs(dx).lt.range) then
               dy=yi-y(j)
               dy=dy-anint(dy)
               if (abs(dy).lt.range) then
                  dz=zi-z(j)
                  dz=dz-anint(dz)
                  if (abs(dz).lt.range) then
                     rsq=dx*dx+dy*dy+dz*dz
                     if (rsq.lt.rasq) then
                        nay(i)=nay(i)+1
                        nay(j)=nay(j)+1
                        nab(i,nay(i))=j
                        nab(j,nay(j))=i
                     end if
                  end if
               end if
            end if
         end do
      end do
      
      maxnay=0
      do i=1,N
         if (nay(i).gt.maxnay) maxnay=nay(i)
      end do
!     write(*,*) 'range maxnay'
!     write(*,*)range,maxnay
      return
      end
      
!     **********************************************
      subroutine printout
      USE my_param
      IMPLICIT NONE

      INTEGER :: k
      REAL*8, DIMENSION(10) :: s
      REAL*8 :: t1, m, L,n1,n2,pvort


      


      t1=timesum+tim
!.................................................
      if (t1.ne.0.0) then
         pVoRT=1.0e0+psum/(3.0e0*Real(N,kind(fract))*Temp*t1)
      end if
      M=runfor/4.0e0
      if (mod(Ncount,M).eq. 0) then
         write(*,*)'   N    Ncount   Nstep   Maxnay  fract pV/RT   '
      End if
      write(*,12)N,INT(Ncount),INT(Nstep),INT(Maxnay),&
     sqrt(2.)*zed*PI/6.0e0, pVoRT
 12   format(I7,2x,I10,2x,I8,4x,I4,2x,f6.4,2x,f20.3)


      if (equilibrated.gt.0.0e0) then
!         write(23,*) INT(Ncount/N), pVorT
         prifreq=prifreq+1
      end if


!.......collect data
      if (equilibrated.gt.0.0e0) then
         Nrecord=Nrecord+1.0e0
         L=Nrecord
         ss(INT(L),1)=t1
         ss(INT(L),2)=Psum
         
         if (Ncount.gt.Nmax-2) then
!.........................averages at end of run
            do k=1,10
               s(k)=0.0e0
            end do
            do k=1,3
               ss(0,k)=0.0e0
            end do
            
            do k=1,INT(Nrecord)
               PVoRT=1.0e0+(ss(k,2)-ss(k-1,2))/&
              (3.0e0*real(N,kind(fract))*temp*&
              (ss(k,1)-ss(k-1,1)))
               s(3)=s(3)+PVoRT
               s(4)=s(4)+PVoRT*PVoRT
            end do
            do k=3,4
               s(k)=s(k)/Nrecord
            end do
            s(4)=sqrt(s(4)-s(3)**2.0e0)
            N1=Nequil/N
            N2=Nmax/N
            
!.....................................................
!            open(unit=9,status='unknown',file='data.dat')

            
          
!            write(9,22) N,INT(N1),INT(N2),INT(Nstep),
!     :           zmax, zed*sqrt(2.)*PI/6.,s(3),s(4),ds
 22         format(4I5,4F10.5,e5.3)
            
!......................
!            close(9)
            Ncount=Nmax+5.0e0 
!     to ensure restart at the next density
            runfor = 0.0e0
!     Calculate contact value of g(r)
            gcon = (s(3)-1)/(4.0e0*fract)
         end if  
         
      end if 
!     equilibrated
      return
      end


      

!     #######################################
      subroutine saveconfig
      USE my_param 
      IMPLICIT NONE


      CHARACTER (LEN=6), PARAMETER :: fname='config'
      REAL*8 ::  dens,side,sinv,cideal,const,rlower,rupper
      REAL*8 :: RXIJ,RYIJ,RZIJ,RIJRT,RIJSQ,msd
      INTEGER :: ibin, i, j, N1shell
      INTEGER :: N2shell,N3shell,N4shell,N5shell,N6shell,N7shell
      REAL*8 :: dlat,tolr
      REAL*8 :: dlat2,dlat3,dlat4,dlat5,dlat6,dlat7
      REAL*8 :: cN1,cN2,cN3,cN4,cN5,cN6,cN7
      REAL*8 :: ccN1,ccN2,ccN3,ccN4,ccN5,ccN6,ccN7
      REAL*8 :: rNfcc, rNid, rNsys
      if(ncfg.eq.0) then
         if(istart.eq.1) then       
            istart=0
            open(unit=11,status='unknown',file=fname)
         else
            open(unit=11,position='append',file=fname)
         end if
      else
         open(unit=11,status='unknown',file=fname)
      end if
      write(11,*) n
      write(11,*) fract      
!      write(19,*) n     
!      write(19,*) fract      
      do i=1,N
         
         write(11,*)x(i)-anint(x(i)),y(i)-anint(y(i))&
        ,z(i)-anint(z(i))
!         write(19,*)x(i)-anint(x(i)),y(i)-anint(y(i))
!     :        ,z(i)-anint(z(i))
      end do
!      do i=1,N
!         write(11,*)vx(i),vy(i),vz(i)
!      end do
      close(11)
      

      



      
!.......preliminaries for g(r)
      dens=sqrt(2.0e0)*zed
      
      side = (REAL(N,kind(fract))/dens)**(1.0e0/3.0e0)
      sINV = 1.0e0/side

      
      msd=0.0e0
      do i=1,N
         msd=msd+ (x(i)-xd(i))**2.0e0+(y(i)-yd(i))**2.0e0&
        +(z(i)-zd(i))**2.0e0
      end do

!      write(33,*)INT(Ncount/real(N,kind(msd))),
!     :     msd/real(N,kind(zed))*side*side
      
      
!.......prevents rediculously large r range

      DO iBIN=1,MAXBIN
         HIST(iBIN)=0  
      END DO
      

!     dlat is nearest neighbor distance for open fcc

      dlat=(.7404805e0/fract)**(1./3.)
      dlat2=1.414214e0*dlat
      dlat3=1.73205e0*dlat
      dlat4=2.0e0*dlat
      dlat5=2.23607*dlat
      dlat6=2.4495*dlat
      dlat7=2.64575*dlat

!      write(*,*)dlat,dlat2,dlat3,dlat4,dlat5
      tolr=.098e0*dlat
      
      transord = 0.0e0
      N1shell = 0
      N2shell = 0
      N3shell = 0
      N4shell = 0
      N5shell = 0
      N6shell = 0
      N7shell = 0

!     ****BUILD HISTOGRAM
      
      DO 100 I=1,N-1
         DO 99 J = I+1,N
            RXIJ = X(I)-X(J)
            RYIJ = Y(I)-Y(J)
            RZIJ = Z(I)-Z(J)
            RXIJ = RXIJ - ANINT(RXIJ)
            RYIJ = RYIJ - ANINT(RYIJ) 
            RZIJ = RZIJ - ANINT(RZIJ)
         RIJSQ=RXIJ*RXIJ +RZIJ*RZIJ + RYIJ*RYIJ

!        scale to particle size
         RIJRT=(RIJSQ)**(0.5e0)*side

         if((RIJRT.gt.(dlat-tolr)).and.(RIJRT.lt.(dlat+tolr))) then
            N1shell=N1shell+2
         end if

         if((RIJRT.gt.(dlat2-tolr)).and.&
        (RIJRT.lt.(dlat2+tolr))) then
            N2shell=N2shell+2
         end if

         if((RIJRT.gt.(dlat3-tolr)).and.&
        (RIJRT.lt.(dlat3+tolr))) then
            N3shell=N3shell+2
         end if
         
         if((RIJRT.gt.(dlat4-tolr)).and.&
        (RIJRT.lt.(dlat4+tolr))) then
            N4shell=N4shell+2
         end if

         if((RIJRT.gt.(dlat5-tolr)).and.&
        (RIJRT.lt.(dlat5+tolr))) then
            N5shell=N5shell+2
         end if         

         if((RIJRT.gt.(dlat6-tolr)).and.&
        (RIJRT.lt.(dlat6+tolr))) then
            N6shell=N6shell+2
         end if         

         if((RIJRT.gt.(dlat7-tolr)).and.&
        (RIJRT.lt.(dlat7+tolr))) then
            N7shell=N7shell+2
         end if         

         iBIN = INT(RIJRT/DELR) + 1
         IF (iBIN.LE.MAXBIN) THEN
            HIST(iBIN) = HIST(iBIN) + 2
         END IF
         
 99   CONTINUE
 100  CONTINUE
      
        CONST = 4.0e0 * PI * dens/3.0e0
        
        DO 110 iBIN=1,MAXBIN
           RLOWER= REAL(iBIN-1)*DELR
           RUPPER= RLOWER + DELR
           cIDEAL = CONST*(RUPPER**3-RLOWER**3)
           GRA(iBIN) = REAL(HIST(iBIN))/(REAL(N)*cIDEAL)
           grat(ibin)=grat(ibin) + gra(ibin)   
 110    CONTINUE
        
        cN1=const*((dlat+tolr)**3.-(dlat-tolr)**3.)
        cN2=const*((dlat2+tolr)**3.-(dlat2-tolr)**3.)
        cN3=const*((dlat3+tolr)**3.-(dlat3-tolr)**3.)
        cN4=const*((dlat4+tolr)**3.-(dlat4-tolr)**3.)
        cN5=const*((dlat5+tolr)**3.-(dlat5-tolr)**3.)
        cN6=const*((dlat6+tolr)**3.-(dlat6-tolr)**3.)
        cN7=const*((dlat7+tolr)**3.-(dlat7-tolr)**3.)

        rNfcc = real(134.,kind(fract))
        rNid  = cN1+cN2+cN3+cN4+cN5+cN6+cN7
        rNsys = real(N1shell+N2shell+N3shell+N4shell&
       +N5shell+N6shell+N7shell,kind(fract))&
       /real(n,kind(fract))


!        ccN1 = abs((real(N1shell,kind(fract))/
!     :       real(N,kind(fract))-cN1)/(12.0e0-cN1))
!
!        ccN2 = abs((real(N2shell,kind(fract))/
!     :       real(N,kind(fract))-cN2)/(6.0e0-cN2))
!
!        ccN3 = abs((real(N3shell,kind(fract))/
!     :       real(N,kind(fract))-cN3)/(24.0e0-cN3))
!        
!        ccN4 = abs((real(N4shell,kind(fract))/
!     :       real(N,kind(fract))-cN4)/(12.0e0-cN4))
!        
!        ccN5 = abs((real(N5shell,kind(fract))/
!     :       real(N,kind(fract))-cN5)/(24.0e0-cN5))        

        transord = (rNsys-rNid)/(rNfcc-rNid)

!        write(*,*)transord,N1shell/real(n),N2shell/real(N)
!     :       ,N3shell/real(n),n4shell/real(n),n5shell/real(n),
!     :       N6shell/real(n),n7shell/real(n)
        
        return
        end
      
!     *****************************************************

 	SUBROUTINE q6r
	
        USE my_param

        IMPLICIT NONE
	integer, DIMENSION(N,200) :: aneib
        INTEGER, DIMENSION(N) :: ancount
        real*8, DIMENSION(N) :: q6
	real*8, DIMENSION(13,N) :: qreal,qimag
        real*8, DIMENSION(0:12) :: factl
	real*8, DIMENSION(13) ::  qr,qi
        REAL*8  :: RXI, RYI, RZI, dx, dy,dz,R2
        REAL*8  :: SIGSQ,binwidth,rcutsq,qlm2
        INTEGER  :: I, J,mm,k,ibondsum,mabs,ibin,ii,jj
        INTEGER  :: nbond
        REAL*8 :: d2,dist,rqcutsq
!     ***********************************************************


        binwidth = (1.0e0/100.0e0)
        

        rqcutsq =  RQ6NEIGHB* RQ6NEIGHB
!       **create relevant factorials**   

	factl(0)=1.0e0
        do mm = 1, 12
	   factl(mm)=factl(mm-1)*real(mm,KIND(rxi))
          
	enddo
       
       
        
!     initialize aneib

        aneib=0
        
        

        qreal= 0.0e0
        qimag= 0.0e0
        
        

!     Make up Neighbor lists and count...
        
        ancount=0
        

        SIGSQ  = diaset * diaset

        
!     ** LOOPS OVER MOLECULES **
        
        DO 100 I = 1, N
           
           RXI = X(I)
           RYI = Y(I)
           RZI = Z(I)
           
           DO 99 J=1,N
              if(i.ne.j)then
                 DX= X(J)-RXI
                 DY= Y(J)-RYI
                 DZ= Z(J)-RZI
                 
                 DX= DX -ANINT(DX)
                 DY= DY -ANINT(DY)
                 DZ= DZ -ANINT(DZ)
                 
                 dist = dx*dx+dy*dy+dz*dz
                 
                 if(dist.lt.rqcutsq)then
                    
                    ancount(i)=ancount(i) +1
                    aneib(i,ancount(i))=j
                 end if
	      end if
              
 99        CONTINUE
 100    CONTINUE
        
        

!        write(*,*)'made neighbor list in q6'

!	**calcaulate total # of bonds

	ibondsum=0

	do i=1,n
	  ibondsum=ibondsum + ancount(i)	
	end do

	ibondsum=int(real(ibondsum)/2.0e0) 
	
	
!       **Begin loop over central atoms
        nbond=0
        q6=0.0e0
        do i = 1,N
          
!          write(*,*) 'for i = ',i,' ancount = ',ancount(i)
        
           do j = 1,ancount(i)
!     calc phi,theta
              
              dx=x(aneib(i,j))-x(i)
	      dy=y(aneib(i,j))-y(i)
	      dz=z(aneib(i,j))-z(i)
              
	      dx=dx-anint(dx)
	      dy=dy-anint(dy)
              dz=dz-anint(dz)
              
              
	      d2=dx*dx +dy*dy +dz*dz
              
	      call legendre(dx,dy,dz,d2,qreal,qimag,nbond,i,N)
      	   enddo
           


	   do k=1,13
	      qlm2=qreal(k,i)*qreal(k,i)+qimag(k,i)*qimag(k,i)
	      mm=k-7
	      mabs=abs(mm)
	      q6(i)=q6(i)+qlm2*factl(6-mabs)/factl(6+mabs)
           enddo
	enddo
	
        
	qr=0.0e0
	qi=0.0e0
        do k=1,13
           do i=1,N
              qr(k)=qr(k)+qreal(k,i)
              qi(k)=qi(k)+qimag(k,i)
           end do
        end do
	  
	gq6=0.0e0
        
	do k=1,13
	   qlm2=qr(k)*qr(k)+qi(k)*qi(k)
!     write(*,*)qlm2
	   mm=k-7
	   mabs=abs(mm)
	   gq6=gq6+qlm2*factl(6-mabs)/factl(6+mabs)
	enddo
        
	if (nbond .eq. 0) then
	   gq6=0.0e0
	else
	   gq6=sqrt(gq6)/real(nbond,KIND(rxi))
	endif
        
	
	do i=1,N
	   if (ancount(i).eq.0) then
              q6(i)=0.0
              iq6dist(1)=iq6dist(1)+1
	   else
              q6(i)=sqrt(q6(i))/real(ancount(i),kind(rxi))
              ibin= int(q6(i)/real(binwidth,kind(rxi)))
              iq6dist(ibin)=iq6dist(ibin)+1
	   endif
           
!     write(*,*)i,q6(i)
	end do
!     write(*,*)'global q6=',gq6
!     write(*,*)'# of neighbs/molecule',real(nbond)/real(N)
        avgbond = avgbond+real(nbond,kind(rxi))/real(N,kind(rxi))
	end
      
!     *****************************************************
      
      
!     ****************************************************
!     **calculate Qlm between neighbors                 **
!     **determine polar and azimuthal angles of vector r**
!     **with respect to x,y,z axis (center of box)    
!     **z=cos(theta)                                    **
!     ****************************************************
      
      subroutine legendre(rxij,ryij,rzij,rijsq,qreal&
     ,qimag,nbond,i,N)
      IMPLICIT NONE
      INTEGER :: i,n,nbond
      real*8, DIMENSION(13,N) :: qreal,qimag
      real*8 :: phi, rij,sinth,x,y,z,rxij,ryij,rzij
      real*8 :: p61,p62,p63,p64,p65,p66,argc,args,rijsq


      rij=sqrt(rijsq)
      x = rxij/rij
      y = ryij/rij
      z = rzij/rij
      
      
      sinth = real(sqrt(1.0e0 - z*z),kind(rxij))
      if (sinth .eq. 0.0e0) then
         phi=0.0e0
      else
         if (abs(x/sinth) .gt. 1.0e0) then
            phi=0.0e0
         else
            phi=atan2(x,y)
         endif
      endif
      
      
      nbond=nbond+1
      
      qreal(7,i)=qreal(7,i)+(231.0e0*z**6.0e0-315.0e0*z**4.0e0+&
     105.0e0*z**2.0e0-5.0e0)/16.0e0
      
      p61=21.0e0*sinth*(33.0e0*z**5.0e0-30.0e0*z**3.0e0+5.0e0*z)&
     /8.0e0
      argc=cos(phi)
      args=sin(phi)		
      qreal(6,i)=qreal(6,i)-p61*argc
      qimag(6,i)=qimag(6,i)-p61*args
      
      qreal(8,i)=qreal(8,i)+p61*argc
      qimag(8,i)=qimag(8,i)-p61*args
      
      p62=105.0e0*sinth*sinth*(33.0e0*z**4.0e0&
     -18.0e0*z**2.0e0+1.0e0)/8.0e0
      argc=cos(2.0e0*phi)
      args=sin(2.0e0*phi)		
      qreal(5,i)=qreal(5,i)+p62*argc
      qimag(5,i)=qimag(5,i)+p62*args
      
      qreal(9,i)=qreal(9,i)+p62*argc
      qimag(9,i)=qimag(9,i)-p62*args
      
      p63=105.0e0*(sinth**3.0e0)*(33.0e0*z**3.0e0-9.0e0*z)/2.0e0
      argc=cos(3.0e0*phi)
      args=sin(3.0e0*phi)		
      qreal(4,i)=qreal(4,i)-p63*argc
      qimag(4,i)=qimag(4,i)-p63*args
      qreal(10,i)=qreal(10,i)+p63*argc
      qimag(10,i)=qimag(10,i)-p63*args
      
      p64=945.0e0*(sinth**4.0e0)*(11.0e0*z**2.0e0-1.0e0)/2.0e0
      argc=cos(4.0e0*phi)
      args=sin(4.0e0*phi)		
      qreal(3,i)=qreal(3,i)+p64*argc
      qimag(3,i)=qimag(3,i)+p64*args
      qreal(11,i)=qreal(11,i)+p64*argc
      qimag(11,i)=qimag(11,i)-p64*args
      
      p65=10395.0e0*z*sinth**5.0e0
      argc=cos(5.0e0*phi)
      args=sin(5.0e0*phi)		
      qreal(2,i)=qreal(2,i)-p65*argc
      qimag(2,i)=qimag(2,i)-p65*args
      qreal(12,i)=qreal(12,i)+p65*argc
      qimag(12,i)=qimag(12,i)-p65*args
      
      p66=10395.0e0*sinth**6.0e0
      argc=cos(6.0e0*phi)
      args=sin(6.0e0*phi)		
      qreal(1,i)=qreal(1,i)+p66*argc
      qimag(1,i)=qimag(1,i)+p66*args
      qreal(13,i)=qreal(13,i)+p66*argc
      qimag(13,i)=qimag(13,i)-p66*args
      
      return
      end
      

      SUBROUTINE TORDER 
      USE my_param
      IMPLICIT NONE

      
      REAL*8  ::    KLATX, KLATY, KLATZ
      
      
      INTEGER ::    I
      REAL*8   ::     SINSUM, COSSUM,CONST
      
      
!     *************
      
      SINSUM = 0.0e0
      COSSUM = 0.0e0
      CONST=PI*(2*REAL(N,KIND(KLATX)))**(1./3.)
      
      KLATX = -1.0e0*CONST
      KLATY = CONST
      KLATZ = -1.0e0*CONST
      
      DO 100 I = 1, N
         
         COSSUM = COSSUM + COS (  KLATX * X(I)&
        + KLATY * Y(I)&
        + KLATZ * Z(I) )
         SINSUM = SINSUM + SIN (  KLATX * X(I)&
        + KLATY * Y(I)&
        + KLATZ * Z(I) )
         
 100  CONTINUE

      COSSUM = COSSUM / REAL ( N,KIND(KLATX) )
      SINSUM = SINSUM / REAL ( N,KIND(KLATX) )
      PRAM  = SQRT ( COSSUM ** 2 + SINSUM ** 2 )
      
      RETURN
      END
      

