!!****m* BrownianDynamics/PositionParticles
!!
!! NAME
!!  PositionParticles
!!
!! SYNOPSIS
!!  This subroutine positions all the particles within the simulation
!!  by calling the spheres module to generate a configuration of spheres
!!  in a cubic simulation cell and then replicating that cell periodically
!!  to fill the actual simulation cell completely.  Additionally, the
!!  externally forced particles are positioned in a user specified location.
!!
!! NOTES
!!  Located in file gen.f90
!!
!! CREATION DATE
!!    02/28/05
!!
!! AUTHOR 
!!    James Swan
!!
!! COPYRIGHT
!!    28-Feb-2005 - James William Swan
!!******
!subroutine PositionParticles
!  use Globals
!  implicit none
!  
!  ! Loop counters.
!  Integer :: nPartGen, i, j, k, l, c
!  !*** For the case phi_bat and phi_motor
!  Integer :: nPartGen1, nPartGen2
!  ! The cubic configuration of spheres.
!  type ( Vector ), allocatable, dimension( : ) :: smallConfig
!  ! Holder vector.
!  type ( Vector ) ::  hold, dvecs, MinImage
!  ! Cubic cell dimensions.
!  Real( 8 ) :: volume, hcell, area
!  Real, allocatable, dimension( : ) :: u
!  
!  Integer( 8 ), dimension( 3 ) :: time
!  Real( 8 ) :: seed0
!  Real( 8 ) :: radID1m, radID2m, radfract1, radfract2
!  
!  gLCells = gNCells( 1 ) * gNCells( 2 ) * gNCells( 3 )   
!  
!  if ( g2Dproblem .eq. 0 ) then
!     
!     !***For three dimension cell
!     
!     !   	nPartGen = gNPart / gCells
!     !    	volume = gPi43 * dble( nPartGen ) / gPhi
!     !   	hcell = volume ** ( 1.0 / 3.0 )
!     
!     nPartGen1 = gNMove / gLCells
!     nPartGen2 = (gNPart-gNMove*gPhip) / gLCells
!     nPartGen = nPartGen1*gPhip + nPartGen2 
!     volume = (gPi43/(1+gPhip))*(gPhip*dble(nPartGen1)/gPhim + dble(nPartGen2)/(gPhi*(1-gPhip)+gPhib*gPhip)) 
!     hcell = volume ** ( 1.0 / 3.0 )
!     
!     !        write(*,*) nPartGen1, nPartGen2, nPartGen
!     !        write(*,*) hcell
!     !        stop
!     
!     allocate( smallConfig( nPartGen ) )
!     
!     if ( useconfig2 .eq. 0 ) then  ! Use equlibrium configuration  0 - NO , 1 - YES
!        
!        ! Generate the cubic cell configuration
!        call spheres90( nPartgen, gPhi )
!        !print*, "After call spheres90"
!        
!        !        write(*,*)'Hasta aqui esta bien'
!        !        stop
!        
!        ! Read in the cubic cell configuration
!        open( unit = 11, file = 'config', status = 'old' )
!        
!        
!        read( 11, * )
!        read( 11, * )
!        
!        do i = 1, nPartGen
!           read( 11, * ) hold % x, hold % y, hold % z
!           smallConfig( i ) = hcell * ( hold + Vector( 0.5, 0.5, 0.5 ) )
!           gPart( i ) % pos = smallConfig( i )
!        end do
!        
!        !        write(*,*)'Hasta aqui esta bien'
!        !        stop
!        
!        close( 11 )
!        
!        c = 0
!        
!        ! Replicate the cubic configuration in the simulation cell.
!        do i = 1, gNCells( 1 )
!           do j = 1, gNCells( 2 )
!              do k = 1, gNCells( 3 )
!                 if ( ( i .ne. 1 ) .or. ( j .ne. 1 ) .or. ( k .ne. 1 ) ) then
!                    c = c + 1
!                    do l = 1, nPartGen
!                       gPart( l + c * nPartGen ) % pos % x = smallConfig( l ) % x &
!                            + ( i - 1 ) * hcell 
!                       gPart( l + c * nPartGen ) % pos % y = smallConfig( l ) % y &
!                            + ( j - 1 ) * hcell 
!                       gPart( l + c * nPartGen ) % pos % z = smallConfig( l ) % z &
!                            + ( k - 1 ) * hcell 
!                    end do
!                 end if
!              end do
!           end do
!        end do
!        
!     else
!        
!        !*******************************************************************************************************
!        open( unit = 1600, file = 'config2') 
!        
!        do i = 1, gNPart
!           
!           read( 1600, * ) gPart( i ) % pos % x, gPart( i ) % pos % y, gPart( i ) % pos % z
!           
!        end do
!        !*******************************************************************************************************
!        
!     end if
!     
!     
!  else     ! This part executes when  g2DProblem = 1
!     
!     !****** For two dimension cell
!     
!     !    	volume = gPi43 * gNPart / gPhi
!     
!     volume = (gPi43/(1+gPhip))*(gPhip*gNMove/gPhim+(gNPart-gNMove*gPhip)/(gPhi*(1-gPhip)+gPhib*gPhip))
!     
!     area = DBLE(gNpart) * gPi / gPhi
!     
!     !hcell = volume ** ( 1.d0 / 3.d0 )
!      hcell = dsqrt( area  )
!     
!     
!     
!     !volume = volume / gLCells
!     
!     area = area / gLCells
!     
!     !gCDim( 1 ) = gNCells( 1 ) * volume ** ( 1.d0 / 3.d0 )
!     !gCDim( 2 ) = gNCells( 2 ) * volume ** ( 1.d0 / 3.d0 )
!     !gCDim( 3 ) = gNCells( 3 ) * volume ** ( 1.d0 / 3.d0 )
!     
!     gCDim( 1 ) = gNCells( 1 ) * area ** ( 1.D0 / 2.D0 )
!     gCDim( 2 ) = gNCells( 2 ) * area ** ( 1.D0 / 2.D0 )
!     gCDim( 3 ) = gNCells( 3 ) * area ** ( 1.D0 / 2.D0 )    
!
!
!     call RandomVectors
!     allocate( u( 3 * gNPart ) )
!     call iTime( time )
!     seed0 = 10000 * time( 3 ) + 100 * time( 2 ) + time( 1 )
!     seed0 = seed0 + time( 1 ) * time( 2 ) * ( time( 3 ) + time( 1 ) ) + time( 3 )
!     call random_number( harvest = seed0 )
!     call random_number( u )
!     
!     do i = 1, gNPart
!        gPart( i ) % pos % x = u( 3 * ( i - 1 ) + 1 ) * gCDim( 1 )
!        gPart( i ) % pos % y = u( 3 * ( i - 1 ) + 2 ) * gCDim( 2 )
!        gPart( i ) % pos % z = 5.D0 ! was 3.0
!     end do
!     
!  end if
!  
!  !************************************************************    
!  ! Adjust the cell volume to account for larger particles.
!  ! For the case where a / b > 1
!  !***********************************************************
!  !    volume = gNPart
!  
!  !    do i = 1, gNMove
!  !        volume = volume + ( gPart( i ) % rat ** 3. ) - 1. 
!  !    end do
!  
!  !    volume = volume * gPi43 / gPhi
!  
!  !volume = (gNPart-gNMove*gPhip)/((gPhib-1)*gPhip+1)
!  
!  !do i = 1, gNMove
!  !   volume = volume + ( gPart(i)%rat**3.D0)/((gPhim-1)*gPhip+1)-1.D0*(1-gPhip)
!  !end do
!  
!  !volume = volume * ( gPi43 / (gPhi*(1-gPhip)+2*gPhip) )  
!  
!  !gCDim( 1 ) = hcell * gNCells( 1 )
!  !gCDim( 2 ) = hcell * gNCells( 2 )
!  !gCDim( 3 ) = hcell * gNCells( 3 )
!  
!  !volume = ( volume / ( gCDim( 1 ) * gCDim( 2 ) * gCDim( 3 ) ) ) ** ( 1.D0 / 3.D0 )
!  
!  !gCDim( 1 ) = gCDim( 1 ) * volume
!  !gCDim( 2 ) = gCDim( 2 ) * volume
!  !gCDim( 3 ) = gCDim( 3 ) * volume
!  
!  !****Recalculate the volume fraction of the system
!  
!  !gPhis = gNPart
!  
!  !do i = 1, gNMove
!  !   gPhis = gPhis + ( gPart(i) % rat ** 3.D0 ) - 1.D0
!  !end do
!  
!  !gPhis = gPhis * gPi43 / hcell**3.D0
!  !****
!  
!  !***Recalculate the correct values of volume fractions
!  
!  !gPhinm = 0
!  !do i = 1, gNMove
!  !   gPhinm = gPhinm+(gPart(i)%rat**3)
!  !end do
!  
!  !gPhinm = (gPi43/hcell**3)*gPhinm 
!  
!  !gPhinb = (gPi43/hcell**3)*(gNPart-gNMove)
!  !!***************************************************** 
!  
!  gCDimH( 1 ) = gCDim( 1 ) / 2.D0
!  gCDimH( 2 ) = gCDim( 2 ) / 2.D0
!  gCDimH( 3 ) = gCDim( 3 ) / 2.D0
!  
!  
!  gMinDim = dsqrt( gCDim( 1 ) ** 2 + gCDim( 2 ) ** 2 + gCDim( 3 ) ** 2 ) / 2.D0
!  
!  gDPXFixed = gNCells( 1 ) * gCDim( 1 ) / ( 1. + gNCells( 1 ) )
!  gDPDim( 1 ) = gDPBins( 1 ) / gCDim( 1 )
!  gDPDim( 2 ) = gDPBins( 2 ) / gCDim( 2 )
!  gDPDim( 3 ) = gDPBins( 3 ) / gCDim( 3 )
!  
!  
!  ! Position the externally forced particles.
!  
!  if ( useconfig2 .eq. 0 ) then
!     
!     do i = 1, gNMove
!        gPart( i ) % pos % x = gPart( i ) % init % x * gCDim( 1 )
!        gPart( i ) % pos % y = gPart( i ) % init % y * gCDim( 2 )
!        gPart( i ) % pos % z = gPart( i ) % init % z * gCDim( 3 )
!        
!        if ( g2Dproblem .eq. 1 ) then
!           gPart( i ) % pos % z = 5.D0
!        end if
!        
!     end do
!     
!  end if
!  
!  ! Store the absolute initial position of particle 1 for later.
!  
!  do i = 1, gNMove
!     gPart( i ) % init % x = gPart( i ) % init % x * gCDim( 1 )
!     gPart( i ) % init % y = gPart( i ) % init % y * gCDim( 2 )
!     gPart( i ) % init % z = gPart( i ) % init % z * gCDim( 3 )
!     
!     if ( g2Dproblem .eq. 1 ) then
!        gPart( i ) % init % z = 5.D0
!     end if
!     
!  end do
!  
!  radID1m = gPart( 1 ) % rat
!  radID2m = gPart( 2 ) % rat
!  
!  radfract1 = 1.D0 / ( 1.0 + ( radID2m / radID1m ) ** 3.D0 )
!  radfract2 = 1.D0 / ( 1.0 + ( radID1m / radID2m ) ** 3.D0 )
!  
!  
!  dvec( 1 ) = Vector( ( ( 1.0 - radfract1 ) * gPart( 1 ) % pos % x - radfract2 * gPart( 2 ) % pos % x ), &
!       ( ( 1.0 - radfract1 ) * gPart( 1 ) % pos % y - radfract2 * gPart( 2 ) % pos % y ), &
!       ( ( 1.0 - radfract1 ) * gPart( 1 ) % pos % z - radfract2 * gPart( 2 ) % pos % z ) )
!  
!  dvec( 2 ) = Vector( ( ( 1.0 - radfract2 ) * gPart( 2 ) % pos % x - radfract1 * gPart( 1 ) % pos % x ), &
!       ( ( 1.0 - radfract2 ) * gPart( 2 ) % pos % y - radfract1 * gPart( 1 ) % pos % y ), &
!       ( ( 1.0 - radfract2 ) * gPart( 2 ) % pos % z - radfract1 * gPart( 1 ) % pos % z ) )
!  
!  
!  dvecs = MinImage( dvec( 1 ) - dvec( 2 ) )
!  
!  dfixed = dvecs * dvecs
!  
!  
!  ! Calculate initial center of mass for same size motors
!  
!  dcmfixed = Vector( ( radfract1 * gPart( 1 ) % pos % x + radfract2 * gPart( 2 ) % pos % x ), &
!       ( radfract1 * gPart( 1 ) % pos % y + radfract2 * gPart( 2 ) % pos % y ), &
!       ( radfract1 * gPart( 1 ) % pos % z + radfract2 * gPart( 2 ) % pos % z ) )
!  
!  
!  !write( *, * ) dcmfixed
!  
!  if ( gNMove .eq. 2 ) then
!     write( *, * ) 'Distance between probes: ', dsqrt( dvecs * dvecs )
!  end if
!  
!  
!end subroutine PositionParticles


!!****f* BrownianDynamics/ChangeForceMode
!!
!! NAME
!!  UpdateParticles
!!
!! SYNOPSIS
!!  This subroutine changes the force mode governing the hardsphere
!!  interaction between the externally forced particles and the
!!  Brownian ones.  This is used on the fly to produce certain
!!  equilibrium externally forced particle configurations.
!!
!! USAGE
!!  subroutine ChangeForceMode( i )
!!
!! INPUTS
!!  i - The new force mode.
!!
!! NOTES
!!  Located in file gen.f90
!!
!! CREATION DATE
!!    02/28/05
!!
!! AUTHOR 
!!    James Swan
!!
!! COPYRIGHT
!!    28-Feb-2005 - James William Swan
!!******
subroutine ChangeForceMode( i )
  use Globals
  implicit none
  
  ! The new force mode
  Integer, intent( in ) :: i
  
  gFMTemp = gForceMode
  gForceMode = i
  
end subroutine ChangeForceMode

!!****f* BrownianDynamics/InitPosHex2D
!!
!! NAME
!!  InitPosHex2D
!!
!! SYNOPSIS
!!  This subroutine positions the particle in a hexagonal
!!  two-dimension lattice. It calulates the distance "a"
!!  between particles and the size of the simulation box
!!  based on the number of particles N and the number den-
!!  sity "n".   
!!
!! USAGE
!!  call InitPosHex2D
!!
!! INPUTS
!!  N/A
!!
!! NOTES
!!  Located in file gen.f90
!!
!! CREATION DATE
!!  01/28/14
!!
!! AUTHOR 
!!    Luis Y. Rivera-Rivera
!!
!! COPYRIGHT
!!    28-Jan-2014 - Luis Rivera
!!******
subroutine InitPosHex2D
   Use Globals
   implicit none

   integer :: KX, KY, K, P, NP, IFACE, I, J
   real(8) :: XL, YL, NDENS, RCOFF
   real(8) :: A, AX, AY, C1
   real(8) :: RX0, RY0, RXI, RYI

  NDENS = gPhi / gPi 
   
  A = ( ( 2.D0 / 3.D0 ** 0.5 ) / NDENS ) ** 0.5
 
  P = IDNINT( DSQRT(DBLE( gNPart / 4)  ))
  
  XL = 3.D0**0.5*A*DBLE(P)
  YL = 2.D0 * A * DBLE(P) 
  
  AX = 3.D0**0.5*A
  AY = 2.D0 * A
  
  KX = P
  KY = P
  C1 = 0.01D0

  K = 0
  
  DO IFACE = 1, 4
     IF ( IFACE .EQ. 1 ) THEN
        RX0 = C1
        RY0 = C1
     ELSE IF (IFACE .EQ. 2) THEN
        RX0 = C1
        RY0 = C1 + A
     ELSE IF (IFACE .EQ. 3) THEN
        RX0 = AX /2.D0 + C1
        RY0 = A/2.D0 + C1
     ELSE IF (IFACE .EQ. 4) THEN
        RX0 = AX /2.D0 + C1
        RY0 = A * 3.D0 / 2.D0 + C1
     END IF
   
     DO J = 0, KY - 1
        RYI = DBLE(J)*AY + RY0
        IF ( RYI .GE. YL ) CYCLE
        DO I = 0, KX - 1
           RXI = DBLE(I) * AX + RX0
           IF ( RXI .GE. XL ) CYCLE
           K = K + 1
           
           gPart(K) % pos % x = RXI
           gPart(K) % pos % y = RYI
           gPart(K) % pos  % z = 0.D0 !YL / 2.D0 
        END DO
     END DO
   END DO
   
   gCDim(1) = XL
   gCDim(2) = YL
   gCDim(3) = YL  
   
  gCDimH( 1 ) = gCDim( 1 ) / 2.D0
  gCDimH( 2 ) = gCDim( 2 ) / 2.D0
  gCDimH( 3 ) = gCDim( 3 ) / 2.D0

  
end subroutine InitPosHex2D





!!****f* BrownianDynamics/InitPosSquare2D
!!
!! NAME
!!  InitPosSquare2D
!!
!! SYNOPSIS
!!  This subroutine positions the particle in a squared
!!  two-dimension lattice. It calulates the distance "a"
!!  between particles and the size of the simulation box
!!  based on the number of particles N and the number den-
!!  sity "n".   
!!
!! USAGE
!!  call InitPosSquare2D
!!
!! INPUTS
!!  N/A
!!
!! NOTES
!!  Located in file gen.f90
!!
!! CREATION DATE
!!  01/28/14
!!
!! AUTHOR 
!!    Luis Y. Rivera-Rivera
!!
!! COPYRIGHT
!!    28-Aug-2014 - Luis Rivera
!!******
subroutine InitPosSquare2D
   Use Globals
   implicit none

   integer :: KX, KY, K, P, NP, IFACE, I, J
   real(8) :: XL, YL, NDENS, RCOFF
   real(8) :: A, AX, AY, C1
   real(8) :: RX0, RY0, RXI, RYI

  NDENS = gPhi / gPi

  A =  1.D0 / DSQRT(NDENS)  ! A = ( ( 2.D0 / 3.D0 ** 0.5 ) / NDENS ) ** 0.5

  P = IDNINT( DSQRT(DBLE( gNPart)  ))

  !XL = 3.D0**0.5*A*DBLE(P)
  !YL = 2.D0 * A * DBLE(P)
   
  XL = A * DBLE(P)
  YL = A * DBLE(P)

  AX = A   !DSQRT(3) * A
  AY = A   !2.D0 * A

  KX = P
  KY = P
  C1 = 0.01D0

  K = 0

  DO IFACE = 1, 4
     IF ( IFACE .EQ. 1 ) THEN
        RX0 = C1
        RY0 = C1
     ELSE IF (IFACE .EQ. 2) THEN
        RX0 = C1
        RY0 = C1 + A
     ELSE IF (IFACE .EQ. 3) THEN
        RX0 = AX + C1
        RY0 = A  + C1
     ELSE IF (IFACE .EQ. 4) THEN
        RX0 = AX + C1
        RY0 = A  + C1
     END IF

     DO J = 0, KY - 1
        RYI = DBLE(J)*AY + RY0
        IF ( RYI .GE. YL ) CYCLE
        DO I = 0, KX - 1
           RXI = DBLE(I) * AX + RX0
           IF ( RXI .GE. XL ) CYCLE
           K = K + 1

           gPart(K) % pos % x = RXI
           gPart(K) % pos % y = RYI
           gPart(K) % pos % z = 0.D0 !YL / 2.D0 
        END DO
     END DO
   END DO

   gCDim(1)    = XL
   gCDim(2)    = YL
   gCDim(3)    = YL
   gCDimH( 1 ) = gCDim( 1 ) / 2.D0
   gCDimH( 2 ) = gCDim( 2 ) / 2.D0
   gCDimH( 3 ) = gCDim( 3 ) / 2.D0
  
  gCutOff = gCDim(1) / 2.D0
  gCutOff2 = gCutOff ** 2.D0 
  
  gReciprocalVec % x = gPi2 / A
  gReciprocalVec % y = gPi2 / A
  gReciprocalVec % z = gPi2 / A

END SUBROUTINE InitPosSquare2D
