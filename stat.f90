!!****m* BrownianDynamics/Statistics
!!
!! NAME
!!  Statistics
!!
!! SYNOPSIS
!!  This subroutine calls all the statistical routines.
!!
!! NOTES
!!  Located in file stat.f90
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
subroutine Statistics
use Globals
implicit none

    Integer :: time ! Counter.

    time = gT * gDTInv


    if ( gT .ge. gAvgTiming ) then

        call OsmPress
        call PartVel
!        call PairDist
!        call RadDist
        
        if ( gForceMode .eq. 3 ) then
            write( *, * ) 'Force mode changed'
            call ChangeForceMode( 1 )
        end if

        if ( mod( time, gDenTiming ) .eq. 0 ) then
            call DenProfile
        end if
        
    else
    
    	gHSF( : ) = Vector( 0.d0, 0.d0, 0.d0 );
    	gHST( : ) = Vector( 0.d0, 0.d0, 0.d0 );
    	gHSFRh( : ) = Vector( 0.d0, 0.d0, 0.d0 );
    	gHSFnRh( : ) = Vector( 0.d0, 0.d0, 0.d0 );
    	
        gOPMotor( : ) = 0.0;
        gFCM = Vector( 0.d0, 0.d0, 0.d0 );
        gTensileF( : ) = Vector( 0.d0, 0.d0, 0.d0 );

    end if

	! if you add the following, the HS force is reseted at each time step and it mess up the averaging
    ! gHSF( : ) = Vector( 0.d0, 0.d0, 0.d0 )

    call ZeroMeasures

end subroutine Statistics

!!****m* BrownianDynamics/OsmPress
!!
!! NAME
!!  OsmPress
!!
!! SYNOPSIS
!!  This subroutine records the osmotic pressure in the fluid and updates
!!  the average.
!!
!! NOTES
!!  Located in file stat.f90
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
subroutine OsmPress
use Globals
implicit none

    gOP = 1. + gOP / ( 3.d0 * gNPart * gDT )
    

    gAP = gAP + gOP

end subroutine OsmPress

!!*********************************************************************************!!
!!****m* BrownianDynamics/PartVel
!!
!! NAME
!!  PartVel
!!
!! SYNOPSIS
!!  This subroutine records the velocity of all the externally forced
!!  particles and updates the average velocity.
!!
!! NOTES
!!  Located in file stat.f90
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
subroutine PartVel
use Globals
implicit none

    Integer :: i, j, iu1, iu2, iu3 ! Counter and unit of files
    type ( Vector ) diff, MinImage, Contain, disp ! Displacements.
    ! Holder variables.
    Real( 8 ) :: tScale2
    
    Real( 8 ) :: radID1m, radID2m, radfract1, radfract2
    type ( Vector ) :: dcenter
    type ( Vector ), allocatable, dimension( : ) :: reldist
    Real( 8 ), allocatable, dimension( : ) :: gVelx, gVely, gVelz, rdx, rdy, rdz
    Real( 8 ), allocatable, dimension( : ) :: magnr, magnrsq
    
    allocate( reldist( 3*gNPart ) ) !change gNMove by gNPart
    allocate( gVelx( 3*gNPart ) )
    allocate( gVely( 3*gNPart ) )
    allocate( gVelz( 3*gNPart ) )
    allocate( rdx( 3*gNPart ) )
    allocate( rdy( 3*gNPart ) )
    allocate( rdz( 3*gNPart ) )
    allocate( magnr( 3*gNPart ) )
    allocate( magnrsq( 3*gNPart ) )
    
    tScale2 = 1.d0 / ( ( gT - gAvgTiming ) * gDTInv )
    
    
    do i = 1, gNMove
    
        diff = gPart( i ) % pos - gOldPos( i )
        
        if ( gForceMode .eq. 2 ) then
        	gVel( i ) =  diff * gDTInv
        else
        	gVel( i ) =  MinImage( diff ) * gDTInv
        end if
        
        gAV( i ) = gAV( i ) + gVel( i )
 
        
       	! Free rotor: Angular velocity calculation
       	
       	if ( gFreeRotor .eq. 1 ) then
        
        	radID1m = gPart( 1 ) % rat
		radID2m = gPart( 2 ) % rat
		
		radfract1 = 1.0 / ( 1.0 + ( radID2m / radID1m ) ** 3.0 )
		radfract2 = 1.0 / ( 1.0 + ( radID1m / radID2m ) ** 3.0 )
		
		dcenter = Vector( ( radfract1 * gPart( 1 ) % pos % x + radfract2 * gPart( 2 ) % pos % x ), &
             	( radfract1 * gPart( 1 ) % pos % y + radfract2 * gPart( 2 ) % pos % y ), &
             	( radfract1 * gPart( 1 ) % pos % z + radfract2 * gPart( 2 ) % pos % z ) )
             	
        
        	reldist( i ) = gPart( i ) % pos - dcenter
        
       		gVelx( i ) = gVel( i ) % x
        	gVely( i ) = gVel( i ) % y
        	gVelz( i ) = gVel( i ) % z
        
        	rdx( i ) = reldist( i ) % x
        	rdy( i ) = reldist( i ) % y
        	rdz( i ) = reldist( i ) % z
        
        	magnr( i ) = sqrt( rdx( i ) ** 2.0 + rdy( i ) ** 2.0 + rdz( i ) ** 2.0 )
        	magnrsq( i ) = magnr( i ) * magnr( i )
        
        	gAngVel( i ) = 1.0 / magnrsq( i ) * Vector( ( rdy( i ) * gVelz( i ) - rdz( i ) * gVely( i ) ), &
        				-1.0 * ( rdx( i ) * gVelz( i ) - rdz( i ) * gVelx( i ) ), &
        				( rdx( i ) * gVely( i ) - rdy( i ) * gVelx( i ) ) )
        
        	gAAngVel( i ) = gAAngVel( i ) + gAngVel( i )
        
        end if
        
        
        if ( mod( int( gT / gDT ), 100 ) .eq. 0 ) then
           
           if ( useconfig2 .eq. 1 ) then
           
              write( *, * ) ( gAV( i ) % x ) * tScale2
              
           end if
           
        end if
        
        
    end do
    
    if ( gSpecifiedFluxMode .eq. 1 ) then
        AKo( 1 ) = AKo( 1 ) + Xfinal
    end if

     
end subroutine PartVel


!subroutine PairDist
!use Globals
!implicit none


!end subroutine PairDist


!!****m* BrownianDynamics/RadDist
!!
!! NAME
!!  RadDist
!!
!! SYNOPSIS
!!  This subroutine measures the radial distribution between the
!!  externally forced particles by storing the distance between
!!  those particles in statistical bins.
!!
!! NOTES
!!  Located in file stat.f90
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
!subroutine RadDist
!use Globals
!implicit none
!
!    Integer :: i, j, k ! Loop counters.
!    type ( Vector ) :: dist, MinImage ! Displacements.
!    Real( 8 ) :: rad ! Distance between particles.
!
!    k = 1
!
!    do i = 1, gNMove
!        do j = i + 1, gNMove
!
!            dist = gPart( i ) % pos - gPart( j ) % pos
!            dist = MinImage( dist )
!            rad = sqrt( dist * dist )
!
!            gRD( k, int( gRDBins * rad / gMinDim ) ) = &
!                            gRD( k, int( gRDBins * rad / gMinDim ) ) + 1
!            k = k + 1
!        end do
!    end do
!
!end subroutine RadDist


!!****m* BrownianDynamics/DenProfile
!!
!! NAME
!!  DenProfile
!!
!! SYNOPSIS
!!  This subroutine stores the ID and location of the particles
!!  at different time steps for use in generating the average
!!  density profile of in the simulation cell.
!!
!! NOTES
!!  Located in file stat.f90
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
subroutine DenProfile
use Globals
implicit none

    type ( Vector ) :: disp, hold, Contain, MinImage
    Integer :: i, hX, hY, hZ
    Real( 8 ) :: factor1, factor2, r, theta, phi, phi2

    if ( gForceMode .eq. 1 ) then
        disp = gPart( 1 ) % pos
       
        if ( gFreeRotor .eq. 1 ) then
 
            do i = 2, gNMove
               factor1 = i - 1
               factor2 = 1.d0 / dble( i )

                hold = gPart( i ) % pos - disp
                hold = MinImage( hold )
                hold = hold + disp

                disp = factor1 * disp + hold 
                disp = disp * factor2
            end do

        end if
        
    else if ( gForceMode .eq. 2 ) then
        disp = gPart( 1 ) % pos

        do i = 2, gNMove
           factor1 = i - 1
           factor2 = 1.d0 / dble( i )

            hold = gPart( i ) % pos - disp
            hold = MinImage( hold )
            hold = hold + disp

            disp = factor1 * disp + hold 
            disp = disp * factor2
        end do
        
    end if


    do i = 1, gNMove
    
    
        hold = MinImage( gPart( i ) % pos - disp )

		if ( i .eq. 1 ) then 
        ! transformation to spherical coordinates
		r = sqrt( ( hold % x ) **2 + ( hold % y ) **2 + ( hold % z ) ** 2 )
		theta = atan( ( hold % y ) / ( hold % x ) )
		phi = acos( ( hold % z ) / r )
		
		end if

		
		! Locate two motors
		if ( i .eq. 1 ) then
		
			hold % x = 0
			hold % y = r
			
		end if
		
		if ( i .eq. gNMove ) then
		
			hold % x = 0
			hold % y = -r
			
		end if
		

                hold % x = hold % x + gDPXFixed
                hold % y = hold % y + gCDimH( 2 ) 
		hold % z = hold % z + gCDimH( 3 )
        
        hold = Contain( hold )


        hX = int( hold % x * gDPDim( 1 ) ) + 1
        hY = int( hold % y * gDPDim( 2 ) ) + 1
        hZ = int( hold % z * gDPDim( 3 ) ) + 1

        if ( hX .lt. 1 ) then
            hX = 1
        else if ( hX .gt. gDPBins( 1 ) ) then
            hX = gDPBins( 1 )
        end if

        if ( hY .lt. 1 ) then
            hY = 1
        else if ( hY .gt. gDPBins( 2 ) ) then
            hY = gDPBins( 2 )
        end if

        if ( hZ .lt. 1 ) then
            hZ = 1
        else if ( hZ .gt. gDPBins( 3 ) ) then
            hZ = gDPBins( 3 )
        end if

        gDP( hX, hY, hZ, i + 1 ) = gDP( hX, hY, hZ, i + 1 ) + 1.d0
    end do

    do i = gNMove + 1, gNPart
    	hold = gPart( i ) % pos - disp 
        !hold = MinImage( gPart( i ) % pos - disp )
        
        ! transformation to spherical coordinates
        
        if ( gFreeRotor .eq. 1 ) then
			r = sqrt( ( hold % x ) **2 + ( hold % y ) **2 + ( hold % z ) ** 2 )
			theta = atan( ( hold % y ) / ( hold % x ) )
			phi2 = acos( ( hold % z ) / r )
		


			hold % x = r * cos( theta ) * sin( phi2 - phi )
			hold % y = r * sin( theta ) * sin( phi2 - phi )
			hold % z = r * cos( phi2 - phi )
		
			hold = hold + disp

		end if
		
		!hold = hold + disp
		
        hold % x = hold % x + gDPXFixed
        hold % y = hold % y + gCDimH( 2 )
        hold % z = hold % z + gCDimH( 3 )


        hold = Contain( hold )

        hX = int( hold % x * gDPDim( 1 ) ) + 1
        hY = int( hold % y * gDPDim( 2 ) ) + 1
        hZ = int( hold % z * gDPDim( 3 ) ) + 1

        if ( hX .lt. 1 ) then
            hX = 1
        else if ( hX .gt. gDPBins( 1 ) ) then
            hX = gDPBins( 1 )
        end if

        if ( hY .lt. 1 ) then
            hY = 1
        else if ( hY .gt. gDPBins( 2 ) ) then
            hY = gDPBins( 2 )
        end if

        if ( hZ .lt. 1 ) then
            hZ = 1
        else if ( hZ .gt. gDPBins( 3 ) ) then
            hZ = gDPBins( 3 )
        end if

        gDP( hX, hY, hZ, 1 ) = gDP( hX, hY, hZ, 1 ) + 1.d0
        gPD3D( hX, hY, hZ ) = gPD3D( hX, hY, hZ ) + 1
    end do

end subroutine DenProfile

!!****m* BrownianDynamics/ZeroMeasures
!!
!! NAME
!!  ZeroMeasures
!!
!! SYNOPSIS
!!  This subroutine calls zeros out the statistical measures.
!!
!! NOTES
!!  Located in file stat.f90
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
subroutine ZeroMeasures
use Globals
implicit none

	Integer :: i
	
    gOP = 0
    
    !do i = 1, gNMove
        !gOPMotor( i ) = 0.0
    !end do

end subroutine ZeroMeasures


!!****m* BrownianDynamics/StatFinal
!!
!! NAME
!!  StatFinal
!!
!! SYNOPSIS
!!  This subroutine performs all the final statistical calculations.
!!
!! NOTES
!!  Located in file stat.f90
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
subroutine StatFinal
use Globals
implicit none

    ! Loop counters.
    Integer :: i, j, k, l, x, y, z, hX, hY, hZ
    ! Holder variables.
    Real( 8 ) :: tScale, avg, dist, tScale2
    ! The density profile.
    Integer, allocatable, dimension( :, :, : ) :: den
    ! Holder variables.
    Real( 8 ), dimension( gNMove + 1, 3 ) :: rad2Inv
    Integer, dimension( gNMove + 1, 3 ) :: rad

    ! Calculate the average meausures (velocity, osmotic pressure).
    tScale = 1.d0 / ( ( gTf - gAvgTiming ) * gDTInv )
    tScale2 = 1.d0 / ( ( gT - gAvgTiming ) * gDTInv )

    gAP = gAP * tScale 
    
    if ( gFreeRotor .eq. 1 ) then
    	gFCM = gFCM * tScale
    end if
    
    if ( gSpecifiedFluxMode .eq. 1 ) then
        	AKo( 1 ) = AKo( 1 ) * tScale
    end if
    
    do i = 1, gNMove
    
        gAV( i ) = gAV( i ) * tScale
        
        
        if ( gFreeRotor .eq. 1 ) then
        	gAAngVel( i ) = gAAngVel( i ) * tScale
        	gTensileF( i ) = gTensileF( i ) * tScale
        end if
        
        write( *, * ) 'Avg. x-velo:', gAV( i ) % x
        
        gHSF( i ) = gHSF( i ) * tScale
        gHST( i ) = gHST( i ) * tScale
        gHSFRh( i ) = gHSFRh( i ) * tScale
        gHSFnRh( i ) = gHSFnRh( i ) * tScale
        gOPMotor( i ) = gOPMotor( i ) * tScale

    end do
    
    do i = 1, gNMove * ( gNMove - 1 ) / 2

        avg = 0.d0

        do j = 1, gRDBins
            dist = gMinDim * j / gRDBins
            avg = avg + dist * dist * gRD( i, j )
        end do

        do j = 1, gRDBins
            dist = gMinDim * j / gRDBins
            gRD( i, j ) = dist * dist * gRD( i, j ) / avg
        end do

        avg = 0.d0

        !do j = 1, gPDBins( 1 )
            !do k = 1, gPDBins( 2 )
                !do l = 1, gPDBins( 3 )
                    !avg = avg + gPD( i, j, k, l )
                !end do
            !end do
        !end do

        !do j = 1, gPDBins( 1 )
            !do k = 1, gPDBins( 2 )
                !do l = 1, gPDBins( 3 )
                    !gPD( i, j, k, l ) = gPD( i, j, k, l ) / avg
                !end do
            !end do
        !end do

    end do

    !if ( gNMove .le. 1 ) then

        !avg = 0.d0

        !do j = 1, gPDBins( 1 )
            !do k = 1, gPDBins( 2 )
                !do l = 1, gPDBins( 3 )
                    !avg = avg + gPD( 1, j, k, l )
                !end do
            !end do
        !end do

        !do j = 1, gPDBins( 1 )
            !do k = 1, gPDBins( 2 )
                !do l = 1, gPDBins( 3 )
                    !gPD( 1, j, k, l ) = gPD( 1, j, k, l ) / avg
                !end do
            !end do
        !end do

    !end if

	!Normalize g at contact
	!do i = 1, gNMove
		!avg = 0
		!do j = 1, gThethaBins
			!avg = avg + gPDC( j, i )
		!end do
		
		!do j = 1, gThethaBins
			!gPDC( j, i ) = gPDC( j, i ) / avg
		!end do
	!end do
	

    ! Calculate the density profile (this could take a while).
    avg = 0.d0
    allocate( den( gDPBins( 1 ), gDPBins( 2 ), gDPBins( 3 ) ) )

    do j = 1, 3
        rad( 1, j ) = int( gDPBins( j ) / gCDim( j ) )
        rad2Inv( 1, j ) = ( gCDim( j ) / gDPBins( j ) ) ** 2
    end do

    do i = 2, gNMove + 1
        do j = 1, 3
            rad( i, j ) = int( gPart( i - 1 ) % rat * gDPBins( j ) / gCDim( j ) )
            rad2Inv( i, j ) = ( gCDim( j ) & 
                                  / ( gPart( i - 1 ) % rat * gDPBins( j ) ) ) ** 2
        end do
    end do

    do i = 1, gDPBins( 1 )
        do j = 1, gDPBins( 2 )
            do k = 1, gDPBins( 3 )
                do l = 1, gNMove + 1
                    if ( gDP( i, j, k, l ) .gt. 0 ) then
                        do x = 0, 2 * rad( l, 1 )
                            do y = 0, 2 * rad( l, 2 )
                                do z = 0, 2 * rad( l, 3 )
                                    dist = rad2Inv( l, 1 ) * ( x & 
                                           - rad( l, 1 ) ) ** 2 &
                                       + rad2Inv( l, 2 ) * ( y &
                                           - rad( l, 2 ) ) ** 2 &
                                       + rad2Inv( l, 3 ) * ( z &
                                           - rad( l, 3 ) ) ** 2

                                    if ( dist .le. 1.0 ) then
                                        hX = i + x - rad( l, 1 )
                                        hY = j + y - rad( l, 2 )
                                        hZ = k + z - rad( l, 3 )

                                        if ( hX .gt. gDPBins( 1 ) ) then
                                            hX = hX - gDPBins( 1 )
                                        else if ( hX .lt. 1 ) then
                                            hX = hX + gDPBins( 1 )
                                        end if

                                        if ( hY .gt. gDPBins( 2 ) ) then
                                            hY = hY - gDPBins( 2 )
                                        else if ( hY .lt. 1 ) then
                                            hY = hY + gDPBins( 2 )
                                        end if

                                        if ( hZ .gt. gDPBins( 3 ) ) then
                                            hZ = hZ - gDPBins( 3 )
                                        else if ( hZ .lt. 1 ) then
                                            hZ = hZ + gDPBins( 3 )
                                        end if

                                        den( hX, hY, hZ ) = den( hX, hY, hZ ) &
                                                        + gDP( i, j, k, l )
                                        avg = avg + gDP( i, j, k, l )
                                    end if 
                                end do
                            end do
                        end do
                    end if
                end do
            end do
        end do
    end do

    do i = 1, gDPBins( 1 )
        do j = 1, gDPBins( 2 )
            do k = 1, gDPBins( 3 )
                gDP( i, j, k, 1 ) = den( i, j, k ) / avg
            end do
        end do
    end do

	

  end subroutine StatFinal
  
  
  SUBROUTINE OrderParameter
    USE GLOBALS
    USE VECTORCLASS
    IMPLICIT NONE
    
    INTEGER :: I
    REAL(8) :: DOTPROD, SUMROT
    REAL(8) :: SUMTRANS
    SUMROT   = 0.D0
    SUMTRANS = 0.D0
    
    DO I = 1, gNPart
       DOTPROD = gInitOrient(i) * orient(i)
       SUMROT = SUMROT + DOTPROD
       SUMTRANS = SUMTRANS + DCOS( gReciprocalVec * gPart(i) % pos )
    END DO
    
    gRotOrderPar = SUMROT / DFLOAT(gNPart)
    gTransOrderPar = SUMTRANS / DFLOAT(gNPart)
    WRITE(1800,101) gT, gRotOrderPar
    WRITE(1900,101) gT, gTransOrderPar
    101 format(2(1X, E24.15E3))
  END SUBROUTINE OrderParameter
  
