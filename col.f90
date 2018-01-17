!!***m* BrownianDynamics/CollisionDetection
!!
!! NAME
!!  CollisionDetection
!!
!! SYNOPSIS
!!  This SUBROUTINE sets up the collision detection algorithm by building
!!  an outer loop that ensures all collisions are resolved and generating
!!  the linked cells for making the collision detection faster.  The
!!  collision detection itself is handled in a separate routine.  This
!!  routine is only CALLed from the main program loop.
!!
!! NOTES
!!  Located in file col.f90
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
!!******************************************************************************
!SUBROUTINE CollisionDetection
!USE Globals
!IMPLICIT NONE
!
!    type ( IntLL ), target :: hold ! A holder for the linked lists.
!    type ( IntLL ), pointer :: ptr ! A pointer to the linked lists.
!    ! The IDs of particles in two neighboring linked cells.
!    Integer, allocatable, dimension( : ) :: sVals, tVals
!    ! Loop counters.
!    Integer :: i, j, k, l, m, n, o, ml, mm, mn
!    
!    
!    ! Vectors USEd for measuring the displacement between particles.
!    type ( Vector ) :: diff, disp1, disp2, motDisp, motDisp2, MinImage, ranVec
!    type ( Vector ) :: w, dcenter
!
!    ! Scale factors USEd for moving particles the correct amount.
!    Real( 8 ) :: radID1, radID2, radSum2, diff2, radfract1, radfract2
!
!    ! Variables USEd for free rotor mode
!    Real( 8 ) :: dvecs2, radID1m, radID2m, radSum2m, Rx, Ry, Rz
!    Real( 8 ) :: dvecsx, dvecsy, dvecsz, dnewx, dnewy, dnewz, orientnorm
!
!    type ( Vector ) :: dispm1, dispm2, dnew, dvecs
!
!
!    gCollisionFlag = .true. ! Loop through the routine ones
!
!    ! Free rotor mode is defined like two linked motors
!
!!!!!*********************************************************************
!    if ( gFreeRotor .eq. 1 ) then
!                ! Calculate relative distance between motors
!
!                dvecs = gPart( 1 ) % pos - gPart( 2 ) % pos
!                dvecs = MinImage( dvecs )
!
!                dvecs2 = dvecs * dvecs
!
!                ! write( *, *) 'relative distance:', dsqrt( dvecs2 )
!
!
!                ! Checking if relative distance is not what we want
!                if ( abs( dfixed - dvecs2 ) .gt. gMinSepMotors ) then
!
!                        radID1m = gPart( 1 ) % rat
!                        radID2m = gPart( 2 ) % rat
!
!
!                        radSum2m = ( radID1m + radID2m ) * ( radID1m + radID2m )
!
!                        ! Include these variables
!                        radfract1 = 1.0 / ( 1.0 + ( radID2m / radID1m ) ** 3.0 )
!                        radfract2 = 1.0 / ( 1.0 + ( radID1m / radID2m ) ** 3.0 )
!
!
!                        Rx = orient( 1 ) % x  !
!                        Ry = orient( 1 ) % y  ! Look!!!
!                        Rz = orient( 1 ) % z  !
!
!
!                        dvecsx = dvecs % x
!                        dvecsy = dvecs % y
!                        dvecsz = dvecs % z
!
!                        dvecs = dvecs * ( sqrt( dfixed / dvecs2 ) - 1. )
!
!                        dispm1 = dvecs * ( radID2m / ( radID1m + radID2m ) )
!                        dispm2 = dvecs * ( radID1m / ( radID1m + radID2m ) )
!
!                        gPart( 1 ) % pos = gPart( 1 ) % pos + dispm1
!                        gPart( 2 ) % pos = gPart( 2 ) % pos - dispm2
!
!                       ! Free rotor mode - Fixed center of mass
!                          if ( gFixedCenterMass .eq. 1 ) then
!             
!                        dcenter = Vector( ( radfract1 * gPart( 1 ) % pos % x + radfract2 * gPart( 2 ) % pos % x ), &
!                                ( radfract1 * gPart( 1 ) % pos % y + radfract2 * gPart( 2 ) % pos % y ), &
!                                ( radfract1 * gPart( 1 ) % pos % z + radfract2 * gPart( 2 ) % pos % z ) )
!
!                       !dcenter = MinImage( dcenter )
!
!                                dcenter = dcenter - dcmfixed
!
!                                gPart( 1 ) % pos = gPart( 1 ) % pos - dcenter
!                                gPart( 2 ) % pos = gPart( 2 ) % pos - dcenter
!
!                               ! gFCM is the force required to bring center of mass back to its initial and imposed position
!                               gFCM = gFCM + dcenter
!
!                           END if
!
!                       ! gTF calculates the force (tensile force) to bring the two motors back to initial and imposed relative distance
!                        gTensileF( 1 ) = gTensileF( 1 ) + dispm1
!                        gTensileF( 2 ) = gTensileF( 2 ) - dispm2
!
!                        dnew = MinImage( gPart( 1 ) % pos - gPart( 2 ) % pos )
!
!                        dnewx = dnew % x
!                        dnewy = dnew % y
!                        dnewz = dnew % z
!
!
!                        ! Calculate new direction of reaction vector: Rnew = ( dvecs x R ) x dnew
!                        ! R is perpENDicular to dvecs
!
!                        orientnorm = sqrt( ( dvecsy * dnewy * Rx + dvecsz * dnewz * Rx - dvecsx * dnewy * Ry - dvecsx * dnewz * Rz) ** 2 &
!                                                + ( -1.d0 * dnewx * dvecsy * Rx + dvecsx * dnewx * Ry + dvecsz * dnewz * Ry - dvecsy * dnewz * Rz ) ** 2 &
!                                                + ( -1.d0 * dnewx * dvecsz * Rx - dnewy * dvecsz * Ry + dvecsx * dnewx * Rz + dvecsy * dnewy * Rz ) ** 2 )
!
!                        orient( 1 ) = 1.0 / orientnorm * Vector( dvecsy * dnewy * Rx + dvecsz * dnewz * Rx - dvecsx * dnewy * Ry - dvecsx * dnewz * Rz, &
!                                                 -1.d0 * dnewx * dvecsy * Rx + dvecsx * dnewx * Ry + dvecsz * dnewz * Ry - dvecsy * dnewz * Rz, &
!                                                 -1.d0 * dnewx * dvecsz * Rx - dnewy * dvecsz * Ry + dvecsx * dnewx * Rz + dvecsy * dnewy * Rz )
!
!                        orient( 2 ) = 1.0 / orientnorm * Vector( -1.d0 * dvecsy * dnewy * Rx - dvecsz * dnewz * Rx + dvecsx * dnewy * Ry + dvecsx * dnewz * Rz, &
!						dnewx * dvecsy * Rx - dvecsx * dnewx * Ry - dvecsz * dnewz * Ry + dvecsy * dnewz * Rz, &
!						dnewx * dvecsz * Rx + dnewy * dvecsz * Ry - dvecsx * dnewx * Rz - dvecsy * dnewy * Rz )
!
!
!                END if
!       END if
!!!!!!!!!*****************************************************************************************************************************************************************    
!       
!        if ( gNoBathColFlag .eq. 1 ) then        !!Like ideal gas, gNoBathColFlag = 1
!
!                allocate( sVals( gNMove ) )              !!Number of motor particles
!                allocate( tVals( gNPart - gNMove ) )     !!Number of bath particles
!        
!                do while ( gCollisionFlag )
!
!                   gCollisionFlag = .false. ! Make sure we don't get stuck
!
!                   do i = 1, gNMove
!
!                      sVals( i ) = i
!
!                   END do
!
!                   do i = 1, gNPart - gNMove
!
!                      tVals( i ) = i + gNMove
!
!                   END do
!
!                    CALL CheckCollision( sVals, tVals, gNMove, gNPart - gNMove )
!
!                END do
! 
!        else
!	
!    	do while ( gCollisionFlag )
!               CALL FillLCells ! Set up the linked cells
!
!        	gCollisionFlag = .false. ! Make sure we don't get stuck
!
!        	! Scan through all the cells and check neighboring ones for collisions
!        	do i = 1, gNLCells( 1 )
!            	do j = 1, gNLCells( 2 )
!                	do k = 1, gNLCells( 3 )
!
!                    	! These are the particles in one cell.
!                    	allocate( sVals( gLCell( i, j, k ) % sz ) )
!
!                    	do l = 1, gLCell( i, j, k ) % sz
!                        	sVals( l ) = gLCell( i, j, k ) % part( l ) 
!                    	END do
!
!                    	! Choose the other neighboring cells.
!                    	do l = gLLS( 1 ), 1
!                        	do m = gLLS( 2 ), 1
!                            	do n = gLLS( 3 ), 1
!                                	ml = i + l
!                                	mm = j + m
!                                	mn = k + n
!
!                                	if ( ml .eq. 0 ) then
!                                    	ml = gNLCells( 1 )
!                                	else if ( ml .gt. gNLCells( 1 ) ) then
!                                    	ml = 1
!                                	END if
!
!                                	if ( mm .eq. 0 ) then
!                                    	mm = gNLCells( 2 )
!                                	else if ( mm .gt. gNLCells( 2 ) ) then
!                                    	mm = 1
!                                	END if
!
!                                	if ( mn .eq. 0 ) then
!                                    	mn = gNLCells( 3 )
!                                	else if ( mn .gt. gNLCells( 3 ) ) then
!                                    	mn = 1
!                                	END if
!
!
!                                	! This is the other neightboring cell.
!                                	allocate( tVals( gLCell( ml, mm, mn ) % sz ) )
!
!                                	do o = 1, gLCell( ml, mm, mn ) % sz
!                                    	tVals( o ) = gLCell( ml, mm, mn ) % part( o )
!                                	END do
!         
!                                	! Check for collision between particles in each cell.
!                                	CALL CheckCollision( sVals, tVals, &
!                                 	gLCell( i, j, k ) % sz, gLCell( ml, mm, mn ) % sz )
!
!                                	deallocate( tVals )
!
!                            	END do
!                        	END do
!                    	END do 
!
!                    	deallocate( sVals )
!
!                	END do
!            	END do
!        	END do
!    	END do
!
!	END if
!	
!    if ( gOsmMotorCol ) then
!        gOsmMotorCount = 1
!        gOsmMotorCol = .false.
!    else
!        gOsmMotorCount = gOsmMotorCount + 1
!    END if
!
!
!END SUBROUTINE CollisionDetection
!!******************************************************************************************************


!!****m* BrownianDynamics/CheckCollision
!!
!! NAME
!!  CheckCollision
!!
!! SYNOPSIS
!!  This SUBROUTINE handles all of the particle-particle collision detection.  
!!  The parameters passed to it are simply arrays of particle IDs for each of 
!!  the neighboring cells.  Collision is checked by computing the center to center 
!!  distance between particles and seeing if its is less than the sum of the 
!!  particles' radii.  Then the particles are moved into contact along the line
!!  connecting their centers and the osmotic pressure is calculated based upon this
!!  distance.
!!
!! USAGE
!!  SUBROUTINE CheckCollision( sVals, tVals, sSz, tSz )
!!
!! INPUTS
!!  sVals, tVals - arrays containing the IDs of particles to be checked for
!!  for collision.
!!
!!  sSz, tSz - the size of arrays sVals and tVals respectively.
!!
!! NOTES
!!  Located in file col.f90
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
!SUBROUTINE CheckCollision( sVals, tVals, sSz, tSz )
!USE Globals
!IMPLICIT NONE
!
!    ! The size of the ID arrays passed to the SUBROUTINE.
!    Integer :: sSz, tSz
!
!    ! Loop counters.
!    Integer :: i, j, id1, id2, diffid,  id3, id4
!
!    ! An array of particle IDs to check for collision.
!    Integer, dimension( sSz ), intent( in ) :: sVals
!
!    ! An array of particle IDS to check for collision.
!    Integer, dimension( tSz ), intent( in ) :: tVals
!
!    ! Vectors USEd for measuring the displacement between particles.
!    type ( Vector ) :: diff, disp1, disp2, motDisp, motDisp2, MinImage, ranVec
!    type ( Vector ) :: w, dcenter
!
!    ! Scale factors USEd for moving particles the correct amount.
!    Real( 8 ) :: radID1, radID2, radSum2, diff2
!
!    ! Variables USEd for free rotor mode
!    Real( 8 ) :: dvecs2, radID1m, radID2m, radSum2m, Rx, Ry, Rz
!    Real( 8 ) :: dvecsx, dvecsy, dvecsz, dnewx, dnewy, dnewz, orientnorm
!    type ( Vector ) :: dispm1, dispm2, dnew, dvecs
!    Real, dimension ( 4 ) :: random
!
!    ! Activate only when using uniformly random distribution
!    !Real, allocatable, dimension ( : ) :: random
!
!    ! Difference between particle 1 and bath particles, reaction probability
!    Real :: diffx, ran0, diffz, ratiof, Daeff, Df, diffy
!    Real( 8 ) :: xx, yy, zz, rr, theta, phiangle, mu_th, hthetha 
!    Real( 8 ) :: magdiff, magR, costhetarxn
!    Real( 8 ) :: tScale2, ratfact1, ratfact2, numerator, denominator, stericfactor
!
!
!    do i = 1, sSz
!        do j = 1, tSz
!
!            id1 = sVals( i )    !For motor particles
!            id2 = tVals( j )    !For bath particles
!
!            if ( ( gNoBathColFlag .eq. 0 ) .or. & 
!               ( ( id1 .le. gNMove ) .or. ( id2 .le. gNMove ) ) ) then
!            if ( id1 .ne. id2 ) then
!                diff = gPart( id1 ) % pos - gPart( id2 ) % pos
!                
!                ! Set up the proper scale for collision.
!                if ( ( id1 .le. gNMove ) .or. ( id2 .le. gNMove ) ) then
!                    radID1 = gPart( id1 ) % rat ! 
!                    radID2 = gPart( id2 ) % rat ! 
!                    radSum2 = ( radID1 + radID2 ) * ( radID1 + radID2 )
!                else
!                    radID1 = 1.
!                    radID2 = 1.
!                    radSum2 = 4.
!                END if
!
!                ! Find the minimum image of the collision.
!                diff = MinImage( diff )
!
!                diff2 = diff * diff
!
!                if ( ( radSum2 - diff2 ) .gt. gMinSep ) then
! 
!                    diff = diff * ( sqrt( radSum2 / diff2 ) - 1. )
!                    
!                    ! Make sure the particles move the appropriate amount.
!                    if ( ( id1 .le. gNMove ) .or. ( id2 .le. gNMove ) ) then
!                        disp1 = diff * ( radID2 / ( radID1 + radID2 ) )
!                        disp2 = diff * ( radID1 / ( radID1 + radID2 ) )
!
!                        ! Make sure the particles move as specified
!                        if ( gForceMode .gt. 1 ) then  ! In our case gForceMode = 1, therefore this sentence do not applied
!                            if ( ( id1 .le. gNMove ) .and. ( id2 .le. gNMove ) ) then
!                                disp1 = Vector( 0, 0, 0 )
!                                disp2 = Vector( 0, 0, 0 )
!                            else if ( id1 .le. gNMove ) then
!                                disp1 = Vector( 0, 0, 0 )
!                                disp2 = diff
!                                w = diff * ( 1.0 / sqrt( diff2 ) ) * radID1
!                                w = w + dvec( id1 )
!                                sgn( id1 ) = ( 2.0 * mod( id1, 2 ) - 1.0 )
!                                gHSF( id1 ) = gHSF( id1 ) + diff
!                                gHST( id1 ) % x = gHST( id1 ) % x + ( w % y ) * ( diff % z ) - ( w % z ) * ( diff % y )
!                                gHST( id1 ) % y = gHST( id1 ) % y - ( ( w % x ) * ( diff % z ) - ( w % z ) * ( diff % x ) )
!                                gHST( id1 ) % z = gHST( id1 ) % z + ( w % x ) * sgn( id1 ) * ( diff % y ) - ( w % y ) * ( diff % x )
!                            else if ( id2 .le. gNMove ) then
!                                disp1 = diff
!                                disp2 = Vector( 0, 0, 0 )
!                                w = diff * ( 1.0 / sqrt( diff2 ) ) * radID2
!                                w = w + dvec( id2 )
!                                sgn( id2 ) = ( 2.0 * mod( id2, 2 ) - 1.0 )
!                                gHSF( id2 ) = gHSF( id2 ) - diff
!                                gHST( id2 ) % x = gHST( id2 ) % x + ( w % y ) *  ( -1.0 * diff % z ) - ( w % z ) * ( -1.0 * diff % y )
!                                gHST( id2 ) % y = gHST( id2 ) % y - ( ( w % x ) * ( -1.0 * diff % z ) - ( w % z ) * ( -1.0 * diff % x ) )
!                                gHST( id2 ) % z = gHST( id2 ) % z + ( w % x ) * ( -1.0 * diff % y ) - ( w % y ) * ( -1.0 * diff % x )
!
!                            END if
!                        END if
!                    else
!                        disp1 = diff * 5D-1 
!                        disp2 = diff * 5D-1
!                    END if
!
!                    ! Move the particles.
!                    gPart( id1 ) % pos = gPart( id1 ) % pos + disp1
!                    gPart( id2 ) % pos = gPart( id2 ) % pos - disp2
!
!                    
!                    !gPart( id1 ) % abspos = gPart( id1 ) % abspos + disp1
!                    !gPart( id2 ) % abspos = gPart( id2 ) % abspos - disp2
!
!
!                    gHSF( id1 ) = gHSF( id1 ) + disp1 
!                    gHSF( id2 ) = gHSF( id2 ) - disp2
!
!
!					! Measure the force on each side of the motor: reactive and non reactive
!					! Only works properly for fixed motor problem
!					!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!					if ( gForceMode .gt. 1 ) then    ! In our case gForceMode = 1
!					
!					    if ( ( id1 .le. gNMove ) .or. ( id2 .le. gNMove ) ) then
!					
!						    if ( IEOR( (id1 .le. gNMove ), (id2 .le. gNMove ) ) ) then
!						
!  							    if ( id1 .le. gNMove ) then
!       							    diffid = id1
!   							    else
!								    diffid = id2
!							    END if
!							
!   						    END if
!
!						    if ( ( diff * orient(diffid) ) * abs( id2 - id1 ) / ( id2 - id1 ) < 0 ) then
!						
!							    if ( id1 .le. gNMove ) then
!								    gHSFRh( id1 ) = gHSFRh( id1 ) + disp2
!							    else if ( id2 .le. gNMove ) then
!								    gHSFRh( id2 ) = gHSFRh( id2 ) - disp1
!							    END if
!							
!						    else
!
!   							    if ( id1 .le. gNMove ) then
!   								    gHSFnRh( id1 ) = gHSFnRh( id1 ) + disp2
!							    else if ( id2 .le. gNMove ) then
!								    gHSFnRh( id2 ) = gHSFnRh( id2 ) - disp1
!							    END if
!							
!						    END if
!
!					    END if
!					
!					END if
!                                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    ! Measure the osmotice pressure. Q: Shoud this be after rotor problem code?
!                    diff = MinImage( gPart( id1 ) % pos - gPart( id2 ) % pos )
!                    gOP = gOP + diff * disp1
!
!					! Measure the stress of the motor
!					if ( ( id1 .le. gNMove ) .or. ( id2 .le. gNMove ) ) then
!						if ( id1 .le. gNMove ) then
!							gOPMotor( id1 ) = gOPMotor( id1 ) + diff * disp2
!						else if ( id2 .le. gNMove ) then
!							gOPMotor( id2 ) = gOPMotor( id2 ) + diff * disp1
!						END if
!
!					END if
!					
!                    gCollisionFlag = .true.
!                    
!
!                    ! Surface reaction emulator at half sphere
!                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                    if ( IEOR( (id1 .le. gNMove ), (id2 .le. gNMove ) ) ) then
!                       if ( id1 .le. gNMove ) then
!                           diffid  =  id1
!                           motDisp = disp1
!                           motDisp2 = disp2
!                       else
!                           diffid = id2
!                           motDisp = Vector( -1.d0 * disp2 % x, -1.d0 * disp2 % y, &
!                                   -1.d0 * disp2 % z )
!                           motDisp2 = Vector( -1.d0 * disp1 % x, -1.d0 * disp1 % y, &
!                                    -1.d0 * disp1 % z ) 
!                       END if
!
!
!                       diffx = diff % x
!                       diffy = diff % y
!                       diffz = diff % z
!                       magdiff = ( diffx ** 2.0 + diffy ** 2.0 + diffz ** 2.0 ) ** ( 1.0 / 2.0 )
!                       Rx = orient( diffid ) % x
!		       Ry = orient( diffid ) % y
!		       Rz = orient( diffid ) % z
!                       magR = ( Rx ** 2.0 + Ry ** 2.0 + Rz ** 2.0 ) ** ( 1.0 / 2.0 )RANDOMPOSITIONING
!                       
!                       ! Collision with the reactive side
!                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
!                       if ( ( diff * orient(diffid) * ( 1.0 / magdiff ) * ( 1.0 / magR ) ) * &
!                          abs( id2 - id1 ) / ( id2 - id1 ) < -1.0 * gCosThetaRxn( diffid ) ) then
!                        
!                          gCountCol = gCountCol + 1.
!                          
!                          !Following computations are for Gaussian random distribution
!                          random( 1 ) = ran0( seed2 )
!                          random( 2 ) = ran0( seed2 )
!                          random( 3 ) = ran0( seed2 )
!                          random( 4 ) = ran0( seed2 )
!
!                                   
!                           ! Deleted part (see col.f90 original file)
!
!                          
!                          gOsmMotorCount = 1.d0
!
!                           
!
!                          gOsmMotorCol = .true.
!                          coltime( id2 ) = 0
!  
!                          if (( rxnp(diffid) .gt. random(1)) )  then
!                              if ( id1 .le. gNMove ) then
!                                  id3 = id2
!                                  id4 = id1
!                              else
!                                  id3 = id1
!                                  id4 = id2
!                              END if
!
!                              
!                              gCountColReact = gCountColReact + 1.
!                              
!                              if ( gRandLoc .eq. 1 ) then    ! gRandLoc : Reaction probability/ Mean Specific flux and relocation mode
!                                                             !            Relocate particles: 1 -randomly, 0 -other side of motor
!                                                             !            In our case, gRandLoc = 1
!                                  random( : ) = 2.d0 * random( : ) - 1.d0
!                                  ranVec = Vector( random( 2 ) * gCDimH( 1 ), random( 3 ) * gCDimH( 2 ), random( 4 ) * gCDimH( 3 ) )
!                                  
!                                  if ( g2Dproblem .eq. 1 ) then
!                                  	ranVec = Vector( random( 2 ) * gCDimH( 1 ), random( 3 ) * gCDimH( 2 ), 0.0 )
!                                  END if
!                                  
!                                  ! xb = xm + n*(a+b) ranVec/sqrt(ranVec^2) == at contact if n =1. Play with n 
!                                  gPart( id3 ) % pos = gPart( id4 ) % pos + ranVec * ( 1.d0 + 2.d0 / sqrt( ranVec * ranVec ) ) 
!
!                                  
!                                  !Following computations are for uniformly random distribution
!                                  !gPart( id3 ) % pos % x = random(2) * gCDim(1)
!                                  !gPart( id3 ) % pos % y = random(3) * gCDim(2)
!                                  !gPart( id3 ) % pos % z = random(4) * gCDim(3)
!                                  
!                              else 
!                                  rr = 1.d0 + gPart( id4 ) % rat
!                                  theta = random(2) * gPi
!                                  phiangle = random(3) * gPi + gPi/2.
!                                  xx = rr * sin(theta) * cos(phiangle)
!                                  yy = rr * sin(theta) * sin(phiangle)
!                                  zz = rr * cos(theta)
!                                  gPart( id3 ) % pos = gPart( id4 ) % pos + Vector(xx, yy, zz)
!                             END if
!
!                              ! Move particle back to where it was (no HS step when bath particle reacts)
!                              if ( hsMode .eq. 1 ) then
!                                   gPart( diffid ) % pos = gPart( diffid ) % pos - motDisp
!                                   gHSF( diffid ) = gHSF( diffid ) - motDisp2
!                              END if
!                              
!                          END if
!                       
!                       ! Collision with no rective side
!                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!                       else
!                          
!                          gCountColNoReactive = gCountColNoReactive + 1
!                          
!                       END if
!                    
!                    END if
!                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    ! Make sure the particles are inside the simulation cell.
!                    CALL ContainParticle( id1 )
!                    CALL ContainParticle( id2 )
!                    
!                    !Pair-distribution function at contact
!                    !do i = 1, gNMove
!                    	!diffx = diff % x
!                    	!mu_th = diffx / sqrt( diff * diff )
!                    	!hthetha = int( ( mu_th + 1.d0 ) / 2.d0 * ( gThetaBins - 1 ) ) + 1
!                   	!gPDC( hthetha, i ) = gPDC( hthetha, i ) + 1
!                    !END do
!
!                END if
!            END if
!            END if
!        END do
!    END do
!
!END SUBROUTINE CheckCollision


!!****f* BrownianDynamics/MinImage
!!
!! NAME
!!  CheckCollision
!!
!! SYNOPSIS
!!  This function finds the minimum image of a vector measuring the 
!!  difference between particle positions and returns it.
!!
!! USAGE
!!  type ( Vector ) function MinImage( diff )
!!
!! INPUTS
!!  diff - any vector type object.
!!
!! RESULT
!!  A vector with size and orientation corresponding to the minimum image
!!  of the input vector.
!!
!! NOTES
!!  Located in file col.f90
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
!type ( Vector ) function MinImage( diff )
!USE Globals
!IMPLICIT NONE

! The vector for whom the minimum image is desired
!    type ( Vector ), intent( in ) :: diff
!
!    if ( diff % x .gt. gCDimH( 1 ) ) then
!        MinImage % x = diff % x - gCDim( 1 )
!!    else if ( diff % x .lt. -1 * gCDimH( 1 ) ) then
!        MinImage % x = diff % x + gCDim( 1 )
!    else
!        MinImage % x = diff % x
!    END if
!
!    if ( diff % y .gt. gCDimH( 2 ) ) then
!        MinImage % y = diff % y - gCDim( 2 )
!    else if ( diff % y .lt. -1 * gCDimH( 2 ) ) then
!        MinImage % y = diff % y + gCDim( 2 )
!    else
!        MinImage % y = diff % y
!    END if
!
!    if ( diff % z .gt. gCDimH( 3 ) ) then
!        MinImage % z = diff % z - gCDim( 3 )
!    else if ( diff % z .lt. -1 * gCDimH( 3 ) ) then
!        MinImage % z = diff % z + gCDim( 3 )
!    else
!        MinImage % z = diff % z
!    END if
!
!END function MinImage

!REAL(8), DIMENSION(gND) function MinImage( diff )
!  USE Globals
!  IMPLICIT NONE
!  INTEGER :: i  
!  REAL(8), DIMENSION(gND), intent( in ) :: diff
    
!  DO i = 1, gND
!     IF ( diff(i) .GT. gCDimH(i) ) THEN
!        MinImage(i) = diff(i) - gCDim( i )
!     ELSE IF  ( diff(i) .LT. -1.0 * gCDimH(i) ) THEN
!        MinImage(i) = diff(i) + gCDim( i )
!     ELSE
!        MinImage(i) = diff(i) 
!  END DO

!END function MinImage



!REAL(8) FUNCTION vdot(v1,v2)
!  USE Globals
!  IMPLICIT NONE
!  REAL(8), DIMENSION(gND) :: v1
!  REAL(8), DIMENSION(gND) :: v2
!  
!  vdot = SUM( v1(:) * v2(:) )  
!  
!END FUNCTION vdot


!!****m* BrownianDynamics/InitLCells
!!
!! NAME
!!  InitLCells
!!
!! SYNOPSIS
!!  This SUBROUTINE initializes the linked lists and cells by determining 
!!  how many linked cells are needed and allocating the memory for the lists.  
!!  The number of cells is determined by the largest distance over which two 
!!  particles can collide.  This is USEd exclusively for collision detection.
!!
!! NOTES
!!  Located in file col.f90
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
SUBROUTINE InitLCells
USE Globals

    gNLCells( 1 ) = int( 0.95 * gCDim( 1 ) / ( gLargestRad + gSecondLargestRad ) )
    gNLCells( 2 ) = int( 0.95 * gCDim( 2 ) / ( gLargestRad + gSecondLargestRad ) )
    gNLCells( 3 ) = int( 0.95 * gCDim( 3 ) / ( gLargestRad + gSecondLargestRad ) )

    gLCells = gNLCells( 1 ) * gNLCells( 2 ) * gNLCells( 3 )
    
    if ( gNLCells( 1 ) .eq. 1 ) then
        gLLS( 1 ) = 1
    else if ( gNLCells( 1 ) .eq. 2 ) then
        gLLS( 1 ) = 0
    else 
        gLLS( 1 ) = -1
    END if

    if ( gNLCells( 2 ) .eq. 1 ) then
        gLLS( 2 ) = 1
    else if ( gNLCells( 2 ) .eq. 2 ) then
        gLLS( 2 ) = 0
    else 
        gLLS( 2 ) = -1
    END if

    if ( gNLCells( 3 ) .eq. 1 ) then
        gLLS( 3 ) = 1
    else if ( gNLCells( 3 ) .eq. 2 ) then
        gLLS( 3 ) = 0
    else 
        gLLS( 3 ) = -1
    END if

    allocate( gLCell( gNLCells( 1 ), gNLCells( 2 ), gNLCells( 3 ) ) )
 
END SUBROUTINE InitLCells


!!****m* BrownianDynamics/FillLCells
!!
!! NAME
!!  FillLCells
!!
!! SYNOPSIS
!!  This SUBROUTINE fills the linked cells by first filling the linked lists and 
!!  then allocating memory for the arrays storing the IDs of the particles in each
!!  cell.  This is USEd exclusively for collision detection.
!!
!! NOTES
!!  Located in file col.f90
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
SUBROUTINE FillLCells
USE Globals

    ! Counters, array sizes and bin positions.
    Integer :: i, x, y, z, sz
    ! A pointer to linked list data.
    type ( IntLL ), pointer :: ptr

    CALL ClearLCells

    !print*, "BEFORE FIRST DO LOOP"
    do i = 1, gNPart
        x = int( gNLCells( 1 ) * gPart( i ) % pos % x / gCDim( 1 ) ) + 1
        y = int( gNLCells( 2 ) * gPart( i ) % pos % y / gCDim( 2 ) ) + 1
        z = int( gNLCells( 3 ) * gPart( i ) % pos % z / gCDim( 3 ) ) + 1

        sz = gLCell( x, y, z ) % sz

        gLCell( x, y, z ) % sz = sz + 1

        ! Fill the linked lists.
        if ( sz .eq. 0 ) then
            gLCell( x, y, z ) % ll % val = i
        else   
            ptr => gLCell( x, y, z ) % ll

            CALL IntAddValue( ptr, i )
        END if

    END do
  
    !print*, "AFTER FIRST DO LOOP"

    ! Empty the linked lists and fill the arrays.
    do x = 1, gNLCells( 1 )
        do y = 1, gNLCells( 2 )
            do z = 1, gNLCells( 3 )
                ptr => gLCell( x, y, z ) % ll

                allocate( gLCell( x, y, z ) % part( gLCell( x, y, z ) % sz ) )

                do i = 1, gLCell( x, y, z ) % sz
                    gLCell( x, y, z ) % part( i ) = IntGetValue( ptr )
                END do

                deallocate( ptr )
            END do
        END do
    END do

END SUBROUTINE FillLCells


!!****m* BrownianDynamics/ClearLCells
!!
!! NAME
!!  ClearLCells
!!
!! SYNOPSIS
!!  This SUBROUTINE clears the linked cells and frees up memory USEd by the
!!  linked lists and arbitrary sized arrays.  This is USEd exclusively for 
!!  collision detection.
!!
!! NOTES
!!  Located in file col.f90
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
SUBROUTINE ClearLCells
USE Globals

    ! Loop counters.
    Integer :: i, j, k

    do i = 1, gNLCells( 1 )
        do j = 1, gNLCells( 2 )
            do k = 1, gNLCells( 3 )
                gLCell( i, j, k ) % sz = 0

                ! Free up old memory.
                if ( allocated( gLCell( i, j, k ) % part ) ) then
                    deallocate( gLCell( i, j, k ) % part )
                END if

                ! Create new linked lists.
                nullify( gLCell( i, j, k ) % ll )
                allocate( gLCell( i, j, k ) % ll )
                nullify( gLCell( i, j, k ) % ll % next )
            END do
        END do
    END do

END SUBROUTINE ClearLCells






