!!************************************************************************
!!****m* BrownianDynamics/Initialize
!!
!! NAME
!!  Initialize
!!
!! SYNOPSIS
!!  This SUBROUTINE CALLs all of the initialization routines.
!!
!! NOTES
!!  Located in file init.f90
!!
!! CREATION DATE
!!    02/28/05
!!
!! AUTHOR 
!!    James Swan
!!
!! COPYRIGHT
!!    28-Feb-2005 - James William Swan
!!************************************************************************
!!************************************************************************
!!************************************************************************
SUBROUTINE Initialize
  USE Globals
  USE VectorFunctions
  IMPLICIT NONE
  
  open( unit = 101, file = 'conf.in', status = 'old' )
  
  CALL StartDataRead 
  
  CALL AllocateArrays
  
  CALL FinishDataRead
  
  CALL InitializeValues
   
  CALL RandomPositioning
  
  CALL OverlapCheck

  CALL InitQuat  !quat.f90

  !CALL DipolePosition

  CALL ZeroValues
  
 ! CALL RandomizePos

END SUBROUTINE Initialize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

!!****m* BrownianDynamics/Restart
!!
!! NAME
!!  Restart
!!
!! SYNOPSIS
!!  This SUBROUTINE save data for restarting when the code is stopped.
!!  The code restarts using this data and continues to the final.
!!
!! NOTES
!!  Located in file init.f90
!!
!! CREATION DATE
!!    04/22/2013
!!
!! AUTHOR
!!    Ronal De La Cruz
!!
!! COPYRIGHT
!!    22-Apr-2013 - Ronal De La Cruz
!!**********************************************
SUBROUTINE Restart
  USE Globals
  IMPLICIT NONE
  
  Integer :: i
  Integer :: j
  !Real(8), dimension(2,2) :: readtime
  
  !type ( Vector ) :: disp, Contain  !, tocenter ! Displacement vectors.
  
  open( unit = 3010, file = 'Restart_posit.out', status = 'old' )
  open( unit = 3110, file = 'Restart_orient.out', status = 'old' )
  open( unit = 3210, file = 'Restart_point_data.out', status = 'old' )
  open( unit = 3310, file = 'Restart_AbsPosit.out', status = 'old' )
  
  DO i=1,gNPart
     read(3010,3011) gPart(i) % pos
     read(3310,3011) gPartAbs(i) % pos
     read(3110,3111) gT, orient( i ) % x ,orient( i ) % y ,orient( i ) % z   !in the three directions
3111 format( 4(1X, E24.15E3))
3011 format( 3(1X, E24.15E3))
  ENDDO
  !DO j=2,2
  ! DO i=1,2
  read(3210, *) gT !readtime(j,i)
  !ENDDO
  !ENDDO
  !gT = readtime(2,1)
  print*, 'Restarting from Time(gT) =',gT
  
END SUBROUTINE Restart


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!!************************************************************************
!!****m* BrownianDynamics/OpenFiles
!!
!! NAME
!!  OpenFiles
!!
!! SYNOPSIS
!!  This SUBROUTINE opens all the files USEd in the simulation with
!!  the exception of the input file since that one contains the path to
!!  the output files.
!!
!!  Files:
!!   sum.out - a summary of the input and output data from the simulation.
!!
!!   trac.out - the positions of all the externally forced particles at 1000
!!   evenly spaced intervals throughout the simulation.
!!   
!!   bath.out - the positions of all the Brownian particles at 1000 evenly
!!   spaced intervals throughout the simulation.
!!
!!   rd.out - the radial distribution of distances between different externally
!!   forced particles.
!!
!!   dpxy.out - the average density profile projected in the X-Y plane.
!!
!! NOTES
!!  Located in file init.f90
!!
!! CREATION DATE
!!    02/28/05
!!
!! AUTHOR 
!!    James Swan
!!
!! COPYRIGHT
!!    28-Feb-2005 - James William Swan
!!***************************************************************************************
!!***************************************************************************************
!!***************************************************************************************
SUBROUTINE OpenFiles
  USE Globals
  IMPLICIT NONE
  
  INTEGER :: i ! Index counter.
  INTEGER :: j, iu1, iu2, iu3 ! Osmotic motor counter and unit for new files
  character ( Len = 3 ):: jstring ! File identification for each osmotic motor
  
  gNFiles = 6
  
  i = index( gFileLoc, ' ' ) - 1
  
  open( unit = 102, file = gFileLoc( 1 : i ) // '/sum'//TRIM(gsim)//'.out', status = 'unknown' )
  
  IF (gMovies .eq. 1) THEN
     open( unit = 103, file = gFileLoc( 1 : i ) // '/trac.out', status = 'unknown' )
     open( unit = 104, file = gFileLoc( 1 : i ) // '/bath.out', status = 'unknown' )
  END IF
  
  !open(  unit = 1800, file = gFileLoc( 1 : i ) // '/rotorderpar.out', status = 'unknown' )
  !open(  unit = 1900, file = gFileLoc( 1 : i ) // '/transorderpar.out', status = 'unknown' ) 
  !open( unit = 1000, file = gFileLoc( 1 : i ) // '/newabspos.out', status = 'unknown' )
  open( unit = 2000, file = gFileLoc( 1 : i ) // '/Restart_pos_orient.out', status = "unknown" )  !RONAL-LUIS RESTARTING
 ! open( unit = 1600, file = gFileLoc( 1 : i ) // '/restartdata.out', status = "unknown" )  !
 ! open( unit = 1601, file = gFileLoc( 1 : i ) // '/restarttime.out', status = "unknown" )  !
  open( unit = 1400, file = gFileLoc( 1 : i ) // '/periodic_pos_orient.out', status = 'unknown' )  !

   OPEN( UNIT = 3, FILE = gFileLoc( 1 : i ) // '/TransClustStats.out', STATUS = 'REPLACE', ACTION = 'WRITE' )
   OPEN( UNIT = 4, FILE = gFileLoc( 1 : i ) // '/StdyStClusterStats.out', STATUS = 'REPLACE', ACTION = 'WRITE')
   OPEN( UNIT = 5, FILE = gFileLoc( 1 : i ) // '/ClusterDist.out', STATUS = 'REPLACE', ACTION = 'WRITE')
   OPEN( UNIT = 6, FILE = gFileLoc( 1 : i ) // '/Fractal.out', STATUS = 'REPLACE', ACTION = 'WRITE')
!   OPEN( UNIT = 7, FILE = gFileLoc( 1 : i ) // '/ProbDist.out', STATUS = 'REPLACE', ACTION = 'WRITE')
   OPEN( unit = 8, FILE = gFileLoc( 1 : i ) // '/RadDist.out', status = 'unknown' )
   OPEN( UNIT = 9, FILE = gFileLoc( 1 : i ) // '/TransPopulat.out', STATUS = 'REPLACE', ACTION = 'WRITE' )
  ! OPEN( unit = 1010, file = gFileLoc( 1 : i ) // '/initconf.xyz', status = 'unknown' )
  !    open( unit = 106, file = gFileLoc( 1 : i ) // '/dpxy.out', status = 'unknown' )
  !    open( unit = 1500, file = gFileLoc( 1 : i ) // '/orient.out', status = 'unknown' )
  !    open( unit = 1610, file = gFileLoc( 1 : i ) // '/numcol.out', status = 'unknown' )
  !    open( unit = 1612, file = gFileLoc( 1 : i ) // '/avgvel.out', status = 'unknown' )
  !    open( unit = 1613, file = gFileLoc( 1 : i ) // '/dpyz.out', status = 'unknown' )
  !    open( unit = 1614, file = gFileLoc( 1 : i ) // '/gpxy.out', status = 'unknown' )
  !    open( unit = 1615, file = gFileLoc( 1 : i ) // '/gpyz.out', status = 'unknown' )
  !    open( unit = 1616, file = gFileLoc( 1 : i ) // '/aangvel.out', status = 'unknown' ) 
  !    open( unit = 1617, file = gFileLoc( 1 : i ) // '/cm.out', status = 'unknown' ) 
  !    open( unit = 1618, file = gFileLoc( 1 : i ) // '/motorforce.out', status = 'unknown' ) 
  !    open( unit = 1619, file = gFileLoc( 1 : i ) // '/position'//TRIM(gsim)//'.plo', status = 'unknown' ) 
  !    open( unit = 1620, file = gFileLoc( 1 : i ) // '/bath_pos.plo', status = 'unknown' )
  !    open( unit = 1621, file = gFileLoc( 1 : i ) // '/bath_angles.plo', status = 'unknown' )
  !    open( unit = 1622, file = gFileLoc( 1 : i ) // '/angles'//TRIM(gsim)//'.plo', status = 'unknown' )
  !    open( unit = 1623, file = gFileLoc( 1 : i ) // '/velolang'//TRIM(gsim)//'.plo', status = 'unknown' )
 END SUBROUTINE OpenFiles
!!***************************************************************************************
!!***************************************************************************************



!!***************************************************************************************
!!***************************************************************************************
!!****m* BrownianDynamics/StartDataRead
!!
!! NAME
!!  StartDataRead
!!
!! SYNOPSIS
!!  This SUBROUTINE begins reading in data from the file conf.in
!!  This data includes the number of particles and the sizes of
!!  the various arrays to be allocated.  After allocation, the
!!  rest of the data can be read.
!!
!! NOTES
!!  Located in file init.f90
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
SUBROUTINE StartDataRead
  USE Globals
  IMPLICIT NONE
  
  CALL ReadJunk                     
  
  read( 101, * ) gFileLoc           ! exp
  
  CALL ReadJunk
  

  read( 101, * ) gNrun              ! Total Number of runs 
  read( 101, * ) gRestart           ! Restarting Flag
  read( 101, * ) gNPart             ! Total Number of particles
  read( 101, * ) gNMove             ! Number of motor particles
  read( 101, * ) gPhi               ! Volume fraction of particles  
  read( 101, * ) gPhip              ! Volume fraction parameter
  read( 101, * ) gPhim              ! Volume fraction of the motor
  read( 101, * ) gPhib              ! Volume fraction of the bath
  read( 101, * ) gForceMode         ! Simulation mode ( 1: force, 2: velocity, 3: mixed )
  read( 101, * ) gNDimensionality   ! Dimensionality of the system: q2D or 3D
  read( 101, * ) gFreeRotor         ! Free rotor mode: 0 - NO, 1 - Yes
  read( 101, * ) gFixedCenterMass   ! Fixed center of mass - rotor: 0 - No, 1 - Yes
  read( 101, * ) gSpecifiedFluxMode ! Specified flux problem: 0 - No, 1 - Yes
  read( 101, * ) gParaMag           ! Type of magnetism: TRUE for paramagnetic, FALSE for ferromagnetic
  read( 101, * ) gH(1)              ! Magnetic field in x direction: 0 - no magnetic field, 1 - magnetic field
  read( 101, * ) gH(2)              ! Magnetic field in y direction: 0 - no magnetic field, 1 - magnetic field
  read( 101, * ) gH(3)              ! Magnetic field in z direction: 0 - no magnetic field, 1 - magnetic field
  read( 101, * ) gLag               ! Langevin Parameter
  read( 101, * ) gDDinter           ! Dipole-Dipole interaction factor
  read( 101, * ) gRepulinter        ! Repulsive interaction factor
  read( 101, * ) gLambda            ! Dipolar coupling constant parameter
  read( 101, * ) gPe                ! Peclet number: F / (kT/a)    By Luis
  read( 101, * ) gSh                ! Shifted dipole displacement from center of mass
  read( 101, * ) gSigma             ! Lennard Jones Characteristic Length  
  read( 101, * ) gEps               ! Lennard Jones to strength  parameter  
  read( 101, * ) gBrownianFlag      ! Brownian Flag for BM - by Luis
  read( 101, * ) gUSEAV             ! Avarage velocity
  
  CALL ReadJunk
  
  read( 101, * ) gTf                ! End time 
  read( 101, * ) gDT                ! Step size
  read( 101, * ) gNDT               ! Number of time steps to record - Luis
  read( 101, * ) gNT                ! Time calculation
  read( 101, * ) gRestTiming        ! Time to write for restarting   
  read( 101, * ) gNDTCLUS           ! Cluster analysis 

  CALL ReadJunk
  
  read( 101, * ) gMovies
  read( 101, * ) gNFrames
  
  CALL ReadJunk
  
  read( 101, * ) gNCells( 1 )
  read( 101, * ) gNCells( 2 )
  read( 101, * ) gNCells( 3 )
  
  CALL ReadJunk
  
  read( 101, * ) gDPBins( 1 ) 
  read( 101, * ) gDPBins( 2 ) 
  read( 101, * ) gDPBins( 3 )
  
  CALL ReadJunk
  
  read( 101, * ) gPDBins( 1 ) 
  read( 101, * ) gPDBins( 2 ) 
  read( 101, * ) gPDBins( 3 ) 
  read( 101, * ) gRDBins
  read( 101, * ) gThetaBins
  
END SUBROUTINE StartDataRead
!!***************************************************************************************
!!***************************************************************************************

!!***************************************************************************************
!!***************************************************************************************
!!****m* BrownianDynamics/FinishDataRead
!!
!! NAME
!!  FinishDataRead
!!
!! SYNOPSIS
!!  This SUBROUTINE finishes reading in the data fromt the conf.in file.
!!  Now that the memory is allocated, specific details about the type of
!!  simulation including the size of externally forced particles, and how
!!  hard those particles are being forced can be read in.
!!
!! NOTES
!!  Located in file init.f90
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
SUBROUTINE FinishDataRead
  USE Globals
  IMPLICIT NONE
  
  Integer :: i ! Loop counter.
  type ( Vector ) :: dvecs, MinImage
  
  CALL ReadJunk
  
  read ( 101, * ) gRandLoc        ! Relocate particles: 1 - ranDOmly, 0 - other side of motor
  
  CALL ReadJunk
  
  read ( 101, * ) gRotateFlag     ! Which mode: 0 - no rotation, 1 - run and tumble, 2 - diffusive
  read ( 101, * ) gRotateFactor   ! Rotation parameter 0 - no rotation, 1 - full rotation
  read ( 101, * ) gRotateFactorb  ! Rotation of the bath particles: 0 - no rotation, 1 - full rotation
  read ( 101, * ) gDDinterRot     ! Rotation due to DD interactions -----------   Luis  -------------
  read ( 101, * ) gNoBathColFlag  ! Which collision mode: 0 - non-ideal gas, 1 - ideal gas
  read ( 101, * ) hsMode          ! Which HS motor mode: 0 - Move the motor, 1 - Don't move motor
  read ( 101, * ) USEconfig2      ! Use equilibrium configuration: 0 - No, 1 - Yes

END SUBROUTINE FinishDataRead


SUBROUTINE ParticleOrientation
  USE Globals
  USE VectorFunctions
  
  IMPLICIT NONE
  INTEGER :: i
  REAL(8) :: magu, magn
  REAL(8), DIMENSION(gND) :: n, u
  CALL ranDOm_seed()
  

  DO i = 1, gNPart
     CALL ranDOm_number(n)
     gN(i,:) = 2.D0 * n(:) - 1.D0
     magn = DSQRT(vdot(gN(i,:),gN(i,:)))
     gN(i,:) = gN(i,:) / magn
  END DO
  
END SUBROUTINE ParticleOrientation


!!***************************************************************************************
!!***************************************************************************************
!!****m* BrownianDynamics/ReadJunk
!!
!! NAME
!!  ReadJunk
!!
!! SYNOPSIS
!!  This SUBROUTINE reads the three blank lines that separate the
!!  different sections of data in the conf.in input file.
!!
!! NOTES
!!  Located in file init.f90
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
SUBROUTINE ReadJunk
  
  read( 101, * )
  read( 101, * )
  read( 101, * )
  
END SUBROUTINE ReadJunk
!!***************************************************************************************
!!***************************************************************************************

!!***************************************************************************************
!!***************************************************************************************
SUBROUTINE AllocateArrays
  USE Globals
  IMPLICIT NONE
  
  Integer :: mult ! The number of particle-particle displacements among
  ! externally forced particles.
  
  mult = gNPart * ( gNPart - 1 ) / 2
  
  allocate( gPart( gNPart ) )
  allocate( gOldPos( gNPart ) )
  allocate( gVel( gNPart ) )
  allocate( gAV( gNPart ) )
  allocate( gDP( gDPBins( 1 ), gDPBins( 2 ), gDPBins( 3 ), gNMove + 1 ) )
  allocate( gPD3D( gDPBins( 1 ), gDPBins( 2 ), gDPBins( 3 ) ) )
  allocate( gPDC( gThetaBins, gNMove ) )
  allocate( gPD( gPDBins( 1 ), gPDBins( 2 ), gPDBins( 3 ) ) )
  !allocate( gPD( mult, gPDBins( 1 ), gPDBins( 2 ), gPDBins( 3 ) ) )
  !allocate( gRD( mult, gRDBins ) )
  ALLOCATE( gRD(gRDBins) )
  allocate( gBVec( gNPart ) )
  
  allocate( orient( gNPart ) ) 
  allocate( phian( gNPart ) )
  allocate( position( gNPart ) )
  allocate( gPartold( gNPart ) )
  allocate( orientold( gNPart ) )
  allocate(gPartAbs(gNPart) )    !----- Luis ----
  
  allocate( dvec( gNPart ) )
  allocate( coltime( gNPart) )
  allocate( sgn( gNPart) )
  
  allocate( gOPMotor( gNPart ) )
  
  allocate( rxnp( gNPart ) )
  allocate( gCosThetaRxn( gNPart ) )
  
  allocate( Ko( gNPart ) )
  allocate( AKo( gNPart ) )          
  allocate( gProductColor( gNPart ) )
  
  allocate( gHSF( gNPart ) )
  allocate( gHSFRh( gNPart ) )
  allocate( gHSFnRh( gNPart ) )
  allocate( gHST( gNPart ) )
  
  allocate( gInitOrient( gNPart ) )
  allocate( gAngVel( gNPart ) )
  allocate( gAAngVel( gNPart ) )
  
  allocate( gTensileF( gNPart ) )
  
  ALLOCATE( gCDimH(gND) )
  ALLOCATE( gCDim(gND) )
  ALLOCATE( gR(gNPart, gND) )
  ALLOCATE( gN(gNPart, gND) )
  ALLOCATE( gNdip(gNPart, gND) )  !Ronal
  ALLOCATE( gFor(gNPart,gND) )
  ALLOCATE( gTor(gNPart,gND) )
  ALLOCATE( RG2(gNPart) )
  ALLOCATE( RCM(gNPart, gND) )
  ALLOCATE( CLUSTER_INDEX_LIST(gNPart) )
  ALLOCATE( CLUSTDIST(gNPart) )
  ALLOCATE( PROB(gNPart) ) 
  ALLOCATE( NpC( gNPart ) )
  ALLOCATE( REFF( gNPart) )
  !ALLOCATE( gRDipole(gNPart,gND) )
  ALLOCATE( gShVec(gNPart,gND) )
  ALLOCATE( gQuat(gNPart,4) )
  ALLOCATE( gNRB(gND) )
  ALLOCATE( gNdipRB(gND) )
  ALLOCATE( gShVecRB(gND) )
  ALLOCATE( gHRB(gNPart,gND) )
  ALLOCATE( gTorRB(gNPart,gND) )
END SUBROUTINE AllocateArrays
!!***************************************************************************************
!!***************************************************************************************

!!***************************************************************************************
!!***************************************************************************************
!!****m* BrownianDynamics/InitializeValues
!!
!! NAME
!!  InitializeValues
!!
!! SYNOPSIS
!!  This SUBROUTINE initializes most of the values that are USEd
!!  in calculations later in the program.  This simply speeds
!!  up the program since these things DO not need to be calculated
!!  upon every loop through the program.
!!
!! NOTES
!!  Located in file init.f90
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
SUBROUTINE InitializeValues
  USE Globals
  USE VectorClass
  USE VectorFunctions
  IMPLICIT NONE
  
  Integer :: time( 3 )
  Integer :: i
  REAL(8) :: LANGFUNSQ, NDENS, L
 
  gNRB(1) = 0.D0  !rigid body orientation (mobile coordinate system)
  gNRB(2) = 0.D0
  gNRB(3) = 1.D0

  gNdipRB(1)=1.d0
  gNdipRB(2)=0.d0
  gNdipRB(3)=0.d0
 ! gShVecRB(:) = gSh * gNRB(:) !this or the below definition
 
   gShVecRB(1) = 0.d0 
   gShVecRB(2) = 0.d0
   gShVecRB(3) = gSh   
  
  gLargestRad = 1.D0                               ! The largest radius of a particle in the simulation.
  gSecondLargestRad = 1.D0                         ! The second largest radius of a particle in the simulation.
  gEQLCOUNT = 0
  
  IF ( gParaMag .EQ. .TRUE. ) gLambda = gLambda * LANGFUNSQ( gLag )

  gTVInv = 1.D0 / gTV
  gStericTVInv = gTVInv * gSteric

  IF (gPe .GT. 1.D0) gDT = gDT / gPe
  
  gPeDT = gPe * gDT

  RMN = 2.D0 / ( 1.D0 + gTV )
  RMN2 = RMN ** 2
  
  gLambda4 =   4.D0 * gLambda 
  gLambda12 = 12.D0 * gLambda
  
  gNDTRAND = 10000
  RB2 = ( 1.2D0 * 2.D0 ) ** 2  ! Include in conf.in ?
  gPi = 4.D0 * atan( 1.D0 )
  gPi43 = gPi * 4.D0 / 3.D0
  gPi2 = gPi * 2.D0
  gSQRT2DT = DSQRT( gBrownianFlag * 2.D0 * gDT )
  gDTInv = 1.D0 / ( gDT )
  
  gBrowDT   = DSQRT( (3.D0 / 2.D0) * gRotateFactor * gDT )
  gCells = gNCells( 1 ) * gNCells( 2 ) * gNCells( 3 ) 
  gDT34 = (3.D0 / 4.D0) * gDT
  gMinSep = 1d-12         ! The minimum overlap measurable by the collision detection scheme.
  gMinSepMotors = gMinSep
  gReseed = 10000         ! The number of time steps to wait before reseeding the ranDOm number generator.
  !gAvgTiming = 0.2D0 * gTf  ! Time to wait before calc. statistics
  gAvgTiming = 0.8D0 * gTF
  gAvgTiming2 = 1.D0 * gTf ! Time to wait to start using react prob for free motor
   
  gRestCount = 0     
  gOutCount = 0           ! The timing for outputting particle position data.
  gOutClusCounter = 0
  gOutTiming = gTf / ( gNDT * gDT ) ! For stats USE gOutTiming = 1  ! The timing for outputting particle position data.
  gOutClust  = gTf / ( gNDTCLUS * gDT )
  gOutFrames = gTf / ( gNFrames * gDT )
  gFramesCounter = 0
  gDenTiming = 10         !  The time for measuring the density profile.
  
  gOsmMotorCol = .false.
  gOsmMotorCount = 1
   

  IF ( gNDimensionality .EQ. 2 ) THEN !! Surface fraction when system in q2D
     gNdensity = gPhi / gPi
     L = DSQRT( gNPart * gPi / gPhi )
  ELSE IF ( gNDimensionality .EQ. 3 ) THEN !! Volume fraction when system in 3D
     gNdensity = ( 3.D0 * gPhi ) / ( 4.D0 * gPi )
     L = ( ( 4.D0 * gPi * gNPart ) / ( 3.D0 * gPhi ) ) ** (1.D0 / 3.D0)
  END IF
 

  gCDim(:)    =  L
  gCDimH(:)   = gCDim(:) / 2.D0
  gSigma2     = gSigma ** 2.D0
  gEps24      = 24.D0 * gEps 
  gCutOffWCA  = 2**(1.D0/6.D0)*gSigma
  gCutOffWCA2 = gCutOffWCA ** 2.D0
  gCutOff     = gCDimH(1)
  gCutOff2    = gCutOff ** 2.D0
  
  !write(*,*) gSigma2, gEps24, gCutOffWCA, gCutOffWCA2
  !stop

  gAreaVolumeNoDim = gCDim(1)**(dble(gNDimensionality)) ! measures the area or the volume considering gCDim(1)= gCDim(3)= gCDim(2)
  
 print*, 'gAreaVolumeNoDim', gAreaVolumeNoDim
 print*, 'gNdensity',gNdensity
!print*, gNPart/(gCDim(2)*gCDim(1))  
  
  gDR = MIN(gCDimH(1),gCDimH(2)) / DFLOAT(gRDBins)
  gRho = DFLOAT(gNPart) / ( PRODUCT(gCDim) )
  CALL iTime( time )
  seed2 = 10000 * time( 3 ) + 100 * time( 2 ) + time( 1 )
  seed2 = seed2 + time( 1 ) * time( 2 ) * ( time( 3 ) + time( 1 ) ) + time( 3 ) 
  
  
  !jmax = nint( 0.1 * ( gNPart - gNMove ) ) 
  jmax = 15
  stdvFlux = ( jmax / 2.0 ) ** ( 1.0 / 2.0 )
  
  ! Print data to make external statistics with particle trayectories: 0-NO, 1-Yes
  gTracStat = 0
  !        gTracStat = 1
  
  ! Print data to make external statistics with bath particle trayectories: 0-NO, 1-Yes
  gBathStat = 0
  
  ! Print data to make external statistics with osmotic rotor: 0-NO, 1-Yes
  gRotorStat = 0
  
  
END SUBROUTINE InitializeValues
!!***************************************************************************************

!!***************************************************************************************
!!****m* BrownianDynamics/ZeroValues
!!
!! NAME
!!  ZeroValues
!!
!! SYNOPSIS
!!  This SUBROUTINE zeros out all of the arrays USEd in the program.
!!
!! NOTES
!!  Located in file init.f90
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
SUBROUTINE ZeroValues
USE Globals
use omp_lib
IMPLICIT NONE

    Integer :: i, j, k, l ! Loop counters.
!$omp parallel 

    !$omp do
    DO i = 1, gDPBins( 1 )
        DO j = 1, gDPBins( 2 )
            DO k = 1, gDPBins( 3 )
                DO l = 1, gNMove + 1
                    gDP( i, j, k, l ) = 0
                END DO
            END DO
        END DO
    END DO
    !$omp end do
    
    !$omp do
    DO i = 1, gPDBins( 1 )
        DO j = 1, gPDBins( 2 )
            DO k = 1, gPDBins( 3 )
                gPD( i, j, k ) = 0
            END DO
        END DO
    END DO
    !$omp end do
    
    !DO i = 1, gPDBins( 1 )
        !DO j = 1, gPDBins( 2 )
            !DO k = 1, gPDBins( 3 )
                !gPD( l, i, j, k ) = 0
            !END DO
        !END DO
    !END DO
    !$omp do
    DO i = 1, gRDBins
        !DO j = 1, gNPart * ( gNPart - 1 ) / 2
            gRD( i ) = 0.D0
        !END DO
    END DO
    !$omp end do

    !$omp do
    DO i = 1, gNMove
        gAV( i ) = Vector( 0.0, 0.0, 0.0 ) 
        AKo( i ) = 0.0
    END DO
    !$omp end do
!$omp end parallel

    gAP = 0
    
   ! DO i = 1, gNPart
   !     coltime( i ) = 0
   !     gPartAbs(i) % pos % x = 0.   !----- Luis ----
   !     gPartAbs(i) % pos % y = 0.
   !     gPartAbs(i) % pos % z = 0.
   ! END DO 



  ! DO i = 1, gNPart
  !      gPart(i) % poscounter % x = 0.d0
  !      gPart(i) % poscounter % y = 0.d0 
  !      gPart(i) % poscounter % z = 0.d0   
  !         
  ! END DO 

END SUBROUTINE ZeroValues
!!****************************************************************************



!!****m* BrownianDynamics/CloseFiles
!!
!! NAME
!!  CloseFiles
!!
!! SYNOPSIS
!!  This SUBROUTINE closes all the files USEd in the program.
!!
!! NOTES
!!  Located in file init.f90
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
SUBROUTINE CloseFiles
  USE Globals
  IMPLICIT NONE
  
  Integer :: i ! Loop counters.
  
  DO i = 1, gNFiles
     close( 100 + i )
  END DO
  
  !    close(1500)
END SUBROUTINE CloseFiles
!!****************************************************************************
!!****************************************************************************


!!****************************************************************************
!!****************************************************************************
recursive SUBROUTINE generate_X_value( Xint )
  USE Globals
  USE Knuth_RanDOm
  IMPLICIT NONE
  Integer, INTENT(OUT) :: Xint
  Integer :: Kint, factorial, i
  Real( 8 ) :: PX, PK, FX, FXX, fctx, c, Xbath
  Real, dimension ( 4 ) :: randnum
  Real( 8 ) :: ran0
  Real, allocatable, dimension( : ) :: u 
  allocate( u( 3 * gNPart ) )
  
  CALL rand( u, 3 * gNPart ) ! mean = 0, stdv = 1 
  c =  u( 1 )
  Xbath = stdvFlux * c + Ko( 1 )
  Xint = nint( Xbath )
  
  !write(*,*) 'Check point 1', Xbath, Xint, jmax
  !paUSE
  
  if ( ( Xint .lt. 0 ) .or. ( Xint .gt. jmax ) ) then
     CALL generate_X_value( Xint )
  else
     PX = Ko( 1 ) ** Xint * exp( -1.0 * Ko( 1 ) ) / factorial( Xint )
     Kint = nint( Ko( 1 ) )
     PK = Ko( 1 ) ** Ko( 1 ) * exp( -1.0 * Ko( 1 ) ) / factorial( Kint )
     
     fctx = 1.10 ! Slightly greater than 1 
     FX = fctx * PK * exp( -1.0 * ( Xint - Ko( 1 ) ) ** 2 / jmax )
     !FX = PK
     
     randnum( 1 ) = ran0( seed2 )
     !write(*,*) 'ranDOm', randnum( 1 )
     !paUSE
     FXX = FX * randnum( 1 )
     
     if ( FXX .gt. PX ) then
        CALL generate_X_value( Xint )
     END if
  END if
  
  !write(*,*) 'Check point 2', Xint, FXX, PX, PK, factorial( Xint )
  !paUSE
  
END SUBROUTINE generate_X_value
!!****************************************************************************

!!****************************************************************************
recursive function factorial( n ) result (fact)
  IMPLICIT NONE
  Integer :: fact
  INTEGER, intent( IN ) :: n
  
  if ( n .le. 0 ) then
     fact = 1
  else
     fact = n * factorial( n - 1 )
  END if
  
END function factorial
!!****************************************************************************

REAl(8) FUNCTION LANGFUNSQ(ALPHA)
  USE GLOBALS
  IMPLICIT NONE
  REAL(8) :: ALPHA, COTH
  
  IF ( ALPHA .EQ. 0.D0 ) THEN
     LANGFUNSQ = 0.D0			!! Udated from 0 to 1 on July 15, 2015 by Angel González !cambié a cero porque regula ferromag y paramag-- chequear bien -- Ronal
  ELSE IF ( ALPHA .GT. 0.D0 ) THEN
     COTH = 1.D0 / TANH( ALPHA )
     LANGFUNSQ = COTH - 1.D0 / ALPHA
     LANGFUNSQ = LANGFUNSQ * LANGFUNSQ
  END IF
  
END FUNCTION LANGFUNSQ
!!************************************************************************
!!************************************************************************
!!************************************************************************


SUBROUTINE RandomPositioning
  USE Globals
  USE VectorFunctions
  IMPLICIT NONE
  INTEGER :: i
  REAL(8), DIMENSION(gND) :: ran
  
  CALL RANDOM_SEED()
  
  DO i = 1, gNPart
     CALL RANDOM_NUMBER(ran)
     gR(i,:) = 2.D0 * ran - 1.D0

     IF ( gNDimensionality .EQ. 2 ) THEN !! Added to define dimensionality on August 15, 2015
	gR(i,3) = 0.0
     END IF

     gR(i,:) = gR(i,:) * gCDimH(:)
  END DO

END SUBROUTINE RandomPositioning
!!************************************************************************
!!************************************************************************

SUBROUTINE OverlapCheck
  USE GLOBALS
  USE VectorFunctions
  IMPLICIT NONE
  
  LOGICAL :: OVERLAP_FLAG
  INTEGER :: I, J
  REAL(8) :: TOL, RMAG, RMAG2, A, D, D2
  REAL(8), DIMENSION(gND) :: RIJ, u
  
  A = 1.D0  ! particle radius
  OVERLAP_FLAG = .TRUE.
  TOL = 1.0D-10
  
  DO WHILE ( OVERLAP_FLAG )
     OVERLAP_FLAG = .FALSE.
     DO I = 1, gNPart - 1
        DO J = I + 1 , gNPart
           RIJ(:) = gR(i,:) - gR(j,:)
           RIJ  = MinImage( RIJ )
           RMAG2 = vdot( RIJ, RIJ)
           RMAG  = DSQRT(RMAG2)
           D  = 2.D0*A - RMAG
           D2 = D / 2.D0
           IF ( D .GT. 0.D0 + TOL ) THEN
              OVERLAP_FLAG = .TRUE.
              u(:) = RIJ(:) / RMAG
              gR(i,:) = gR(i,:) + D2 * u(:)
              gR(j,:) = gR(j,:) - D2 * u(:)
              CALL ContainParticle( I )
              CALL ContainParticle( J )
           END IF
        END DO
     END DO
  END DO
END SUBROUTINE OverlapCheck


