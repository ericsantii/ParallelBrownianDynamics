!!****h* BrownianDynamics/Globals
!! NAME
!!  Globals
!!
!! SYNOPSIS
!!  This module contains all of the global variables for running
!!  the simulation with the exception of randomization variables
!!  and initial particle distribution variables. Those are
!!  contained in separate modules.  This module is called by virtually
!!  every subroutine for access to this data.
!!
!! NOTES
!!  There are other variables in this module not cited in the
!!  documentation.  These are only values calculated in advance to
!!  spped the codes' execution.
!!
!!  Located in file glob.f90
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
module Globals
use CellClass
use ParticleClass
implicit none
save




REAL(8) :: gRotOrderPar
REAL(8) :: gTransOrderPar
REAL(8) :: gTransSUM
REAL(8) :: gRotSUM

TYPE(VECTOR), DIMENSION(:), ALLOCATABLE :: gInitOrient
TYPE(VECTOR) :: gReciprocalVec
    
    !!****v* Globals/gND
    !! NAME
    !!  gND
    !! NOTES
    !!  Dimensionality: 2 for 2D & 3 for 3D
    !!******
    INTEGER ::  gND

    !!****v* Globals/gNDimensionality
    !! NAME
    !!  gND
    !! NOTES
    !!  Dimensionality: 0 for q2D or 1 for 3D
    !!******
    INTEGER ::  gNDimensionality
    
    !!****v* Globals/gFileLoc
    !! NAME
    !!  gFileLoc
    !! NOTES
    !!  The path to where the results are stored.
    !!******
    Character( Len = 50 ) :: gFileLoc
    
    !!****v* Globals/gPeDT
    !! NAME
    !!  gPeDT
    !! NOTES
    !!  Peclet number multiplied by time step
    !!******
    REAL(8) :: gPeDT

    
    !!****v* Globals/gTScaleChange
    !! NAME
    !!  gTScaleChange
    !! NOTES
    !!  Scaling for Peclet values
    !!******
    REAL(8) :: gTScaleChange

    !!****v* Globals/gRestart
    !! NAME
    !!  gRestart
    !! NOTES
    !!  Restarting run flag
    !!******
    INTEGER :: gRestart

    !!****v* Globals/gNPart
    !! NAME
    !!  gNPart
    !! NOTES
    !!  The number of particles in the simulation.
    !!******
    Integer :: gNPart

    !!****v* Globals/gNMove
    !! NAME
    !!  gNMove
    !! NOTES
    !!  The number of particles experiencing other non-brownian forces.
    !!******
    Integer :: gNMove

    !!****v* Globals/gPart
    !! NAME
    !!  gPart
    !! NOTES
    !!  The data for all of the particles in the simulation.
    !!******
    type ( Particle ), ALLOCATABLE, DIMENSION( : ) :: gPart

    !!****v* Globals/gPartAbs
    !! NAME
    !! gPartAbs
    !! NOTES
    !! The data for absolute position of particles
    !!******
    type ( Particle ), ALLOCATABLE, DIMENSION( : ) :: gPartAbs

    !!****v* Globals/gLCell
    !! NAME
    !!  gLCell
    !! NOTES
    !!  The 3-D linked cells for doing faster collision detection.
    !!******
    type ( Cell ), ALLOCATABLE, DIMENSION( :, :, : ) :: gLCell

    !!****v* Globals/gPhi
    !! NAME
    !!  gPhi
    !! NOTES
    !!  The total volume fraction of particles in the periodic cell.
    !!******
    REAL( 8 ) :: gPhi

    !!****v* Globals/gT
    !! NAME
    !!  gT
    !! NOTES
    !!  The current simulation time.
    !!******
    REAL( 8 ) :: gT
    
    !!****v* Globals/gDT
    !! NAME
    !!  gDT
    !! NOTES
    !!  The simulation time step.
    !!******
    REAL( 8 ) :: gDT

    !!****v* Globals/gNT
    !! NAME
    !!  gNT
    !! NOTES
    !! Number of time steps for which to calculate
    !!******
    REAL(8) :: gNT
    
    !!****v* Globals/gTf
    !! NAME
    !!  gTf
    !! NOTES
    !!  The end time of the simulation.
    !!******
    REAL( 8 ) :: gTf

    !!****v* Globals/gR
    !! NAME
    !!  gR
    !! NOTES
    !!  The position vector of the particles
    !!******
    REAL( 8 ), ALLOCATABLE, DIMENSION(:,:) :: gR
    
    !!****v* Globals/gN
    !! NAME
    !!  gN
    !! NOTES
    !!  The orientation of the particles
    !!******
    REAL( 8 ), ALLOCATABLE, DIMENSION(:,:) :: gN

    !!****v* Globals/gNdip
    !! NAME
    !!  gNdip
    !! NOTES
    !!  The orientation of the magnetic dipole of the particles
    !!******
    REAL( 8 ), ALLOCATABLE, DIMENSION(:,:) :: gNdip

    !!****v* Globals/gCDim
    !! NAME
    !!  gCDim
    !! NOTES
    !!  The DIMENSIONs of the periodic simulation cell.
    !!******
    REAL( 8 ), ALLOCATABLE, DIMENSION(:) :: gCDim

    !!****v* Globals/gAreaVolumeNoDim
    !! NAME
    !!  gAreaVolumeNoDim
    !! NOTES
    !!  The dimensionless area or volume depends if the simulation is in 2D or 3D
    !!******
    REAL( 8 ) :: gAreaVolumeNoDim



    !!****v* Globals/gNLCells
    !! NAME
    !!  gNLCells
    !! NOTES
    !!  The number of linked cells for fast collision detection.
    !!******
    Integer, DIMENSION( 3 ) :: gNLCells
    
    !!****v* Globals/gLLS
    !! NAME
    !!  gLLS
    !! NOTES
    !!  The starting point for looping through linked cells.
    !!******
    Integer, DIMENSION( 3 ) :: gLLS
    
    !!****v* Globals/gNCells
    !! NAME
    !!  gNCells
    !! NOTES
    !!  The number of cubic cells composing the overall simulation cell.
    !!******
    Integer, DIMENSION( 3 ) :: gNCells
    
    !!****v* Globals/gDPBins
    !! NAME
    !!  gDPBins
    !! NOTES
    !!  The number of bins for generating the density profile.
    !!******
    Integer, DIMENSION( 3 ) :: gDPBins
    
    !!****v* Globals/gPDBins
    !! NAME
    !!  gPDBins
    !! NOTES
    !!  The number of bins for generating the pair distribution.
    !!****** 
    Integer, DIMENSION( 3 ) :: gPDBins
    
    !!****v* Globals/gRDBins
    !! NAME
    !!  gRDBins
    !! NOTES
    !!  The number of bins for generating the radial distribution.
    !!******
    Integer :: gRDBins
    
    !!****v* Globals/gOP
    !! NAME
    !!  gOP
    !! NOTES
    !!  The osmotic pressure for the simulation cell in one time step.
    !!******
    REAL( 8 ) :: gOP   
    
    !!****v* Globals/gOPMotor
    !! NAME
    !!  gOP
    !! NOTES
    !!  The osmotic pressure for the motor in one time step.
    !!******
    REAL( 8 ), ALLOCATABLE, DIMENSION( : ) :: gOPMotor 
    
    !!****v* Globals/gAP
    !! NAME
    !!  gAP
    !! NOTES
    !!  The average osmotic pressure.
    !!******
    REAL( 8 ) :: gAP 
    
    !!****v* Globals/gOldPos
    !! NAME
    !!  gOldPos
    !! NOTES
    !!  The last position of the externally forced particles.
    !!******
    type ( Vector ), ALLOCATABLE, DIMENSION( : ) :: gOldPos    
    
    !!****v* Globals/gVel
    !! NAME
    !!  gVel
    !! NOTES
    !!  The velocity of the externally forced particles for one time step.
    !!******
    type ( Vector ), ALLOCATABLE, DIMENSION( : ) :: gVel  
    
    !!****v* Globals/gAV
    !! NAME
    !!  gAV
    !! NOTES
    !!  The average velocity of the externally forced particles.
    !!******
    type ( Vector ), ALLOCATABLE, DIMENSION( : ) :: gAV
    
    !!****v* Globals/gDP
    !! NAME
    !!  gDP
    !! NOTES
    !!  The density profile bins.
    !!******
    
    REAL( 8 ), ALLOCATABLE, DIMENSION( :, :, :, : ) :: gDP    
    !!****v* Globals/gPD
    !! NAME
    !!  gPD
    !! NOTES
    !!  The pair distribution bins.
    !!******
    REAL( 8 ), ALLOCATABLE, DIMENSION( :, :, : ) :: gPD 
    !REAL( 8 ), ALLOCATABLE, DIMENSION( :, :, :, : ) :: gPD 
    !! Pair distribution function in 3D
    REAL( 8 ), ALLOCATABLE, DIMENSION( :, :, : ) :: gPD3D   
    !! Pair distribution function at contact
    REAL( 8 ), ALLOCATABLE, DIMENSION( :, : ) :: gPDC
    
    !!****v* Globals/gRD
    !! NAME
    !!  gRD
    !! NOTES
    !!  The radial distribution bins.
    !!******
    REAL( 8 ), ALLOCATABLE, DIMENSION( : ) :: gRD
    
    !!****v* Globals/gBVec
    !! NAME
    !!  gBVec
    !! NOTES
    !!  The Brownian displacement vectors for all the particles. 
    !!******
    type ( Vector ), ALLOCATABLE, DIMENSION( : ) :: gBVec
    
    !!****v* Globals/gNFiles
    !! NAME
    !!  gNFiles
    !! NOTES
    !!  The number of open files for the simulation.
    !!******
    Integer :: gNFiles    
    
    !!****v* Globals/gReseed
    !! NAME
    !!  gReseed
    !! NOTES
    !!  The number of time steps to wait before reseeding the random number generator.
    !!******
    Integer :: gReseed    
    
    !!****v* Globals/gForceMode
    !! NAME
    !!  gForceMode
    !! NOTES
    !!  The type of force being appied to the driven particles.
    !!******
    Integer :: gForceMode   
    
    !!****v* Globals/gOutCount
    !! NAME
    !!  gOutCount
    !! NOTES
    !!  The timing for outputting particle position data.
    !!******
    Integer :: gOutCount  

    !!****v* Globals/gRestCount
    !! NAME
    !!  gRestCount
    !! NOTES
    !!  The timing for outputting particle position and orientation data
    !!  to make restarting.
    !!******
    Integer :: gRestCount

    !!****v* Globals/gRestTiming
    !! NAME
    !!  gRestTiming
    !! NOTES
    !!  The timing for outputting a poit of restarting.
    !!******
    Integer :: gRestTiming  
    
    !!****v* Globals/gOutTiming
    !! NAME
    !!  gOutTiming
    !! NOTES
    !!  The timing for outputting particle position data.
    !!******
    Integer :: gOutTiming  
    INTEGER :: gOutFrames
    INTEGER :: gFramesCounter
    !!**********************************************************************************************************
    
    !!****v* Globals/gNDT
    !! NAME
    !!  gNDT
    !! NOTES
    !!  How many time step between calculations
    !!******
    Integer :: gNDT  
    
    !!**********************************************************************************************************
    
    !!****v* Globals/gDenTiming
    !! NAME
    !!  gDenTiming
    !! NOTES
    !!  The time for measuring the density profile.
    !!******
    Integer :: gDenTiming 
    
    !!****v* Globals/gMinSep
    !! NAME
    !!  gMinSep
    !! NOTES
    !!  The minimum overlap measurable by the collision detection scheme.
    !!******
    REAL( 8 ) :: gMinSep    
    
    !!****v* Globals/gLargestRad
    !! NAME
    !!  gLargestRad
    !! NOTES
    !!  The largest radius of a particle in the simulation.
    !!******
    REAL( 8 ) :: gLargestRad    
    
    !!****v* Globals/gSecondLargestRad
    !! NAME
    !!  gSecondLargestRad
    !! NOTES
    !!  The second largest radius of a particle in the simulation.
    !!******
    REAL( 8 ) :: gSecondLargestRad  
    
    !!****v* Globals/gAvgTiming
    !! NAME
    !!  gAvgTiming
    !! NOTES
    !!  The amount of time to wait before starting to measure the average properties.
    !!******
    REAL( 8 ) :: gAvgTiming 
    
    !!****v* Globals/gAvgTiming2
    !! NAME
    !!  gAvgTiming
    !! NOTES
    !!  The amount of time to wait before starting to use the reaction prob for free motor.
    !!******
    REAL( 8 ) :: gAvgTiming2
    
    !!****v* Globals/gMinDim
    !! NAME
    !!  gMinDim
    !! NOTES
    !!  Half the longest corner to corner distance in the simulation cell.
    !!******
    REAL( 8 ) :: gMinDim 
    
    !!****v* Globals/gDa
    !! NAME
    !! gDa
    !! NOTES
    !! Damkholer number, ratio between kinetic velocity and diffusion velocity
    !!****
    REAL( 8 ) :: gDa
    
    !!****v* Globals/gPe
    !! NAME
    !!  gPe
    !! NOTES
    !!  Peclet number, ratio between driving force and thermal energy
    !!****
    REAL( 8 ) :: gPe  

    !!***v* Globals/gLag
    !! NAME
    !! gLag
    !! NOTES
    !! Langevin parameter
    !!***
    REAL( 8 ) :: gLag
    
    !!***v* Globals/gH
    !! NAME
    !!  gH 
    !! NOTES
    !!  Magnetic field vector
    !!***
    REAL( 8 ), DIMENSION(3) :: gH

    !!***v* Globals/guseAV
    !! NAME
    !! guseAV
    !! NOTES
    !! Use Average velocity data
    !!***
    REAL( 8 ) :: guseAV 
   
    !!***v* Globals/gPhip
    !! NAME
    !! gPhip 
    !! NOTES
    !! Volume fraction parameter 
    !!***
    Integer :: gPhip 

    !!***v* Globals/gPhim
    !! NAME
    !! gPhim
    !! NOTES
    !! Motor volume fraction 
    !!***
    REAL( 8 ) :: gPhim

    !!***v* Globals/gPhib
    !! NAME
    !! gPhib
    !! NOTES
    !! Bath volume fraction 
    !!***
    REAL( 8 ) :: gPhib

    !!***v* Globals/gPhis
    !! NAME
    !! gPhib
    !! NOTES
    !! New volume fraction of the system
    !!***
    REAL( 8 ) :: gPhis

    !!***v* Globals/gPhinm
    !! NAME
    !! gPhinm
    !! NOTES
    !! New volume fraction of the motor 
    !!***
    REAL( 8 ) :: gPhinm

    !!***v* Globals/gPhinb
    !! NAME
    !! gPhinb
    !! NOTES
    !! New volume fraction of the bath 
    !!***
    REAL( 8 ) :: gPhinb

    !!***v* Globals/gNdensity
    !! NAME
    !! gNdensity
    !! NOTES
    !! Number density (in 2D or in 3D) 
    !!***
    REAL( 8 ) :: gNdensity


    !!***v* Globals/gsim
    !! NAME
    !! gsim
    !! NOTES
    !! String for number of runs 
    !!***
    character ( len = 4 ) gsim 
    
    !!***v* Globals/gNrun
    !! NAME
    !! gNrun
    !! NOTES
    !! Number of run 
    !!***
    Integer :: gNrun

    !!***v* Globals/gDDinter
    !! NAME
    !! gDDinter
    !! NOTES
    !! A factor to asignate dipole-dipole magnetic interaction 
    !!***
    REAL(8) :: gDDinter   


 !!***v* Globals/gDDinterRot
    !! NAME
    !! gDDinterRot
    !! NOTES
    !! A factor to asignate dipole-dipole magnetic rotation
    !!***
    REAL(8) :: gDDinterRot 
   
    !!***v* Globals/eps3
    !! NAME
    !!  eps3
    !! NOTES
    !!  Levi-Civita symbol of rank 3
    !!***
    REAL(8), DIMENSION(3,3,3) :: eps3


    !!***v* Globals/eps2
    !! NAME
    !!  eps2
    !! NOTES
    !!  Levi-Civita symbol of rank 3
    !!***
    REAL(8), DIMENSION(2,2) :: eps2

    !!***v* Globals/gCutOff
    !! NAME
    !! gCutOff
    !! NOTES
    !! Cut off radius for dipole-dipole magnetic interaction
    !!***
    REAL(8) ::  gCutOff 
    
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: gFor
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: gTor
    
    !!***v* Globals/gCutOff
    !! NAME
    !! gCutOff
    !! NOTES
    !! Cut off radius squared for dipole-dipole magnetic interaction
    !!***
    REAL(8) ::  gCutOff2  

    !!***v* Globals/gRepulinter
    !! NAME
    !! gRepulinter
    !! NOTES
    !! A factor to asignate repulsive interaction 
    !!***
    REAL(8) :: gRepulinter  
  
    !!***v* Globals/gLambda
    !! NAME
    !! gLambda
    !! NOTES
    !! Dipolar coupling constant 
    !!***
    REAL( 8 ) :: gLambda

    !!***v* Globals/gRotateFactorb
    !! NAME
    !! gRotateFactorb
    !! NOTES
    !! This parameter asigned or not rotation to the bath particles
    !!***
    REAL( 8 ) :: gRotateFactorb

    type ( Vector ), ALLOCATABLE, DIMENSION( : ) :: orient      
    type ( Vector ), ALLOCATABLE, DIMENSION( : ) :: phian       
    type ( Vector ), ALLOCATABLE, DIMENSION( : ) :: position    
    type ( Vector ), ALLOCATABLE, DIMENSION( : ) :: gPartold
    type ( Vector ), ALLOCATABLE, DIMENSION( : ) :: orientold

    INTEGER :: gNDTRAND
    type ( Vector ), ALLOCATABLE, DIMENSION( : ) :: dvec        
    type ( Vector ) :: dcmfixed
    Integer :: gRotateFlag
    REAL( 8 ) :: gRotateFactor
    REAL( 8 ) :: gBrownianFlag     ! Brownian motion flag: 1 for BM, 0 for no BM
    Integer :: gNoBathColFlag
    Integer :: gRandLoc
    Integer :: gCountColNoReactive
    Integer :: gCountCol
    Integer :: gCountColReact
    Integer :: hsMode
    Integer :: useconfig2
    Integer :: gOldNumCol
    Integer :: gNumCol
    Integer :: gNumColReactSide
    Integer :: gOldNumColReactSide
    Integer :: gNumColNoReactSide
    Integer :: gOldNumColNoReactSide
    Integer :: gFreeRotor
    Integer :: g2Dproblem
    Integer :: gFixedCenterMass

    Logical :: gOsmMotorCol
    REAL( 8 ) :: gOsmMotorCount
    REAL( 8 ), ALLOCATABLE, DIMENSION( : ) :: coltime
    REAL( 8 ), ALLOCATABLE, DIMENSION( : ) :: sgn
    REAL( 8 ) :: gMinSepMotors
    REAL( 8 ) :: dfixed
    REAL( 8 ) :: gtol
   
    
    REAL, ALLOCATABLE, DIMENSION( : ) :: rxnp
    REAL, ALLOCATABLE, DIMENSION( : ) :: gCosThetaRxn
    
    REAL, ALLOCATABLE, DIMENSION( : ) :: Ko
    REAL, ALLOCATABLE, DIMENSION( : ) :: AKo
    Integer, ALLOCATABLE, DIMENSION( : ) :: gProductColor
    REAL( 8 ) :: stdvFlux
    Integer :: jmax
    Integer :: Xfinal
    Integer :: gSpecifiedFluxMode
    
    !!****v* Globals/gCollisionFlag
    !! NAME
    !!  gCollisionFlag
    !! NOTES
    !!  A flag for detecting whether or not there's a collision.
    !!******
    Logical :: gCollisionFlag

    type ( Vector ), DIMENSION( : ), ALLOCATABLE :: gHSF
    type ( Vector ), DIMENSION( : ), ALLOCATABLE :: gHSFnRh
    type ( Vector ), DIMENSION( : ), ALLOCATABLE :: gHSFRh
    type ( Vector ), DIMENSION( : ), ALLOCATABLE :: gHST
    
    ! For rotor
    type ( Vector ), DIMENSION( : ), ALLOCATABLE :: gAngVel
    type ( Vector ), DIMENSION( : ), ALLOCATABLE :: gAAngVel
    type ( Vector ) :: gFCM
    type ( Vector ), ALLOCATABLE, DIMENSION( : ) :: gTensileF

    ! Factors of pi used in various calculations.
    REAL( 8 ) :: gPi, gPi43, gPi2

    ! ( 2 * DT ) ^ 1/2, used for Brownian steps.
    REAL( 8 ) :: gSqrt2DT

    ! 1 / DT, used for statistics
    REAL( 8 ) :: gDTInv
    
    ! The x-direction reference point for fixing the origin of the density profile.
    REAL( 8 ) :: gDPXFixed
    ! Half the DIMENSIONs of the simulation cell.
    REAL( 8 ), ALLOCATABLE, DIMENSION( : ) :: gCDimH
    ! The DIMENSIONs of the density profile bins.
    REAL( 8 ), DIMENSION( 3 ) :: gDPDim
    ! The total number of cubic cells in the simulation cell.
    Integer :: gCells
    ! The total number of linked cells for fast collision detection.
    Integer :: gLCells
    
    ! A temporary variable for storing the last force mode used.
    Integer :: gFMTemp

    REAL( 8 ) :: RMN, RMN2
    REAL( 8 ) :: gLagDT34
    Real( 8 ) :: gSigma
    Real( 8 ) :: gSigma2
    Real( 8 ) :: gEps24
    Real( 8 ) :: gEps
    Real( 8 ) :: gCutOffWCA
    Real( 8 ) :: gCutOffWCA2
    
    !!****v* Globals/gParaMag
    !! NAME
    !!  gParaMag
    !! NOTES
    !!  A flag for specifying the type of magnetism of the particles:
    !!  TRUE for paramagnetic partices, FALSE for ferromagnetic
    !!******
    LOGICAL :: gParaMag
 
    !!****v* Globals/gSteric
    !! NAME
    !!  gSteric
    !! NOTES
    !!  Steric repulsion DIMENSIONless parameterr defined as \lambda_v (Satoh A.)
    !!******
    REAL(8) :: gSteric

    !!****v* Globals/gTV
    !! NAME
    !!  gTV
    !! NOTES
    !!  Ratio of thickness of steric layer \delta to radius of particle (Satoh A.)
    !!******
    REAL(8) :: gTV

    REAL(8) :: gTVInv        ! 1.d0 / gTV
    REAL(8) :: gStericTVInv  ! \lambda_v / gTV
    REAL(8) :: gBrowDT 
    REAL(8) :: gDT34 
    REAL(8) :: gLambda12
    REAL(8) :: gLambda4

    Integer :: seed2  !-->Important to Random-number generation
 
    Integer :: gThetaBins


    Integer :: gMovies
    Integer :: gTracStat
    Integer :: gBathStat
    Integer :: gRotorStat
    INTEGER :: gNFrames

!******************************************************************
!**************** Cluster Statistics Variables ********************
!******************************************************************

REAL(8)                              :: RB2
REAL(8)                              :: EPS
REAL(8)                              :: THETA
REAL(8)                              :: S_MEAN
REAL(8)                              :: L_MEAN
REAL(8)                              :: NpC_MEAN
REAL(8)                              :: RG2_MEAN
REAL(8)                              :: REFF_MEAN
REAL(8)                              :: NPC_EQL
REAL(8)                              :: RG2_EQL
REAL(8)                              :: REFF_EQL
REAL(8)                              :: S_EQL
REAL(8)                              :: L_EQL
REAL(8)                              :: THETA_EQL
REAL(8)                              :: EPS_EQL
REAL(8)                              :: nSINGLETS, nDOUBLETS, nTRIPLETS, nQUADRUPLETS, nPENTUPLETS
INTEGER                              :: SINGLETS, DOUBLETS, TRIPLETS, QUADRUPLETS, PENTUPLETS
INTEGER                              :: NC
INTEGER                              :: gEQLCOUNT
INTEGER                              :: gNDTCLUS
INTEGER                              :: gOutClusCounter
INTEGER                              :: gOutClust
INTEGER, ALLOCATABLE, DIMENSION(:)   :: NpC
REAL(8), ALLOCATABLE, DIMENSION(:)   :: RG2
REAL(8), ALLOCATABLE, DIMENSION(:)   :: REFF
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: RCM
INTEGER, ALLOCATABLE, DIMENSION(:)   :: CLUSTER_INDEX_LIST
INTEGER, ALLOCATABLE, DIMENSION(:)   :: CLUSTDIST
REAL(8), ALLOCATABLE, DIMENSION(:)   :: PROB
REAL(8)                              :: gDR
REAL(8)                              :: gRho

    REAL(8) :: gSh
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: gQuat
    !REAL(8), ALLOCATABLE, DIMENSION(:,:) :: gRDipole
    REAL(8), ALLOCATABLE, DIMENSION(:)   :: gNRB
    REAL(8), ALLOCATABLE, DIMENSION(:)   :: gNdipRB
    REAL(8), ALLOCATABLE, DIMENSION(:)   :: gShVecRB
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: gShVec
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: gHRB
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: gTorRB   
END MODULE Globals




