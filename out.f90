!!****m* BrownianDynamics/OutputInit
!!
!! NAME
!!  OutputInit
!!
!! SYNOPSIS
!!  This SUBROUTINE echos all of the USEr specIFied data and inFORMATion
!!  about the simulation cell calculated during the initialization steps in
!!  a FORMATted output file.
!!
!! NOTES
!!  Located in file out.f90
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
SUBROUTINE OutputInit
  USE Globals
  IMPLICIT NONE
  
  Integer :: i, iu1, iu2 ! The following loop WRITE initial inFORMATion into output files
  
  !IF ( gTracStat .eq. 1 ) THEN
  !END IF
  
  WRITE( 102, 5001 ) 'Simulation Summary' 
  
  WRITE( 102, 5002 ) 'Particle data'
  WRITE( 102, 501 )
  
  WRITE( 102, * ) 'Total number of particles:', gNPart
  WRITE( 102, * ) 'Peclet Number:' , gPe
  
  IF ( gND .eq. 2 ) THEN
     WRITE( 102, * ) 'This simulation is in 2D'
     WRITE( 102, * ) 'Surface fraction of particles:', gPi * gNPart / ( gCDim( 1 ) * gCDim( 2 ) )
     WRITE( 102, * ) 'Area density:', gNPart / ( gCDim( 1 ) * gCDim( 2 ) )
  else
     WRITE( 102, * ) 'Volume fraction of particles:', gPhi
     !  WRITE( 102, * ) 'Volume fraction of the motor particle:', gPhim
     !  WRITE( 102, * ) 'Volume fraction of the bath particle:', gPhib
     WRITE( 102, * ) 'Calculate new volume fraction of the system:', gPhis
     WRITE( 102, * ) '********************************************'
     !  WRITE( 102, * ) 'Correct value of volume fraction of the motor:', gPhinm
     !  WRITE( 102, * ) 'Correct value of volume fraction of the bath:', gPhinb
     WRITE( 102, * ) '********************************************'
     WRITE( 102, * ) 'Langevin parameter:', gLag
     WRITE( 102, * ) 'Dipolar coupling constant:', gLambda
     WRITE( 102, * ) 'Num. density:', gNPart / ( gCDim( 1 ) * gCDim( 2 ) * gCDim( 3 ) )
  END IF
  
  WRITE( 102, 5004 ) 'Rotation factor:', gRotateFactor
  WRITE( 102, 5003 ) 'Simulation mode:', gForceMode
  WRITE( 102, 5003 ) 'HS mode:', hsMode 
  WRITE( 102, 5003 ) 'Bath particles cols?:', gNoBathColFlag
  WRITE( 102, 501 )
  
  !DO i = 1, gNMove
  !WRITE( 102, * ) 'Imposed Ps:', rxnp( i )
  !WRITE( 102, * ) 'Cosine of react patch:', gCosThetaRxn( i )
  !END DO
  
  WRITE( 102, 501 )
  !WRITE( 102, 5005 ) 'Number of moving particles: ', gNMove
  WRITE( 102, 501 )
  
  !WRITE( 102, 5006 ) 'ID', 'Radius', 'Peclet'
  !DO i = 1, gNMove
  !    WRITE( 102, 5007 ) i, gPart( i ) % rat, gPart( i ) % pec
  !END DO
  WRITE( 102, 500 )
  
  WRITE( 102, 5008 ) 'Simulation timing data'
  
  WRITE( 102, * ) 'Total time: ', gTf
  WRITE( 102, * ) 'Time step: ', gDT
  WRITE( 102, 500 )
  
  ! Include the timing of statistics measurements and outputs
  
  WRITE( 102, 5011 ) 'Simulation cell data'
  WRITE( 102, 501 )
  
  WRITE( 102, 5012 ) 'Cell dimensions:'
  WRITE( 102, 5013 ) 'X:', gCDim( 1 )
  WRITE( 102, 5013 ) 'Y:', gCDim( 2 )
  WRITE( 102, 5013 ) 'Z:', gCDim( 3 )
  WRITE( 102, 501 ) 
  
  WRITE( 102, 5012 ) 'Number of cubic cells:'
  WRITE( 102, 5014 ) 'X:', gNCells( 1 )
  WRITE( 102, 5014 ) 'Y:', gNCells( 2 )
  WRITE( 102, 5014 ) 'Z:', gNCells( 3 )
  WRITE( 102, 501 )
  
  WRITE( 102, 5012 ) 'Number of linked list cells:'
  WRITE( 102, 5014 ) 'X:', gNLCells( 1 )
  WRITE( 102, 5014 ) 'Y:', gNLCells( 2 )
  WRITE( 102, 5014 ) 'Z:', gNLCells( 3 )
  WRITE( 102, 501 )
  
  WRITE( 102, 5012 ) 'Number of density profile bins:'
  WRITE( 102, 5014 ) 'X:', gDPBins( 1 )
  WRITE( 102, 5014 ) 'Y:', gDPBins( 2 )
  WRITE( 102, 5014 ) 'Z:', gDPBins( 3 )
  WRITE( 102, 501 )
  
  WRITE( 102, 5012 ) 'Number of pair dist. bins:'
  WRITE( 102, 5014 ) 'X:', gPDBins( 1 )
  WRITE( 102, 5014 ) 'Y:', gPDBins( 2 )
  WRITE( 102, 5014 ) 'Z:', gPDBins( 3 )
  WRITE( 102, 500 )
  
500 FORMAT( 1x, '|', t6, 41('-'), t51, '|', /, 1x, '|', t51, '|' )
501 FORMAT( 1x, '|', t6, '|', t46, '|', t51, '|' )
5001 FORMAT( 1x, 14('-'), 2x, A, 2x, 14('-'), /, 1x, '|', t51, '|' )
5002 FORMAT( 1x, '|', t6, 12('-'), 2x ,A, 2x, 12('-'), t51, '|' )
5003 FORMAT( 1x, '|', t6, '|', 2x, A, 4x, I5, t46, '|', t51, '|' )
5004 FORMAT( 1x, '|', t6, '|', 2x, A, F6.3, t46, '|', t51, '|' )
6004 FORMAT( 1x, '|', t6, '|', 2x, A, F10.3, t46, '|', t51, '|' )
5005 FORMAT( 1x, '|', t6, '|', 2x, A, 4x, I3, t46, '|', t51, '|' )
5006 FORMAT( 1x, '|', t6, '|', 4x, A, 2x ,'---', 2x, A, 2x, '---', 2x, A, &
          t46, '|', t51, '|' )
5007 FORMAT( 1x, '|', t6, '|', 3x, I3, 7x, F6.1, 7x, F6.1, t46, '|', t51, '|' )
5008 FORMAT( 1x, '|', t6, 8('-'), 2x, A, 2x, 8('-'), t51, '|', &
          /, 1x, '|', t6, '|', t46, '|', t51, '|' )
5009 FORMAT( 1x, '|', t6, '|', 2x, A, F6.2, t46, '|', t51, '|' )
5010 FORMAT( 1x, '|', t6, '|', 2x, A, F7.3, t46, '|', t51, '|' )
5011 FORMAT( 1x, '|', t6, 8('-'), 2x, A, 2x, 9('-'), t51, '|' )
5012 FORMAT( 1x, '|', t6, '|', 2x, A, t46, '|', t51, '|' )
5013 FORMAT( 1x, '|', t6, '|', 4x, A, F9.3, t46, '|', t51, '|' )  
5014 FORMAT( 1x, '|', t6, '|', 4x, A, I4, t46, '|', t51, '|' )  
  
END SUBROUTINE OutputInit


!!****m* BrownianDynamics/OutputInter
!!
!! NAME
!!  OutputInter
!!
!! SYNOPSIS
!!  This SUBROUTINE prints out verious data while the simulation is running.
!!
!! NOTES
!!  Located in file out.f90
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

SUBROUTINE OutputInter
  USE Globals
  USE VectorFunctions
  IMPLICIT NONE
  
  INTEGER :: i 
  REAL(8), DIMENSION(gND) :: disp
  Real(8) :: tScale2
  
  tScale2 = 1.d0 / ( ( gT - gAvgTiming ) * gDTInv )
  
  IF ( mod( int( gT / gDT ), 100 ) .eq. 0 ) THEN
     WRITE( *, 300 ) 100 * gT / gTf, ' % complete'
  END IF
  
  IF ( mod( gOutCount, gOutTiming ) .eq. 0 ) THEN
     DO i = 1, gNPart
        WRITE( 1400, 1203 ) gT, Contain( gR(i,:) ), gN(i,:), gNDip(i,:), gQuat(i,:)
         ! 1203 FORMAT(3(1X, E24.15E3) )
         1203 FORMAT(1X, E15.7E3, 3(1X, E24.15E3), 6(1X, E24.15E3),4(1X, E24.15E3) )  
     END DO
  END IF
  gOutCount = gOutCount + 1
  
  IF ( gMovies .eq. 1 ) THEN
     IF ( MOD( gFramesCounter, gOutFrames ) .EQ. 0 ) THEN
        DO i = 1, gNMove
           WRITE( 103, * ) Contain( gR(i,:) )
           IF ( gSpecIFiedFluxMode .eq. 0 ) THEN
              WRITE( 103, * ) gN(i,:)
           END IF
        END DO
        WRITE( 103, * ) '--------'
        DO i = gNMove + 1, gNPart
           IF ( gSpecIFiedFluxMode .eq. 1 ) THEN
              WRITE( 104, 1219 ) Contain( gR(i,:) ) 
1219          FORMAT( 3(1X, F6.3), I3 )
           else
              WRITE( 104, 5555 ) Contain( gR(i,:) )  , gN(i,:)
5555          FORMAT(3(1X, E24.15E3), 3(1X, E24.15E3))
           END IF
        END DO
        WRITE( 104, * ) '------'
     END IF
  END IF
  gFramesCounter = gFramesCounter + 1
  
300 FORMAT( 1x, F6.3, A )
111 FORMAT( 1x, 5f30.16, 5f30.16, 5f30.16, 5f30.16 )
END SUBROUTINE OutputInter




!!****m* BrownianDynamics/OutputRestartData
!!
!! NAME
!!  OutputRestartData
!!
!! SYNOPSIS
!!  This SUBROUTINE save data for restarting when the code is stopped.
!!  The code restarts using this data and continues to the final.
!!
!! NOTES
!!  Located in file out.f90
!!
!! CREATION DATE
!!    04/22/2013
!!
!! AUTHOR
!!    Ronal De La Cruz
!!
!! COPYRIGHT
!!    22-Apr-2013 - Ronal De La Cruz
!!******
!SUBROUTINE OutputRestartData
!  USE Globals
!  IMPLICIT NONE
!
!  Integer :: i, j
!
!  type ( Vector ) :: disp, Contain
!
!  IF ( mod( gRestCount, gRestTiming ) .eq. 0 ) THEN
!
!      DO i = 1, gNPart
!        WRITE( 1600, 3001 ) orient(i)
!        3001 FORMAT(3(1X, E24.9E3))
!      END DO
!
!     WRITE(1601,*) gT
!
! END IF
!
!    gRestCount = gRestCount + 1
!
!END SUBROUTINE OutputRestartData



!!****m* BrownianDynamics/OutputFinal
!!
!! NAME
!!  OutputFinal
!!
!! SYNOPSIS
!!  This SUBROUTINE wraps up the program by writing all the
!!  final statistical inFORMATion to the appropriate files.
!!
!! NOTES
!!  Located in file out.f90
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
!SUBROUTINE OutputFinal
!  USE Globals
!  IMPLICIT NONE
!  
!  Real( 8 ) :: csOP, avg, dp, dp2, geq, inveffvis ! Average properties.
!  Integer :: i, j, k, l, thick, thick2 ! Counters.
!  type ( Vector ) :: dvecs, MinImage
!  
!  
!  WRITE( 102, 5015 ) 'Simulation results'
!  WRITE( 102, 501 )
!  
!  geq = ( 1. - 0.5 * gPhi ) / ( 1. - gPhi ) ** 3.
!  
!  csOP = 1. + 4. * gPhi  * geq
!  
!  inveffvis = 1. / ( 1. + 2. * gPhi * geq )
!  
!  WRITE( 102, * ) 'Avg. osmotic pressure:', gAP 
!  WRITE( 102, 5016 ) 'C-S pressure:', csOP
!  
!  DO i = 1, gNMove
!     !        WRITE( 1623, * ) gAV( i ) % x
!     !        WRITE( 1623, * ) gAV( i ) % y
!     !        WRITE( 1623, * ) gAV( i ) % z
!     WRITE( 102, 501 )
!     WRITE( 102, 5017 ) 'Moving particle:', i
!     WRITE( 102, * ) 'Avg. x-velo:', gAV( i ) % x
!     WRITE( 102, * ) 'Avg. y-velo:', gAV( i ) % y
!     WRITE( 102, * ) 'Avg. z-velo:', gAV( i ) % z
!     WRITE( 102, 5018 ) 'Avg. speed:', sqrt( gAV( i ) * gAV( i ) )
!     WRITE( 102, * ) 'HS x-force:', ( gHSF( i ) % x ) / gDT
!     WRITE( 102, * ) 'HS y-force:', ( gHSF( i ) % y ) / gDT
!     WRITE( 102, * ) 'HS z-force:', ( gHSF( i ) % z ) / gDT
!     
!     IF ( gForceMode .eq. 2 ) THEN
!        WRITE( 102, * ) 'HS x-force rxn side:', ( gHSFRh( i ) % x ) / gDT
!        WRITE( 102, * ) 'HS y-force rxn side:', ( gHSFRh( i ) % y ) / gDT
!        WRITE( 102, * ) 'HS z-force rxn side:', ( gHSFRh( i ) % z ) / gDT
!        WRITE( 102, * ) 'HS x-force inert side:', ( gHSFnRh( i ) % x ) / gDT
!        WRITE( 102, * ) 'HS y-force inert side:', ( gHSFnRh( i ) % y ) / gDT
!        WRITE( 102, * ) 'HS z-force inert side:', ( gHSFnRh( i ) % z ) / gDT
!        
!        WRITE( 102, * ) 'Avg. stress on motor:', gOPMotor( i ) / gDT
!        
!     END IF
!     
!     IF ( gNMove .eq. 2 ) THEN
!        WRITE( 102, * ) 'Avg x-torque:', ( gHST( i ) % x ) / gDT
!        WRITE( 102, * ) 'Avg y-torque:', ( gHST( i ) % y ) / gDT
!        WRITE( 102, * ) 'Avg z-torque:', ( gHST( i ) % z ) / gDT
!        dvecs = MinImage( dvec( 1 ) - dvec( 2 ) )
!        WRITE( 102, * ) 'Avg z-torque per dist:', ( gHST( i ) % z ) / gDT / dsqrt( dvecs * dvecs )
!     END IF
!     
!     IF ( gFreeRotor .eq. 1 ) THEN
!        WRITE( 102, * ) 'Angular x-velo:', gAAngVel( i ) % x
!        WRITE( 102, * ) 'Angular y-velo:', gAAngVel( i ) % y
!        WRITE( 102, * ) 'Angular z-velo:', gAAngVel( i ) % z
!        WRITE( 102, * ) 'Tensile x-force:', ( gTensileF( i ) % x ) / gDT
!        WRITE( 102, * ) 'Tensile y-force:', ( gTensileF( i ) % y ) / gDT
!        WRITE( 102, * ) 'Tensile z-force:', ( gTensileF( i ) % z ) / gDT
!        WRITE( 102, * ) 'CM x-force:', ( gFCM % x ) / gDT
!        WRITE( 102, * ) 'CM y-force:', ( gFCM % y ) / gDT
!        WRITE( 102, * ) 'CM z-force:', ( gFCM % z ) / gDT
!     END IF
!     
!  END DO
!  
!  WRITE( 102, 5017 ) 'BEFORE HS SCHEME'
!  WRITE( 102, * ) 'Cols NON-RX side', gCountColNoReactive
!  WRITE( 102, * ) 'Cols RX side', gCountCol
!  WRITE( 102, * ) 'Total num coll', gCountColNoReactive + gCountCol
!  WRITE( 102, 5017 ) 'Reaction Statistics'
!  WRITE( 102, * ) 'Consumed particles', gCountColReact
!  WRITE( 102, * ) 'Reaction probability', dble(gCountColReact) / dble(gCountCol)
!  WRITE( 102, 5017) 'AFTER HS SCHEME'
!  WRITE( 102, * ) 'Cols NON-RX side 2', gOldNumColNoReactSide
!  WRITE( 102, * ) 'Cols RX side 2', gOldNumColReactSide
!  WRITE( 102, * ) 'Total num cols 2', dble(gOldNumCol)
!  WRITE( 102, * ) 'Fraction:', dble(gCountColNoReactive + gCountCol) / dble(gOldNumCol)
!  WRITE( 102, 5017 ) 'ADDITIONAL DATA'
!  
!  IF ( gSpecIFiedFluxMode .eq. 1 ) THEN
!     WRITE( 102, * ) 'Avg. Flux per time step', AKo( 1 )
!  END IF
!  
!  IF ( gNMove .eq. 2 ) THEN
!     WRITE( 102, * ) 'Distance between motors: ', dsqrt( dvecs * dvecs )
!  END IF
!  
!  
!  l = 1
!  
!  DO i = 1, gNMove
!     DO j = i + 1, gNMove
!        
!        !            WRITE( 105, * ) 'Particle pair: ', i, j
!        
!        avg = 0
!        
!        DO k = 1, gRDBins
!           
!           !                WRITE( 105, * ) k * gMinDim / gRDBins, gRD( l, k ) 
!           
!           avg = k * gMinDim / gRDBins * gRD( l, k ) + avg
!           
!        END DO
!        
!        !            WRITE( 105, * ) 'Average: ', avg
!        !            WRITE( 105, * ) 'Avg / Half Diag: ', avg / gMinDim
!        
!        l = l + 1
!     END DO
!  END DO
!  
!  thick = 2 * int( gPart( i ) % rat * gDPBins( 3 ) / gCDim( 3 ) )
!  ! thick = int( gDPBins( 3 ) / 2.5)
!  
!  !    WRITE( 106, * ) gDPBins( 1 )
!  !    WRITE( 106, * ) gDPBins( 2 )
!  !    WRITE( 1614, * ) gDPBins( 1 )
!  !    WRITE( 1614, * ) gDPBins( 2 )
!  
!  DO i = 1, gDPBins( 1 )
!     DO j = 1, gDPBins( 2 )
!        dp = 0
!        dp2 = 0
!        
!        DO k = -1 * thick, thick
!           dp = dp + gDP( i, j, k + int( gDPBins( 3 ) / 2 ), 1 )
!           dp2 = dp2 + gPD3D( i, j, k + int( gDPBins( 3 ) / 2 ) )
!        END DO
!        
!        !            WRITE( 106, * ) dp / ( 2. * real( thick ) + 1. )
!        !            WRITE( 1614, * ) dp2 / ( 2. * real( thick ) + 1. )
!     END DO
!     
!  END DO
!  
!  
!  thick2 = 2 * int( gPart( i ) % rat * gDPBins( 1 ) / gCDim( 1 ) )
!  ! thick2 = int( gDPBins( 1 ) / 2.5 )
!  !    WRITE( 1613, * ) gDPBins( 2 )
!  !    WRITE( 1613, * ) gDPBins( 3 )
!  !    WRITE( 1615, * ) gDPBins( 1 )
!  !    WRITE( 1615, * ) gDPBins( 2 )
!  
!  DO i = 1, gDPBins( 2 )
!     DO j = 1, gDPBins( 3 )
!        dp = 0
!        dp2 = 0
!        
!        DO k = -1 * thick2, thick2
!           dp = dp + gDP( k + int( gDPBins( 1 ) / 2 ), i, j, 1 )
!           dp2 = dp2 + gPD3D( k + int( gDPBins( 1 ) / 2 ), i, j )
!        END DO
!        
!     END DO
!     
!  END DO
!  
!500 FORMAT( 1x, '|', t6, 41('-'), t51, '|', /, 1x, '|', t51, '|' )
!501 FORMAT( 1x, '|', t6, '|', t46, '|', t51, '|' )
!5015 FORMAT( 1x, '|', t6, 9('-'), 2x, A, 2x, 10('-'), t51, '|' )
!5016 FORMAT( 1x, '|', t6, '|', 2x, A, F7.3, t46, '|', t51, '|' )
!5017 FORMAT( 1x, '|', t6, '|', 2x, A, I3, t46, '|', t51, '|' )
!5018 FORMAT( 1x, '|', t6, '|', 4x, A, F12.3, t46, '|', t51, '|' )
!  
!END SUBROUTINE OutputFinal



SUBROUTINE RESTARTRUN
  USE Globals
  USE VectorClass
  IMPLICIT NONE
  INTEGER :: nlines 
  INTEGER :: nlcut
  INTEGER :: i
  CHARACTER(LEN=60) :: cmd
  CHARACTER(LEN=10) :: nlcut_str
  
  nlines = 0
  
  OPEN ( unit = 1400, file='/Data/periodic_pos.out', status="unknown")
  
  DO
     READ (1400,*, END=10)
     nlines = nlines + 1
  END DO
10 CLOSE(1400)
  
  nlcut = mod(nlines,gNPart)
  nlcut = nlines - nlcut
  
  cmd = 'head -'//trim(ADJUSTL(nlcut_str))//' Data/periodic_pos.out > Data/periodic_pos.out_tmp'
  
  CALL SYSTEM(cmd)
  
  cmd = 'rm Data/periodic_pos.out ; mv Data/periodic_pos.out_tmp Data/periodic_pos.out'
  
  CALL SYSTEM(cmd)
  
  OPEN ( unit = 1400, file='/Data/periodic_pos.out', status="unknown")
  
  DO i = 1, nlines - gNPart
     READ(1400,*) 
  END DO
  
  DO i = 1, gNPart
     READ(1400,*) gPart(i) % pos
  END DO
  
  CLOSE(1400)
  
  OPEN ( unit = 1400, file='Data/periodic_pos.out' , status="unknown", access='appEND')
  
  
END SUBROUTINE RESTARTRUN





