
Subroutine CalculateFT
  USE Globals
  USE VectorFunctions
  IMPLICIT NONE

  INTEGER :: i, j
  REAL(8) :: RijM, Rij2, Rij3, Rij4, RijC, RdipIJ2, RdipIJM, RdipIJ3, RdipIJ4, OrientFunc
  REAL(8), DIMENSION(gND) :: Rij
  REAL(8), DIMENSION(gND) :: RdipIJ
  REAL(8), DIMENSION(gND) :: t
  REAL(8), DIMENSION(gND) :: tdip
  REAL(8), DIMENSION(gND) :: Nj
  REAL(8), DIMENSION(gND) :: Ni
  REAL(8), DIMENSION(gND) :: Ndipj  !Ronal
  REAL(8), DIMENSION(gND) :: Ndipi  !Ronal
  REAL(8), DIMENSION(gND) :: Fcmi, FcmIJ, FcmJI, FdipJI
  REAL(8), DIMENSION(gND) :: Tcmi, TcmIJ, TcmJI
  REAL(8), DIMENSION(gND) :: FdipIJ, FrepIJ
  REAL(8), DIMENSION(gND) :: TdipIJ, TshIJ, TdipJI, TshJI
  REAL(8), DIMENSION(gND) :: NixH
  REAL(8), DIMENSION(gND) :: NixFij
  REAL(8), DIMENSION(gND) :: B1
  REAL(8), DIMENSION(gND) :: B2

  gFor(:,:) = 0.D0
  gTor(:,:) = 0.D0
  Fcmi(:)   = 0.D0
  Tcmi(:)   = 0.D0
  FcmIJ(:)  = 0.D0
  TcmIJ(:)  = 0.D0
  TcmJI(:)  = 0.D0
  FcmJI(:)  = 0.D0
  FdipIJ(:) = 0.D0
  FdipJI(:) = 0.D0
  FrepIJ(:) = 0.D0
  TdipIJ(:) = 0.D0
  TshIJ(:)  = 0.D0
  TdipJI(:) = 0.D0
  TshJI(:)  = 0.D0

  DO i = 1, gNPart - 1

     Ndipi(:)= gNdip(i,:)
     Fcmi(:) = gFor(i,:)
     Tcmi(:) = gTor(i,:)

     DO j = i + 1, gNpart

        FcmIJ(:)  = 0.D0
        TcmIJ(:)  = 0.D0
        FdipIJ(:) = 0.D0
        TdipIJ(:) = 0.D0
        FrepIJ(:) = 0.D0
        TshIJ(:)  = 0.D0
        TcmJI(:)  = 0.D0
        FcmJI(:)  = 0.D0
        FdipJI(:) = 0.D0
        TshJI(:)  = 0.D0
        TcmJI(:)  = 0.D0
        FcmJI(:)  = 0.D0
        FdipJI(:) = 0.D0
        TshJI(:)  = 0.D0

        CALL distdip(i,j,Rij,RIJ2,RdipIJ,RdipIJ2) !!************************* !vfun.f90

        IF ( RIJ2 .GE. gCutOff2) CYCLE

        !*****************************************!
        !      Interactions on Virtual Site       !
        !*****************************************!
        RdipIJM = DSQRT(RdipIJ2)
        tdip(:) = RdipIJ(:) / RdipIJM
        Ndipj(:)= gNdip(j,:)

        FdipIJ(:) = MagneticForce(RdipIJM,tdip,Ndipi,Ndipj)
        TdipIJ(:) = MagneticTorque(RdipIJM,tdip,Ndipi,Ndipj)
        TshIJ(:)  = vcross( gShVec(i,:), FdipIJ )

        !*****************************************!
        !      Interactions on Center of Mass     !
        !*****************************************!
        RIJM = DSQRT(RIJ2)
        t(:)    = RIJ(:) / RIJM
        IF ( RIJ2 .LE. gCutOffWCA2 ) FrepIJ(:) = WCAForce(Rij2,t)

        !******************************************!
        ! Transform interactions to Center of Mass !
        !******************************************!
        FcmIJ(:) = FdipIJ(:) + FrepIJ(:)
        TcmIJ(:) = TdipIJ(:) + TshIJ(:)

        !**********************************************************!
        !  Calculate Interactions on particle i due to particle j  !
        !**********************************************************!
        Fcmi(:)  = Fcmi(:) + FcmIJ(:)
        Tcmi(:)  = Tcmi(:) + TcmIJ(:)

        !*****************************************************!
        ! Recalculate Torques on Particle j due to particle i !
        !*****************************************************!
        FdipJI(:) = - FdipIJ(:)
      ! TdipJI(:) =   MagneticTorque(RdipIJM,tdip,Nj,Ni)       !Ronal
        TdipJI(:) =   MagneticTorque(RdipIJM,tdip,Ndipj,Ndipi) !Ronal
        TshJI(:)  =   vcross( gShVec(j,:), FdipJI(:))
        FcmJI(:)  = - FcmIJ(:)  !!!!changed from - FcmJI to - FcmIJ  positive --RONAL---
        TcmJI(:)  =   TdipJI(:) + TshJI(:)

        gFor(j,:) = gFor(j,:) + FcmJI(:)
        gTor(j,:) = gTor(j,:) + TcmJI(:)

     END DO

     gFor(i,:) = Fcmi(:)
     gTor(i,:) = Tcmi(:)

  END DO

  !***************************************************!
  !  Add External Field Contribution to Total Torque  !
  !***************************************************!
  DO i = 1, gNPart
   !  NixH      = vcross(gN(i,:),gH(:))  !(:) for gH
     NixH      = vcross(gNdip(i,:),gH(:))   !Ronal
     gTor(i,:) = gTor(i,:) + gLag * NixH(:)
  END DO

END SUBROUTINE CalculateFT


SUBROUTINE UpdateParticles
  USE Globals
  USE VectorClass
  USE VectorFunctions

  IMPLICIT NONE
  INTEGER :: i
  REAL(8) magu, magn
  REAL(8), DIMENSION(gNPart,gND) :: Xi
  REAL(8), DIMENSION(gNPart,gND) :: TxN
  REAL(8), DIMENSION(gNPart,gND) :: XixN
  REAL(8), DIMENSION(gNPart,gND) :: Omega

  CALL RandomVectors

  DO i = 1, gNPart
     Xi(i,1) = gBVec(i) % x
     Xi(i,2) = gBVec(i) % y
     Xi(i,3) = gBVec(i) % z

     IF ( gNDimensionality .EQ. 2 ) THEN !! Added to define dimensionality on August 15, 2015 (q2D)
        Xi(i,3) = 0.0
     END IF

     !gR(i,:) = gR(i,:) + gSQRT2DT * Xi(i,:) !Luis había comentado la parte de BM (last term)
     gR(i,:) = gR(i,:) + gPeDT * gN(i,:) + gDT * gFor(i,:) + gSQRT2DT * Xi(i,:) !Luis había comentado la parte de BM (last term)

     IF ( gNDimensionality .EQ. 2 ) THEN !! Added to define dimensionality on August 15, 2015 (q2D)
	gR(i,3) = 0.0
     END IF
                                                                                !orientation of the self propulsion in the gN direction
     CALL ContainParticle( i )
    ! print*, '  gPeDT=',   gPeDT
  END DO

  CALL UpdateQuat               !quat.f90
  CALL UpdateOrientationQuat    !quat.f90
  !CALL DipolePosition

END SUBROUTINE UpdateParticles



 !!****m* BrownianDynamics/ContainParticle
 !!
 !! NAME
 !!  ContainParticle
 !!
 !! SYNOPSIS
 !!  This SUBROUTINE is a wrapper for the FUNCTION Contain and makes
 !!  sure that a specified particle is inside the simulation cell.
 !!
 !! USAGE
 !!  SUBROUTINE ContainParticle( i )
 !!
 !! INPUTS
 !!  i - the ID of a particle to be checked for containment within the
 !!  simulation cell.
 !!
 !! NOTES
 !!  Located in file dyn.f90
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
 !SUBROUTINE ContainParticle( i )
 !  USE Globals
 !  IMPLICIT NONE
 !
 !  Integer :: i
 !  type ( Vector ) :: Contain
 !
 !  gPart( i ) % pos = Contain( gPart( i ) % pos )

 !END SUBROUTINE ContainParticle

 SUBROUTINE ContainParticle( i )
   USE Globals
   USE VectorFunctions
   IMPLICIT NONE

   INTEGER :: i

   gR(i,:) = Contain( gR(i,:) )

 END SUBROUTINE ContainParticle

 !!****f* BrownianDynamics/Contain
 !!
 !! NAME
 !!  Contain
 !!
 !! SYNOPSIS
 !!  This FUNCTION checks that a position vector lies within the
 !!  simulation cell.  If it isn't, it implements a periodic boundary
 !!  condition and returns a new vector corresponding to the periodic
 !!  image which is inside the box.
 !!
 !! USAGE
 !!  type ( Vector ) FUNCTION Contain( v )
 !!
 !! INPUTS
 !!  v - a vector to be checked for containment within the simulation cell.
 !!
 !! RESULT
 !!  The periodic image of the input vector that lies with the simulation cell.
 !!
 !! NOTES
 !!  Located in file dyn.f90
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
 type ( Vector ) FUNCTION Contain( v )
   USE Globals
   IMPLICIT NONE

   type ( Vector ) :: v

   if ( v % x .lt. 0 ) then
      Contain % x = v % x + gCDim( 1 )
   else if ( v % x .gt. gCDim( 1 ) ) then
      Contain % x = v % x - gCDim( 1 )
   else
      Contain % x = v % x
   END if

   if ( v % y .lt. 0 ) then
      Contain % y = v % y + gCDim( 2 )
   else if ( v % y .gt. gCDim( 2 ) ) then
      Contain % y = v % y - gCDim( 2 )
   else
      Contain % y = v % y
   END if

   if ( v % z .lt. 0 ) then
      Contain % z = v % z + gCDim( 3 )
   else if ( v % z .gt. gCDim( 3 ) ) then
      Contain % z = v % z - gCDim( 3 )
   else
      Contain % z = v % z
   END if

 END FUNCTION Contain


! REAL(8) DIMENSION(gND) FUNCTION Contain( v )
!   USE Globals
!   IMPLICIT NONE
!
!   INTEGER :: i
!   REAL(8), DIMENSION(gND) :: v
!
!   DO i = 1, gND
!      IF ( v(i) .LT. 0.0 ) THEN
!         Contain(i) = v(i) + gCDim(i)
!      ELSE IF ( v(i) .GT. gCDim(i) ) THEN
!         Contain(i) = v(i) - gCDim(i)
!      ELSE
!         Contain(i) = v(i)
!      END IF
!   END DO
!
! END FUNCTION Contain




SUBROUTINE RandomizePos
  USE Globals
  USE VectorFunctions
  IMPLICIT NONE
  REAL(8), DIMENSION(gNPart,gND) :: Xi
  INTEGER :: i, j, k

  WRITE(*,*) "--------------------------------------------------------"
  WRITE(*,*) "----------------- BEGIN RANDOMISATION ------------------"
  WRITE(*,*) "--------------------------------------------------------"

  K = 0

  DO WHILE ( gT .LE. gTF )

     CALL RandomVectors

     DO i = 1, gNPart
        Xi(i,1) = gBVec(i) % x
        Xi(i,2) = gBVec(i) % y
        Xi(i,3) = gBVec(i) % z
        gR(i,:) = gR(i,:) + gSqrt2DT * Xi(i,:)
        CALL ContainParticle( i )
     END DO

     CALL OverlapCheck

     WRITE(*,*) REAL(K) / REAL(gNDTRAND) * 100.D0, "%"

     IF ( K .GE. gNDTRAND ) THEN
        gT = 0.D0
        EXIT
     END IF

     gT = gT + gDT
     K  = K  + 1

  END DO

  WRITE(*,*) "--------------------------------------------------------"
  WRITE(*,*) "--------------RANDOMISATION IS NOW COMPLETE-------------"
  WRITE(*,*) "--------------------------------------------------------"
  WRITE(*,*) "----------SETTING PARAMETERS TO ORIGINAL VALUES---------"
  WRITE(*,*) "--------------------------------------------------------"
  WRITE(*,*) "----------------SIMULATION WILL NOW BEGIN---------------"
  WRITE(*,*) "--------------------------------------------------------"

END SUBROUTINE RandomizePos



! SUBROUTINE RandomizePos
!   uSE GLOBALS
!   IMPLICIT NONE
!   Integer :: I, J, K
!   REAL(8) :: LambTemp4, LambTemp12, LangTemp, DtTemp, StericTVTemp
!
!   WRITE(*,*) "--------------------------------------------------------"
!   WRITE(*,*) "----------------- BEGIN RANDOMISATION ------------------"
!   WRITE(*,*) "--------------------------------------------------------"
!
!   K = 0
!
!   LambTemp4    = gLambda4
!   LambTemp12   = gLambda12
!   LangTemp     = gLag
!   DtTemp       = gDT
!   StericTVTemp = gStericTVInv
!
!   gLambda4     = 0.D0
!   gLambda12    = 0.D0
!   gLag         = 0.D0
!   gDT          = 0.001D0
!   gStericTVInv = 50.D0 / gTV
!
!   DO WHILE ( gT .LE. gTF )
!
!      CALL UPDATEPARTICLES
!      !CALL ORDERPARAMETER
!
!      WRITE(*,*) REAL(K) / REAL(gNDTRAND) * 100.D0, "%"
!
!      IF ( K .GE. gNDTRAND ) THEN
!
!         gLambda4     = LambTemp4
!         gLambda12    = LambTemp12
!         gLag         = LangTemp
!         gDT          = DtTemp
!         gT           = 0.D0
!         gStericTVInv = StericTVTemp
!
!         WRITE(*,*) "--------------------------------------------------------"
!         WRITE(*,*) "--------------RANDOMISATION IS NOW COMPLETE-------------"
!         WRITE(*,*) "--------------------------------------------------------"
!         WRITE(*,*) "----------SETTING PARAMETERS TO ORIGINAL VALUES---------"
!         WRITE(*,*) "--------------------------------------------------------"
!         WRITE(*,*) "----------------SIMULATION WILL NOW BEGIN---------------"
!         WRITE(*,*) "--------------------------------------------------------"
!         EXIT
!
!      END IF
!
!      gT = gT + gDT
!      K  = K  + 1
!
!   END DO
!
! END SUBROUTINE RandomizePos


!REAL(8) DIMENSION(3) FUNCTION vcross(u,v)
!  USE Globals
!  IMPLICIT NONE
!  INTEGER :: i, j
!  REAL(8), DIMENSION(3) :: u
!  REAL(8), DIMENSION(3) :: v
!
!  vcross(:) = 0.D0
!
!  DO i =1, gND
!    DO j=1, gND
!         vcross(:) = vcross(:) + u(i) * v(j) * eps3(i,j,:)
!    END DO
! END DO
!
! RETURN
!
!END FUNCTION
