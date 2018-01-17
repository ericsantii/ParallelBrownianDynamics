!!****************************************************************************!!
!!****m* BrownianDynamics/InitQuat
!!
!!  NAME
!!   
!!
!!  SYNOPSIS
!!
!!  NOTES
!!   Located in file quat.f90
!!
!!  AUTHOR
!!   Luis Y. Rivera-Rivera
!!
!! COPYRIGHT
!!    4-March-2014
!!*****
!! UPDATE
!!  
!!******        
!!
!!****************************************************************************!!
SUBROUTINE InitQuat
  USE GLOBALS
  USE VectorFunctions
  IMPLICIT NONE 
  INTEGER :: I
  REAL(8) :: Mag, x, y ,z
  Real(8), dimension( 3, 3 ) :: A
  Real(8), dimension( 3, 3 ) :: AT

  call random_seed()

  DO I = 1, gNPart
 
    call random_number(x) 
    gQuat(I,1) = 2.D0*x - 1.D0
  
    call random_number(y) 
    gQuat(I,3) = 2.D0*y - 1.D0

    call random_number(z)
    gQuat(I,2) = 2.D0*z - 1.D0

    gQuat(I,4) = 1.D0 - gQuat(I,1) ** 2 - gQuat(I,2) ** 2 - gQuat(I,3) ** 2
   
   !debugging  
   !print*, gQuat(I,1), gQuat(I,2),gQuat(I,3),gQuat(I,4), 'qi'
   !pause
  
    Mag = DSQRT( vdot(gQuat(i,:),gQuat(i,:)) )

    gQuat(I,:) = gQuat(I,:) / Mag
   ! debbuging  -- the squared sum is 1 --it is correct! --Ronal
     !print*, gQuat(I,1), gQuat(I,2),gQuat(I,3),gQuat(I,4), 'new qi'
     !Mag = DSQRT( vdot(gQuat(i,:),gQuat(i,:)) )
     !print*, Mag, 'Mag'
     !pause
   
    CALL MATRICES(3, gQuat(i,1), gQuat(i,2), gQuat(i,3), gQuat(i,4), A, AT)
    
    gN(i,:)     = TRANSFORM(AT,gNRB) !
    gNdip(i,:)  = TRANSFORM(AT,gNdipRB) ! Ronal
    gShVec(i,:) = TRANSFORM(AT,gShVecRB) !
  END DO

END SUBROUTINE InitQuat


!!****************************************************************************!!
!!****m* BrownianDynamics/UpdateOrientationQuat
!!
!!  NAME
!!   UpdateOrientationQuat
!!
!!  SYNOPSIS
!!    The orientations of the particles is updated in the Lab Frame
!!    using the recently updated quaternion parameters
!!     
!!  NOTES
!!   Located in file quat.f90
!!
!!  AUTHOR
!!   Luis Y. Rivera-Rivera
!!
!! COPYRIGHT
!!    5-March-2014
!!*****
!! UPDATE
!!  
!!******        
!!
!!****************************************************************************!!
SUBROUTINE UpdateOrientationQuat
  USE GLOBALS
  USE VectorFunctions
  IMPLICIT NONE
  INTEGER :: i
  Real(8) :: NormN, NormDip  !Ronal
  Real( 8 ), dimension( 3, 3 ) :: A
  Real( 8 ), dimension( 3, 3 ) :: AT

  DO i = 1, gNPart

     call Matrices ( 3, gQuat(I, 1), gQuat(I,2), gQuat(I,3), gQuat(I,4), A, AT )
     
     !**************************************************!
     !  Calculating the new orientations and Shifting   !
     !   vector of the particles based on the new       !
     !               quaternions parameters             !
     !**************************************************!

     gN(i,:) = TRANSFORM(AT,gNRB)
     gNdip(i,:)= TRANSFORM(AT,gNdipRB)  !Ronal
     NormN    = DSQRT( vdot(gN(i,:), gN(i,:)))  !Ronal
     NormDip    = DSQRT( vdot(gNdip(i,:), gNdip(i,:)))  !Ronal
     
     gN(i,:) = gN(i,:) / NormN  !Ronal
     gNdip(i,:) = gNdip(i,:) / NormDip  !Ronal
     
     gShVec(i,:) = TRANSFORM(AT,gShVecRB)

  END DO
END SUBROUTINE UpdateOrientationQuat




SUBROUTINE UpdateQuat
  USE GLOBALS
  USE VectorFunctions
  IMPLICIT NONE
  INTEGER :: I
  REAL(8) :: MagQuat
  REAL(8), DIMENSION(gNPart,gND) :: Xi
  REAL(8), DIMENSION(3) :: DeltaOmega
  REAL(8), DIMENSION(4) :: DeltaQuat
  Real(8), dimension(3,3) :: A
  Real(8), dimension(3,3) :: AT

  CALL RandomVectors
  DO I = 1, gNPart

     call Matrices ( 3, gQuat(I, 1), gQuat(I,2), gQuat(I,3), gQuat(I,4), A, AT ) 

     gTorRB(i,:) = TRANSFORM(A,gTor(i,:))
  
     	Xi(i:,1)  = gBVec(i) % x
 	   Xi(i:,2)  = gBVec(i) % y
      Xi(i:,3)  = gBVec(i) % z


     
     DeltaOmega(:) = gDT34 * gTorRB(i,:) + gBrowDT * Xi(i,:) !no veo contribucion Browniana -acabo de implementarla-- it is correct?
   
  
     DeltaQuat(1) = 0.5D0 *( - gQuat( i , 2 ) * DeltaOmega(1) - gQuat( i , 3 ) * DeltaOmega(2) - gQuat( i , 4 ) * DeltaOmega(3) )
     DeltaQuat(2) = 0.5D0 *(   gQuat( i , 1 ) * DeltaOmega(1) - gQuat( i , 4 ) * DeltaOmega(2) + gQuat( i , 3 ) * DeltaOmega(3) )
     DeltaQuat(3) = 0.5D0 *(   gQuat( i , 4 ) * DeltaOmega(1) + gQuat( i , 1 ) * DeltaOmega(2) - gQuat( i , 2 ) * DeltaOmega(3) )
     DeltaQuat(4) = 0.5D0 *( - gQuat( i , 3 ) * DeltaOmega(1) + gQuat( i , 2 ) * DeltaOmega(2) + gQuat( i , 1 ) * DeltaOmega(3) )
     
     gQuat(i,:) = gQuat(i,:) + DeltaQuat(:)     
     MagQuat    = DSQRT( vdot(gQuat(i,:),gQuat(i,:)) )
     gQuat(i,:) = gQuat(i,:) / MagQuat

  END DO


END SUBROUTINE UpdateQuat



!!****************************************************************************!!
!!****m* BrownianDynamics/TransformHN
!!
!!  NAME
!!   TransformHN
!!
!!  SYNOPSIS
!!   Transform H -> H' :  H' = A.H
!!   Transform N'-> N  :  N = AT.N'
!!  NOTES
!!   Located in file quat.f90
!!
!!  AUTHOR
!!   Luis Y. Rivera-Rivera
!!
!! COPYRIGHT
!!    4-March-2014
!!*****
!! UPDATE
!!  
!!******        
!!
!!****************************************************************************!!
!SUBROUTINE TransformHN
!USE GLOBALS
!USE VECTORCLASS

!IMPLICIT NONE
!INTEGER :: I
!Real( 8 ), dimension( 3, 3 ) :: A
!Real( 8 ), dimension( 3, 3 ) :: AT

!DO I = 1, gNPart
!
!CALL MATRICES(3, gQuat(I,1), gQuat(I,2), gQuat(I,3), gQuat(I,4), A, AT)

!gHRB(I) % x = A(1,3)
!gHRB(I) % y = A(2,3)
!gHRB(I) % z = A(3,3)
!
!gN(I) % x = A(3,1)
!gN(I) % y = A(3,2)
!gN(I) % z = A(3,3)
!
!gU(I) % x = AT(1,2) * gSinPsi + AT(1,3) * gCosPsi
!gU(I) % y = AT(2,2) * gSinPsi + AT(2,3) * gCosPsi
!gU(I) % z = AT(3,2) * gSinPsi + AT(3,3) * gCosPsi
!
!orient(I) % x = gN(I) % x
!orient(I) % y = gN(I) % y
!orient(I) % z = gN(I) % z
!
!END DO
!
!END SUBROUTINE TransformHN

