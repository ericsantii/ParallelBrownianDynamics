MODULE VectorFunctions
  IMPLICIT NONE
  SAVE

  CONTAINS

  FUNCTION vcross(u, v)
    REAL(8), DIMENSION(3) :: u, v
    REAL(8), DIMENSION(3) :: vcross
    vcross(1) = u(2) * v(3) - u(3) * v(2)
    vcross(2) = u(3) * v(1) - u(1) * v(3)
    vcross(3) = u(1) * v(2) - u(2) * v(1)
  END FUNCTION vcross


  FUNCTION sqrdist(i,j)
    USE Globals
    INTEGER, INTENT(IN)     :: i, j
    REAL(8), DIMENSION(gND) :: drij
    REAL(8) :: sqrdist

    drij(:) = gR(i,:) - gR(j,:) 
    sqrdist = SUM( drij(:) * drij(:) )

  END FUNCTION sqrdist

  SUBROUTINE distdip(i,j,rij,d2,rdipij,ddip2) !!*************************
     USE Globals
     IMPLICIT NONE
     INTEGER :: i,j
     REAL(8) :: d2, ddip2
     REAL(8), DIMENSION(gND) :: rij, rdipij, gShVeci, gShVecj

     rij    = gR(i,:) - gR(j,:)
     rij    = MinImage( rij )
     gShVeci = gShVec(i,:)
     gShVecj = gShVec(j,:)  !this was the error --Ronal
     rdipij = rij(:) + (gShVeci(:) - gShVecj(:)) 	         !! *Update July 14, 2015
     d2     = SUM(rij(:) * rij(:))
     ddip2  = SUM(rdipij(:) * rdipij(:))
     RETURN

  END SUBROUTINE
  
  SUBROUTINE dist(i,j,rij,d2)
    USE Globals
    IMPLICIT NONE
    INTEGER :: i,j
    REAL(8) :: d2, ddip2
    REAL(8), DIMENSION(gND) :: rij, rdipij

    rij    = gR(i,:) - gR(j,:)
    rij    = MinImage( rij )
    d2     = SUM(rij(:) * rij(:))
    RETURN

  END SUBROUTINE 

 
  FUNCTION vdot(a,b)
    REAL(8) :: vdot
    REAL(8), DIMENSION(:), INTENT(IN) :: a
    REAL(8), DIMENSION(:), INTENT(IN) :: b
    vdot = SUM( a(:) * b(:) )
  END FUNCTION vdot

 
 
  FUNCTION MinImage( diff )
    USE Globals
    IMPLICIT NONE
    INTEGER :: i
    REAL(8), DIMENSION(gND), INTENT(IN)  :: diff
    REAL(8), DIMENSION(gND)              :: MinImage
    
    DO i = 1, gND
       MinImage(i) = diff(i) - gCDim(i) * ANINT( diff(i) / gCDim(i) )
    END DO

  END FUNCTION MinImage



  FUNCTION Contain(v) 
    USE Globals
    IMPLICIT NONE
    INTEGER :: i
    REAL(8), DIMENSION(gND), INTENT(IN)  :: v
    REAL(8), DIMENSION(gND)              :: Contain
  
    DO i = 1, gND
       Contain(i) = v(i) - gCDim(i) * ANINT( v(i) / gCDim(i) )
    END DO

  END FUNCTION Contain

! FUNCTION MagneticForce(RijM,Rij4,t,Ni,Nj)
!  USE Globals
!  IMPLICIT NONE
!  REAL(8) :: C1, C2, C3, C4
!  REAL(8) :: Lamb12Rij4
!  REAL(8) :: RijC, RijM, Rij4
!  REAL(8) :: CORR1, CORR2
!  REAL(8), DIMENSION(gND) :: t, Ni, Nj
!  REAL(8), DIMENSION(gND) :: MagneticForce
!
!  Lamb12Rij4 = gLambda12 / Rij4
!  C1 = vdot(Ni,Nj)
!  C2 = vdot(Ni,t )
!  C3 = vdot(Nj,t )
!  C4 = C1 - 3.D0 * C2 * C3
!  RijC  = RijM / gCutOff
!  CORR1 = 1.D0 - ( RijC ) ** 4
!  CORR2 = 1.D0 - 4.D0 * ( RijC ) ** 3 + 3.D0 * ( RijC ) ** 4
!
!  MagneticForce(:) = Lamb12Rij4 * ( ( CORR1 * C4 * t(:) ) + CORR2 * ( C3 * Ni(:) + C2 * Nj(:) - 2.D0 * C2 * C3 * t(:) ) )
!
!END FUNCTION MagneticForce

 FUNCTION MagneticForce(RijM,t,Ni,Nj)
  USE Globals
  IMPLICIT NONE
  REAL(8) :: C1, C2, C3, C4
  REAL(8) :: Lamb12Rij4
  REAL(8) :: RijM, Rij4
  REAL(8), DIMENSION(gND) :: t, Ni, Nj
  REAL(8), DIMENSION(gND) :: MagneticForce
  
  RIJ4 = RIJM ** 4
  Lamb12Rij4 = gLambda12 / Rij4
  C1 = vdot(Ni,Nj)
  C2 = vdot(Ni,t )
  C3 = vdot(Nj,t )
  C4 = C1 - 5.D0 * C2 * C3

  MagneticForce(:) = Lamb12Rij4 * ( C4 * t(:)  +  C3 * Ni(:) + C2 * Nj(:) )

END FUNCTION MagneticForce


!FUNCTION StericForce(Fij,RijM,t)
!  USE Globals
!  IMPLICIT NONE
!  REAL(8), DIMENSION(gND) :: StericForce, Fij, t
!  REAL(8)                 :: RijM
!  
!  StericForce(:) = Fij(:) + gStericTVInv * DLOG( 2.D0 / RijM) * t(:)
!
!END FUNCTION StericForce

FUNCTION StericForce(RijM,t)
  USE Globals
  IMPLICIT NONE
  REAL(8), DIMENSION(gND) :: StericForce, Fij, t
  REAL(8)                 :: RijM

  StericForce(:) = gStericTVInv * DLOG( 2.D0 / RijM) * t(:)

END FUNCTION StericForce


FUNCTION WCAForce(Rij2,t)
  USE Globals
  IMPLICIT NONE
  REAL(8), DIMENSION(gND) :: WCAForce, Fij, t
  REAL(8)                 :: Rij2, RIJM, SR2, SR6, SR12, Eps24R
  
  
  RIJM = DSQRT(RIJ2)    ! Papers de sofia o buttinoni
  SR2  = gSigma2 / RIJ2
  SR6  = SR2 ** 3.D0
  SR12 = SR6 ** 2.D0
  Eps24R = gEps24 / RIJM
  
  WCAForce(:) = Eps24R * (2.D0*SR12 - SR6) * t(:)

END FUNCTION WCAForce

!FUNCTION MagneticTorque(RijM,Rij3,t,Ni,Nj)
!  USE Globals
!  IMPLICIT NONE
!  REAL(8), DIMENSION(gND) :: MagneticTorque
!  REAL(8), DIMENSION(gND) :: Ni, Nj, t
!  REAL(8), DIMENSION(gND) :: B1, B2, C3
!  REAL(8)                 :: TCORR, TCORR3
!  REAL(8)                 :: RijM, Rij3, RijC
!
!  B1 = vcross(Ni,Nj)
!  B2 = vcross(Ni,t)
!  C3 = vdot(Nj,t )
!  RijC  = RijM / gCutOff
!  TCORR  = 1.D0 - 4.D0 * ( RijC ) ** 3 + 3.D0 * ( RijC ) ** 4
!  TCORR3 = TCORR / Rij3
!
!  MagneticTorque(:) = - TCORR3 * ( B1(:) - 3.D0 * C3 * B2(:) )
!  MagneticTorque(:) = gLambda4 * MagneticTorque(:)
!END FUNCTION MagneticTorque

FUNCTION MagneticTorque(RIJM,t,Ni,Nj)
  USE Globals
  IMPLICIT NONE
  REAL(8), DIMENSION(gND) :: MagneticTorque
  REAL(8), DIMENSION(gND) :: Ni, Nj, t
  REAL(8), DIMENSION(gND) :: B1, B2, C3
  REAL(8)                 :: TCORR, TCORR3
  REAL(8)                 :: RijM, Rij3, RijC

  RIJ3 = RIJM ** 3
  B1 = vcross(Ni,Nj)
  B2 = vcross(Ni,t)
  C3 = vdot(Nj,t )
 
  MagneticTorque(:) = - ( gLambda4 / RIJ3 ) * ( B1(:) - 3.D0 * C3 * B2(:) )

END FUNCTION MagneticTorque

!FUNCTION ShiftedDipoleInducedTorque(Ni,FMij)
!  USE GLOBALS
!  IMPLICIT NONE
!  REAL(8), DIMENSION(gND) :: ShiftedDipoleInducedTorque
!  REAL(8), DIMENSION(gND) :: Ni
!  REAL(8), DIMENSION(gND) :: FMij
!
!  ShiftedDipoleInducedTorque(:) = gSh * vcross(Ni,FMij)
!  
!END FUNCTION ShiftedDipoleInducedTorque

FUNCTION TRANSFORM(A,V)
  USE Globals
  IMPLICIT NONE
  INTEGER                     :: i
  REAL(8), DIMENSION(gND)     :: TRANSFORM
  REAL(8), DIMENSION(gND,gND) :: A
  REAL(8), DIMENSION(gND)     :: V
  
  DO i = 1, gND
     TRANSFORM(i) = SUM( A(i,:) * V(:) )
  END DO

END FUNCTION TRANSFORM

END MODULE VectorFunctions
