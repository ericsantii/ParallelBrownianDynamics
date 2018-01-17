!!****************************************************************************!!
!!****m* BrownianDynamics/Matrices
!!
!!  NAME
!!   Orientation
!!
!!  SYNOPSIS
!!   This subroutine builts the A matrix transformation for X-CONVENTION 
!!   using the Euler parameters (quaternions) e0, e1, e2 and e3.
!!     
!!      E0, E1, E2, E3      : Represent the quaternion parameters
!!      AMATRIX             : Represent the matrix transformation A(3x3)
!!      AMATRIXTRANS        : Represent the A( i, j ) transpose
!!      N                   : Represent the dimension of the matrices 
!!
!!  NOTES
!!   Located in file matrices.f90
!!
!!  AUTHOR
!!   Glenn C. Vidal-Urquiza
!!
!! COPYRIGHT
!!    31-August-2012 - Glenn C. Vidal
!!*****
!! UPDATE
!!    11-July-2013 - Glenn C. Vidal-Urquiza
!!******	
!!
!!****************************************************************************!!
subroutine Matrices ( N, E0, E1, E2, E3, AMATRIX, AMATRIXTRANS )
    use Globals
    IMPLICIT NONE
    REAL( 8 ) :: E0, E1, E2, E3
    REAL( 8 ), DIMENSION ( N , N ) :: AMATRIX
    REAL( 8 ), DIMENSION ( N, N ) :: AMATRIXTRANS
    INTEGER :: N
    INTEGER :: J, K

!** Transformation matrix A( i, j ) - (X-CONVENTION)              

    AMATRIX( 1, 1 ) =  E0 ** 2 + E1 ** 2 - E2 ** 2 - E3 ** 2                    !(+)
    AMATRIX( 1, 2 ) = 2.0D0 * ( E1 * E2 + E0 * E3 )                             !(+)
    AMATRIX( 1, 3 ) =  2.0D0 * ( E1 * E3 - E0 * E2 )                            !(+)
    
    AMATRIX( 2, 1 ) = 2.0D0 * ( E1 * E2 - E0 * E3 )                             !(+)
    AMATRIX( 2, 2 ) =  E0 ** 2 - E1 ** 2 + E2 ** 2 - E3 ** 2                    !(+)
    AMATRIX( 2, 3 ) =  2.0D0 * ( E2 * E3 + E0 * E1 )                            !(+)

    AMATRIX( 3, 1 ) = 2.0D0 * ( E1 * E3 + E0 * E2 )                             !(+)
    AMATRIX( 3, 2 ) = 2.0D0 * ( E2 * E3 - E0 * E1 )                             !(+)
    AMATRIX( 3, 3 ) =  E0 ** 2 - E1 ** 2 - E2 ** 2 + E3 ** 2                    !(+)

!** Transpose of A( i, j )

    DO J = 1, 3
       DO K = 1, 3
          AMATRIXTRANS( J, K ) = AMATRIX( K, J )       
       END DO
    END DO
       
RETURN
end subroutine Matrices
