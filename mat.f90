program matrixtest
  implicit none
  real(8), dimension(3,3) :: A
  real(8), dimension(3)   :: V, B, TRANSFORM
  
  A(1,1) = 1.d0
  A(1,2) = 0.d0
  A(1,3) = -2.d0
  
  A(2,1) = 3.d0
  A(2,2) = -1.d0
  A(2,3) = 1.d0
  
  A(3,1) = 5.d0
  A(3,2) = 1.d0
  A(3,3) = 3.d0

  V(1) = 1.d0
  V(2) = -4.d0
  V(3) = 2.d0
  
  !B = TRANSFORM(A,V) 
  !B(1:3) = SUM( A(1:3,:) * V(:) )
  B(1) = SUM( A(1,:) * V(:) )
  B(2) = SUM( A(2,:) * V(:) )
  B(3) = SUM( A(3,:) * V(:) )
  write(*,*) B
end program matrixtest

FUNCTION TRANSFORM(A,V)
  IMPLICIT NONE
  INTEGER                 :: i
  REAL(8), DIMENSION(3)   :: TRANSFORM
  REAL(8), DIMENSION(3,3) :: A
  REAL(8), DIMENSION(3)   :: V

  DO i = 1, 3
     TRANSFORM(i) = SUM( A(i,:) * V(:) )
  END DO

END FUNCTION TRANSFORM

