program dotproduct
use vector_functions

implicit none
integer, parameter :: gND = 3
real(8), dimension(2,2) :: eps2
real(8), dimension(3) :: u
real(8), dimension(3) :: v
real(8), dimension(3) :: w
real(8) :: a
!real(8), dimension(3) :: cross

u = [1.0, 3.0, -1.0]
v = [2.0, 1.0, 0.0]

!write(*,*) sum( u(:) * v(:) )


  eps2(1,:)=(/ 0.D0 , 1.D0 /)
  eps2(2,:)=(/ -1.D0, 0.D0 /)

  w = cross(u,v)
  a = vdot(u,v)
 ! write(*,*) w
 ! write(*,*) a
   
end program dotproduct

!function cross(a, b)
! Parameters.
!real(8), dimension(3) :: a, b
! vectors a and b
! Return value.
!real(8), dimension(3) :: cross
! cross = a x b
! Compute cross product.
!cross(1) = a(2) * b(3) - a(3) * b(2)
!cross(2) = a(3) * b(1) - a(1) * b(3)
!cross(3) = a(1) * b(2) - a(2) * b(1)
!end function cross

!FUNCTION vcross(a,b) result(c)
!  IMPLICIT NONE
!  INTEGER :: i, j
!  REAL(8), INTENT(IN),  DIMENSION(3) :: a
!  REAL(8), INTENT(IN),  DIMENSION(3) :: b
!  REAL(8), INTENT(OUT), DIMENSION(3) :: c
!  real(8), dimension(3,3,3) :: eps3
!
!  eps3(1,:,1)=(/ 0.D0, 0.D0, 0.D0 /)
!  eps3(2,:,1)=(/ 0.D0, 0.D0, 1.D0 /)
!  eps3(3,:,1)=(/ 0.D0, -1.D0, 0.D0/)
!
!  eps3(1,:,2)=(/ 0.D0, 0.D0, -1.D0/)
!  eps3(2,:,2)=(/ 0.D0, 0.D0, 0.D0 /)
!  eps3(3,:,2)=(/ 1.D0, 0.D0, 0.D0 /)
!
!  eps3(1,:,3)=(/ 0.D0, 1.D0, 0.D0 /)
!  eps3(2,:,3)=(/ -1.D0, 0.D0, 0.D0/)
!  eps3(3,:,3)=(/ 0.D0, 0.D0, 0.D0 !/)
!
!  c(:) = 0.D0
!  DO i =1, 3
!    DO j=1, 3
!         c(:) = c(:) + a(i) * b(j) * eps3(i,j,:)
!    END DO
!  END DO
!
!END FUNCTION vcross
