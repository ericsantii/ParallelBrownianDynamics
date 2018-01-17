module vector_functions
  implicit none
  
contains
  
  function cross(a, b)
    ! Parameters.
    real(8), dimension(3) :: a, b
    ! vectors a and b
    ! Return value.
    real(8), dimension(3) :: cross
    ! cross = a x b
    ! Compute cross product.
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross
   
   FUNCTION vdot(a,b)
    REAL(8) :: vdot
    !INTEGER :: nel = size(a)
    REAL(8), DIMENSION(:) :: a
    REAL(8), DIMENSION(:) :: b

    vdot = SUM( a(:) * b(:) )

   END FUNCTION vdot  
  
end module vector_functions
