!!****h* BrownianDynamics/VectorClass
!!
!! NAME
!!  VectorClass
!!
!! SYNOPSIS
!!  This module comprises all the data and methods for implementing
!!  a vector type object which makes the geometrical manipulations
!!  in the program faster and easier to read.
!!
!! NOTES
!!  Located in file vec.f90
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
module VectorClass
implicit none

    !!****t* VectorClass/Vector
    !!
    !! NAME
    !!  VectorClass
    !!
    !! SYNOPSIS
    !!  This derived type stores data for three orthonormal vectors.
    !!
    !! NOTES
    !!  Located in file vec.f90
    !!******
    type :: Vector
        Real( 8 ) :: x
        Real( 8 ) :: y
        Real( 8 ) :: z
    end type Vector

    ! Allows the + operator to be used between vectors.
    interface operator ( + )
        module procedure VectorAdd
    end interface

    ! Allows the - operator to be used between vectors.
    interface operator ( - )
        module procedure VectorSub
    end interface

    ! Allows the * operator to be used between vectors (dot product)
    ! and between a vector and a scalar.
    interface operator ( * )
        module procedure VectorDot
        module procedure ScalarMult1
        module procedure ScalarMult2
    end interface

   ! Allows the x operator to be used between vectors (cross product)
   interface operator ( / )
        module procedure VectorCross
    end interface

contains

    !!****f* VectorClass/VectorAdd
    !!
    !! NAME
    !!  VectorAdd
    !!
    !! SYNOPSIS
    !!  This subroutine adds two vectors.
    !!
    !! USAGE
    !!  type ( Vector ) Function VectorAdd( v1, v2 )
    !!
    !! INPUTS
    !!  v1, v2 - two vector objects.
    !!
    !! RESULT
    !!  The sum of the two vector inputs.
    !!
    !! NOTES
    !!  Located in file vec.f90
    !!******
    type ( Vector ) Function VectorAdd( v1, v2 )
    implicit none
        type ( Vector ), intent( in ) :: v1
        type ( Vector ), intent( in ) :: v2

        VectorAdd%x = v1%x + v2%x
        VectorAdd%y = v1%y + v2%y
        VectorAdd%z = v1%z + v2%z

    end function VectorAdd

    !!****f* VectorClass/VectorSub
    !!
    !! NAME
    !!  VectorSub
    !!
    !! SYNOPSIS
    !!  This subroutine subtracts two vectors.
    !!
    !! USAGE
    !!  type ( Vector ) Function VectorSub( v1, v2 )
    !!
    !! INPUTS
    !!  v1, v2 - two vector objects.
    !!
    !! RESULT
    !!  The difference of the two vector inputs.
    !!
    !! NOTES
    !!  Located in file vec.f90
    !!******
    type ( Vector ) function VectorSub( v1, v2 )
    implicit none
        type ( Vector ), intent( in ) :: v1
        type ( Vector ), intent( in ) :: v2

        VectorSub%x = v1%x - v2%x
        VectorSub%y = v1%y - v2%y
        VectorSub%z = v1%z - v2%z

    end function VectorSub


    !!****f* VectorClass/VectorDot
    !!
    !! NAME
    !!  VectorDot
    !!
    !! SYNOPSIS
    !!  This subroutine gets the dot product of two vectors.
    !!
    !! USAGE
    !!  Real( 8 ) Function VectorDot( v1, v2 )
    !!
    !! INPUTS
    !!  v1, v2 - two vector objects.
    !!
    !! RESULT
    !!  The dot product of the two vector inputs.
    !!
    !! NOTES
    !!  Located in file vec.f90
    !!******
    Real( 8 ) function VectorDot( v1, v2 )
        type ( Vector ), intent( in ) :: v1
        type ( Vector ), intent( in ) :: v2

        VectorDot = v1 % x * v2 % x + v1 % y * v2 % y + v1 % z * v2 % z
 
    end function VectorDot


    !!****f* VectorClass/ScalarMult1
    !!
    !! NAME
    !!  ScalarMult1
    !!
    !! SYNOPSIS
    !!  This subroutine gets the dot product of two vectors.
    !!
    !! USAGE
    !!  type ( Vector ) function ScalarMult1( v, s )
    !!
    !! INPUTS
    !!  v - a vector object.
    !!
    !!  s - a real scalar.
    !!
    !! RESULT
    !!  The input vector rescaled by the input scalar.
    !!
    !! NOTES
    !!  Located in file vec.f90
    !!******
    type ( Vector ) function ScalarMult1( v, s )
        type ( Vector ), intent( in ) :: v
        Real( 8 ), intent( in ) :: s

        ScalarMult1 % x = v % x * s
        ScalarMult1 % y = v % y * s
        ScalarMult1 % z = v % z * s

    end function ScalarMult1


    !!****f* VectorClass/ScalarMult2
    !!
    !! NAME
    !!  ScalarMult2
    !!
    !! SYNOPSIS
    !!  This subroutine gets the dot product of two vectors.  Same
    !!  as above, but with the scalar/vector order reversed.
    !!
    !! USAGE
    !!  type ( Vector ) function ScalarMult2( s, v )
    !!
    !! INPUTS
    !!  v - a vector object.
    !!
    !!  s - a real scalar.
    !!
    !! RESULT
    !!  The input vector rescaled by the input scalar.
    !!
    !! NOTES
    !!  Located in file vec.f90
    !!******
    type ( Vector ) function ScalarMult2( s, v )
        Real( 8 ), intent( in ) :: s
        type ( Vector ), intent( in ) :: v

        ScalarMult2 % x = v % x * s
        ScalarMult2 % y = v % y * s
        ScalarMult2 % z = v % z * s

    end function ScalarMult2



!------------ Cross product implemented by Luis------------                                                                          
   type ( Vector ) Function VectorCross( v1, v2 )
    implicit none
        type ( Vector ), intent( in ) :: v1
        type ( Vector ), intent( in ) :: v2

        VectorCross%x = v1%y * v2%z - v2%y * v1%z
        VectorCross%y = v2%x * v1%z - v1%x * v2%z
        VectorCross%z = v1%x * v2%y - v1%y * v2%x

    end function VectorCross

end module VectorClass
