!!****h* BrownianDynamics/IntLLClass
!!
!! NAME
!!  IntLLClass
!!
!! SYNOPSIS
!!  This module is used to track groups of particles which are close enough
!!  to one another to collide.  It uses linked lists to track the position
!!  of particles within these smaller cells which are sized so that particles
!!  can only collide with other particles within the same cell or within a
!!  neighboring cell.
!!
!! NOTES
!!  Located in file ll.f90
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
module IntLLClass
implicit none

    !!****t* IntLLClass/IntLL
    !!
    !! NAME
    !!  IntLL
    !!
    !! SYNOPSIS
    !!  This derived type contains all the data necessary for building the
    !!  the linked cells.  That data includes the number of particles in the cell
    !!  and an array storing the ID of those particles.  This data is updated
    !!  several times per time step while checking for particle collisions.
    !!
    !! NOTES
    !!  Located in file ll.f90
    !!******
    type :: IntLL
    Integer :: val
    type ( IntLL ), pointer :: next 
    end type IntLL

contains

    !!****m* IntLLClass/IntAddValue
    !!
    !! NAME
    !!  IntAddValue
    !!
    !! SYNOPSIS
    !!  This subroutine adds a specified value to a specified linked list.
    !!
    !! USAGE
    !!  subroutine IntAddValue( list, value )
    !!
    !! INPUTS
    !!  list - an initialized linked list.
    !!
    !!  value - the value to be added to the linked list.
    !!
    !! NOTES
    !!  Located in file ll.f90
    !!******
    subroutine IntAddValue( list, value )
    implicit none

        type ( IntLL ), pointer :: list
        type ( IntLL ), pointer :: hold
        Integer :: value

        allocate( hold )
        
        hold % val = list % val
        hold % next => list % next

        list % val = value

        list % next => hold

    end subroutine IntAddValue


    !!****f* IntLLClass/IntGetValue
    !!
    !! NAME
    !!  IntGetValue
    !!
    !! SYNOPSIS
    !!  This subroutine gets the tail value of a specified linked list.
    !!
    !! USAGE
    !!  Integer function IntGetValue( list )
    !!
    !! INPUTS
    !!  list - an initialized linked list.
    !!
    !! RESULTS
    !!  The integer value composing the tail of the linked list.
    !!
    !! NOTES
    !!  Located in file ll.f90
    !!******
    Integer function IntGetValue( list )
    implicit none

        type ( IntLL ), pointer :: list
        type ( IntLL ), pointer :: hold

        IntGetValue = list % val

        hold => list % next

        deallocate( list )

        list => hold

    end function IntGetValue

end module IntLLClass
