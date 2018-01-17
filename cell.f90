!!***** BrownianDynamics/CellClass
!!
!! NAME
!!  CellClass
!!
!! SYNOPSIS
!!  This module is used to track groups of particles which are close enough
!!  to one another to collide. It uses linked lists to track the position
!!  of particles within these smaller cells which are sized so that particles
!!  can only collide with other particles within the same cell or within a
!!  neighboring cell.
!!
!! NOTES
!!  Located in file cell.f90
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
module CellClass
use IntLLClass
implicit none

    !!****** CellClass/Cell
    !!
    !! NAME
    !!  CellClass
    !!
    !! SYNOPSIS
    !!  This derived type contains all the data necessary for building the
    !!  the linked cells.  That data includes the number of particles in the cell
    !!  and an array storing the ID of those particles.  This data is updated
    !!  several times per time step while checking for particle collisions.
    !!
    !! NOTES
    !!  Located in file cell.f90
    !!******
    type :: Cell
        ! The number of particles in the cell.
        Integer :: sz
        ! A linked list temporarily storing the IDs of the particles within the cell.
        type ( IntLL ), pointer :: ll

        ! An array storing the IDs of all the particles within the cell.
        Integer, allocatable, dimension( : ) :: part 

    end type Cell

end module CellClass
