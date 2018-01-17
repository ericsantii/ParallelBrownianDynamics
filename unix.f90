PROGRAM unixsystem
IMPLICIT NONE
CHARACTER(LEN=30) :: cmd ! string to store the Unix command
cmd="wc -l Data/periodic_pos.out"
! As an example print the working directory
CALL SYSTEM(cmd)
END PROGRAM unixsystem
