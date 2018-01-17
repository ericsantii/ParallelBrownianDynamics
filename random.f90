program test_random_number
REAL :: r(5,5)
CALL random_seed()         ! see example of RANDOM_SEED
CALL RANDOM_NUMBER(r)

write(*,*) r
end program
