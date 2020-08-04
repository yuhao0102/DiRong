module parse_namelist_mod
    integer, public :: time_step, coupling_freq, decomp_type_id
    integer, public :: num_point
    real(kind=4), public :: overlapping_rate
    logical, public      :: read_ncfile
    character(len=30),public :: grid_file_name
contains
    subroutine parse_namelist
        implicit none
        namelist /atm_demo_nml/ time_step ,decomp_type_id  , &
            coupling_freq, num_point, overlapping_rate, read_ncfile, &
            grid_file_name
        open(10, file="./atm_demo.nml")
        read(10, nml=atm_demo_nml)
    end subroutine parse_namelist
end module
