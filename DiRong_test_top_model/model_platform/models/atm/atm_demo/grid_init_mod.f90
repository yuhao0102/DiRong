module grid_init_mod

    real(kind=4), public, allocatable :: lat(:), lon(:)
    integer, public           :: latlen, lonlen

contains


    subroutine grid_init

        use spmd_init_mod, only:masterproc, mpicomm, ier
        use parse_namelist_mod, only: num_point, read_ncfile, grid_file_name
        use mpi
        implicit none
        include "netcdf.inc"

        character*1024 :: input_data_dir, input_file_name
        character*1024 :: input_file_dir_name
        integer        :: ncid_input, ret
        integer        :: latid, lonid
        integer        :: latdimid, londimid
        integer        :: ii, i
        
        if (.not. read_ncfile) then
          ii = sqrt(real(num_point))
          do i=1, ii
              if( (mod(num_point, i).eq.0) .and. (i .ne. num_point/i)) then
                  if (i .ge. num_point/i) then
                      lonlen = i
                      latlen = num_point/i
                  else
                      latlen = i
                      lonlen = num_point/i
                  endif
              endif
          enddo

          allocate(lon(lonlen),lat(latlen))
          do i=1,latlen
              lat(i) = -89.0 + 179.0 * i / latlen
          enddo
          do i=1,lonlen
              lon(i) = 0.0 + 359.0 * i / lonlen
          enddo
          if (masterproc) then
              print *, "atm use auto-generated lon-lat grid" 
          endif
      else
          if (masterproc) then
              open(1, file=grid_file_name, status="old")
              read(1, *) latlen
              lonlen = latlen
            
              allocate(lat(latlen), lon(lonlen))
              do i = 1, latlen
                  read(1, *) lon(i), lat(i)
              enddo
          end if

          call mpi_bcast(latlen, 1, mpi_integer, 0, mpicomm, ier)
          call mpi_bcast(lonlen, 1, mpi_integer, 0, mpicomm, ier)

          if (masterproc == .false.) then
              allocate(lat(latlen), lon(lonlen))
          endif

          call mpi_bcast(lat, latlen, mpi_real8, 0, mpicomm, ier)
          call mpi_bcast(lon, lonlen, mpi_real8, 0, mpicomm, ier)
        endif
     end subroutine grid_init

end module grid_init_mod
