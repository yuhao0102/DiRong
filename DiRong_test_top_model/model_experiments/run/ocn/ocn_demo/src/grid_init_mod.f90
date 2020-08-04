module grid_init_mod

    real(kind=4), public, allocatable :: lat(:), lon(:)
    integer, public           :: latlen, lonlen

contains


    subroutine grid_init

        use spmd_init_mod, only:masterproc, mpicomm, ier
        use parse_namelist_mod
        use mpi
        implicit none
        include "netcdf.inc"

        character*1024 :: input_data_dir, input_file_name
        character*1024 :: input_file_dir_name
        integer        :: ncid_input, ret
        integer        :: latid, lonid
        integer        :: latdimid, londimid
        integer        :: ii, i

        if(.not. read_ncfile) then
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
        else
          if (masterproc) then
              ret = nf_open (grid_file_name, nf_nowrite, ncid_input)
            
              ret = nf_inq_dimid (ncid_input, "grid_size", londimid)
              ret = nf_inq_dimlen (ncid_input, londimid, latlen)
              lonlen = latlen 
              allocate(lat(latlen), lon(lonlen))
              ret = nf_inq_varid (ncid_input, "grid_center_lat", latid)
              ret = nf_inq_varid (ncid_input, "grid_center_lon", lonid)
              ret = nf_get_var_real (ncid_input, latid, lat)
              ret = nf_get_var_real (ncid_input, lonid, lon)
          end if

          call mpi_bcast(latlen, 1, mpi_integer, 0, mpicomm, ier)
          call mpi_bcast(lonlen, 1, mpi_integer, 0, mpicomm, ier)

          if (masterproc == .false.) then
              allocate(lat(latlen), lon(lonlen))
          endif
 
          call mpi_bcast(lat, latlen, mpi_real4, 0, mpicomm, ier)
          call mpi_bcast(lon, lonlen, mpi_real4, 0, mpicomm, ier)
        endif
    end subroutine grid_init

end module grid_init_mod
