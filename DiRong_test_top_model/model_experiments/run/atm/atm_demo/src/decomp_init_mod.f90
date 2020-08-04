module decomp_init_mod
    integer, public              :: decomp_size
    integer, public, allocatable :: local_grid_cell_index(:)
contains
    subroutine decomp_init
        use spmd_init_mod, only: mpicomm, mytask_id, npes, ier
        use mpi
        use grid_init_mod, only:latlen, lonlen
        use parse_namelist_mod, only : decomp_type_id, overlapping_rate, read_ncfile
        implicit none

        integer :: i, ii,j, k, t, additional_size, num_total_point,max_num_point
        integer :: x_size, y_size, task_id, displ
        integer :: xdisp,xcount,ydisp,ycount
        integer :: displs(1:npes)
        integer :: sndcnts(1:npes)
        integer :: max_num_point_proc_num,num_remaining_point 
        if(.not. read_ncfile) then
            num_total_point = latlen*lonlen
        else
            num_total_point = latlen
        endif
        decomp_size = num_total_point/npes
        if(mod(mytask_id, npes).lt.mod(num_total_point, npes)) then
            decomp_size = decomp_size+1
        endif

        if (decomp_type_id == 1) then
            allocate(local_grid_cell_index(decomp_size))
            do i = 1, decomp_size
                local_grid_cell_index(i) = mytask_id+1+(i-1)*npes
                if (local_grid_cell_index(i) .gt. num_total_point) then
                    local_grid_cell_index(i) = local_grid_cell_index(i) - num_total_point
                endif
            end do
        else if (decomp_type_id == 2) then
            additional_size = int(real(decomp_size)*(overlapping_rate - 1.0))
            ! 这个参数大于1的话是重叠，重叠的个数这个，如果参数小于1，则是缺失了，缺失的话addition是负的
            allocate(local_grid_cell_index(decomp_size+additional_size))

            decomp_size = num_total_point/npes
            if(mytask_id .ge. mod(num_total_point, npes)) then
                displ = decomp_size*mytask_id + mod(num_total_point, npes)
            else
                displ = decomp_size*mytask_id + mytask_id
                decomp_size = decomp_size + 1
            endif
            do i = 1, decomp_size + additional_size
                local_grid_cell_index(i) = i+displ
                if (local_grid_cell_index(i) .gt. num_total_point) then
                    local_grid_cell_index(i) = local_grid_cell_index(i) - num_total_point
                endif
            end do
            decomp_size = decomp_size + additional_size
        else if (decomp_type_id == 3) then
            ii = int(sqrt(real(npes)))
            do i=ii, 1, -1
                if( (mod(npes, i).eq.0) ) then
                    y_size = i
                    x_size = npes/i
                    exit
                endif
            enddo
            xcount = lonlen / x_size
            ycount = latlen / y_size
            if(mod(mytask_id, x_size).lt.mod(lonlen, x_size)) then
                xcount = xcount + 1
            endif
            if(((mytask_id)/x_size).lt.mod(latlen, y_size)) then
                ycount = ycount + 1
            endif              

            if(mod(mytask_id,x_size) .ge. mod(lonlen, x_size)) then
                xdisp = (lonlen/x_size)*mod(mytask_id,x_size)+mod(lonlen, x_size)
            else
                xdisp = (lonlen/x_size)*mod(mytask_id,x_size)+mod(mytask_id, x_size)
            endif
            if (mytask_id/x_size .ge. mod(latlen, y_size)) then
               ydisp = (latlen/y_size)*(mytask_id/x_size)+mod(latlen, y_size)
            else
               ydisp = (latlen/y_size)*(mytask_id/x_size)+mytask_id/x_size
            endif
            allocate(local_grid_cell_index(xcount*ycount))
            decomp_size = xcount*ycount
            k=1
            do j=ydisp, ydisp+ycount-1
                do i=xdisp, xdisp+xcount-1
                    local_grid_cell_index(k)= i + 1 + lonlen*j
                    k=k+1
                enddo
            enddo
        else if (decomp_type_id == 4) then
            ! 负载不平衡
            if ( overlapping_rate .lt. 1.0) then
                print *, "ERROR! The load balancing rate must greater than 1"
                call mpi_finalize(ier)
            endif
            
            max_num_point = int(real(num_total_point/npes)*3.5)
            max_num_point_proc_num = num_total_point/2/max_num_point
            do i=1, max_num_point_proc_num
                sndcnts(i) = max_num_point
            enddo
            num_remaining_point=num_total_point-sum(sndcnts(1:max_num_point_proc_num))
            do i=max_num_point_proc_num+1, npes
                sndcnts(i)=num_remaining_point/(npes-max_num_point_proc_num)
            enddo
            num_remaining_point=num_total_point-sum(sndcnts)
            if (num_remaining_point .gt. 0) then
                do i=1,num_remaining_point
                    sndcnts(i) = sndcnts(i)+1
                enddo
            else
                do i=1,-num_remaining_point
                    sndcnts(i) = sndcnts(i)-1
                enddo
            endif
            displs(1) = 0
            do i=2, npes
                displs(i) = displs(i-1) + sndcnts(i-1)
            end do
            if(mytask_id==0) then
                print *, displs
                print *, sndcnts
            endif
            decomp_size=sndcnts(mytask_id+1)
            allocate(local_grid_cell_index(decomp_size))
            do i = 1, decomp_size
                local_grid_cell_index(i) = i+displs(mytask_id+1)
            end do
        endif
        if(mytask_id.eq.0) then
            open(1, file="record_time", position="append")
            write(1,*), decomp_type_id, num_total_point, npes
            close(1)
        endif
!        print *, "the decomp_size is ", decomp_size
    end subroutine decomp_init

end module decomp_init_mod
