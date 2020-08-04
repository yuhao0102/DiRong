module coupling_ocn_model_mod

    use m_MCTWorld,only: MCTWorld_init => init
    use m_MCTWorld,only: MCTWorld_clean => clean
    use m_MCTWorld,only: ThisMCTWorld
    use m_MCTWorld,only: ComponentToWorldRank
!---Domain Decomposition Descriptor DataType and associated methods
    use m_GlobalSegMap,only: GlobalSegMap
    use m_GlobalSegMap,only: GlobalSegMap_init => init
    use m_GlobalSegMap,only: GlobalSegMap_lsize => lsize
    use m_GlobalSegMap,only: GlobalSegMap_clean => clean
    use m_GlobalSegMap,only: GlobalSegMap_Ordpnts => OrderedPoints
!  Field storage data types and associated methods.
    USE m_AttrVect, ONLY : AttrVect
    USE m_AttrVect, ONLY : AttrVect_init => init
    USE m_AttrVect, ONLY : AttrVect_zero => zero
    USE m_AttrVect, ONLY : AttrVect_clean => clean
    USE m_AttrVect, ONLY : AttrVect_indxR => indexRA
    USE m_AttrVect, ONLY : AttrVect_importRAttr => importRAttr
    USE m_AttrVect, ONLY : AttrVect_exportRAttr => exportRAttr
!  Intercomponent communitcations scheduler.
    USE m_Router, ONLY : Router
    USE m_Router, ONLY : Router_init => init
    USE m_Router, ONLY : Router_clean => clean
!  Intercomponent transfer.
    USE m_Transfer, ONLY: MCT_send => send
    USE m_Transfer, ONLY: MCT_recv => recv
    USE m_Transfer, ONLY: MCT_isend => isend
    USE m_Transfer, ONLY: MCT_irecv => irecv
    USE m_Transfer, ONLY: MCT_waitr => waitrecv
    USE m_Transfer, ONLY: MCT_waits => waitsend
 
    USE m_SparseMatrix, ONLY : SparseMatrix
    USE m_SparseMatrix, ONLY : SparseMatrix_init => init
    USE m_SparseMatrix, ONLY : SparseMatrix_importGRowInd => importGlobalRowIndices
    USE m_SparseMatrix, ONLY : SparseMatrix_importGColInd => importGlobalColumnIndices
    USE m_SparseMatrix, only : SparseMatrix_importMatrixElts => importMatrixElements

    USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus
    USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus_init => init
    USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus_clean => clean
    USE m_SparseMatrixPlus, ONLY : Xonly ! Decompose matrix by row
    USE m_SparseMatrixPlus, ONLY : time_in_SparseMatrixPlus

    USE m_MatAttrVectMul, ONLY : MCT_MatVecMul => sMatAvMult
    use mpi
    implicit none
    include "netcdf.inc"

    integer, public        :: atm_demo_comp_id
    integer, public        :: ncomps = 2 ! total component number
    integer, public        :: AtmID = 1, OcnID = 2  ! pick an id for the atmosphere
    TYPE(GlobalSegMap)     :: GSMap_A, GSMapInterp_O
    TYPE(AttrVect)         :: atm2ocn_A, ocn2atm_A, atm2ocn_A2, ocn2atm_A2
    TYPE(SparseMatrix)     :: sMat_A, sMat_O
    TYPE(SparseMatrixPlus) :: A2OMatPlus, O2AMatPlus
    TYPE(Router)           :: Rout
    INTEGER, DIMENSION(:), POINTER :: rows, columns, rows2, columns2
    REAL, DIMENSION(:), POINTER    :: weights, weights2
contains

    subroutine check_is_nf_success(istatus)
        use netcdf
        implicit none
        integer, intent(in) :: istatus
        if(istatus .ne. nf90_noerr) then
           print *,trim(nf90_strerror(istatus))
        endif
    end subroutine

    subroutine register_atm_demo_component(comm)
        integer, intent(inout) :: comm
        integer                :: ierr

        call MPI_INIT(ierr)
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,AtmID,0,comm,ierr)

        call MCTWorld_init(ncomps,MPI_COMM_WORLD,comm,AtmID)

    end subroutine register_atm_demo_component

    subroutine register_component_coupling_configuration
       use decomp_init_mod
       use spmd_init_mod
       use grid_init_mod
       use parse_namelist_mod
       integer                     :: i, j, k, temp_length, ierr, another_root, another_grid_size22
       integer                     :: block_num, begin_idx, Isize, Asize
       integer, allocatable        :: start(:), length(:)
       integer                     :: num_elements, nRows, nColumns
       real, dimension(:), pointer :: A, A2
       integer                     :: my_grid_size, another_grid_size, ncid_input, ret
       integer                     :: row_id, col_id, s_id, src_id, dst_id, element_id
       real                        :: time1(2), time2(2), time3(2), time4(2), time5(2), time6(2), time7(2), time8(2)
       real                        :: local_time1, local_time2, max_time1, max_time2, local_time3, local_time4, max_time3, max_time4, local_time5, local_time6, max_time5, max_time6, local_time7, max_time7, local_time8, max_time8
       integer status(MPI_STATUS_SIZE) 
    
        another_root = ComponentToWorldRank(0, OcnID, ThisMCTWorld)        
        if(mytask_id.eq.0) call mpi_recv(another_grid_size22, 1, MPI_INTEGER, another_root, 0, MPI_COMM_WORLD, status, ierr)
        call mpi_bcast(another_grid_size22, 1, MPI_INTEGER, 0, mpicomm, ierr)
        
        call mpi_barrier(mpicomm, ierr)
        call cpu_time(time1(1))

        block_num = 1
        begin_idx = local_grid_cell_index(1)
        do i = 2, decomp_size 
            if(begin_idx+1 .ne. local_grid_cell_index(i)) then
                block_num = block_num + 1
            endif
            begin_idx = local_grid_cell_index(i)
        enddo

        allocate(start(block_num), length(block_num))
        k = 1
        temp_length = 1
        start(k) = local_grid_cell_index(1) 
        begin_idx = local_grid_cell_index(1)
        do i = 2, decomp_size
            if(begin_idx+1 .ne. local_grid_cell_index(i)) then
                length(k) = temp_length
                temp_length = 0
                k = k + 1
                start(k) = local_grid_cell_index(i)
            endif
            begin_idx = local_grid_cell_index(i)
           temp_length = temp_length + 1
        enddo
        length(block_num) = temp_length

        call cpu_time(time3(1))
        call GlobalSegMap_init (GSMap_A, start, length, 0, mpicomm, AtmID)
        call cpu_time(time3(2))
        
        deallocate(start, length)
   
        call cpu_time(time4(1))
        if(mytask_id .eq. 0) then
            call check_is_nf_success(nf_open ("weight/weight_ocn_to_atm_"//weight_file_index//".nc", nf_nowrite, ncid_input))
            call check_is_nf_success(nf_inq_dimid (ncid_input, "n_a", src_id))
            call check_is_nf_success(nf_inq_dimid (ncid_input, "n_b", dst_id))
            call check_is_nf_success(nf_inq_dimid (ncid_input, "n_s", element_id))
            call check_is_nf_success(nf_inq_dimlen (ncid_input, src_id, my_grid_size))
            call check_is_nf_success(nf_inq_dimlen (ncid_input, dst_id, another_grid_size))
            call check_is_nf_success(nf_inq_dimlen (ncid_input, element_id, num_elements))

            if(read_ncfile .and. (another_grid_size .ne. latlen)) then
                print *, "some thing wrong in weight file"
                stop
            endif
            allocate(rows(num_elements), columns(num_elements),weights(num_elements), stat=ierr)

            call check_is_nf_success(nf_inq_varid (ncid_input, "row", row_id))
            call check_is_nf_success(nf_inq_varid (ncid_input, "col", col_id))
            call check_is_nf_success(nf_inq_varid (ncid_input, "S", s_id))
            call check_is_nf_success(nf_get_var_int (ncid_input, row_id, rows))
            call check_is_nf_success(nf_get_var_int (ncid_input, col_id, columns))
            call check_is_nf_success(nf_get_var_real (ncid_input, s_id, weights))

            nRows = another_grid_size
            nColumns = my_grid_size
            call SparseMatrix_init(sMat_O, nRows, nColumns, num_elements)
            call SparseMatrix_importGRowInd(sMat_O, rows, num_elements)
            call SparseMatrix_importGColInd(sMat_O, columns, num_elements)
            call SparseMatrix_importMatrixElts(sMat_O, weights, num_elements)
            deallocate(rows, columns, weights, stat=ierr)

            call check_is_nf_success(nf_open ("weight/weight_atm_to_ocn_"//weight_file_index//".nc", nf_nowrite, ncid_input))
            call check_is_nf_success(nf_inq_dimid (ncid_input, "n_a", src_id))
            call check_is_nf_success(nf_inq_dimid (ncid_input, "n_b", dst_id))
            call check_is_nf_success(nf_inq_dimid (ncid_input, "n_s", element_id))
            call check_is_nf_success(nf_inq_dimlen (ncid_input, src_id, my_grid_size))
            call check_is_nf_success(nf_inq_dimlen (ncid_input, dst_id, another_grid_size))
            call check_is_nf_success(nf_inq_dimlen (ncid_input, element_id, num_elements))

            if(read_ncfile .and. (my_grid_size .ne. latlen)) then
                print *, "some thing wrong in weight file"
                stop
            endif
                        
            allocate(rows(num_elements), columns(num_elements), weights(num_elements), stat=ierr)
            
            call check_is_nf_success(nf_inq_varid (ncid_input, "row", row_id))
            call check_is_nf_success(nf_inq_varid (ncid_input, "col", col_id))
            call check_is_nf_success(nf_inq_varid (ncid_input, "S", s_id))
            call check_is_nf_success(nf_get_var_int (ncid_input, row_id, rows))
            call check_is_nf_success(nf_get_var_int (ncid_input, col_id, columns))
            call check_is_nf_success(nf_get_var_real (ncid_input, s_id, weights))
            
            nRows = another_grid_size
            nColumns = my_grid_size
            call SparseMatrix_init(sMat_A, nRows, nColumns, num_elements)
            call SparseMatrix_importGRowInd(sMat_A, rows, num_elements)
            call SparseMatrix_importGColInd(sMat_A, columns, num_elements)
            call SparseMatrix_importMatrixElts(sMat_A, weights, num_elements)
            deallocate(rows, columns, weights, stat=ierr)
        endif
        call cpu_time(time4(2))
        
!  copy from COAWST, create a GSMap for ocean, but the start and length arrays seem differ with the ocean decomposition
! more like a intermediate decomposition of sequential decomposition
! because dst_grid_dims(1)==latlen, dst_grid_dims(2)==lonlen, so 
!        call cpu_time(temo(1))
!        call mpi_bcast(another_grid_size, 1, MPI_INTEGER, 0, mpicomm, ierr)
!        call cpu_time(temo(2))
!        print *, "wowo", temo(2)-temo(1)

        Isize = INT(another_grid_size22 / npes)
        if (mytask_id .eq. npes-1) then
            Isize = another_grid_size22 - Isize * (npes - 1)
        endif
  
        allocate (start(1), length(1))
        start(1) = mytask_id * (INT(another_grid_size22 / npes)) + 1
        length(1) = Isize

        call cpu_time(time5(1))
        CALL GlobalSegMap_init (GSMapInterp_O, start, length, 0, mpicomm, AtmID)
        call cpu_time(time5(2))
        
        if (allocated(start)) deallocate(start)
        if (allocated(length)) deallocate(length)

! Create ATM sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
        call cpu_time(time8(1))
        call SparseMatrixPlus_init(A2OMatPlus, sMat_A, GSMap_A, GSMapInterp_O, Xonly, 0, mpicomm, AtmID)
        local_time6 = time_in_SparseMatrixPlus(1)
        local_time7 = time_in_SparseMatrixPlus(2)
        call SparseMatrixPlus_init(O2AMatPlus, sMat_O, GSMapInterp_O, GSMap_A, Xonly, 0, mpicomm, AtmID)
        local_time6 = local_time6 + time_in_SparseMatrixPlus(1)
        local_time7 = local_time7 + time_in_SparseMatrixPlus(2)
        call cpu_time(time8(2))
     
! Create Ocean sparse matrix plus for interpolation.
! Specify matrix decomposition to be by row.
        Asize=GlobalSegMap_lsize(GSMap_A, mpicomm)
        CALL AttrVect_init(atm2ocn_A, rList="field1:field2:field3:field4:field5", lsize=Asize)
        CALL AttrVect_zero(atm2ocn_A)
        CALL AttrVect_init(ocn2atm_A, rList="field6:field7:field8:field9:field10", lsize=Asize)
        CALL AttrVect_zero(ocn2atm_A)

        Asize=GlobalSegMap_lsize(GSMapInterp_O, mpicomm)
        CALL AttrVect_init(atm2ocn_A2, rList="field1:field2:field3:field4:field5", lsize=Asize)
        CALL AttrVect_zero(atm2ocn_A2)
        CALL AttrVect_init(ocn2atm_A2, rList="field6:field7:field8:field9:field10", lsize=Asize)
        CALL AttrVect_zero(ocn2atm_A2)

        call cpu_time(time2(1))
        CALL Router_init(OcnID, GSMapInterp_O, mpicomm, Rout)
        call cpu_time(time2(2))
    
        call cpu_time(time1(2))
        local_time1 = (time1(2)-time1(1))
        local_time2 = (time2(2)-time2(1))
        local_time3 = (time3(2)-time3(1))
        local_time4 = (time4(2)-time4(1))
        local_time5 = (time5(2)-time5(1))
        local_time8 = (time8(2)-time8(1))
        call mpi_reduce(local_time1, max_time1, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        call mpi_reduce(local_time2, max_time2, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        call mpi_reduce(local_time3, max_time3, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        call mpi_reduce(local_time4, max_time4, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        call mpi_reduce(local_time5, max_time5, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        call mpi_reduce(local_time6, max_time6, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        call mpi_reduce(local_time7, max_time7, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        call mpi_reduce(local_time8, max_time8, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
       
        if(mytask_id.eq.0) then
            open(1, file="record_time_atm", position="append")
            write(1,*), max_time1, " ",max_time2, " ", max_time3
            write(1,*), max_time4, " ",max_time5, " ", max_time6
            write(1,*), max_time7, " ",max_time8
            close(1)
        endif

!!!!! try to send and recv data
        Asize=GlobalSegMap_lsize(GSMap_A, mpicomm)
        allocate(A(Asize))
        A = -2222.0
!        CALL AttrVect_importRAttr(atm2ocn_A, "field1", A, Asize)
!        CALL AttrVect_importRAttr(atm2ocn_A, "field2", A, Asize)
!        CALL MCT_MatVecMul(atm2ocn_A, A2OMatPlus, atm2ocn_A2)
!        CALL MCT_isend (atm2ocn_A2, Rout, 1)
!        CALL MCT_waits(Rout)
        
        Asize=GlobalSegMap_lsize(GSMap_A, mpicomm)
        allocate(A2(Asize))
!        CALL MCT_irecv(ocn2atm_A2, Rout, 2)
!        CALL MCT_waitr(ocn2atm_A2, Rout)
!        CALL MCT_MatVecMul(ocn2atm_A2, O2AMatPlus, ocn2atm_A)
!        CALL AttrVect_exportRAttr(ocn2atm_A, "field6", A2, Asize)

    end subroutine register_component_coupling_configuration

end module coupling_ocn_model_mod
