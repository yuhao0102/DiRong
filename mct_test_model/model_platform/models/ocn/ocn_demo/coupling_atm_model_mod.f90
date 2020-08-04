module coupling_atm_model_mod

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
    implicit none

    integer, public        :: ocn_demo_comp_id
    integer, public        :: ncomps = 2 ! total component number
    integer, public        :: AtmID = 1, OcnID = 2  ! pick an id for the ocean
    integer                :: Asize
    TYPE(GlobalSegMap)     :: GSMap_O
    TYPE(AttrVect)         :: ocn2atm_A, atm2ocn_A
    TYPE(SparseMatrix)     :: sMat_O
    TYPE(SparseMatrixPlus) :: O2AMatPlus
    TYPE(Router)           :: Rout
    INTEGER, DIMENSION(:), POINTER :: rows, columns
    REAL, DIMENSION(:), POINTER    :: weights
    REAL, DIMENSION(:), POINTER    :: A, A2
contains

    subroutine register_ocn_demo_component(comm)
        use mpi
        integer, intent(inout) :: comm
        integer                :: ierr

        call MPI_INIT(ierr)
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,OcnID,0,comm,ierr)
        call MCTWorld_init(ncomps,MPI_COMM_WORLD,comm,OcnID)

    end subroutine register_ocn_demo_component

    subroutine register_component_coupling_configuration
        use mpi
        use decomp_init_mod
        use spmd_init_mod
        use grid_init_mod
        integer                  :: i, j, k, temp_length, ierr, another_root, total_grid_size
        integer                  :: block_num, begin_idx
        integer, allocatable     :: start(:), length(:)
        integer                  :: num_elements, nRows, nColumns
        real                     :: time1(2), time2(2), time3(2), time4(2), time5(2), time6(2)
        real                     :: local_time1, local_time2, max_time1, max_time2, local_time3, local_time4, max_time3, max_time4, local_time5, local_time6, max_time5, max_time6

        another_root = ComponentToWorldRank(0, AtmID, ThisMCTWorld)
        total_grid_size = num_total_point
        if(mytask_id.eq.0) call mpi_send(total_grid_size, 1, MPI_INTEGER, another_root, 0, MPI_COMM_WORLD, ierr)

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
        call GlobalSegMap_init(GSMap_O, start, length, 0, mpicomm, OcnID)
        call cpu_time(time3(2))

        deallocate(start)
        deallocate(length)

        call cpu_time(time4(1))
        call cpu_time(time4(2))

        call cpu_time(time5(1))
        call cpu_time(time5(2))

        Asize=GlobalSegMap_lsize(GSMap_O, mpicomm)
        CALL AttrVect_init(atm2ocn_A, rList="field1:field2:field3:field4:field5", lsize=Asize)
        CALL AttrVect_zero(atm2ocn_A)
        CALL AttrVect_init(ocn2atm_A, rList="field6:field7:field8:field9:field10", lsize=Asize)
        CALL AttrVect_zero(ocn2atm_A)

        call cpu_time(time2(1))        
        CALL Router_init(AtmID, GSMap_O, mpicomm, Rout)
        call cpu_time(time2(2))
        call cpu_time(time1(2))
        
        call cpu_time(time6(1))
        call cpu_time(time6(2))

        local_time1 = (time1(2)-time1(1))
        local_time2 = (time2(2)-time2(1))
        local_time3 = (time3(2)-time3(1))
        local_time4 = (time4(2)-time4(1))
        local_time5 = (time5(2)-time5(1))
        local_time6 = (time6(2)-time6(1))
        call mpi_reduce(local_time1, max_time1, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        call mpi_reduce(local_time2, max_time2, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        call mpi_reduce(local_time3, max_time3, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        call mpi_reduce(local_time4, max_time4, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        call mpi_reduce(local_time5, max_time5, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        call mpi_reduce(local_time6, max_time6, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)

        if(mytask_id.eq.0) then
                open(1, file="record_time_ocn", position="append")
                write(1,*), max_time1, " ",max_time2, " ", max_time3
                write(1,*), max_time4, " ",max_time5, " ", 0
                write(1,*), max_time6, 0
                close(1)
        endif    

!        CALL MCT_irecv (atm2ocn_A, Rout, 1)
!        CALL MCT_waitr (atm2ocn_A, Rout)
!        allocate(A(Asize))
!        CALL AttrVect_exportRAttr(atm2ocn_A, "field1", A, Asize)
!        CALL AttrVect_exportRAttr(atm2ocn_A, "field2", A, Asize)
        
!        allocate(A2(Asize))
!        A2 = 314.0
!        CALL AttrVect_importRAttr(ocn2atm_A, "field6", A2)
!        CALL MCT_isend(ocn2atm_A, Rout, 2)
!        CALL MCT_waits(Rout)

    end subroutine register_component_coupling_configuration

end module coupling_atm_model_mod
