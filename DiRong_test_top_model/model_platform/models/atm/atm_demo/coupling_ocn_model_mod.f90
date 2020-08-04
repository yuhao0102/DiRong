module coupling_ocn_model_mod

    use CCPL_interface_mod

    implicit none

    integer, public             :: atm_demo_comp_id
    integer, private            :: decomp_id, grid_h2d_id
contains

    subroutine register_atm_demo_component(comm)
        use CCPL_interface_mod
        integer, intent(inout) :: comm
        atm_demo_comp_id = CCPL_register_component(-1, "atm_demo", "atm", comm, change_dir=.true., annotation = "register atm model atm_demo")
    end subroutine register_atm_demo_component

    subroutine register_component_coupling_configuration
        use mpi
        use CCPL_interface_mod
        use spmd_init_mod, only: mpicomm, mytask_id, npes
        use parse_namelist_mod,only:time_step, coupling_freq, read_ncfile
        use grid_init_mod, only:latlen, lonlen, lon, lat
        use decomp_init_mod, only:decomp_size, local_grid_cell_index
        use variable_mod, only:pslm, prectm, fldsm, fsdsm, maskm, sst,ssh,shf,mld

        implicit none

        integer          :: export_interface_id, import_interface_id
        integer          :: timer_id, fields_id(5), ierr
        integer          :: field_id_psl, field_id_prect, field_id_flds, field_id_fsds, field_id_flds2
        integer          :: field_id_sst, field_id_ssh, field_id_shf, field_id_mld, field_id_mld2
        real             :: start_time, end_time
        real             :: total_time, max_total_time
        
        call CCPL_set_normal_time_step(atm_demo_comp_id, time_step, annotation="setting the time step for atm_demo")
        if(.not. read_ncfile) then        
            grid_h2d_id = CCPL_register_H2D_grid_via_global_data(atm_demo_comp_id, "atm_demo_H2D_grid", "LON_LAT", "degrees", "cyclic", lonlen, latlen, 0.0, 360.0, -90.0, 90.0, lon, lat, annotation="register atm_demo H2D grid ")
        else
            grid_h2d_id = CCPL_register_H2D_grid_via_global_data(atm_demo_comp_id, "atm_demo_H2D_grid", "LON_LAT", "degrees", "cyclic", lonlen, 0, 0.0, 360.0, -90.0, 90.0, lon, lat, annotation="register atm_demo H2D grid ")
        endif
        decomp_id = CCPL_register_normal_parallel_decomp("decomp_atm_demo_grid", grid_H2D_id, decomp_size, local_grid_cell_index, annotation="allocate decomp for atm_demo grid")
        !------------register field instances to C-Coupler2--------------

        field_id_psl = CCPL_register_field_instance(pslm, "psl", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="Pa", annotation="register field instance of Sea level pressure") 
        field_id_prect = CCPL_register_field_instance(prectm, "prect", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="m/s", annotation="register field instance of precipitation")
        field_id_fsds = CCPL_register_field_instance(fsdsm, "fsds", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="W/m2", annotation="register field instance of Short wave downward flux at surface")
        field_id_flds  = CCPL_register_field_instance(fldsm, "flds", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="W/m2", annotation="register field instance of Long wave downward flux at surface")
        field_id_flds2  = CCPL_register_field_instance(fldsm, "flds2", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="W/m2", annotation="register field instance of Long wave downward flux at surface")
        field_id_sst  = CCPL_register_field_instance(sst, "sst", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="C", annotation="register field instance of Sea surface temperature")
        field_id_shf  = CCPL_register_field_instance(shf, "shf", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="W/m2", annotation="register field instance of Net surface heat flux")
        field_id_ssh = CCPL_register_field_instance(ssh, "ssh", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="m", annotation="register field instance of Sea surface height")
        field_id_mld = CCPL_register_field_instance(mld, "mld", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="m", annotation="register field instance of Mixed layer depth")
        field_id_mld2 = CCPL_register_field_instance(mld, "mld2", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="m", annotation="register field instance of Mixed layer depth")
        !--------register coupling frequency to C-Coupler2-------------
        timer_id = CCPL_define_single_timer(atm_demo_comp_id, "seconds", coupling_freq, 0, 0, annotation="define a single timer for atm_demo")

        !--------register export & import interface to C-Coupler2------
        fields_id(1) = field_id_psl
        fields_id(2) = field_id_prect
        fields_id(3) = field_id_fsds
        fields_id(4) = field_id_flds
        fields_id(5) = field_id_flds2
        export_interface_id = CCPL_register_export_interface("send_data_to_ocn", 5, fields_id, timer_id, annotation="register interface for sending data to atmosphere")

        fields_id(1) = field_id_sst
        fields_id(2) = field_id_shf
        fields_id(3) = field_id_ssh
        fields_id(4) = field_id_mld
        fields_id(5) = field_id_mld2
        import_interface_id = CCPL_register_import_interface("receive_data_from_ocn", 5, fields_id, timer_id, 0, annotation="register interface for receiving data from atmosphere")
        
        call mpi_barrier(mpicomm, ierr)
        call cpu_time(start_time)
        call CCPL_end_coupling_configuration(atm_demo_comp_id, annotation = "component atm_demo ends configuration")       
        call cpu_time(end_time)
        total_time = real(end_time-start_time)
        call mpi_reduce(total_time, max_total_time, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
        
        if(mytask_id .eq. 0) then
            open(1, file="record_time", status='old', position="append")
            write(1, *) "atm coupling time ", max_total_time
            close(1)
        endif
    end subroutine register_component_coupling_configuration

end module coupling_ocn_model_mod
