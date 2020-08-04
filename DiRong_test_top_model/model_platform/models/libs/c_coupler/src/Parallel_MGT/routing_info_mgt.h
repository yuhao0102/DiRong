/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu and then
  *  modified by Dr. Cheng Zhang and Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn or
  *  Dr. Cheng Zhang via zhangc-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef ROUTER_MGT_H
#define ROUTER_MGT_H


#include "common_utils.h" 
#include "decomp_info_mgt.h"
#include "compset_communicators_info_mgt.h"
#include "quick_sort.h"
#include <vector>

void record_time(double, char*);

struct routing_mapping_table_entry
{
    int key;
    int global_index;         
    int local_index;          
    int local_process_id;       
    int remote_local_index;   
    int remote_local_process_id;

    bool operator<(const routing_mapping_table_entry b) const { return this->key < b.key; }
    bool operator>(const routing_mapping_table_entry b) const { return this->key > b.key; }   
    bool operator==(const routing_mapping_table_entry b) const { return this->key == b.key; }
};


struct Routing_info_with_one_process
{
    int remote_proc_global_id;
    int num_elements_transferred;
    int num_local_indx_segments;
    int *local_indx_segment_starts;
    int *local_indx_segment_lengths;
    bool send_or_recv;                           // true is send and false is recv
};


class Routing_info
{
    private:
        int src_comp_id;
        int dst_comp_id;
        char src_comp_full_name[NAME_STR_SIZE];
        char index_dst_comp_full_name[NAME_STR_SIZE];
        char true_dst_comp_full_name[NAME_STR_SIZE];
        char src_decomp_name[NAME_STR_SIZE];
        char dst_decomp_name[NAME_STR_SIZE];
        int src_decomp_size;
        int dst_decomp_size;
        Decomp_info *src_decomp_info;
        Decomp_info *dst_decomp_info;
        Comp_comm_group_mgt_node * src_comp_node;
        Comp_comm_group_mgt_node * dst_comp_node;
        int src_comp_node_id;
        int dst_comp_node_id;
        int current_proc_id_src_comp;
        int current_proc_id_dst_comp;
        char remote_comp_name[NAME_STR_SIZE];
        char local_decomp_name[NAME_STR_SIZE];
        char remote_decomp_name[NAME_STR_SIZE];
        int num_dimensions;
        int total_num_transferred_cells;
        long local_decomp_size;
        long remap_decomp_size;
        int *logical_indx_lookup_table_local;
        int *logical_indx_lookup_table_remote;
        std::vector<Routing_info_with_one_process *> recv_from_remote_procs_routing_info;
        std::vector<Routing_info_with_one_process *> send_to_remote_procs_routing_info;

        routing_mapping_table_entry *src_local_routing_mapping_table_entries;
		routing_mapping_table_entry *dst_local_routing_mapping_table_entries;
        routing_mapping_table_entry *src_remote_routing_mapping_table_entries;
		routing_mapping_table_entry *dst_remote_routing_mapping_table_entries;
        int num_src_local_routing_mapping_table_entries;
		int num_dst_local_routing_mapping_table_entries;
        int num_src_remote_routing_mapping_table_entries;
		int num_dst_remote_routing_mapping_table_entries;
		
public:
        Routing_info(const int, const int, const char*, const char*);
        ~Routing_info();
        Routing_info_with_one_process *get_routing_info(bool, int);
        int get_num_elements_transferred_with_remote_proc(bool is_send, int i) { return get_routing_info(is_send,i)->num_elements_transferred; }
        int *get_local_indx_segment_starts_with_remote_proc(bool is_send, int i) { return get_routing_info(is_send,i)->local_indx_segment_starts; }
        int *get_local_indx_segment_lengths_with_remote_proc(bool is_send, int i) { return get_routing_info(is_send,i)->local_indx_segment_lengths; }
        int get_num_local_indx_segments_with_remote_proc(bool is_send, int i) { return get_routing_info(is_send,i)->num_local_indx_segments; }
        bool match_router(const int, const int, const char*, const char*);
        int get_num_dimensions() { return num_dimensions; }
        long get_local_decomp_size() { return local_decomp_size; }
        long get_remap_decomp_size() { return remap_decomp_size; }
        Comp_comm_group_mgt_node *get_src_comp_node();
        Comp_comm_group_mgt_node *get_dst_comp_node();
        long get_src_decomp_size() { return src_decomp_size; }
        long get_dst_decomp_size() { return dst_decomp_size; }
        
    private:
        void intersect_routing_mapping_tables_between_components(routing_mapping_table_entry**, routing_mapping_table_entry*, int*, int, bool);
        void initialize_routing_mapping_table(routing_mapping_table_entry**, Decomp_info*, Comp_comm_group_mgt_node*, int*, int, int, int );
        void exchange_routing_mapping_tables_between_components(int, int, Decomp_info*, bool);
        void sort_routing_mapping_table_entries(routing_mapping_table_entry**, Decomp_info*, Comp_comm_group_mgt_node*, int*, int, int, int, int, int);
        void output_routing_mapping_table(routing_mapping_table_entry*, char*, int, int);
        void check_routing_mapping_table_entries(routing_mapping_table_entry*, int, int, int, char*);
        void calculate_min_max_grid_index(int, int, int, int, int, int*, int*);	
		void exchange_routing_mapping_tables_between_processes(int, routing_mapping_table_entry **, routing_mapping_table_entry **, int *, int *, int *, int *, MPI_Comm comm, int *);
        int calculate_max_power2(int);
        void merge_routing_mapping_tables(routing_mapping_table_entry**, routing_mapping_table_entry*, int, int, int, int, int*);
        void calculate_routing_mapping_tables();
		void build_router_based_on_routing_mapping_tables();
        void build_2D_router();
        Routing_info_with_one_process *compute_routing_info_between_decomps(int, const int*, int, const int*, int, int, int, bool, routing_mapping_table_entry*, int);		
		Routing_info_with_one_process *generate_routing_info_between_procs(routing_mapping_table_entry *, int, int, bool);
};


class Routing_info_mgt
{
    private:
        std::vector<Routing_info *> routers;
    
    public:
        Routing_info_mgt() {}
        ~Routing_info_mgt();
        Routing_info *search_router(const int, const int, const char*, const char*);
        Routing_info *search_or_add_router(const int, const int, const char*, const char*);
};

#endif

