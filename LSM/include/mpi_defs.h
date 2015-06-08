
#ifdef MPILIB_DOUBLE_UNDERSCORE

#define MPI_SendRecv mpi_sendrecv_
#define MPI_Type_Free mpi_type_free_
#define MPI_Allreduce mpi_allreduce_
#define MPI_Reduce mpi_reduce_
#define MPI_BCAST mpi_bcast_
#define MPI_AllGatherV mpi_allgatherv_
#define MPI_GatherV mpi_gatherv_
#define MPI_Type_vector mpi_type_vector_
#define MPI_Type_extent mpi_type_extent_
#define MPI_Type_struct mpi_type_struct_
#define MPI_Type_Commit mpi_type_commit_
#define MPI_COMM_RANK mpi_comm_rank_
#define mpi_comm_rank mpi_comm_rank_
#define MPI_ALLTOALLV mpi_alltoallv_
#define MPI_ScatterV mpi_scatterv_
#define MPI_BARRIER mpi_barrier_

#define MPI_NULL_COPY_FN mpi_null_copy_fn_
#define MPI_NULL_DELETE_FN mpi_null_delete_fn_
#define MPI_COMM_NULL_COPY_FN mpi_comm_null_copy_fn_
#define MPI_COMM_NULL_DELETE_FN mpi_comm_null_delete_fn_
#define MPI_TYPE_NULL_COPY_FN mpi_type_null_copy_fn_
#define MPI_TYPE_NULL_DELETE_FN mpi_type_null_delete_fn_
#define MPI_WIN_NULL_COPY_FN mpi_win_null_copy_fn_
#define MPI_WIN_NULL_DELETE_FN mpi_win_null_delete_fn_
#define MPI_DUP_FN mpi_dup_fn_
#define MPI_COMM_DUP_FN mpi_comm_dup_fn_
#define MPI_TYPE_DUP_FN mpi_type_dup_fn_
#define MPI_WIN_DUP_FN mpi_win_dup_fn_
#define mpi_abort mpi_abort_
#define mpi_finalize mpi_finalize_
#define MPI_WTIME MPI_WTIME_
#define MPI_WTICK MPI_WTICK_
#define PMPI_WTIME PMPI_WTIME_
#define PMPI_WTICK PMPI_WTICK_

#define MPI_Send mpi_send_
#define MPI_Recv mpi_recv_

#define MPI_Allgather mpi_allgather_
#define MPI_ABORT mpi_abort_

#define mpi_comm_size mpi_comm_size_
#define mpi_gather mpi_gather_
#define mpi_gatherV mpi_gatherv_
#define mpi_GatherV mpi_gatherv_
#define mpi_ScatterV mpi_scatterv_

#define MPI_INIT mpi_init_
#define MPI_COMM_SIZE mpi_comm_size_

#define MPI_Bcast mpi_bcast_
#define MPI_Barrier mpi_barrier_
#define MPI_Finalize mpi_finalize_

#endif

