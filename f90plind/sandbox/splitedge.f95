!! splitedge.f90
!! 
!! A function that splits the edges of a given 
!! array.
!!
!! INPUTS:
!!
!! edge_array: array of edge tuples
!!
!! bad_edge_indices: array of the indices to be split.
!!

subroutine split_edge(edges, bad_edges, n_edges, n_bad, n_pts, new_edges)
! =====================================================
! Returns an edge array which has the split edges. 
! =====================================================
    integer*8, intent(in)   :: n_edges
    integer*8, intent(in)   :: n_bad 
    integer*8, intent(in)   :: n_pts

    integer*8, intent(in)   :: edges(2, n_edges)
    integer*8, intent(in)   :: bad_edges(n_bad)
    integer*8, intent(inout)  :: new_edges(2,n_edges+n_bad)
    integer*8 :: i
    integer*8 :: edge_ind
    integer :: count_new_edge

    edge_ind=1
    count_new_edge=0
    do i = 1, n_edges
        if (ANY(bad_edges == i)) then
           count_new_edge=count_new_edge+ 1 
           
           ! Add the edge connection the old first point to the midpoint
           new_edges(1, edge_ind) = edges(1, i)
           new_edges(2, edge_ind)= n_pts+ count_new_edge
           edge_ind=edge_ind+1

           ! Add the edge connecting the old second point to the midpoint
           new_edges(1, edge_ind) = edges(2, i)
           new_edges(2, edge_ind)= n_pts+ count_new_edge

           edge_ind= edge_ind+1

         else 

         new_edges(1, edge_ind) = edges(1, i)
         new_edges(2, edge_ind) = edges(2, i)

         edge_ind=edge_ind+1
         end if

    end do
    return
end subroutine

