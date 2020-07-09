!! splitedgepoints.f90
!! 
!! A function that splits the edges of a given array and 
!! resets the points of the array accordingly.
!!
!! Function currently runs as O(n_edges) linearly.  
!!
!! INPUTS:
!!
!! edges: original array of the edges
!!
!! bad_edges: array of indices of the edges to be split
!!
!! points: original array of the points
!!
!! new_points: new array of the split points. It should be populated
!!            from 0 to n_pts with existing edges. 
!!
!! new_edges: new array of the split edges. 
!!
!!

subroutine split_edge_points(edges, bad_edges, new_edges, new_points, n_edges, n_bad, n_dim, n_pts)
! =====================================================
! Alters a new edge and a new points array. 
! =====================================================
    implicit none 

    integer*8, intent(in)   :: n_edges
    integer*8, intent(in)   :: n_bad 
    integer*8, intent(in)   :: n_pts
    integer*8, intent(in)   :: n_dim

    integer*8, intent(in)   :: edges(2, n_edges)
    integer*8, intent(in)   :: bad_edges(n_bad)
    integer*8, intent(inout)  :: new_edges(2,n_edges+n_bad)
    complex*16, intent(inout) :: new_points(n_dim, n_pts+n_bad)
    integer*8 :: edge_ind
    integer*8 :: i, j, k
    integer*8 :: count_new_edge
    !complex*8 :: new_midpoint

    edge_ind=1
    count_new_edge=0
    do i = 1, n_edges
        if (ANY(bad_edges == i)) then

           count_new_edge=count_new_edge+ 1
           do j= 1, n_dim
             !new_midpoint= (points(j, edges(1,i))+points(j, edges(2,i))/2
             new_points(j, n_pts+count_new_edge)= (new_points(j, edges(1,i))+new_points(j, edges(2,i)))/2
           enddo
            
           
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

