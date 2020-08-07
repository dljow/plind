!! badedgesort.f95
!! 
!! Decides what are the bad edges that
!! will be corrected this round, based
!! on the simplices encountered. 

subroutine bad_edge_sort(edges, bad_edges, simplices, new_bad_edges, n_edges, n_bad, n_simp, n_dim)
! =====================================================
! Alters the edge array new_bad_edges in case of a duplicated simplex . 
! =====================================================
    implicit none

    integer*8, intent(in)   :: n_edges
    integer*8, intent(in)   :: n_bad
    integer*8, intent(in)   :: n_simp
    integer*8, integer(in)  :: n_dim

    integer*8, intent(in)   :: edges(n_edges,2)
    integer*8, intent(in)   :: simplices(n_simp, n_dim+1) 
    integer*8, intent(in)   :: bad_edges(n_bad)
    integer*8, intent(inout)  :: new_bad_edges(n_bad)
    integer*8 :: i
    integer*8 :: edge_ind
    integer*8 :: simp(n_dim+1)


    do i=1, n_bad
       
       do j= 1, n_simp
            simp=simplices(j)
            if (ANY(simp== edges(1, bad_edges(i))).or.ANY(simp== edges(1, bad_edges(i)))) 

            # To combat this, only the first edge is split, and the second is assumed to be caught
        # by the subsequent time step
        for i, bad_edge in enumerate(bad_edges):
            # Keep track of the simplices associated with an edge
            simplices_tag = np.isin(self.simplices, bad_edge).sum(axis=-1) > 1
            simplices_tag = np.where(simplices_tag)[0]

            if not np.any(np.in1d(simplices_tag, used_simps)):  # Check if we have used this simplex
                # Flag simplices_tag to not reuse simplices
                used_simps = np.append(used_simps, simplices_tag)
                # Add edge to new bad edge array
                uni_bad_edges = np.append(uni_bad_edges, bad_edge)
                # Remove bad edges from the list of edges
                edges_tag = np.isin(self.edges, bad_edge).sum(axis=-1) == 2
                self.edges = self.edges[~(edges_tag)]

                # Add simplice(s) with the proper extras populated
                uni_bad_simps = np.append(uni_bad_simps, self.simplices[simplices_tag], axis=0)



    return
end subroutine

