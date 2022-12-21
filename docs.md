hrep(A::AbstractMatrix, b::AbstractVector, linset::BitSet=BitSet())
Creates an H-representation for the polyhedron defined by the inequalities 
$$\langle A_i, x \rangle = b_i$$ if `i in linset`

and 
$$\langle A_i, x \rangle \le b_i$$
otherwise where $A_i$ is the $i$th row of $A$, i.e. $A[i,:]$ and $$b_i$$ is $b[i]$.
