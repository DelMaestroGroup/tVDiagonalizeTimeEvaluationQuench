"""
Create the sparse translational, reflection and particle-hole symmetry block H_(q=0,R=1,P=1) of the hamiltonian of fermionic 1D chains with PBC/APBC.
"""
function sparse_Block_Diagonal_Hamiltonian_q0R1PH1(basis::AbstractFermionsbasis, Cycles:: Array{Int64,2}, CycleSize:: Vector{Int64}, NumOfCycles::Int64, InvCycles_Id:: Vector{Int64}, InvCycles_order:: Vector{Int64}, t::Float64, V::Float64)
    if basis.K!=2*basis.N
        warn("particle-hole symmetry works only at half-filling,", "  quit")
        quit()
    end
    const N=basis.N

    #Creating the block H_(q=0,R=1,P=1) of the hamiltonian.
    rows = Int64[]
    cols = Int64[]
    elements = Float64[]
    end_site = basis.K
    for CycleId =1: NumOfCycles
        # Diagonal part
        Vsum = 0.0
        bra=basis[Cycles[CycleId,1]]
        for j=1:end_site
            j_next = j % basis.K + 1
            Vsum += bra[j] * (bra[j_next])
        end
        push!(rows, CycleId)
        push!(cols, CycleId)
        push!(elements, Vsum*V)

        # Off-diagonal part

        CycleId_CycleSize=CycleSize[CycleId]
        for j=1:end_site
             j_next = j % basis.K + 1
             # Tunnel right, tunnel left.
             for (site1, site2) in [(j, j_next), (j_next, j)]
                 if bra[site1] > 0
                     ket = copy(bra)
                     ket[site1] -= 1
                     ket[site2] += 1
                     if dot(ket,ket)==N 
                         kId=serial_num(basis, ket)
                         CycleId1 = InvCycles_Id[kId]
                         kIdcy = InvCycles_order[kId]
                         CycleId1_CycleSize=CycleSize[CycleId1]
                         k1=CycleId1
                         factor=-t*sqrt(CycleId_CycleSize/CycleId1_CycleSize)*.5
                         push!(rows, k1)
                         push!(cols, CycleId)
                         push!(elements, factor)
                         push!(rows, CycleId)
                         push!(cols, k1)
                         push!(elements, conj(factor))
                     end
                 end
             end
         end
    end
        return sparse(rows, cols, elements, NumOfCycles, NumOfCycles), NumOfCycles
end