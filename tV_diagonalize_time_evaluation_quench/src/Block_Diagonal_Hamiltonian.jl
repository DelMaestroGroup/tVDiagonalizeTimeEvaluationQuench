"""
Create the translational symmetry block H_(q) of the hamiltonian of fermionic 1D chains with PBC/APBC.
"""
function Block_Diagonal_Hamiltonian(basis::AbstractSzbasis, Cycles:: Array{Int64,2}, CycleSize:: Vector{Int64}, NumOfCycles::Int64, t::Float64, V::Float64, q::Int64)

    InvolvedCycles = Int64[]
    NumberOfInvolvedCycles = Int64
    #Finding the involved cycles for q. 
    NumberOfInvolvedCycles=0
    for i=1: NumOfCycles
        if  q%(basis.K /CycleSize[i])==0
            NumberOfInvolvedCycles+=1
            push!(InvolvedCycles, i)
        end
    end
    #Creating the block H_(q) of the hamiltonian.
    Hq=zeros(Complex128, NumberOfInvolvedCycles, NumberOfInvolvedCycles)
    exp_q=zeros(Complex128, basis.K)
    end_site = basis.K
    for i=1: basis.K
        exp_q[i]=-t*exp((i-1)*(0.0-1.0im)*2*pi*q/basis.K)*0.5
    end
    for (i, CycleId) in enumerate(InvolvedCycles)
        CycleSizeCycleId =CycleSize[CycleId]
        # Diagonal part
        Vsum = 0.0+0.0im

        bra=basis[Cycles[CycleId,1]]

        for j=1:end_site
            j_next = j % basis.K + 1
            Vsum += bra[j] * (bra[j_next])
        end
        Hq[i,i]= Vsum*V
        # Off-diagonal part
        for j=1:end_site
            j_next = j % basis.K + 1
            # Tunnel right, tunnel left.
            for (site1, site2) in [(j, j_next), (j_next, j)]
                if bra[site1] > 0
                    ket = copy(bra)
                    ket[site1] -= 1
                    ket[site2] += 1
                    if ket in basis
                        kId=serial_num(basis, ket)
                        for (k1, CycleId1) in enumerate(InvolvedCycles)
                            kIdcy =findfirst(Cycles[CycleId1,:], kId)
                            if kIdcy>0
                                factor=exp_q[kIdcy]*sqrt(CycleSizeCycleId/CycleSize[CycleId1])
                                Hq[k1, i]+= factor
                                Hq[i, k1]+=conj(factor)

                            end
                        end
                    end
                end
            end
        end
    end
    Hq, InvolvedCycles,NumberOfInvolvedCycles
end