"""
Create the translational-particle-hole symmetry block H_(q,P) of the hamiltonian of  fermionic 1D chains with PBC/APBC.
"""
function Block_Diagonal_Hamiltonian_PH(basis::AbstractSzbasis, Cycles:: Array{Int64,2}, CycleSize:: Vector{Int64}, NumOfCycles::Int64, CPH:: Vector{Bool}, t::Float64, V::Float64, q::Int64, P::Int64)
    if abs(P)!=1
        warn("abs(P)!=1, ", "  quit")
        quit()
    end
    InvolvedCycles = Int64[]
    NumberOfInvolvedCycles = Int64
    if basis.K!=2*basis.N
        warn("particle-hole symmetry works only at half-filling,", "  quit")
        quit()
    end
    exp_q=zeros(Complex128, basis.K)
    for i=1: basis.K
        exp_q[i]=exp((i-1)*(0.0-1.0im)*2*pi*q/basis.K)
    end
    #Finding the involved cycles for q and P. 
    NumberOfInvolvedCycles=0
    flip=false
        for i=1: NumOfCycles
            if flip
                flip=~flip
                continue
            end
            if  q%(basis.K /abs(CycleSize[i]))==0
                if CPH[i]
                    NumberOfInvolvedCycles+=1
                    push!(InvolvedCycles, i)
                    flip=~flip
                else
                    bra=basis[Cycles[i,1]]
                    z= exp_q[findfirst(Cycles[i,:], serial_num(basis,1 .-bra))] 
                    if P==sign(real(z))
                        NumberOfInvolvedCycles+=1
                        push!(InvolvedCycles, i)
                    end 
                end 
            end
        end
    for i=1: basis.K
        exp_q[i]=-t*exp_q[i]*0.5
    end
    #Creating the block H_(q,P) of the hamiltonian.
    Hq=zeros(Complex128, NumberOfInvolvedCycles, NumberOfInvolvedCycles)
    end_site = basis.K
    for (i, CycleId) in enumerate(InvolvedCycles)
        # Diagonal part
        Vsum = 0.0+0.0im
        bra=basis[Cycles[CycleId,1]]
        for j=1:end_site
            j_next = j % basis.K + 1
            Vsum += bra[j] * (bra[j_next])
        end
        Hq[i,i]= Vsum*V
        # Off-diagonal part
        CycleId_CycleSize=CycleSize[CycleId]
        factor_CycleId = CPH[CycleId] ? 1/sqrt(2)  : 1.0
        flip=true
        flip1= CPH[CycleId]
        x=1
        y=P
        while flip
            flip= ~flip
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
                                    CycleId1_CycleSize=CycleSize[CycleId1]
                                    factor_CycleId1 = CPH[CycleId1] ? 1/sqrt(2)  : 1.0
                                    factor=exp_q[kIdcy]*x*sqrt(CycleId_CycleSize/CycleId1_CycleSize)* factor_CycleId* factor_CycleId1
                                    Hq[k1, i]+= factor
                                    Hq[i, k1]+=conj(factor)
                                elseif CPH[CycleId1] 
                                    kIdcy =findfirst(Cycles[CycleId1+1,:], kId)
                                    if kIdcy>0
                                        CycleId1_CycleSize=CycleSize[CycleId1]
                                        factor_CycleId1 = CPH[CycleId1] ? 1/sqrt(2)  : 1.0
                                        factor =exp_q[kIdcy]*y*sqrt(CycleId_CycleSize/CycleId1_CycleSize)* factor_CycleId* factor_CycleId1
                                        Hq[k1, i]+= factor
                                        Hq[i, k1]+=conj(factor)
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if flip1
                bra=basis[Cycles[CycleId+1,1]]
                flip=~flip
                flip1=~flip1
                x=P
                y=1
            end
        end
    end
    Hq, InvolvedCycles,NumberOfInvolvedCycles
end