"""
Create a list of occupation basis for each translational symmetry cycle for fermionic 1D chains with PBC/APBC.
"""
function Block_Diagonal_Hamiltonian_ParticleHoleSym(basis::AbstractSzbasis, Cycles:: Array{Int64,2}, CycleSize:: Vector{Int64}, NumOfCycles::Int64, t::Float64, V::Float64, q::Int64, P::Int64)
    if abs(P)!=1
        warn("abs(P)!=1, ", "  quit")
        quit()
    end
    InvolvedCycles = Int64[]
    NumberOfInvolvedCycles = Int64
    exp_q=zeros(Complex128, basis.K)
    for i=1: basis.K
        exp_q[i]=exp((i-1)*(0.0-1.0im)*2*pi*q/basis.K)
    end
    if ((q<0) | (q> basis.K-1))
        warn("q does not satisfy the condition: ",0,"<q<", basis.K-1 ,",", "  quit")
        quit()
    end
#    if q*(q-basis.K/2)!=0

        NumberOfInvolvedCycles=0
        flip=false
        for i=1: NumOfCycles
            if  q%(basis.K /abs(CycleSize[i]))==0
                if flip
                    flip=~flip
                    continue
                end
                if CycleSize[i]<0
                    NumberOfInvolvedCycles+=1
                    push!(InvolvedCycles, i)
                    flip=~flip
                elseif P==(-1)^(q/(basis.K /abs(CycleSize[i])))
                    NumberOfInvolvedCycles+=1
                    push!(InvolvedCycles, i)
                end 
            end
        end

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
            CycleId_CycleSize = CycleId_CycleSize>0 ? CycleId_CycleSize : 2* CycleId_CycleSize
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
                                CycleId1_CycleSize=CycleSize[CycleId1]
                                CycleId1_CycleSize = CycleId1_CycleSize>0 ? CycleId1_CycleSize : 2* CycleId1_CycleSize
                                if kIdcy>0
                                    phase=exp_q[kIdcy]
                                    Hq[k1, i]+=-t* phase*sqrt(abs(CycleId_CycleSize)/abs(CycleId1_CycleSize))/2.0
                                    Hq[i, k1]+=conj(-t* phase*sqrt(abs(CycleId_CycleSize)/abs(CycleId1_CycleSize)))/2.0
                                elseif CycleId1_CycleSize<0
                                    kIdcy =findfirst(Cycles[CycleId1+1,:], kId)
                                    if kIdcy>0
                                        phase=exp_q[kIdcy]*P
                                        Hq[k1, i]+=-t* phase*sqrt(abs(CycleId_CycleSize)/abs(CycleId1_CycleSize))/2.0
                                        Hq[i, k1]+=conj(-t* phase*sqrt(abs(CycleId_CycleSize)/abs(CycleId1_CycleSize)))/2.0
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
"""
    elseif q==0
        NumberOfInvolvedCycles=0
        flip=false
        for i=1: NumOfCycles
            if flip
                flip=~flip
                continue
            end
            if CycleSize[i]<0
                NumberOfInvolvedCycles+=1
                push!(InvolvedCycles, i)
                flip=~flip
            elseif P==1
                NumberOfInvolvedCycles+=1
                push!(InvolvedCycles, i)
            end 
        end
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
            CycleId_CycleSize = CycleId_CycleSize>0 ? CycleId_CycleSize : 2* CycleId_CycleSize
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
                                CycleId1_CycleSize=CycleSize[CycleId1]
                                CycleId1_CycleSize = CycleId1_CycleSize>0 ? CycleId1_CycleSize : 2* CycleId1_CycleSize
                                if kIdcy>0
                                    phase=exp_q[kIdcy]
                                    Hq[k1, i]+=-t* phase*sqrt(abs(CycleId_CycleSize)/abs(CycleId1_CycleSize))/2.0
                                    Hq[i, k1]+=conj(-t* phase*sqrt(abs(CycleId_CycleSize)/abs(CycleId1_CycleSize)))/2.0
                                elseif CycleId1_CycleSize<0
                                    kIdcy =findfirst(Cycles[CycleId1+1,:], kId)
                                    if kIdcy>0
                                        phase=exp_q[kIdcy]*P
                                        Hq[k1, i]+=-t* phase*sqrt(abs(CycleId_CycleSize)/abs(CycleId1_CycleSize))/2.0
                                        Hq[i, k1]+=conj(-t* phase*sqrt(abs(CycleId_CycleSize)/abs(CycleId1_CycleSize)))/2.0
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    elseif (q-basis.K /2) ==0
        NumberOfInvolvedCycles=0
        flip=false
        for i=1: NumOfCycles
            if flip
                flip=~flip
                continue
            end
            if  q%(basis.K /abs(CycleSize[i]))==0
                if CycleSize[i]<0
                    NumberOfInvolvedCycles+=1
                    push!(InvolvedCycles, i)
                    flip=~flip
                else
                    bra=basis[Cycles[i,1]]
                    if P==(-1)^(findfirst(Cycles[i,:], serial_num(basis,reverse(basis[Cycles[i,1]])))+1)
                        NumberOfInvolvedCycles+=1
                        push!(InvolvedCycles, i)
                    end 
                end 
            end
        end
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
            CycleId_CycleSize = CycleId_CycleSize>0 ? CycleId_CycleSize : 2* CycleId_CycleSize
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
                                CycleId1_CycleSize=CycleSize[CycleId1]
                                CycleId1_CycleSize = CycleId1_CycleSize>0 ? CycleId1_CycleSize : 2* CycleId1_CycleSize
                                if kIdcy>0
                                    phase=exp_q[kIdcy]
                                    Hq[k1, i]+=-t* phase*sqrt(abs(CycleId_CycleSize)/abs(CycleId1_CycleSize))/2.0
                                    Hq[i, k1]+=conj(-t* phase*sqrt(abs(CycleId_CycleSize)/abs(CycleId1_CycleSize)))/2.0
                                elseif CycleId1_CycleSize<0
                                    kIdcy =findfirst(Cycles[CycleId1+1,:], kId)
                                    if kIdcy>0
                                        phase=exp_q[kIdcy]*P
                                        Hq[k1, i]+=-t* phase*sqrt(abs(CycleId_CycleSize)/abs(CycleId1_CycleSize))/2.0
                                        Hq[i, k1]+=conj(-t* phase*sqrt(abs(CycleId_CycleSize)/abs(CycleId1_CycleSize)))/2.0
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
"""
    Hq, InvolvedCycles,NumberOfInvolvedCycles
end