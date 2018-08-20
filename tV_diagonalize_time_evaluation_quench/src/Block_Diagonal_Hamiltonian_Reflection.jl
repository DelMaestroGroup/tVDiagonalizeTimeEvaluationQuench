"""
Create the translational-reflection symmetry block H_(q,R) of the hamiltonian of fermionic 1D chains with PBC/APBC.
"""
function Block_Diagonal_Hamiltonian_Reflection(basis::AbstractSzbasis, Cycles:: Array{Int64,2}, CycleSize:: Vector{Int64}, NumOfCycles::Int64, CReflection:: Vector{Bool}, t::Float64, V::Float64, q::Int64, R::Int64)
    if abs(R)!=1
        warn("abs(R)!=1, ", "  quit")
        quit()
    end
    InvolvedCycles = Int64[]
    NumberOfInvolvedCycles = Int64
    if q*(q-basis.K/2)!=0
        warn("q does not satisfy the condition: q= 0 OR q= ", basis.K/2 ,",", "  quit")
        quit()
    end
    exp_q=zeros(Complex128, basis.K)
    for i=1: basis.K
        exp_q[i]=-t*exp((i-1)*(0.0-1.0im)*2*pi*q/basis.K)*0.5
    end
    #Finding the involved cycles for q and R. 
    NumberOfInvolvedCycles=0
    flip=false
    if q==0
        for i=1: NumOfCycles
            if flip
                flip=~flip
                continue
            end
            if CReflection[i]
                NumberOfInvolvedCycles+=1
                push!(InvolvedCycles, i)
                flip=~flip
            elseif R==1
                NumberOfInvolvedCycles+=1
                push!(InvolvedCycles, i)
            end 
        end
    elseif (q-basis.K /2) ==0
        for i=1: NumOfCycles
            if flip
                flip=~flip
                continue
            end
            if  q%(basis.K /abs(CycleSize[i]))==0
                if CReflection[i]
                    NumberOfInvolvedCycles+=1
                    push!(InvolvedCycles, i)
                    flip=~flip
                else
                    bra=basis[Cycles[i,1]]
                    if R==(-1)^(findfirst(Cycles[i,:], serial_num(basis,reverse(basis[Cycles[i,1]])))+1)
                        NumberOfInvolvedCycles+=1
                        push!(InvolvedCycles, i)
                    end 
                end 
            end
        end
    end
    #Creating the block H_(q,R) of the hamiltonian.
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
        factor_CycleId = CReflection[CycleId] ? 1/sqrt(2)  : 1.0
        flip=true
        flip1=CReflection[CycleId]
        x=1
        y=R
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
                                    factor_CycleId1 = CReflection[CycleId1] ? 1/sqrt(2)  : 1.0
                                    factor=exp_q[kIdcy]*x*sqrt(CycleId_CycleSize/CycleId1_CycleSize)* factor_CycleId* factor_CycleId1
                                    Hq[k1, i]+= factor
                                    Hq[i, k1]+=conj(factor)
                                elseif CReflection[CycleId1] 
                                    kIdcy =findfirst(Cycles[CycleId1+1,:], kId)
                                    if kIdcy>0
                                        CycleId1_CycleSize=CycleSize[CycleId1]
                                        factor_CycleId1 = CReflection[CycleId1] ? 1/sqrt(2)  : 1.0
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
                x=R
                y=1
            end
        end
    end
    Hq, InvolvedCycles,NumberOfInvolvedCycles
end