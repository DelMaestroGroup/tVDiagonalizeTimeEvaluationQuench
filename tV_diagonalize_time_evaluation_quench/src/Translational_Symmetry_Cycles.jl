"""
Create a list of occupation basis for each translational symmetry cycle for fermionic 1D chains with PBC/APBC For odd/even N.
"""
function Translational_Symmetry_Cycles(basis::AbstractSzbasis)

    for x = [:NumOfCycles,:Num_cycles_max,:MemberID,:i_next]
       @eval $x = Int64
    end

    IdStatus = Bool
    Num_cycles_max=round(Int64,basis.D/basis.K*1.5)
    Cycles= zeros(Int64, Num_cycles_max, basis.K)
    CycleSize = zeros(Int64, Num_cycles_max)
    Status  = trues(basis.D)
    NumOfCycles=0
    for (i, bra) in enumerate(basis)
        if Status[i]
            NumOfCycles+=1
            MemberID=0
            IdStatus=true
            i_next=i
            while IdStatus
                Status[i_next]=false
                MemberID+=1
                Cycles[NumOfCycles, MemberID]= i_next
                i_next=serial_num(basis, circshift(basis[i_next],1))
                IdStatus = Status[i_next]
            end
            CycleSize[NumOfCycles]= MemberID
        end
    end
    Cycles, CycleSize, NumOfCycles
end
