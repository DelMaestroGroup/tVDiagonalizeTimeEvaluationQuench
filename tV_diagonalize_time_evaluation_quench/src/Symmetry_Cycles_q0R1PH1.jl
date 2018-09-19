"""
Create a list of occupation basis for each translational, reflection and particle-hole symmetry cycle for fermionic 1D chains.
"""
function Symmetry_Cycles_q0R1PH1(basis::AbstractSzbasis)

    if basis.K!=2*basis.N
        warn("particle-hole symmetry works only at half-filling,", "  quit")
        quit()
    end

    for x = [:NumOfCycles,:Num_cycles_max,:MemberID,:i_next]
       @eval $x = Int64
    end
    ll=length(basis)
    IdStatus = Bool
    
    #Num_cycles_max=round(Int64,basis.D/basis.K*1.5/4) #This is a rough upper bound on the number of cycles. It could me smaller.
    Num_cycles_max=round(Int64,basis.D/basis.K*1.5) #This is a rough upper bound on the number of cycles. It could me smaller.

    Cycles= zeros(Int64, Num_cycles_max, basis.K*4)
    CycleSize = zeros(Int64, Num_cycles_max)
    CReflection= zeros(Bool, Num_cycles_max)
    InvCycles_Id = zeros(Int64, ll)
    InvCycles_order = zeros(Int64, ll)

    Status  = trues(basis.D)
    NumOfCycles=0
    MemberID=0
    for (i, bra) in enumerate(basis)
        i_next=i
        j=4
        while j >0
            if Status[i_next]
                if j==4
                   NumOfCycles+=1
                   MemberID=0
                end
                IdStatus=true
                while IdStatus
                    Status[i_next]=false
                    MemberID+=1
                    Cycles[NumOfCycles, MemberID]= i_next
                    InvCycles_Id[i_next]= NumOfCycles
                    InvCycles_order[i_next]= MemberID
                    i_next=serial_num(basis, circshift(basis[i_next],1))
                    IdStatus = Status[i_next]
                end
            else
                j=0
            end
            j-=1
            if j==3
                i_next=serial_num(basis,reverse(bra))
                if ~Status[i_next]
                    i_next=serial_num(basis,1 .-bra) 
                    if ~Status[i_next]
                        j=0
                    else
                        j=1
                    end
                end                
            end
            if j==2
                i_next=serial_num(basis,1 .-bra) 
                if ~Status[i_next]
                    j=0
                else
                    j=2
                end                
            elseif j==1
                i_next=serial_num(basis,reverse(1 .-bra))
            end
            if j==0
                CycleSize[NumOfCycles]= MemberID
            end
        end
    end                
    return Cycles, CycleSize, NumOfCycles, InvCycles_Id, InvCycles_order
end
