"""
Calculate the particle entanglement entropy for an eigenstate of the one-site-translation operator with eigenvalue q=0.
"""
function particle_entropy_Ts(basis::AbstractSzbasis, Asize::Int, d::Vector{Complex128}, MaxOccupation::Int, measure_obdm::Bool)

    SRen = Array{Float64}(3)
    DimA=Int64
    DimB=Int64
    DimAdA=Int64
    for x = [:fN,:Wa,:Wb,:Flips,:SumFlips]
       @eval $x = Int64
    end
    norm= Float64
    facto= Float64
    AdA_elem= Complex128
    L = basis.K
    N = basis.N
    Bsize = N - Asize
    if Asize>Bsize
        Asize,Bsize=Bsize,Asize
    end

    # Dimensions of partition Hilbert spaces
    facto=1.0
    for i=1:Asize
    	facto *=((1.0*L-i+1)/(1.0*i))
    end
    DimA = Int(round(facto))
    facto=1.0
    for i=1:Bsize
    	facto *=((1.0*L-i+1)/(1.0*i))
    end
    DimB = Int(round(facto))

    DimAdA= DimA
    const basisA = RestrictedSzbasis(L, Asize, MaxOccupation)
    const basisB = RestrictedSzbasis(L, Bsize, MaxOccupation)

    CyclesA, CycleSizeA, NumOfCyclesA =Translational_Symmetry_Cycles(basisA)
    CyclesB, CycleSizeB, NumOfCyclesB =Translational_Symmetry_Cycles(basisB)
    # using Int32 instead of Int64 in A matrix Structure works up to L=32 and N=16.
    AmatrixStructure=zeros(Int64,NumOfCyclesA, NumOfCyclesB,L)

    λ=Array{Float64}(NumOfCyclesA*L)

    braA = Array{Int}(L)
    braB = Array{Int}(L)
    bra = Array{Int}(L)

    # Weight factors
    Wa=factorial(Asize)
    Wb=factorial(Bsize)

    # Normalization coefficient
    norm=sqrt(Wa*Wb/factorial(basis.N))
    Aparity= Asize%2
    Bparity= Bsize%2
    element= Complex128

   # construct the AmatrixStructure
    OcupationOverlap =MaxOccupation+1
    for i=1: NumOfCyclesB
        braB=basisB[CyclesB[i,1]] 
        for j=1: NumOfCyclesA 
            minSize=CycleSizeA[j]
            maxSize=CycleSizeB[i]
            if CycleSizeA[j]>CycleSizeB[i]
                minSize,maxSize=maxSize, minSize
            end
            phase0=1
            for k=1: minSize             
                braA= basisA[CyclesA[j,k]]
                bra=braA+braB
                    # absorbing phase changes due to a particle crossing the system boundary.
                    if Bparity==1 && j>1
                        phase0*= 1-2*braA[1]
                    end
                Overlap =findfirst(bra, OcupationOverlap)
                if Overlap < 1 #(Overlap==0)
                    Flips=0
                    SumFlips=0
                    for Index=1:L
                       Flips = Flips +(1-2* Flips)* braA[Index]
                       SumFlips += Flips* braB[Index]
                    end
                    phase=(-1)^SumFlips*phase0
                    AmatrixStructure[j,i,k]=serial_num(basis, bra)* phase
                end
            end
        end
    end

    #find the spectrum of the reduced density matrix
    for q=0:L-1
        Amatrixq=zeros(Complex128, NumOfCyclesA, NumOfCyclesB)
        for i=1: NumOfCyclesB
            for j=1: NumOfCyclesA 
                minSize=CycleSizeA[j]
                maxSize=CycleSizeB[i]
                minparity= Aparity
                maxparity= Bparity
                minCycleparity= div(Asize* CycleSizeA[j],L)%2 
                maxCycleparity= div(Bsize* CycleSizeB[i],L)%2 
                if CycleSizeA[j]>CycleSizeB[i]
                    minSize,maxSize=maxSize, minSize
                    minparity, maxparity= maxparity ,minparity
                    minCycleparity, maxCycleparity= maxCycleparity, minCycleparity
                end
                shiftmin=div(div(L,minSize)-minparity,2)* maxparity* minCycleparity
                shiftmax=div(div(L,maxSize)-maxparity,2)* minparity* maxCycleparity
                element=0.0+0.0im
               if  (q-shiftmin)%div(L,minSize)==0 && (q-shiftmax)%div(L,maxSize)==0
                    halfe=0.5*maxparity * minCycleparity
                    factor1=norm*sqrt(maxSize)/sqrt(minSize)
                    phasefactor=(0.0+1.0im)*2*pi/minSize
                    qValue=(q-shiftmin)/div(L,minSize)
                    for k=1:minSize
                        serialnum= AmatrixStructure[j,i,k]
                        if serialnum != 0 
                            phase= sign(serialnum)
                            serialnum=abs(serialnum)
                            element += factor1*phase*d[serialnum]*exp((1.0+0.0im)*(k-1)*(qValue +halfe)* phasefactor)
                        end 
                    end
                end
                Amatrixq[j,i]= element
          end
       end
       S = svdvals!(Amatrixq)
       S.=S.^2
       λ[q*NumOfCyclesA+1:q*NumOfCyclesA+NumOfCyclesA]=S[1:NumOfCyclesA]
    end


   # construct the spatial OBDM
   if measure_obdm && Asize == 1
       obdm = zeros(Float64, DimA)
       phase=(0.0+1.0im)*2*pi/L
       shift=(0.5* Bparity)
       j_index=Int(L/2)
       for i=1:L
           AdA_elem=0.0+0.0im
           for k=1:L
               AdA_elem += exp(phase*(i-j_index)*(k-1+shift))*λ[k]/L
           end
           obdm[i]=real(AdA_elem)
       end
   end

    err = abs(sum(λ.^1) - 1.0)
    if err > 1e-12
        warn("RDM eigenvalue error ", err)
    end
#  println("  λ ", λ )

    LogNn=log(factorial(basis.N)/factorial(Bsize)/factorial(Asize))
    SRen[1]=0
    for k=1: NumOfCyclesA*L

        if λ[k]>0
            SRen[1] += λ[k]*log(λ[k])
	end
    end
    SRen[1]=-SRen[1]-LogNn
    SRen[2]=-log(sum(λ.^2))-LogNn
    SRen[3]=-log(sum(λ.^3))/2-LogNn

    if measure_obdm && Asize==1
        return SRen,obdm
    else
        return SRen
    end
end
