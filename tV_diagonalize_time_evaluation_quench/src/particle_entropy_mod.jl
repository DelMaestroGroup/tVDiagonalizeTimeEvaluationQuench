"""
Calculate the particle entanglement entropy for a subset A, using the SVD on (A^*)A First test.
"""
function particle_entropy_mod(basis::AbstractSzbasis, Asize::Int, d::Vector{Complex128}, MaxOccupation::Int)
    SRen = Array(Float64, 3)
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

    DimAdA= min(DimA,DimB)
    const basisA = RestrictedSzbasis(L, Asize, MaxOccupation)
    const basisB = RestrictedSzbasis(L, Bsize, MaxOccupation)

    braA = Array(Int, L)
    braB = Array(Int, L)
    bra = Array(Int, L)
    # Matrix A
    Amatrix = zeros(Complex128, DimA, DimB)
    # Weight factors
    Wa=factorial(Asize)
    Wb=factorial(Bsize)
    # Normalization coefficient
    norm=sqrt(Wa*Wb/factorial(basis.N))


    for (i, braA) in enumerate(basisA)
        for (j, braB) in enumerate(basisB)
	    bra=braA+braB
	    if bra in basis
		Flips=0
		SumFlips=0
		for k=1:L
		    Flips = Flips +(1-2* Flips)* braA[k]
		    SumFlips += Flips* braB[k]
		end
		Amatrix[i, j] = (-1)^ SumFlips*norm*d[serial_num(basis, bra)]
	    end
	end
    end

    AdAmatrix = zeros(Complex128, DimAdA, DimAdA)
#print("sts\n")  
#    for i=1:DimA
#        for j=i:DimA
#	    AdA_elem=0.0
#            for k=1:DimB
#               AdA_elem += Amatrix[i,k]* conj(Amatrix[j,k])
#	    end
#	    AdAmatrix[i,j]= AdA_elem
#	    AdAmatrix[j,i]= conj(AdA_elem)
#        end
#    end
#print("end\n")  
#open( "Corr_pbcL22N11.dat", "w") do f
#   write(f, "# i-j Corr_cal  Corr_th\n")
#    for i=1:DimA
#        Corr=0 
#        for j=1:Int((L-2)/4)
#	    Corr+=cos(2*pi*j/L*(i-Int(L/2)))
#        end
#        Corr=(2*Corr+1)*2/L/L
#        write(f, "$(i-Int(L/2)) $(AdAmatrix[i,Int(L/2)]) $(Corr)\n")
#    end
#flush(f)
#end
    #S = svdvals!(AdAmatrix)
    S = svdvals!(Amatrix)
    S.=S.^2
    #print(S,"\n")
    err = abs(sum(S.^1) - 1.0)
    if err > 1e-12
        warn("RDM eigenvalue error ", err)
    end
    LogNn=log(factorial(basis.N)/factorial(Bsize)/factorial(Asize))
    SRen[1]=0
    for k=1:DimAdA
        if S[k]>0
            SRen[1] +=S[k]*log(S[k])
	end
    end
    SRen[1]=-SRen[1]-LogNn
    SRen[2]=-log(sum(S.^2))-LogNn
    SRen[3]=-log(sum(S.^3))/2-LogNn
    #warn(" AdA eigenvalue are", S)
    return SRen
end
