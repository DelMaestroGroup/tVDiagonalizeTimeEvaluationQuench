# Renyi entanglement entropy of Bose-Hubbard chains in 1D.

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "src"))
#using BoseHubbardDiagonalize
using tVDiagonalize

using ArgParse
using JeszenszkiBasis

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
    "M"
        help = "number of sites"
        arg_type = Int
        required = true
    "N"
        help = "number of particles"
        arg_type = Int
        required = true
    "--out"
        metavar = "FILE"
        help = "path to output file"
        required = true
    "--site-max"
        metavar = "N"
        help = "site occupation restriction"
        arg_type = Int
end
add_arg_group(s, "boundary conditions")
@add_arg_table s begin
    "--pbc"
        help = "periodic boundary conditions (default)"
        arg_type = BdryCond
        action = :store_const
        dest_name = "boundary"
        constant = PBC
        default = PBC
    "--obc"
        help = "open boundary conditions"
        arg_type = BdryCond
        action = :store_const
        dest_name = "boundary"
        constant = OBC
        default = PBC
end
add_arg_group(s, "tV parameters")
@add_arg_table s begin
    "--V0"
        metavar = "V0"
        help = "initial V"
        arg_type = Float64
        default = 20.0
    "--V"
        metavar = "V"
        help = "final V"
        arg_type = Float64
        default = 0.0
    "--t"
        metavar = "t"
        help = "t value"
        arg_type = Float64
        default = 1.0
end
add_arg_group(s, "time parameters")
@add_arg_table s begin
    "--time-min"
        metavar = "time"
        help = "minimum time"
        arg_type = Float64
        default = 1.0
    "--time-max"
        metavar = "time"
        help = "maximum time"
        arg_type = Float64
        default = 20.0
    "--time-step"
        metavar = "time"
        help = "time step"
        arg_type = Float64
    "--time-num"
        metavar = "N"
        help = "number of time"
        arg_type = Int
    "--time-log"
        help = "use logarithmic scale for time"
        action = :store_true
end
add_arg_group(s, "entanglement entropy")
@add_arg_table s begin
    "--ee"
        metavar = "XA"
        help = "compute all EEs with partition size XA"
        arg_type = Int
        required = true
end
c = parsed_args = parse_args(ARGS, s, as_symbols=true)

# Number of sites
const M = c[:M]
# Number of particles
const N = c[:N]
# Output file
const output = c[:out]
# Site occupation restriction
const site_max = c[:site_max]
# Boundary conditions
const boundary = c[:boundary]
# Size of region A
const Asize = c[:ee]

# Initial V
const V0 = c[:V0]
# Final V
const V = c[:V]

# Initial time
const time_min = c[:time_min]



if c[:time_log] && c[:time_num] === nothing
    println("--time-log must be used with --time-num")
    exit(1)
end

if c[:time_step] === nothing
    if c[:time_num] === nothing
        time_range = c[:time_min]:0.5:c[:time_max]
    else
        if c[:time_log]
            time_range = logspace(c[:time_min], c[:time_max], c[:time_num])
        else
            time_range = linspace(c[:time_min], c[:time_max], c[:time_num])
        end
    end
else
    if c[:time_num] === nothing
        time_range = c[:time_min]:c[:time_step]:c[:time_max]
    else
        println("--time-step and --time-num may not both be supplied")
        exit(1)
    end
end

if site_max === nothing
    const basis = Szbasis(M, N)
else
    const basis = RestrictedSzbasis(M, N, site_max)
end
#_______________________
ll=length(basis)
mus=zeros(Float64, M)
EigenVectors= zeros(Complex128, ll, ll)
EigenEnergies=zeros(Complex128,ll)
exp_q=zeros(Complex128, basis.K)


Vec0=zeros(Float64, ll)
Vech=zeros(Float64, ll)
Vecl=zeros(Float64, ll)

VRead=zeros(Complex128, ll)

#Initial wave function in terms of the spatial basis
wft0=zeros(Complex128, ll)
#Initial wave function in terms of the symmetry basis
wft0H0=zeros(Complex128, ll)

wftH0=zeros(Complex128, ll)
wft=zeros(Complex128, ll)

Vech.=0.01
Vecl.=0.01
Vec0.=1.0




if boundary==OBC
    num_links= basis.K-1
elseif boundary==PBC
    num_links= basis.K
end

for (i, bra) in enumerate(basis)
    cc=0
    for j=1: num_links
        j_next = j % basis.K + 1
        cc+=bra[j]*bra[j_next]
    end
    if cc== basis.N-1
        Vecl[serial_num(basis, bra)]=1.0
    elseif cc==0
        Vech[serial_num(basis, bra)]=1.0
    end
end


Norm=sqrt(sum(Vec0.^2))
Vec0.=Vec0./Norm

Norm=sqrt(sum(Vech.^2))
Vech.=Vech./Norm

Norm=sqrt(sum(Vecl.^2))
Vecl.=Vecl./Norm

#_______________________
open(output, "w") do f
    if site_max === nothing
        write(f, "# M=$(M), N=$(N),V0=$(V0),V=$(V), $(boundary)\n")
    else
        write(f, "# M=$(M), N=$(N), max=$(site_max),V0=$(V0),V=$(V), $(boundary)\n")
    end
    write(f, "# time S1(n=$(Asize)) S2(n=$(Asize)) S3(n=$(Asize))\n")
    #Create the Hamiltonian
    H = sparse_hamiltonian(basis, c[:t],mus, V0, boundary=boundary)
    print(" sparse_hamiltonian finish\n ")
    #Perform the Lanczos diagonalization to obtain the lowest eigenvector
    # http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos
    if V0/c[:t]>1.5
       d = eigs(H, nev=1, which=:SR,tol=1e-13,v0=Vech)
    elseif V0/c[:t]<-1.5
       d = eigs(H, nev=1, which=:SR,tol=1e-13,v0=Vecl)
    else
       d = eigs(H, nev=1, which=:SR,tol=1e-13,v0=Vec0)
    end

    VRead = vec(d[2][1:ll])
    for il=1:ll
        wft0[il] = VRead[il]
    end

    Norm= dot(wft0, wft0)
    wft0.= wft0./sqrt(Norm)

    
    Cycles, CycleSize, NumOfCycles = Translational_Symmetry_Cycles(basis)


    HRank=0
    for q =0: basis.K-1
       for i=1: basis.K
           exp_q[i]=exp((i-1)*(0.0-1.0im)*2*pi*q/basis.K)
       end
       #Create the Hamiltonian
       Hq,HqBasis,HqRank = Block_Diagonal_Hamiltonian(basis, Cycles, CycleSize, NumOfCycles,c[:t],V,q)
       EigenEnergies_q,VV = eig(Hq)
        
       for i_HqRank =1: HqRank
           for j_HqRank =1: HqRank
               for i_Translation =1: CycleSize[HqBasis[j_HqRank]]
                   EigenVectors[i_HqRank+ HRank,Cycles[HqBasis[j_HqRank],i_Translation]]+=exp_q[i_Translation]* VV[j_HqRank, i_HqRank]/sqrt(CycleSize[HqBasis[j_HqRank]])
               end        
           end
           EigenEnergies[i_HqRank+ HRank]= EigenEnergies_q[i_HqRank]
       end
      
        HRank += HqRank
    end
       warn("eigs finish")

       #print(EigenVectors, "\n ")

        wft0H0=zeros(Complex128, ll)
        #wft0H0.=0.0+0.0im
        for il=1: ll
            for jl=1:ll
                wft0H0[il]=wft0H0[il]+conj(EigenVectors[il,jl])*wft0[jl]
            end
        end
    time0=time_min
    for time in time_range
        Delta=-(0.0+1.0im)*(time-time0)
            time0=time
            # Calculate the entropy
            s3_particle = particle_entropy_mod(basis, Asize, wft0, site_max)
            s2_spatial, s2_operational = spatial_entropy(basis, Asize, wft0)

            for il=1: ll
                wft0H0[il]= wft0H0[il]*exp(Delta*EigenEnergies[il])
            end
            wft0.=0.0+0.0im
            
            for il=1:ll
                for jl=1: ll
                    wft0[il]=wft0[il]+wft0H0[jl]*EigenVectors[jl,il]
                end
            end


            Norm= dot(wft0, wft0)
            wft0.= wft0./sqrt(Norm)

            #print("   ","\n")
            #print(wft0,"\n")
            #print("   ","\n")

            write(f, "$(time) $(s3_particle[1]) $(real(s2_spatial)) $(real(s2_operational))\n")
            flush(f)
    end
end