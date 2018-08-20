# Renyi entanglement entropy of the t-V Model after a quantum quench

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "src"))
using tVDiagonalize
using ArgParse
using JeszenszkiBasis

# ------------------------------------------------------------------------------
function getΨ0_trial(t::Float64, V0::Float64, boundary::BdryCond,
    basis::AbstractSzbasis)

    ll = length(basis)

    if -1.5 < V0/t < 1.5
        Ψ0_trial = ones(Float64,ll)/sqrt(ll)

    else
        Ψ0_trial = 0.01*ones(Float64,ll)

        if boundary==OBC
            num_link = basis.K-1
        elseif boundary==PBC
            num_link = basis.K
        end

        for (i, bra) in enumerate(basis)
            cc=0
            for j=1:num_link
                j_next = j % basis.K + 1
                cc+=bra[j]*bra[j_next]
            end

            if (cc== basis.N-1) & (V0/t < -1.5)
                Ψ0_trial[serial_num(basis, bra)]=1.0
            elseif (cc==0) & (V0/t > 1.5)
                Ψ0_trial[serial_num(basis, bra)]=1.0
            end
        end

        Ψ0_trial.=Ψ0_trial./sqrt(dot(Ψ0_trial,Ψ0_trial))
    end

    Ψ0_trial
end

# ------------------------------------------------------------------------------
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
    "--site-max"
        metavar = "N"
        help = "site occupation restriction"
        arg_type = Int
    "--obdm"
        help = "output the spatial dependence of the OBDM"
        action = :store_true
    "--spatial"
        help = "output the spatial entanglement entropy for ℓ = M/2"
        action = :store_true
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
        default = 0.0
    "--V"
        metavar = "V"
        help = "final V"
        arg_type = Float64
        default = 1.0
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
        default = 0.0
    "--time-max"
        metavar = "time"
        help = "maximum time"
        arg_type = Float64
        default = 5.0
    "--time-step"
        metavar = "time"
        help = "time step"
        arg_type = Float64
        default = 0.1
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
        metavar = "ℓ"
        help = "compute all EEs with partition size ℓ"
        arg_type = Int
        required = true
end
c = parsed_args = parse_args(ARGS, s, as_symbols=true)

# Number of sites
const M = c[:M]
# Number of particles
const N = c[:N]

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
        time_num=length(time_range) 
    else
        if c[:time_log]
            time_range = logspace(c[:time_min], c[:time_max], c[:time_num])
        else
            time_range = linspace(c[:time_min], c[:time_max], c[:time_num])
        end
        time_num=c[:time_num]
    end
else
    if c[:time_num] === nothing
        time_range = c[:time_min]:c[:time_step]:c[:time_max]
        time_num=length(time_range) 
    else
        println("--time-step and --time-num may not both be supplied")
        exit(1)
    end
end

if length(time_range) > 1
    Δt = time_range[2]-time_range[1]
else
    Δt = time_range[1]
end

# Output file
if c[:out] === nothing
     output = @sprintf "partEE_%02d_%02d_%+5.3f_%+5.3f_%6.4f_%06.3f_%06.3f_%1d.dat" M N V0 V Δt time_range[1] time_range[end] Asize
else
     output = c[:out]
end

# output file if we are measuring the spatial entanglement entropy
if c[:spatial]
     spat_output = @sprintf "spatEE_%02d_%02d_%+5.3f_%+5.3f_%6.4f_%06.3f_%06.3f_%1d.dat" M N V0 V Δt time_range[1] time_range[end] Asize
end

# are we restricting the number of particles per site?
if site_max === nothing
    const basis = Szbasis(M, N)
else
    const basis = RestrictedSzbasis(M, N, site_max)
end

#_______________________________________________________________________________
ll=length(basis)
μ=zeros(Float64, M)
exp_q=zeros(Complex128, basis.K)
EigenEnergie=Complex128
EigenVector= zeros(Complex128, ll)
#  Initial wave function in terms of the spatial basis
Ψ=zeros(Complex128, ll)
Ψn=Complex128
# the one body density matrix
obdm=zeros(Float64,M,length(time_range))
#  wave function in terms of the spatial basis at time t
Ψt=zeros(Complex128, ll, time_num)


#_______________________________________________________________________________

# open and prepare files for output
f_part = open(output, "w")
if site_max === nothing
    write(f_part, "# M=$(M), N=$(N),V0=$(V0), V=$(V), $(boundary)\n")
else
    write(f_part, "# M=$(M), N=$(N), max=$(site_max), V0=$(V0), V=$(V), $(boundary)\n")
end
write(f_part,@sprintf "#%11s%24s%24s\n" "time (tJ)" "S₁(n=$(Asize))" "S₂(n=$(Asize))")

if c[:spatial]
    ℓsize = div(M, 2)
    f_spat = open(spat_output, "w")
    if site_max === nothing
        write(f_spat, "# M=$(M), N=$(N),V0=$(V0), V=$(V), $(boundary)\n")
    else
        write(f_spat, "# M=$(M), N=$(N), max=$(site_max), V0=$(V0), V=$(V), $(boundary)\n")
    end
    write(f_spat,@sprintf "#%11s%24s%24s\n" "time (tJ)" "S₁(ℓ=$(ℓsize))" "S₂(ℓ=$(ℓsize))")
end

#_______________________________________________________________________________

# Create the Hamiltonian
H = sparse_hamiltonian(basis, c[:t], μ, V0, boundary=boundary)
#println("size(H) = ",Base.summarysize(H))
print(" sparse_hamiltonian finish\n ")

# H0 = full_hamiltonian(basis, c[:t], V0,boundary=boundary)
# EigenValues, EigenVectors = eig(H0)
# wft0 = EigenVectors[:,1]

#Perform the Lanczos diagonalization to obtain the lowest eigenvector
# http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos

# I don't understand why this copying is necessary, it is a type conversion thing
Ψ = eigs(H, nev=1, which=:SR,tol=1e-13,v0=getΨ0_trial(c[:t],V0,boundary,basis))[2][1:ll].*ones(Complex128,ll)
#println("size(complex) = ", Base.summarysize(Ψ[1]))

# Exploit symmetries of the hamiltonian to perform a bloack diagonalization
Cycles, CycleSize, NumOfCycles = Translational_Symmetry_Cycles(basis)

HRank=0
for q =0: basis.K-1
   for i=1: basis.K
       exp_q[i]=exp((i-1)*(0.0+1.0im)*2*pi*q/basis.K)
   end

   #Create the Hamiltonian
   Hq,HqBasis,HqRank = Block_Diagonal_Hamiltonian(basis, Cycles, CycleSize, NumOfCycles,c[:t],V,q)
   EigenEnergies_q,VV = eig(Hq)

   for i_HqRank =1: HqRank
       for j_HqRank =1: HqRank
           for i_Translation =1: CycleSize[HqBasis[j_HqRank]]
               EigenVector[Cycles[HqBasis[j_HqRank],i_Translation]]+=exp_q[i_Translation]* VV[j_HqRank, i_HqRank]/sqrt(CycleSize[HqBasis[j_HqRank]])
           end
       end
       #println("size(EV) = ", Base.summarysize(EigenVector))
       EigenEnergie= EigenEnergies_q[i_HqRank]
       Ψn=0.0+0.0im
       for α=1:ll
            Ψn=Ψn+conj(EigenVector[α])*Ψ[α]
       end
       for (it, time) in enumerate(time_range)
          TimeEvolutionFactor=exp(-(0.0+1.0im)*time*EigenEnergie)
          for α=1:ll
             Ψt[α,it]=Ψt[α,it]+Ψn*EigenVector[α]*TimeEvolutionFactor
          end        
       end
       EigenVector.=0.0+0.0im
   end 
end   

print(" Block_Diagonal_Hamiltonian finished\n ")

#____________________________________________________________
# H = full_hamiltonian(basis, c[:t], V, boundary=boundary)
# EigenEnergies,EigenVectors = eig(H)
#
# warn("eigs finish")
# #print(EigenVectors, "\n ")
#
# for n=1:ll
#     for α=1:ll
#         wft0H0[n] = wft0H0[n] + conj(EigenVectors[α,n])*wft0[α]
#     end
# end

# for time in time_range
#
#     # Calculate the entropy
#     wft0.=0.0+0.0im
#     for α=1:ll
#         for n=1:ll
#             wft0[α] = wft0[α] + exp(-(0.0+1.0im)*time*EigenEnergies[n])*wft0H0[n]*(EigenVectors[α,n])
#         end
#     end
#
#     s3_particle = particle_entropy_mod(basis, Asize, wft0, site_max)
#     write(f, " $(time)\t $(s3_particle[1])\t $(s3_particle[2]) \n")
#     flush(f)
# end
#__________________________________________________________
# Calculate the entropy if we start from t = 0
if abs(time_range[1]) < 1.0E-12

    if c[:obdm] && Asize == 1
        s_particle,obdm[:,1] = particle_entropy_mod(basis, Asize, Ψ, site_max,c[:obdm])
    else
        s_particle = particle_entropy_mod(basis, Asize, Ψ, site_max,c[:obdm])
    end

    write(f_part, @sprintf "%12.6f%24.12E%24.12E\n" time_range[1] s_particle[1] s_particle[2])
    flush(f_part)
    time_start = 2

    if c[:spatial]
        s_spatial = spatial_entropy(basis, ℓsize, Ψ)
        write(f_spat, @sprintf "%12.6f%24.12E%24.12E\n" time_range[1] s_spatial[1] s_spatial[2])
        flush(f_spat)
    end

else
    time_start = 1
end

# do we start from the 1st or 2nd time step?
it = time_start

for time in time_range[time_start:end] 

    Ψt[:,it].= Ψt[:,it]./sqrt(dot(Ψt[:,it],Ψt[:,it]))

    if c[:spatial]
        s_spatial = spatial_entropy(basis, ℓsize, Ψt[:,it])
        write(f_spat, @sprintf "%12.6f%24.12E%24.12E\n" time s_spatial[1] s_spatial[2])
        flush(f_spat)
    end

    if c[:obdm] && Asize == 1
        s_particle,obdm[:,it] = particle_entropy_mod(basis, Asize, Ψt[:,it], site_max,c[:obdm])
    else
        s_particle = particle_entropy_mod(basis, Asize, Ψt[:,it], site_max,c[:obdm])
    end

    write(f_part, @sprintf "%12.6f%24.12E%24.12E\n" time s_particle[1] s_particle[2])
    flush(f_part)
    it += 1
end

close(f_part)

if c[:spatial]
    close(f_spat)
end

#_______________________________________________________________________________
# output the time dependent OBDM to disk
if c[:obdm] && Asize == 1
    obdm_name = @sprintf "obdm_%02d_%02d_%+5.3f_%+5.3f_%6.4f_%06.3f_%06.3f.dat" M N V0 V Δt time_range[1] time_range[end]
    open(obdm_name, "w") do obdm_f
        write(obdm_f, @sprintf "#%11s" "|i-j|")
        for time in time_range
            write(obdm_f, @sprintf "%16.6f" time)
        end
        write(obdm_f, "\n")
        flush(obdm_f)

        for i = 1:M
            write(obdm_f, @sprintf "%16d" (i-Int(M/2)))
            for (time_index, time) in enumerate(time_range)
                write(obdm_f, @sprintf "%16.6E" obdm[i,time_index])
            end
            write(obdm_f, "\n")
            flush(obdm_f)
        end

     end
end
