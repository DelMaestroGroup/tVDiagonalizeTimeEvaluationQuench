abstract type AbstractFermionsbasis end

"""
Basis of Fermionic occupation vectors.
"""
struct Fermionsbasis <: AbstractFermionsbasis
    "Number of sites."
    K::Int
    "Number of fermions."
    N::Int

    "Number of basis vectors."
    D::Int

    "Occupation vectors (K by D)."
    vectors::Matrix{Int}
end

num_vectors(N::Int, K::Int) = binomial(K, N)
num_vectors(basis::Fermionsbasis, N::Int, K::Int) = num_vectors(N, K)


"""
    Fermionsbasis(K::Int, N::Int)
Create a basis for `K` sites and `N` fermions.
"""
function Fermionsbasis(K::Int, N::Int)
    K >= 1 || throw(DomainError(K, "At least 1 site is required."))
    N >= 0 || throw(DomainError(N, "At least 0 particles are required."))
    N <= K || throw(DomainError(N, "fermions do not fit on the sites."))

    # Basis size.
    D = num_vectors(N, K)
#    dNM = M > 0 ? div(N, M) : 1 # @h dNM=1

    v = zeros(Int, K)
    for j in 1:N
        v[j] = 1
    end
    vectors = Matrix{Int}( K, D)
     #   vectors = Matrix{Int}(undef, K, D)

    vectors[:, 1] .= v

    for i in 2:D
        if v[1] > 0
            j = findfirst(!iszero, v .== 0)

            v[j] += 1
            v[j-1] -= 1 
        else
            j = findfirst(!iszero, v)
            k = j + findfirst(!iszero, @view(v[(j+1):end]) .== 0)

            v[k-j] = v[j] - 1
            v[k] += 1
            for l in 1:(k-j-1)
                v[l] = 1
            end
            # The indices after the first one differ from those in the paper.
            for l in (k-j+1):(k-1)
                v[l] = 0
            end
        end
        vectors[:, i] .= v
    end

    Fermionsbasis(K, N, D, vectors)
end
#--------------------------------------------------------------------------------

"""
    serial_num(K::Int, N::Int, M::Int, v::AbstractVector{Int})
Compute the serial number of occupation vector `v` in a basis with `K` sites
and `N` fermions.
"""
function serial_num(K::Int, N::Int, v::AbstractArray{Int})
    I = 1

    for mu in 1:K
        s = 0
        for nu in (mu+1):K
            s += v[nu]
        end
        for i in 0:(v[mu]-1)
            I += num_vectors(N-s-i, mu-1)
        end
    end

    I
end
serial_num(basis::Fermionsbasis, v::AbstractVector{Int}) = serial_num(basis.K, basis.N, v)

serial_num(basis::Fermionsbasis, K::Int, N::Int, v::AbstractVector{Int}) = serial_num(K, N, v)

"""
    sub_serial_num(basis::AbstractFermionsbasis, v::AbstractVector{Int})
Compute the serial number of the reduced occupation vector `v`, which has only
a subset of the sites present in `basis`.
"""

function sub_serial_num(basis::Fermionsbasis, v::AbstractVector{Int})
    K = length(v)
    N = sum(v)

    # Only one way to have no sites.
    K >= 1 || return 1
    # Only one way to have no particles.
    N >= 1 || return 1

    # Count the zero-particle case.
    I = 1

    for n in 1:(N-1)
        I += num_vectors(n, K)
    end

    I + serial_num(K, N, v)
end
#--------------------------------------------------------------------------------

Base.getindex(basis::AbstractFermionsbasis, i::Int) = @view basis.vectors[:, i]
#--------------------------------------------------------------------------------

mutable struct FermionsbasisIterState
    i::Int
end

function Base.start(basis::AbstractFermionsbasis)
    FermionsbasisIterState(0)
end

function Base.next(basis::AbstractFermionsbasis, state::FermionsbasisIterState)
    state.i += 1

    @view(basis.vectors[:, state.i]), state
end

function Base.done(basis::AbstractFermionsbasis, state::FermionsbasisIterState)
    state.i == basis.D
end

Base.eltype(::Type{AbstractFermionsbasis}) = Vector{Int}
Base.length(basis::AbstractFermionsbasis) = basis.D

function Base.in(v::AbstractVector{Int}, basis::Fermionsbasis)
    length(v) == basis.K && sum(v) == basis.N && maximum(v) <= 1
end