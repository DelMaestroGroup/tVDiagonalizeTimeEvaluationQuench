module tVDiagonalize

using JeszenszkiBasis

export
    BdryCond,
    OBC,
    PBC,

    sparse_hamiltonian,
    spatial_entropy,
    particle_entropy_mod,
    Translational_Symmetry_Cycles,
    Block_Diagonal_Hamiltonian,
    full_hamiltonian
"""
Boundary conditions.
"""
@enum BdryCond PBC OBC
@doc "Periodic boundary conditions." PBC
@doc "Open boundary conditions." OBC

include("sparse_hamiltonian.jl")
include("spatial_entropy.jl")
include("particle_entropy_mod.jl")
include("Translational_Symmetry_Cycles.jl")
include("Block_Diagonal_Hamiltonian.jl")
include("full_hamiltonian.jl")
end
