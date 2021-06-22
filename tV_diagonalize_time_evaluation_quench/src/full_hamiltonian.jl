"""
Number of links for the boundary conditions.
"""
num_links(basis::AbstractFermionsbasis, boundary::BdryCond) = boundary == PBC ? basis.K : basis.K - 1

"""
Create the full Hamiltonian matrix for a PBC/OBC tV chain in 1D.

    H = -\\sum_{<i, j>} t_{i,j} (c_i^\\dagger c_j + c_i c_j^\\dagger) + (V) \\sum_i n_i n_{i + 1}+ (Vp) \\sum_i n_i n_{i + 2} - \\sum_i \\mu_i n_i
"""
function full_hamiltonian(basis::AbstractFermionsbasis, Ts::AbstractVector{Float64}, mus::AbstractVector{Float64}, V::Float64 ,Vp::Float64; boundary::BdryCond=PBC)
    end_site = num_links(basis, boundary)

    length(Ts) == end_site || error("Incorrect number of Ts: $(length(Ts)) != $(end_site)")
    length(mus) == basis.K || error("Incorrect number of mus: $(length(mus)) != $(basis.K)")

    H=zeros(ComplexF64, basis.D, basis.D)

    for (i, bra) in enumerate(basis)
        # Diagonal part
        Vsum = 0
        Vpsum = 0
        musum = 0
        for j=1:end_site
            musum += mus[j] * CheckSite(bra,j)
            j_next = j % basis.K + 1
            j_next_next = j_next % basis.K + 1
            Vsum += CheckSite(bra,j) * CheckSite(bra,j_next)
            Vpsum += CheckSite(bra,j) * CheckSite(bra, j_next_next)
        end
        H[i,i]= V * Vsum +Vp * Vpsum- musum

        # Off-diagonal part
        for j=1:end_site
            j_next = j % basis.K + 1
            # Tunnel right, tunnel left.
            for (site1, site2) in [(j, j_next), (j_next, j)]
                 if CheckSite(bra,site1) == 1
                    ket = copy(bra)
                    if CheckSite(bra,site2) == 0
                        ket =EmptySite(ket,site1)
                        ket =OccupySite(ket,site2)
                        factor = 1
                        if j_next == 1
                            factor = (1)^(basis.N-1)
                        end
                        if serial_num(basis, ket)<i
                            H[serial_num(basis, ket),i]=-Ts[j] * factor
                            H[i,serial_num(basis, ket)]=conj(-Ts[j] * factor)
                        end
                    end
                end
            end
        end
    end

    H
end

function full_hamiltonian(basis::AbstractFermionsbasis, Ts::AbstractVector{Float64}, V::Float64 ,Vp::Float64; boundary::BdryCond=PBC)
    full_hamiltonian(basis, Ts, zeros(basis.K), V, Vp , boundary=boundary)
end

function full_hamiltonian(basis::AbstractFermionsbasis, T::Float64, mus::AbstractVector{Float64}, V::Float64 ,Vp::Float64; boundary::BdryCond=PBC)
    full_hamiltonian(basis, fill(T, num_links(basis, boundary)), mus, V, Vp, boundary=boundary)
end

function full_hamiltonian(basis::AbstractFermionsbasis, T::Float64, V::Float64 ,Vp::Float64; boundary::BdryCond=PBC)
    full_hamiltonian(basis, fill(T, num_links(basis, boundary)), V, Vp, boundary=boundary)
end
