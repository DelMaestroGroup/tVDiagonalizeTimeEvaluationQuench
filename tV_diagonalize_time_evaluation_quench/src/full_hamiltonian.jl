"""
Number of links for the boundary conditions.
"""
num_links(basis::AbstractSzbasis, boundary::BdryCond) = boundary == PBC ? basis.K : basis.K - 1

"""
Create the full Hamiltonian matrix for a PBC/OBC tV chain in 1D.

    H = -\\sum_{<i, j>} t_{i,j} (c_i^\\dagger c_j + c_i c_j^\\dagger) + (V) \\sum_i n_i n_{i + 1} - \\sum_i \\mu_i n_i
"""
function full_hamiltonian(basis::AbstractSzbasis, Ts::AbstractVector{Float64}, mus::AbstractVector{Float64}, U::Float64; boundary::BdryCond=PBC)
    end_site = num_links(basis, boundary)

    length(Ts) == end_site || error("Incorrect number of Ts: $(length(Ts)) != $(end_site)")
    length(mus) == basis.K || error("Incorrect number of mus: $(length(mus)) != $(basis.K)")

    H=zeros(Complex128, basis.D, basis.D)

    for (i, bra) in enumerate(basis)
        # Diagonal part
        Usum = 0
        musum = 0
        for j=1:end_site
            musum += mus[j] * bra[j]
            j_next = j % basis.K + 1
            Usum += bra[j] * (bra[j_next])
        end
        H[i,i]= U * Usum - musum
        # Off-diagonal part
        for j=1:end_site
            j_next = j % basis.K + 1
            # Tunnel right, tunnel left.
            for (site1, site2) in [(j, j_next), (j_next, j)]
                if bra[site1] > 0
                    ket = copy(bra)
                    ket[site1] -= 1
                    ket[site2] += 1
                    if ket in basis
                        factor = 1
                        if j_next == 1
                            factor = (1)^(basis.N-1)
                        end
                        if serial_num(basis, ket)<i
                            H[serial_num(basis, ket),i]=-Ts[j] * sqrt(bra[site1]) * sqrt(bra[site2]+1) * factor
                            H[i,serial_num(basis, ket)]=conj(-Ts[j] * sqrt(bra[site1]) * sqrt(bra[site2]+1) * factor)
                        end
                    end
                end
            end
        end
    end

    H
end

function full_hamiltonian(basis::AbstractSzbasis, Ts::AbstractVector{Float64}, U::Float64; boundary::BdryCond=PBC)
    full_hamiltonian(basis, Ts, zeros(basis.K), U, boundary=boundary)
end

function full_hamiltonian(basis::AbstractSzbasis, T::Float64, mus::AbstractVector{Float64}, U::Float64; boundary::BdryCond=PBC)
    full_hamiltonian(basis, fill(T, num_links(basis, boundary)), mus, U, boundary=boundary)
end

function full_hamiltonian(basis::AbstractSzbasis, T::Float64, U::Float64; boundary::BdryCond=PBC)
    full_hamiltonian(basis, fill(T, num_links(basis, boundary)), U, boundary=boundary)
end
