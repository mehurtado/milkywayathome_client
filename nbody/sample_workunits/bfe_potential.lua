function makePotential()
    local bfe_potential = {
        type = "bfe",
        coeffs_file = "data/my_coeffs.dat",
        nmax = 10,
        lmax = 4,
        mass = 1.0e12,
        scale_radius = 20.0
    }
    return bfe_potential
end
