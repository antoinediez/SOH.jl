"""
    random_init(
        ncellx,ncelly,
        mean_ρ,range_ρ,
        mean_θ,range_θ;
        bcond_x="periodic",bcond_y="periodic"
    )

Sample a random data `(ρ,u,v)` of size `(ncellx+2)*(ncelly+2)`. The density `ρ` is sampled
uniformly around `mean_ρ` in an interval of size `range_ρ`. The velocity component `u` and
`v` are respectively given by `cos(θ)` and `sin(θ)` where `θ` is sampled uniformly around
`mean_θ` in an interval of size `range_θ`. Boundary conditions are then applied.
"""
function random_init(
    ncellx,ncelly,
    mean_ρ,range_ρ,
    mean_θ,range_θ;
    bcond_x="periodic",bcond_y="periodic"
)
    ρ = mean_ρ .+ range_ρ .* (rand(ncellx+2,ncelly+2) .- 0.5)
    θ = mean_θ .+ range_θ .* (rand(ncellx+2,ncelly+2) .- 0.5)
    u = cos.(θ)
    v = sin.(θ)
    boundary_conditions_ρuv!(ρ,u,v,bcond_x,bcond_y)
    return ρ,u,v
end

"""
    quadratic_potential_force(ncellx,ncelly,Δx,Δy,C)

Return the components `Fx` and `Fy` of the force exerted by the quadratic radial
potential V(r) = Cr²/2 centered around the midpoint of the box `[0,Lx]*[0,Ly]`
where `Lx = ncellx * Δx` and `Ly = ncelly * Δy`. For convenience, the outputs `Fx` and `Fy` are
matrices of size `(ncellx+2)*(ncelly+2)` with zeros on the boundary.
"""
function quadratic_potential_force(ncellx,ncelly,Δx,Δy,C)
    Fx = zeros(ncellx+2,ncelly+2)
    Fy = zeros(ncellx+2,ncelly+2)
    Lx = ncellx * Δx
    Ly = ncelly * Δy
    for i in 2:ncellx+1
        for j in 2:ncelly+1
            mx = Lx/2
            my = Ly/2
            xi = (i-2)*Δx + Δx/2
            yj = (j-2)*Δy + Δy/2
            Fx[i,j] = -C * (xi - mx)
            Fy[i,j] = -C * (yj - my)
        end
    end
    return Fx, Fy
end
"""
    flat_quadratic_potential_force(ncellx,ncelly,Δx,Δy,r0,C)

Return the components `Fx` and `Fy` of the force exerted by the quadratic radial
potential defined by V(r) = Cr²/2 for r>r0 and V(r)=0 otherwise,
centered around the midpoint of the box `[0,Lx]*[0,Ly]` where `Lx = ncellx * Δx`
and `Ly = ncelly * Δy`. For convenience, the outputs `Fx` and `Fy` are
matrices of size `(ncellx+2)*(ncelly+2)` with zeros on the boundary.
"""
function flat_quadratic_potential_force(ncellx,ncelly,Δx,Δy,r0,C)
    Fx = zeros(ncellx+2,ncelly+2)
    Fy = zeros(ncellx+2,ncelly+2)
    Lx = ncellx * Δx
    Ly = ncelly * Δy
    for i in 2:ncellx+1
        for j in 2:ncelly+1
            mx = Lx/2
            my = Ly/2
            xi = (i-2)*Δx + Δx/2
            yj = (j-2)*Δy + Δy/2
            rij = sqrt((xi-mx)^2+(yj-my)^2)
            if rij > r0
                Fx[i,j] = -C * (rij - r0) * (xi - mx)/rij
                Fy[i,j] = -C * (rij - r0) * (yj - my)/rij
            end
        end
    end
    return Fx, Fy
end
