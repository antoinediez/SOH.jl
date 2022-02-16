"""
    scheme_iter!(
        ρ,u,v,
        Δx,Δy,Δt,
        c1,c2,λ,
        bcond_x,bcond_y,
        Fx,Fy,
        method
    )

One iteration of the numerical scheme which updates the data `(ρ,u,v)`. The scheme
follows the methodology introduced by

S. Motsch, L. Navoret, *Numerical simulations of a non-conservative hyperbolic
system with geometric constraints describing swarming behavior*, Multiscale Model. Simul.,
Vol. 9, No. 3, pp. 1253-1275, 2011.

The finite volume scheme is based on the following three-step splitting:

1. The conservative part

    ∂ₜρ + c₁∇ₓ⋅(ρΩ) = 0
    ∂ₜ(ρΩ) + c₂∇ₓ⋅(ρΩ⊗Ω) + λ∇ₓρ = 0

2. The relaxation part

    ε∂ₜΩ = (1-|Ω|^2)Ω

3. The source term

    ∂ₜρ = 0
    ∂ₜΩ = λP(Ω)F, F=-∇ₓV

For the first and second steps, a simple dimensional splitting is used. The second
step reduces to a mere normalization. The third step can also be solved explicitly
in dimension 2.
"""
function scheme_iter!(
    ρ,u,v,
    Δx,Δy,Δt,
    c1,c2,λ,
    bcond_x,bcond_y,
    Fx,Fy,
    flux_ρ, flux_u, flux_v,
    method
)
    ncellx = size(ρ)[1] - 2
    ncelly = size(ρ)[2] - 2

    #===============================================================================#
    #======================= Dimensional splitting: x-axis =========================#
    #===============================================================================#

    bad_rho_x = 0   #Count the number of nonpositive ρ (should not happen with HLLE)

    #*********************** Finite volume scheme along the x-axis *****************#
    # U = (ρ,ρu,ρv), Ω = (u,v)
    # 1. (Compute the flux)
    # 2. (Conservative part) U^{n+1/2}_i = U^{n}_i - Δt/Δx(F^n_{i+1/2} - F^n_{i-1/2})
    # 3. (Relaxation) Ω^{n+1}_i = Ω^{n+1/2}_i/|Ω^{n+1/2}_i|
    # 4. (Boundary conditions)
    #*******************************************************************************#

    #----------------------- 1. Compute the flux -----------------------------------#

    @inbounds for j in 1:ncelly
        @simd for i in 1:(ncellx+1)
            F = flux_x(
                ρ[i,j+1],ρ[i+1,j+1],
                u[i,j+1],u[i+1,j+1],
                v[i,j+1],v[i+1,j+1],
                c1,c2,λ,method
            )
            flux_ρ[i,j] = F[1]
            flux_u[i,j] = F[2]
            flux_v[i,j] = F[3]
        end
    end

    #----------------------- Update ------------------------------------------------#

    @inbounds for j in 2:(ncelly+1)
        for i in 2:(ncellx+1)

            #------ 2. Conservative part -------------------------------------------#
            U = (ρ[i,j] - (Δt/Δx) * (flux_ρ[i, j-1] - flux_ρ[i-1, j-1]),
                 ρ[i,j]*u[i,j] - (Δt/Δx) * (flux_u[i, j-1] - flux_u[i-1, j-1]),
                 ρ[i,j]*v[i,j] - (Δt/Δx) * (flux_v[i, j-1] - flux_v[i-1, j-1]))

            #------- 3. Relaxation -------------------------------------------------#
            norm = sqrt(U[2]^2 + U[3]^2)

            if U[1]<0
                bad_rho_x += 1  # Bad ρ<0
                ρ[i,j] = 1e-6   # Set ρ to a small value
            else
                ρ[i,j] = U[1]
            end
            if norm>1e-12   # If the norm is too small, it is treated as zero
                u[i,j] = U[2]/norm
                v[i,j] = U[3]/norm
            else
                normuv = sqrt(u[i,j]^2+v[i,j]^2)
                u[i,j] = u[i,j]/normuv
                v[i,j] = v[i,j]/normuv
            end

        end
    end

    # Warn the user in case of unphysical ρ
    if bad_rho_x>0
        println("Warning! Nonpositive densities in $bad_rho_x cells ($(bad_rho_x/(ncellx*ncelly))%) have been removed.")
    end

    #----------------------- 4. Boundary conditions --------------------------------#
    boundary_conditions_ρuv!(ρ,u,v,bcond_x,bcond_y)


    #===============================================================================#
    #======================= Dimensional splitting: y-axis =========================#
    #===============================================================================#

    bad_rho_y = 0 #Count the number of nonpositive ρ (should not happen with HLLE)

    #----------------------- 1. Compute the flux -----------------------------------#
    # Note: Change basis (u,v) -> (v,-u) in order to use the same function `flux_x`
    #-------------------------------------------------------------------------------#

    @inbounds for j in 1:(ncelly+1)
        @simd for i in 1:ncellx
            F = flux_x(
                ρ[i+1,j],ρ[i+1,j+1],
                v[i+1,j],v[i+1,j+1],
                -u[i+1,j],-u[i+1,j+1],
                c1,c2,λ,method
            )
            flux_ρ[i,j] = F[1]
            flux_u[i,j] = -F[3]
            flux_v[i,j] = F[2]
        end
    end

    #----------------------- Update ------------------------------------------------#
    @inbounds for j in 2:(ncelly+1)
        @simd for i in 2:(ncellx+1)

            #------ 2. Conservative part -------------------------------------------#
            U = (ρ[i,j] - (Δt/Δy) * (flux_ρ[i-1, j] - flux_ρ[i-1, j-1]),
                 ρ[i,j]*u[i,j] - (Δt/Δy) * (flux_u[i-1, j] - flux_u[i-1, j-1]),
                 ρ[i,j]*v[i,j] - (Δt/Δy) * (flux_v[i-1, j] - flux_v[i-1, j-1]))

            #------- 3. Relaxation -------------------------------------------------#
            norm = sqrt(U[2]^2 + U[3]^2)

            if U[1]<0
                bad_rho_y += 1  # Bad ρ<0
                ρ[i,j] = 1e-6   # Set ρ to a small value
            else
                ρ[i,j] = U[1]
            end
            if norm>1e-12 # If the norm is too small, it is treated as zero
                u[i,j] = U[2]/norm
                v[i,j] = U[3]/norm
            else
                normuv = sqrt(u[i,j]^2+v[i,j]^2)
                u[i,j] = u[i,j]/normuv
                v[i,j] = v[i,j]/normuv
            end

        end
    end

    # Warn the user in case of unphysical ρ
    if bad_rho_y>0
        println("Warning! Nonpositive densities in $bad_rho_y cells ($(bad_rho_y/(ncellx*ncelly))%) have been removed.")
    end

    #----------------------- 4. Boundary conditions --------------------------------#
    boundary_conditions_ρuv!(ρ,u,v,bcond_x,bcond_y)


    #===============================================================================#
    #======================= Fractional splitting: exterior force ==================#
    #===============================================================================#
    if !isnothing(Fx) && !isnothing(Fy)
        scheme_potential!(u,v,Fx,Fy,Δt,λ)
        boundary_conditions_ρuv!(ρ,u,v,bcond_x,bcond_y)
    end


end

"""
    scheme_potential!(u,v,Fx,Fy,Δt,λ)

Solve explicitly the source term

∂ₜρ = 0
∂ₜΩ = λP(Ω)F, F=-∇ₓV.

In dimension two, the solution is

Ω(t) = (cos(θ),sin(θ))ᵀ where

θ(t) = ψ + 2atan(C_0 exp(-λ|F|t))

with C_0 = tan((θ(0)-ψ)/2) and F = |F|(cos(ψ),sin(ψ))ᵀ.
"""
function scheme_potential!(u,v,Fx,Fy,Δt,λ)
    ncellx = size(u)[1] - 2
    ncelly = size(u)[2] - 2
    @inbounds for j in 2:ncelly+1
        for i in 2:ncellx+1
            normF = sqrt(Fx[i,j]^2 + Fy[i,j]^2)
            if normF > 1e-9
                θ0 = atan(v[i,j],u[i,j])
                ψ = atan(Fy[i,j],Fx[i,j])
                C0 = tan((θ0-ψ)/2)
                θ = ψ + 2 * atan(C0*exp(-λ*normF*Δt))
                u[i,j] = cos(θ)
                v[i,j] = sin(θ)
            end
        end
    end
end
