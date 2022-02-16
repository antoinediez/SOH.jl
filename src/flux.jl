"""
    eigenvalues_Roe(u,c1,c2,λ)

The eigenvalues of the Roe matrix with a 1D velocity `u` are c₂u, c₂u + √Δ and c₂u - √Δ
with the discriminant Δ = λc₁ - c₂(c₁-c₂)u.
"""
function eigenvalues_Roe(u,c1,c2,λ)
    Δ = (c2^2 - c1*c2) * u^2 + λ*c1
    eval_0 = c2*u
    eval_p = c2*u + sqrt(Δ)
    eval_m = c2*u - sqrt(Δ)
    return eval_p, eval_m, eval_0
end

"""
    abs_Roe_matrix(ρl,ρr,ul,ur,vl,vr,c1,c2,λ)

Return the absolute value of the Roe matrix for the 1D Riemannn problem with respective left and right densities
`(ρl,ul,vl)` and `(ρr,ur,vr)`. The eigenvalues are increased to avoid sonic rarefaction wave (LLF method).
"""
function abs_Roe_matrix(ρl,ρr,ul,ur,vl,vr,c1,c2,λ)

    #### Roe average
    um = (sqrt(ρl)*ul + sqrt(ρr)*ur) / (sqrt(ρl) + sqrt(ρr))
    vm = (sqrt(ρl)*vl + sqrt(ρr)*vr) / (sqrt(ρl) + sqrt(ρr))

    #### Eigenvalues and eigenvectors
    eval_p, eval_m, eval_0 = eigenvalues_Roe(um,c1,c1,λ)
    z_p = (c1*c2*um*vm - c2*vm*eval_p) / (c2*um - eval_p)
    z_m = (c1*c2*um*vm - c2*vm*eval_m) / (c2*um - eval_m)
    detP = c1 * (eval_m - eval_p)

    #### Increase the eigenvalues
    eval_l = eigenvalues_Roe(ul,c1,c2,λ)
    eval_r = eigenvalues_Roe(ur,c1,c2,λ)
    eval_p_fix = max(abs(eval_l[1]),abs(eval_r[1]))
    eval_m_fix = max(abs(eval_l[2]),abs(eval_r[2]))
    eval_0_fix = max(abs(eval_l[3]),abs(eval_r[3]))

    #### Compute Roe matrix A = P*D*P^{-1}
    a11 = c1*(eval_p_fix*eval_m - eval_m_fix*eval_p)/detP
    a12 = c1^2*(eval_m_fix - eval_p_fix)/detP
    a13 = 0.
    a21 = eval_m*eval_p*(eval_p_fix - eval_m_fix)/detP
    a22 = c1*(eval_m_fix*eval_m - eval_p_fix*eval_p)/detP
    a23 = 0.
    a31 = (-eval_m_fix*eval_p*z_m + eval_p_fix*eval_m*z_p + eval_0_fix*(eval_p*z_m - eval_m*z_p))/detP
    a32 = c1*(eval_m_fix*z_m - eval_p_fix*z_p + eval_0_fix*(z_p - z_m))/detP
    a33 = eval_0_fix
    return a11,a12,a13,a21,a22,a23,a31,a32,a33
end

"""
    flux_x(ρl,ρr,ul,ur,vl,vr,c1,c2,λ,method)

Return the flux along the x-axis for the 1D Riemannian problem given by the respective left and right
densities `(ρl,ul,vl)` and `(ρr,ur,vr)`. The chosen method `method` can be either `"Roe"` or `"HLLE"`.
"""
function flux_x(ρl,ρr,ul,ur,vl,vr,c1,c2,λ,method)
    if ρl<1e-9 && ρr<1e-9
        return 0.,0.,0.
    end
    if method == "Roe"
        a11,a12,a13,a21,a22,a23,a31,a32,a33 = abs_Roe_matrix(ρl,ρr,ul,ur,vl,vr,c1,c2,λ)
        avg = 0.5 .* (c1*ρl*ul + c1*ρr*ur,
                    c2*(ρl*ul^2 + ρr*ur^2) + λ*(ρl + ρr),
                    c2*(ρl*ul*vl + ρr*ur*vr))
        Roe_term = 0.5 .* (a11*(ρr - ρl) + a12*(ρr*ur - ρl*ul) + a13*(ρr*vr - ρl*vl),
                        a21*(ρr - ρl) + a22*(ρr*ur - ρl*ul) + a23*(ρr*vr - ρl*vl),
                        a31*(ρr - ρl) + a32*(ρr*ur - ρl*ul) + a33*(ρr*vr - ρl*vl))
        return avg .- Roe_term
    elseif method == "HLLE"
        um = (sqrt(ρl)*ul + sqrt(ρr)*ur) / (sqrt(ρl) + sqrt(ρr))
        vp_l = eigenvalues_Roe(ul,c1,c2,λ)
        vp_r = eigenvalues_Roe(ur,c1,c2,λ)
        vp_Roe = eigenvalues_Roe(um,c1,c2,λ)
        s_l = minimum(min(vp_l,vp_Roe))
        s_r = maximum(max(vp_r,vp_Roe))
        s_lm = min(s_l,0.)
        s_rp = max(s_r,0.)

        f_l = (c1*ρl*ul, c2*ρl*ul^2 + λ*ρl, c2*ρl*ul*vl)
        f_r = (c1*ρr*ur, c2*ρr*ur^2 + λ*ρr, c2*ρr*ur*vr)
        U_l = (ρl,ρl*ul,ρl*vl)
        U_r = (ρr,ρr*ur,ρr*vr)
        return ((s_rp.*f_l .- s_lm.*f_r) .+ (s_rp*s_lm) .* (U_r .- U_l)) ./ (s_rp - s_lm)
    else
        error("Method not defined!")
    end

end
