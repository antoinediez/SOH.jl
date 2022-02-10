module ToolBox
using QuadGK, HCubature
export make_new_dir, coefficients_Vicsek, nice_float2string

"""
Create a new directory `dirname` in the current path. 
Append a number to the name if the directory already exists. 
"""
function make_new_dir(dirname::String)
    if ispath(dirname)
        k=1
        while ispath(dirname*"_$k")
            k+=1
        end
        return mkdir(dirname*"_$k")
    else
        return mkdir(dirname)
    end
end

"""
    nice_number(x,K::Int)

Return a string which is the rounded value of `x` with `K` (decimal) digits and
if necessary, pad `"0"` to the right so that there are exactly `K` decimal places.
"""
function nice_float2string(x,K::Int)
    integer_part = trunc(Int,x)
    y = x - integer_part
    decimal_part = rpad(round(Int,y*10^K),K,"0")
    return "$(integer_part).$decimal_part"
end

"""
Return the coefficients `c1,c2,λ` of the SOH model computed for the `model` `"Fokker-Planck"` (default)
or `"BGK"` with the concentration parameters `κ`. 
These coefficients are given by 

c₁ = ∫\\_[0,π] cos(θ)exp(κcos(θ))dθ / ∫\\_[0,π] exp(κcos(θ))dθ

c₂ = ∫\\_[0,π] cos(θ)sin(θ)g(θ)exp(κcos(θ))dθ / ∫\\_[0,π] sin(θ)g(θ)exp(κcos(θ))dθ

λ = 1/κ

where 

g(θ) = θ/κ - π/κ * ∫\\_[0,θ] exp(-κcos(ψ))dψ / ∫\\_[0,π] exp(-κcos(ψ))dψ in the Fokker-Planck case

and 

g(θ) = 1 in the BGK case. 
"""
function coefficients_Vicsek(κ;model="Fokker-Planck")
    λ = 1/κ
    I1 = quadgk(θ->cos(θ)*exp(κ*cos(θ)),0,pi,rtol=1e-5)[1]
    Z1 = quadgk(θ->exp(κ*cos(θ)),0,pi,rtol=1e-5)[1]
    c1 = I1/Z1
    if model=="BGK"
        I2 = quadgk(θ->cos(θ)*sin(θ)^2*exp(κ*cos(θ)),0,pi,rtol=1e-5)[1]
        Z2 = quadgk(θ->sin(θ)^2*exp(κ*cos(θ)),0,pi,rtol=1e-5)[1]
        c2 = I2/Z2
        return c1,c2,λ
    elseif model=="Fokker-Planck"
        Zinv = quadgk(θ->exp(-κ*cos(θ)),0,pi,rtol=1e-5)[1]
        I21 = quadgk(θ->θ*cos(θ)*sin(θ)*exp(κ*cos(θ)),0,pi,rtol=1e-5)[1]/κ 
        I22 = pi/(κ*Zinv) * hcubature(x->cos(x[1])*sin(x[1])*exp(κ*cos(x[1]))*exp(-κ*cos(x[2]))*(x[2]<x[1]),(0,0),(pi,pi),rtol=1e-3)[1]
        I2 = I21 - I22
        Z21 = quadgk(θ->θ*sin(θ)*exp(κ*cos(θ)),0,pi,rtol=1e-5)[1]/κ 
        Z22 = pi/(κ*Zinv) * hcubature(x->sin(x[1])*exp(κ*cos(x[1]))*exp(-κ*cos(x[2]))*(x[2]<x[1]),(0,0),(pi,pi),rtol=1e-3)[1]
        Z2 = Z21 - Z22
        c2 = I2/Z2
        return c1,c2,λ
    else
        ArgumentError("Model not defined!")
    end
end



end