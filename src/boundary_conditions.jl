"""Contains the functions to apply boundary conditions and define classical the exterior forces """
module BoundaryConditions
export boundary_conditions_xy!, boundary_conditions_ρuv!

"""
    boundary_conditions_x!(density,bcond_x)

Apply the boundary conditions `bcond_x` on the x-axis of the matrix `density`.
Currently, only the boundary conditions `"periodic"`, `"Neumann"` and `"reflecting"` are implemented.

"""
function boundary_conditions_x!(density,bcond_x)
    if bcond_x == "periodic"
        density[1,:] = density[end-1,:]
        density[end,:] = density[2,:]
    elseif bcond_x == "Neumann"
        density[1,:] = density[2,:]
        density[end,:] = density[end-1,:]
    elseif bcond_x == "reflecting"
        density[1,:] = -density[2,:]
        density[end,:] = -density[end-1,:]
    else
        ArgumentError("Boundary condition not defined!")
    end
end

"""
    boundary_conditions_y!(density,bcond_y) 

Apply the boundary conditions `bcond_y` on the y-axis of the matrix `density`.
Currently, only the boundary conditions `"periodic"`, `"Neumann"` and `"reflecting"` are implemented.
"""
function boundary_conditions_y!(density,bcond_y)
    if bcond_y == "periodic"
        density[:,1] = density[:,end-1]
        density[:,end] = density[:,2]
    elseif bcond_y == "Neumann"
        density[:,1] = density[:,2]
        density[:,end] = density[:,end-1]
    elseif bcond_y == "reflecting"
        density[:,1] = -density[:,2]
        density[:,end] = -density[:,end-1]
    else
        ArgumentError("Boundary condition not defined!")
    end
end

"""
    boundary_conditions_xy!(density,bcond_x,bcond_y)

Apply successively the functions `boundary_conditions_x!` and `boundary_conditions_x!` to
the matrix `density`. 
"""
function boundary_conditions_xy!(density,bcond_x,bcond_y)
    boundary_conditions_x!(density,bcond_x)
    boundary_conditions_y!(density,bcond_y)
end

"""
    boundary_conditions_ρuv!(ρ,u,v,bcond_x,bcond_y)

Apply the boundary conditions to the data `(ρ,u,v)`. 
"""
function boundary_conditions_ρuv!(ρ,u,v,bcond_x,bcond_y)
    if bcond_x == "periodic" || bcond_x == "Neumann"
        boundary_conditions_x!(ρ,bcond_x)
        boundary_conditions_x!(u,bcond_x)
        boundary_conditions_x!(v,bcond_x)
    elseif bcond_x == "reflecting"
        boundary_conditions_x!(ρ,"Neumann")
        boundary_conditions_x!(u,"reflecting")
        boundary_conditions_x!(v,"Neumann")
    else
        ArgumentError("Boundary condition not defined!")
    end

    if bcond_y == "periodic" || bcond_y == "Neumann"
        boundary_conditions_y!(ρ,bcond_y)
        boundary_conditions_y!(u,bcond_y)
        boundary_conditions_y!(v,bcond_y)
    elseif bcond_y == "reflecting"
        boundary_conditions_y!(ρ,"Neumann")
        boundary_conditions_y!(u,"Neumann")
        boundary_conditions_y!(v,"reflecting")
    else
        ArgumentError("Boundary condition not defined!")
    end
end

end    #module