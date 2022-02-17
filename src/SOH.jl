module SOH
include("boundary_conditions.jl")
include("flux.jl")
include("init.jl")
include("plot_save.jl")
include("run.jl")
include("scheme.jl")
include("toolbox.jl")

export boundary_conditions_xy!, boundary_conditions_ρuv!
export flux_x
export random_init, quadratic_potential_force, flat_quadratic_potential_force
export plot_rhoUV, update_plot!, save_data!, radial_density
export run!
export scheme_iter!
export make_new_dir, coefficients_Vicsek, nice_float2string



end     #module
