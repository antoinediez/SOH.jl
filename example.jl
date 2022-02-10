#-------- Add the package to the current path if not installed globally -------#
push!(LOAD_PATH,joinpath(pwd(),"src"))
#------------------------------------------------------------------------------#

using SOH   # import the package


#-------- Model parameters ----------------------------------------------------#
# The coefficients are computed using the function `coefficients_Vicsek` in the 
# script `toolbox.jl` for the Fokker-Planck and the BGK models. 

κ = 5.0     # concentration parameter
c1,c2,λ = coefficients_Vicsek(κ)

#-------- Domain parameters ---------------------------------------------------#

# Rectangular domain of size Lx*Ly
Lx = 10.0   
Ly = 10.0

# Boundary conditions (possible choices "periodic", "Neumann", "reflecting")
bcond_x = "reflecting"
bcond_y = "reflecting"


#-------- Numerical parameters ------------------------------------------------#

# Numer of cells and spatial step
ncellx = 400
ncelly = 400
Δx = Lx / ncellx
Δy = Ly / ncelly

# Time step 
Δt = 0.01

# Final time
T = 100.

# Method ("Roe" or "HLLE")
method = "HLLE"


#-------- Exterior force ------------------------------------------------------#

# If no exterior force:
# Fx = nothing
# Fy = nothing

# Otherwise define the x and y components as (ncellx+2)*(ncelly+2) matrices.
# See the examples in the script `init.jl`
Fx, Fy = flat_quadratic_potential_force(ncellx,ncelly,Δx,Δy,Lx/3,5.)

#-------- Saving parameters ---------------------------------------------------#

should_save = false
simu_name = "simu"
should_plot = false
step_plot = 2   # A plot every `step_plot` iterations
save_video = true   

#-------- Initial conditions --------------------------------------------------#

ρ,u,v = random_init(ncellx,ncelly,1.,0.,0.,2*pi,bcond_x=bcond_x,bcond_y=bcond_y)


#-------- Finally run the simulation ------------------------------------------#

run!(
    ρ,u,v;
    Lx=Lx,Ly=Ly,Δt=Δt,
    c1=c1, c2=c2, λ=λ,
    final_time=T,
    bcond_x=bcond_x,bcond_y=bcond_y,
    Fx=Fx,Fy=Fy,
    method=method,
    simu_name="simu",
    should_save=should_save,save_step=1,
    should_plot=should_plot,plot_step=2,
    save_video=save_video,
    range=(0.,5.),resolution=(800,600),theme="dark",fps=60
);




