using JLD2, Dates, ProgressMeter
"""
    run!(
        ρ,u,v;
        Lx,Ly,Δt,
        c1, c2, λ,
        final_time,
        bcond_x="periodic",bcond_y="periodic",
        Fx=nothing,Fy=nothing,
        method="HLLE",
        simu_name="simu",
        should_save=false,save_step=1,
        should_plot=false,plot_step=1,
        save_video=false,fps=40,
        range=(0.,5.),
        kwargs...
    )

Run a simulation of the SOH system starting from the data `(ρ,u,v)`. The domain
grid has the following structure.

```
|------- |----------------|---------------|----------------------|--------|
| ghost  |    ghost       |       ...     |    ghost             | ghost  |
|--------║================|===============|======================║--------|
| ghost  ║  (2,ncelly+1)  |      ...      | (ncellx+1,ncelly+1)  ║ ghost  |
|--------║----------------|---------------|----------------------║--------|
| ghost  ║  (2,ncelly)    |      ...      | (ncellx+1,ncelly)    ║ ghost  |
|--------║----------------|---------------|----------------------║--------|
|        ║                |               |                      ║        |
|--------║----------------|---------------|----------------------║--------|
| ghost  ║    (2,3)       |      ...      |   (ncellx+1,3)       ║ ghost  |
|--------║----------------|---------------|----------------------║--------|
| ghost  ║    (2,2)       |      ...      |   (ncellx+1,2)       ║ ghost  |
---------║================|===============|======================║--------|
| ghost  |    ghost       |      ...      |    ghost             | ghost  |
|--------|----------------|---------------|----------------------|--------|
````

The ghost cells (first and last column and rows) are used for the boundary
conditions.

# Arguments
- `ρ` -- density (matrix of size `ncellx*ncelly`)
- `u` -- x-coordinate of the velocity (matrix of size `ncellx*ncelly`)
- `v` -- y-coordinate of the velocity (matrix of size `ncellx*ncelly`)
- `Lx` -- length of the domain along the x-axis
- `Ly` -- length of the domain along the y-axis
- `Δt` -- time-step
- `c1` -- coefficient c1
- `c2` -- coefficient c1
- `λ` -- coefficient λ
- `final_time` -- final time of the simulation
- `bcond_x` -- (optional, default : `"periodic"`) boundary condition along the x-axis
- `bcond_y` -- (optional, default : `"periodic"`) boundary condition along the y-axis
- `Fx` -- (optional, default : `nothing`) if specified, x-coordinate of the exterior force
- `Fy` -- (optional, default : `nothing`) if specified, y-coordinate of the exterior force
- `method` -- (optional, default : `"HLLE`) can be either `"Roe"` or `"HLLE"`
- `simu` -- (optional, default : `"simu"`) name of the simulation directory
- `should_save` -- (optional, default : `false`) flag for saving the data (ρ,u,v)
- `save_step` -- (optional, default : `1`) number of iterations between two saved data
- `should_plot` -- (optional, default : `false`) flag for plotting and saving the heatmap
- `plot_step` -- (optional, default : `1`) number of iterations between two saved plots
- `save_video` -- (optional, default : `false`) flag for saving a video
- `fps` -- (optional, default : `40`) frame per seconds
- `range` -- (optional, default : `(0.0,5.0)`) range of the density on the heatmap
- `kwargs` -- Keyword arguments to be passed in the function `plot_rhoUV`
"""
function run!(
    ρ,u,v;
    Lx,Ly,Δt,
    c1, c2, λ,
    final_time,
    bcond_x="periodic",bcond_y="periodic",
    Fx=nothing,Fy=nothing,
    method="HLLE",
    simu_name="simu",
    should_save=false,save_step=1,
    should_plot=false,plot_step=1,
    save_video=false,fps=40,
    range=(0.,5.),
    kwargs...
)
    #=======================================================================================#
    #============ Collect and display the parameters =======================================#
    #=======================================================================================#

    ncellx = size(ρ)[1] - 2
    Δx = Lx/ncellx
    ncelly = size(ρ)[2] - 2
    Δy = Ly/ncelly
    parameters = Dict{Symbol,Any}(
        :Lx => Lx, :Ly => Ly, :ncellx => ncellx, :ncelly => ncelly,
        :Δx => Δx, :Δy => Δy, :Δt => Δt,
        :c1 => c1, :c2 => c2, :λ => λ, :CFL => c1 * Δt/Δx * sqrt(λ/(c1-c2)),
        :final_time => final_time,
        :bcond_x => bcond_x, :bcond_y => bcond_y,
        :Fx => Fx, :Fy => Fy,
        :method => method,
        :date => string(Dates.now())
    )

    println("\n************* Model parameters **************\n")
    println("c1 = $(parameters[:c1])")
    println("c2 = $(parameters[:c2])")
    println("λ = $(parameters[:λ])")
    println("\n*********************************************\n")

    println("************* Domain parameters *************\n")
    println("Lx = $(parameters[:Lx])")
    println("Ly = $(parameters[:Ly])")
    println("Boundary condition in x = $(parameters[:bcond_x])")
    println("Boundary condition in y = $(parameters[:bcond_y])")
    println("Exterior force : $(!isnothing(parameters[:Fx]) && !isnothing(parameters[:Fy]))")
    println("\n*********************************************\n")

    println("************* Numerical method **************\n")
    println("Δx = $(parameters[:Δx])")
    println("Δy = $(parameters[:Δy])")
    println("Δt = $(parameters[:Δt])")
    println("CFL = $(parameters[:CFL])")
    println("Scheme : $(parameters[:method])")
    println("Final time : $(parameters[:final_time])")
    println("\n*********************************************\n")

    println("************* Saving parameters *************\n")
    if should_save
        if save_step == 1
            println("Save all data")
        else
            println("Save data every $save_step iterations")
        end
        estimated_size = round(8*ncellx*ncelly*3*floor(Int, final_time / Δt)/save_step*1e-9, digits=2)
        println("Estimated size : $estimated_size GB")
    else
        println("Only the initial and final data will be saved")
    end
    if should_plot
        if plot_step == 1
            println("Plot all data")
        else
            println("A plot will be save every $plot_step iterations")
        end
    else
        println("No plot will be saved")
    end
    if save_video
        println("A video with $fps fps will be saved")
    else
        println("No video will be saved")
    end
    println("\n*********************************************\n")

    #=======================================================================================#
    #============ Create the simulation directory and save the initial data ================#
    #=======================================================================================#

    dir_name = make_new_dir(simu_name)
    data_dir = mkdir(joinpath(dir_name,"data"))
    data = Dict(:ρ => ρ, :u => u, :v => v)
    save(joinpath(data_dir,"data_0.jld2"),"iter_0",data)

    if should_plot || save_video
        print("Initializing plot... ")
        plot_dir = make_new_dir(joinpath(dir_name,"plots"))
        if save_video
            fig, ax, hm, arrows, stream = plot_rhoUV(ρ,u,v,Δx,Δy,range,should_plot,save_video,plot_dir,"0.png";fps=fps,kwargs...)
        else
            fig, ax, hm, arrows = plot_rhoUV(ρ,u,v,Δx,Δy,range,should_plot,save_video,plot_dir,"0.png";kwargs...)
        end
        println("Done.\n")
    end

    #=======================================================================================#
    #============ The big loop =============================================================#
    #=======================================================================================#

    println("Run the simulation...\n")
    ntime = floor(Int, final_time / Δt)
    p = Progress(ntime)
    for itime in 1:ntime
        scheme_iter!(ρ,u,v,Δx,Δy,Δt,c1,c2,λ,bcond_x,bcond_y,Fx,Fy,method)
        if should_save
            if itime%save_step == 0
                save_data!(data,ρ,u,v,data_dir,"data_$itime.jld2",key="iter_$itime")
            end
        end
        if should_plot || save_video
            if itime%plot_step == 0
                round_time = nice_float2string(itime*Δt,2)
                if save_video
                    update_plot!(fig,ax,hm,arrows,ρ,u,v,"time=$round_time",plot_dir,"$itime.png",should_plot,stream)
                else
                    update_plot!(fig,ax,hm,arrows,ρ,u,v,"time=$round_time",plot_dir,"$itime.png")
                end
            end
        end
        next!(p)
    end

    #=======================================================================================#
    #============ Collect the simulation time and save the final data ======================#
    #=======================================================================================#

    parameters[:duration] = p.tlast-p.tinit
    save(joinpath(dir_name,"parameters.jld2"), "parameters", parameters)
    if !should_save
        save_data!(data,ρ,u,v,data_dir,"data_$ntime.jld2",key="iter_$ntime")
    end

    if should_plot || save_video
        if save_video
            save(joinpath(dir_name,"video.mp4"), stream)
            return fig, ax, hm, arrows, stream
        else
            return fig, ax, hm, arrows
        end
    end


end
