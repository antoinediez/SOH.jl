module PlotSave
using CairoMakie, ProgressMeter, Statistics
export plot_rhoUV, update_plot!, save_data!, radial_density

"""
    plot_rhoUV(
        ρ,u,v,
        Δx,Δy,
        range,
        save_plot=true,save_video=false,dir_name=pwd(),file_name="0.png"; 
        theme="dark", resolution=(800,600),arrow_step=20,fps=40,kwargs...
    )

Return the figure, axes, heatmap and arrows associated to the density `ρ` and velocity field
`(u,v)`. The plot is saved in the directory `dir_name` in the file `file_name`. Two themes 
are available, either `"dark"` (default) or `"light"`. 
"""
function plot_rhoUV(
    ρ,u,v,
    Δx,Δy,
    range,
    save_plot=true,save_video=false,dir_name=pwd(),file_name="0.png"; 
    theme="dark", resolution=(800,600),arrow_step=20,fps=40,kwargs...
)
    if theme == "dark"
        set_theme!(theme_black())
        colormap = :linear_kryw_5_100_c67_n256
        arrowcolor = :black
    elseif theme == "light"
        set_theme!()
        colormap = Reverse(:linear_kryw_5_100_c67_n256)
        arrowcolor = :white
    else
        AttributeError("Theme not defined!")
    end
    ncellx = size(ρ)[1] - 2
    ncelly = size(ρ)[2] - 2
    Lx = ncellx*Δx
    Ly = ncelly*Δy
    x = Δx .* collect(0:(ncellx-1)) .+ Δx/2
    y = Δy .* collect(0:(ncelly-1)) .+ Δy/2
    rho = ρ[2:ncellx+1,2:ncelly+1]
    fig = Figure(resolution=resolution)
    ax = Axis(fig[1, 1], title="time=0.00")
    hm = heatmap!(ax,x,y,rho,colorrange=range,colormap=colormap)
    ax.aspect = AxisAspect(1)
    xlims!(ax,0,Lx)
    ylims!(ax,0,Ly)
    Colorbar(fig[1,2],hm)

    extract_x = floor(Int,arrow_step/2):arrow_step:ncellx
    extract_y = floor(Int,arrow_step/2):arrow_step:ncelly
    xx = x[extract_x]
    yy = y[extract_y]
    u_grid = u[2:ncellx+1,2:ncelly+1]
    v_grid = v[2:ncellx+1,2:ncelly+1]
    uu = u_grid[extract_x,extract_y]
    vv = v_grid[extract_x,extract_y]
    arrows = arrows!(ax,xx,yy,uu,vv)
    arrows.color = arrowcolor
    arrows.lengthscale = 2/3*arrow_step*Δx
    arrows.arrowsize = floor(Int,1/80*resolution[1])
    arrows.linewidth = floor(Int,3/800*resolution[1])
    arrows.origin = :center
    if save_plot
        save(joinpath(dir_name,file_name),fig)
    end
    if save_video
        stream = VideoStream(fig,framerate=fps)
        return fig, ax, hm, arrows, stream
    else
        return fig, ax, hm, arrows
    end
end

"""
    update_plot!(fig,ax,hm,arrows,ρ,u,v,title,dir_name,file_name,save_plot=true,stream=nothing)

Update a plot created with the function `plot_rhoUV` and save it. 
"""
function update_plot!(fig,ax,hm,arrows,ρ,u,v,title,dir_name,file_name,save_plot=true,stream=nothing)
    ncellx = size(ρ)[1] - 2
    ncelly = size(ρ)[2] - 2
    rho = ρ[2:ncellx+1,2:ncelly+1]
    hm[3] = rho
    u_grid = u[2:ncellx+1,2:ncelly+1]
    v_grid = v[2:ncellx+1,2:ncelly+1]
    extract_x = 10:20:ncellx
    extract_y = 10:20:ncelly
    uu = u_grid[extract_x,extract_y]
    vv = v_grid[extract_x,extract_y]
    arrows[:directions] = vec(Vec2f.(uu,vv))
    ax.title = title
    if save_plot
        save(joinpath(dir_name,file_name),fig)
    end
    if !isnothing(stream)
        recordframe!(stream)
    end
end

"""
    save_data!(
        data,
        ρ,u,v,
        dir_name,file_name;key="data")

Update the dictionary `data` with the values `(ρ,u,v)` and save it. 
"""
function save_data!(
    data,
    ρ,u,v,
    dir_name,file_name;key="data")
    data[:ρ] = ρ
    data[:u] = u
    data[:v] = v
    save(joinpath(dir_name,file_name),key,data)
end

"""
    radial_density(ρ,Δx,Δy,K=400)

Compute the radial density as the mean of `ρ` in `K` annuli of 
equal lengths and evenly spread radii between 0 and the maximal 
radius in the box. Return the vector of radii and the radial density. 
"""
function radial_density(ρ,Δx,Δy,K=400)
    ncellx = size(ρ)[1] - 2
    ncelly = size(ρ)[2] - 2
    Lx = ncellx*Δx
    Ly = ncelly*Δy
    mx = Lx/2
    my = Ly/2
    rmat = zeros(ncellx,ncelly)
    for i in 1:ncellx
        for j in 1:ncellx
            xij = (i-1)*Δx + Δx/2
            yij = (j-1)*Δy + Δy/2
            rmat[i,j] = sqrt((xij-mx)^2 + (yij-my)^2)
        end
    end
    rvec = vec(rmat)
    rho = vec(ρ[2:ncellx+1,2:ncelly+1])
    r = LinRange(0.,maximum(rvec),K)
    dr = r[2]-r[1]
    rho_r = zeros(K)
    for k in 1:K
        rk = r[k]+dr
        rho_r[k] = mean(rho[abs.(rvec.-rk).<dr/2])
    end
    return r, rho_r
end

end    #module

