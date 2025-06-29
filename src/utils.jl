function plot_bridge_with_supports(bo::BridgeOptions, supports::Vector{SupportElement}; 
                                 mode_shape=nothing, scale_factor=1.0, T=20.0,
                                 title="Bridge Structure", show_nodes=true)
    
    # Bridge nodes (assuming horizontal bridge)
    bridge_x = collect(range(0, bo.L, length=bo.n_nodes))
    bridge_y = zeros(bo.n_nodes)
    
    # If mode shape is provided, add deformation DIRECTLY
    if mode_shape !== nothing
        for i in 1:bo.n_nodes
            if 3*(i-1)+1 <= length(mode_shape) && 3*(i-1)+2 <= length(mode_shape)
                bridge_x[i] += scale_factor * mode_shape[3*(i-1)+1]  # x displacement
                bridge_y[i] += scale_factor * mode_shape[3*(i-1)+2]  # y displacement
            end
        end
    end
    
    # Plot bridge
    p = plot(bridge_x, bridge_y, linewidth=3, color=:blue, label="Bridge", 
             aspect_ratio=:equal, grid=true)
    
    if show_nodes
        scatter!(bridge_x, bridge_y, color=:blue, markersize=4, label="Bridge Nodes")
    end
    
    # Plot each support
    support_dof_maps, total_dofs = create_support_dof_mapping(bo, supports)
    
    for (i, support) in enumerate(supports)
        # Connection point on bridge
        conn_node = support.connection_node
        conn_x_undeformed = bo.L * (conn_node - 1) / (bo.n_nodes - 1)  # Undeformed bridge x-position
        conn_y_undeformed = 0.0  # Undeformed bridge y-position
        
        # Support geometry
        n_support_nodes = support.n_elem + 1
        angle_rad = deg2rad(support.angle)
        
        # UNDEFORMED support positions
        support_x_undeformed = zeros(n_support_nodes)
        support_y_undeformed = zeros(n_support_nodes)
        
        for j in 1:n_support_nodes
            # Distance along support from connection point
            distance = (j-1) * support.L / (n_support_nodes - 1)
            support_x_undeformed[j] = conn_x_undeformed + distance * cos(angle_rad)
            support_y_undeformed[j] = conn_y_undeformed + distance * sin(angle_rad)
        end
        
        # DEFORMED support positions (if mode shape provided)
        support_x_deformed = copy(support_x_undeformed)
        support_y_deformed = copy(support_y_undeformed)
        
        if mode_shape !== nothing
            dof_map = support_dof_maps[i]
            
            for j in 1:n_support_nodes
                local_x_idx = 3*(j-1) + 1
                local_y_idx = 3*(j-1) + 2
                
                if local_x_idx <= length(dof_map) && local_y_idx <= length(dof_map)
                    global_x_dof = dof_map[local_x_idx]
                    global_y_dof = dof_map[local_y_idx]
                    
                    # Add deformation DIRECTLY from global mode shape
                    if global_x_dof <= length(mode_shape)
                        support_x_deformed[j] += scale_factor * mode_shape[global_x_dof]
                    end
                    if global_y_dof <= length(mode_shape)
                        support_y_deformed[j] += scale_factor * mode_shape[global_y_dof]
                    end
                end
            end
        end
        
        # Plot support
        plot!(support_x_deformed, support_y_deformed, linewidth=2, color=:red, 
              label=i==1 ? "Supports" : "")
        
        if show_nodes
            scatter!(support_x_deformed, support_y_deformed, color=:red, markersize=3, 
                    label="")
        end
        
        # Fixed base at UNDEFORMED position (last node)
        scatter!([support_x_undeformed[end]], [support_y_undeformed[end]], 
                color=:black, marker=:square, markersize=8, 
                label=i==1 ? "Fixed Base" : "")
    end
    
    xlabel!("X Position (m)")
    ylabel!("Y Position (m)")
    title!(title)
    
    return p
end

function plot_mode_shape(bo::BridgeOptions, supports::Vector{SupportElement}, 
                        mode_shape::Vector{Float64}, mode_num::Int;
                        scale_factor=50.0, T=20.0)
    
    # Plot undeformed structure
    p1 = plot_bridge_with_supports(bo, supports, T=T, 
                                  title="Undeformed Structure", show_nodes=false)
    plot!(p1, linewidth=2, alpha=0.3, color=:gray)
    
    # Plot deformed structure (mode shape)
    p2 = plot_bridge_with_supports(bo, supports, mode_shape=mode_shape, 
                                  scale_factor=scale_factor, T=T,
                                  title="Mode $mode_num (Scaled by $scale_factor)", 
                                  show_nodes=true)
    
    # Combine plots
    plot(p1, p2, layout=(1,2), size=(1200, 400))
end

function animate_mode(bo::BridgeOptions, supports::Vector{SupportElement}, 
                     mode_shape::Vector{Float64}, mode_num::Int;
                     scale_factor=50.0, T=20.0, n_frames=30, fsize = (800,800))
    
    # Calculate sensible fixed limits
    bridge_extent = bo.L
    support_extent = maximum([s.L for s in supports], init=0.0)
    max_deformation = scale_factor * maximum(abs.(mode_shape), init=0.0)
    
    padding = 0.0
    x_limits = (-padding - support_extent, bridge_extent + padding + support_extent)
    y_limits = (-padding - support_extent - max_deformation, 
                padding + max_deformation)
    
    # Create animation with fixed layout
    anim = Animation()
    
    for i in 1:n_frames
        phase = 2π * (i-1) / n_frames
        current_scale = scale_factor * sin(phase)
        
        p = plot_bridge_with_supports(bo, supports, 
                                    mode_shape=mode_shape, 
                                    scale_factor=current_scale, 
                                    T=T,
                                    title="Mode $mode_num Animation",
                                    show_nodes=true)
        
        # Apply consistent layout
        plot!(p, 
            xlims=x_limits, 
            ylims=y_limits,
            aspect_ratio=1,        # Reduce this value to make plot wider
            size=fsize,        # Increase figure width
            legend=:topright,        
            legend_background_color=:white)
        
        frame(anim)
    end
    
    return anim
end

function animate_dynamic_response(bo::BridgeOptions, supports::Vector{SupportElement}, 
                                 u::Matrix{Float64}, time::Vector{Float64};
                                 scale_factor=1.0, n_frames=100, fps=20, 
                                 filename="bridge_dynamics.gif", show_nodes=false,
                                 title_prefix="Bridge Dynamic Response")
    
    # Calculate fixed limits based on maximum deformation
    bridge_extent = bo.L
    support_extent = maximum([s.L for s in supports], init=0.0)
    
    # Find maximum displacements for proper scaling
    max_x_disp = maximum(abs.(u[1:3:min(bo.n_dofs, size(u,1)), :]), init=0.0)
    max_y_disp = maximum(abs.(u[2:3:min(bo.n_dofs, size(u,1)), :]), init=0.0)
    max_deformation = max(max_x_disp, max_y_disp) * scale_factor
    
    padding = 20.0
    x_limits = (-padding - support_extent - max_deformation, 
                bridge_extent + padding + support_extent + max_deformation)
    y_limits = (-padding - support_extent - max_deformation, 
                padding + max_deformation)
    
    # Sample time indices for animation
    if n_frames > length(time)
        time_indices = 1:length(time)
    else
        time_indices = round.(Int, range(1, length(time), length=n_frames))
    end
    
    anim = Animation()
    
    @info "Creating animation with $(length(time_indices)) frames..."
    
    for (frame_num, i) in enumerate(time_indices)
        t = time[i]
        current_displacement = u[:, i]
        
        # Create plot with current displacement
        p = plot_bridge_with_supports(bo, supports, 
                                    mode_shape=current_displacement, 
                                    scale_factor=scale_factor,
                                    title="$title_prefix (t = $(round(t, digits=2))s)",
                                    show_nodes=show_nodes)
        
        # Apply fixed layout
        plot!(p, 
              xlims=x_limits, 
              ylims=y_limits,
              aspect_ratio=2.0,
              size=(1200, 400),
              legend=:topright,
              legend_background_color=:white,
              left_margin=3Plots.mm,
              right_margin=10Plots.mm,
              top_margin=3Plots.mm,
              bottom_margin=5Plots.mm)
        
        frame(anim)
        
        # Progress indicator
        if frame_num % 10 == 0
            @info "Processed frame $frame_num/$(length(time_indices))"
        end
    end
    
    @info "Saving animation: $filename"
    gif(anim, filename, fps=fps)
    return anim
end

# Convenience function that combines reconstruction and animation
function animate_from_modal_response(bo::BridgeOptions, supports::Vector{SupportElement},
                                   q::Matrix{Float64}, Φ_interp, T_func, time::Vector{Float64};
                                   scale_factor=1000.0, n_frames=100, fps=20,
                                   filename="bridge_modal_dynamics.gif", show_nodes=false)
    
    @info "Reconstructing physical displacements..."
    u, du = reconstruct_physical(bo, q, Φ_interp, T_func, time; supports=supports)
    
    @info "Creating animation..."
    return animate_dynamic_response(bo, supports, u, time;
                                  scale_factor=scale_factor, n_frames=n_frames, 
                                  fps=fps, filename=filename, show_nodes=show_nodes)
end

# Function to animate specific DOFs (e.g., only bridge or only supports)
function animate_bridge_response(bo::BridgeOptions, supports::Vector{SupportElement},
                               u::Matrix{Float64}, time::Vector{Float64}, dof_range::UnitRange{Int};
                               scale_factor=1000.0, n_frames=100, fps=20,
                               filename="bridge_partial_dynamics.gif")
    
    # Extract only the specified DOF range
    u_partial = u[dof_range, :]
    
    # Create a temporary displacement vector with zeros for non-selected DOFs
    u_full_temp = zeros(size(u, 1), size(u, 2))
    u_full_temp[dof_range, :] = u_partial
    
    return animate_dynamic_response(bo, supports, u_full_temp, time;
                                  scale_factor=scale_factor, n_frames=n_frames,
                                  fps=fps, filename=filename,
                                  title_prefix="Bridge Response (Partial DOFs)")
end