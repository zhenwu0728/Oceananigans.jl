module Particles

mutable struct Particle
    x :: Array
    y :: Array
    z :: Array
    t :: Array
end

function Particle()
    x = Float64[]
    y = Float64[]
    z = Float64[]
    t = Float64[]
    return Particle(x,y,z,t)
end

function Base.push!(particle::Particle, x::Number, y::Number, z::Number, t::Number)
    push!(particle.x, x)
    push!(particle.y, y)
    push!(particle.z, z)
    push!(particle.t, t)
    return nothing
end

"""
    double_grid(velᵇ)
Compute double grid of each time step from Oceananigans grid (with halo point)
Change the grids from C grid to A grid, easier for interpolation
'velᵇ' is the velocitiy fields on big grids of current time step
"""
function double_grid_3D(velᵇ)
    uᵇ = velᵇ.u.data.parent
    vᵇ = velᵇ.v.data.parent
    wᵇ = velᵇ.w.data.parent
    Nx, Ny, Nz = size(uᵇ)
    ny2h = Ny*2-1; nx2h = Nx*2-1; nz2h = Nz*2-1;
    vel_sh = (nx2h, ny2h, nz2h);
    u2h = zeros(vel_sh);
    v2h = zeros(vel_sh);
    w2h = zeros(vel_sh);
    # Compute values for new grids
    uX = 0.5*(uᵇ[1:Nx-1,:,:] + uᵇ[2:Nx,:,:]);
    vY = 0.5*(vᵇ[:,1:Ny-1,:] + vᵇ[:,2:Ny,:]);
    wZ = 0.5*(wᵇ[:,:,1:Nz-1] + wᵇ[:,:,2:Nz]);
    u2h[2:2:nx2h,1:2:ny2h,1:2:nz2h] = uX;
    u2h[1:2:nx2h,1:2:ny2h,1:2:nz2h] = uᵇ;
    v2h[1:2:nx2h,2:2:ny2h,1:2:nz2h] = vY;
    v2h[1:2:nx2h,1:2:ny2h,1:2:nz2h] = vᵇ;
    w2h[1:2:nx2h,1:2:ny2h,2:2:nz2h] = wZ;
    w2h[1:2:nx2h,1:2:ny2h,1:2:nz2h] = wᵇ;
    uY = 0.5*(u2h[:,1:2:ny2h-1,:] + u2h[:,3:2:ny2h,:]);
    vX = 0.5*(v2h[1:2:nx2h-1,:,:] + v2h[3:2:nx2h,:,:]);
    wX = 0.5*(w2h[1:2:nx2h-1,:,:] + w2h[3:2:nx2h,:,:]);
    u2h[:,2:2:ny2h,:] = uY;
    v2h[2:2:nx2h,:,:] = vX;
    w2h[2:2:nx2h,:,:] = wX;
    uZ = 0.5*(u2h[:,:,1:2:nz2h-1] + u2h[:,:,3:2:nz2h]);
    vZ = 0.5*(v2h[:,:,1:2:nz2h-1] + v2h[:,:,3:2:nz2h]);
    wY = 0.5*(w2h[:,1:2:ny2h-1,:] + w2h[:,3:2:ny2h,:]);
    u2h[:,:,2:2:nz2h] = uZ;
    v2h[:,:,2:2:nz2h] = vZ;
    w2h[:,2:2:nz2h,:] = wY;
    # Delete the boundaries
    uvel = u2h[3:end,2:end-1,2:end-1];
    vvel = v2h[2:end-1,3:end,2:end-1];
    wvel = w2h[2:end-1,2:end-1,3:end];
    # velᵈ = velocity(uvel, vvel, wvel);
    # Need a function to return a new doubled grid & velocity fields
    return #velᵈ
end

"""
    simple_itpl(x, y, z, a)
Simple interpolation: interpolate according to C grid (velocity on faces), for 1D only
'x', 'y', 'z' are grid indices, 'a' is the w velocity field need to interpolate
"""
function simple_itpl(x, y, z, a)
    x₀, y₀, z₀ = trunc(Int,x), trunc(Int,y), trunc(Int,z)
    zᵈ = z - z₀
    w₋ = a[x₀, y₀, z₀]
    w₊ = a[x₀, y₀, z₀+1]
    vel = w₋ * (1 - zᵈ) + w₊ * zᵈ
    return vel
end

"""
    bilinear_itlp(x, y, z, a)
Bilinear interpolation of horizontal velocities
'x', 'y', 'z' are grid indices, 'a' is the velocity field need to interpolate, e.g. u, v
"""
function bilinear_itpl(x, y, z, a)
    x₀, y₀, z₀ = trunc(Int,x), trunc(Int,y), trunc(Int,z)
    xᵈ = x - x₀
    yᵈ = y - y₀
    vel_00 = a[x₀, y₀, z₀]
    vel_10 = a[x₀+1, y₀, z₀]
    vel_01 = a[x₀, y₀+1, z₀]
    vel_11 = a[x₀+1, y₀+1, z₀]
    vel_0 = vel_00 * (1 - yᵈ) + vel_10 * yᵈ
    vel_1 = vel_01 * (1 - yᵈ) + vel_11 * yᵈ
    vel = vel_0 * (1 - xᵈ) + vel_1 * xᵈ
    return vel
end

"""
    trilinear_itlp(x, y, z, a)
Trilinear interpolation of velocities
'x', 'y', 'z' are doubled grid indices(whether 2D or 3D),
'a' is the velocity field need to interpolate, e.g. u, v, w
"""
function trilinear_itpl(x, y, z, a)
    x₀, y₀, z₀ = trunc(Int,x), trunc(Int,y), trunc(Int,z)
    xᵈ = x - x₀
    yᵈ = y - y₀
    zᵈ = z - z₀
    vel_000 = a[x₀, y₀, z₀]
    vel_100 = a[x₀+1, y₀, z₀]
    vel_001 = a[x₀, y₀, z₀+1]
    vel_010 = a[x₀, y₀+1, z₀]
    vel_110 = a[x₀+1, y₀+1, z₀]
    vel_011 = a[x₀, y₀+1, z₀+1]
    vel_101 = a[x₀+1, y₀, z₀+1]
    vel_111 = a[x₀+1, y₀+1, z₀+1]
    vel_00 = vel_000 * (1 - yᵈ) + vel_100 * yᵈ
    vel_01 = vel_001 * (1 - yᵈ) + vel_101 * yᵈ
    vel_10 = vel_010 * (1 - yᵈ) + vel_110 * yᵈ
    vel_11 = vel_011 * (1 - yᵈ) + vel_111 * yᵈ
    vel_0 = vel_00 * (1 - xᵈ) + vel_10 * xᵈ
    vel_1 = vel_01 * (1 - xᵈ) + vel_11 * xᵈ
    vel = vel_0 * (1-zᵈ) + vel_1 * zᵈ
    return vel
end

"""
    get_vels(x, y, z, velᵈ)
Read velocities at (x,y,z) from velocity fields
'x', 'y', 'z' are original grid indices
the velocity field passed to the function can be boubled grids or original grids
"""
function get_vels(x, y, z, vels)
    u = vels.u.data.parent
    v = vels.v.data.parent
    w = vels.w.data.parent
    g = vels.u.grid
    if (g.Nx > 1) & (g.Ny > 1) & (g.Nz > 1)
        uvel = trilinear_itpl(2*x-1, 2*y-1, 2*z-1, vels.u)
        vvel = trilinear_itpl(2*x-1, 2*y-1, 2*z-1, vels.v)
        wvel = trilinear_itpl(2*x-1, 2*y-1, 2*z-1, vels.w)
    elseif (g.Nx > 1) & (g.Ny > 1) & (g.Nz = 1)
        uvel = bilinear_itpl(2*x-1, 2*y-1, z, vels.u) # unit: m/s, bilinear interpolation
        vvel = bilinear_itpl(2*x-1, 2*y-1, z, vels.v) # unit: m/s, bilinear interpolation
        wvel = 0.0
    elseif (g.Nx = 1) & (g.Ny = 1) & (g.Nz > 1)
        uvel = 0.0; vvel = 0.0
        wvel = simple_itpl(x,y,z, vels.w) # unit: m/s, simple interpolation
    else
        return "only vertical 1D, horizontal 2D and 3D are supported for now"
    end
    return uvel, vvel, wvel
end

function periodic_domain(Nx, x)
    if 1 < x < Nx+1
        return x
    elseif x ≥ Nx+1
        return x - Nx
    elseif x ≤ 1
        return x + Nx
    end
end

"""
    agent_advection(phyts_a,vel,ΔT)
Update grid indices of all the individuals according to velocity fields of each time step
Periodic domain is used
'vel' is 'model.velocities' in Oceananigans
"""
function agent_advection(particle,vels,ΔT::Int64)
    g = vels.u.grid
    x, y, z = copy(particle.x[end]), copy(particle.y[end]), copy(particle.z[end])
    uvel, vvel, wvel = get_vels(x, y, z, vels)
    xi, yi, zi = trunc(Int,x), trunc(Int,y), trunc(Int,z)
    dx = uvel/(g.xF[xi-1]-g.xF[xi])*ΔT # unit: grid/h
    dy = vvel/(g.yF[yi-1]-g.yF[yi])*ΔT # unit: grid/h
    dz = wvel/(g.zF[zi-1]-g.zF[zi])*ΔT # vertical movement, plus sinking, unit: grid/h
    x = x + dx
    y = y + dy
    z = max(1.1,min(g.Nz-0.1,z + dz))
    # periodic domain
    x = periodic_domain(g.Nx, x)
    y = periodic_domain(g.Ny, y)
    return x, y, z
end
"""
    agent_advectionRK4(phyts_a, vel_field, ΔT::Int64)
Require 3D doubled grids.
'vel' is 'model.velocities' in Oceananigans
"""
# not finished yet, RK4 needs 3 time steps of velocity fields
function agent_advectionRK4(particle, vels, ΔT::Int64, grid_type::String)
    g = vels.u.grid
    x, y, z = copy(particle.x[end]), copy(particle.y[end]), copy(particle.z[end])
    u1,v1,w1 = get_vels(x, y, z, g, vel_field[1], grid_type) # velocites at t
    xi1, yi1, zi1 = trunc(Int,x),trunc(Int,y),trunc(Int,z)
    gx1 = periodic_domain(g.Nx, x + u1/g.dxF[xi1,1]*0.5*ΔT)
    gy1 = periodic_domain(g.Ny, y + v1/g.dyF[1,yi1]*0.5*ΔT)
    gz1 = z + w1/g.dzF[zi1]*0.5*ΔT
    gz1 = max(2.0,min(g.Nz-0.1,gz1))
    u2,v2,w2 = get_vels(gx1, gy1, gz1, g, vel_field[2], grid_type) # velocites at t+0.5ΔT
    xi2, yi2, zi2 = trunc(Int,gx1),trunc(Int,gy1),trunc(Int,gz1)
    gx2 = periodic_domain(g.Nx, x + u2/g.dxF[xi2,1]*0.5*ΔT)
    gy2 = periodic_domain(g.Ny, y + v2/g.dyF[1,yi2]*0.5*ΔT)
    gz2 = z + w2/g.dzF[zi2]*0.5*ΔT
    gz2 = max(2.0,min(g.Nz-0.1,gz2))
    u3,v3,w3 = get_vels(gx2, gy2, gz2, g, vel_field[2], grid_type) # velocites at t+0.5ΔT
    xi3, yi3, zi3 = trunc(Int,gx2),trunc(Int,gy2),trunc(Int,gz2)
    gx3 = periodic_domain(g.Nx, x + u3/g.dxF[xi3,1]*0.5*ΔT)
    gy3 = periodic_domain(g.Ny, y + v3/g.dyF[1,yi3]*0.5*ΔT)
    gz3 = z + w3/g.dzF[zi3]*0.5*ΔT
    gz3 = max(2.0,min(g.Nz-0.1,gz3))
    u4,v4,w4 = get_vels(gx3, gy3, gz3, g, vel_field[3], grid_type) # velocites at t+ΔT
    xi4, yi4, zi4 = trunc(Int,gx3),trunc(Int,gy3),trunc(Int,gz3)
    dx = (u1/g.dxF[xi1,1] + 2*u2/g.dxF[xi2,1] + 2*u3/g.dxF[xi3,1] + u4/g.dxF[xi4,1]) / 6 * ΔT
    dy = (v1/g.dyF[1,yi1] + 2*v2/g.dyF[1,yi2] + 2*v3/g.dyF[1,yi3] + v4/g.dyF[1,yi4]) / 6 * ΔT
    dz = (w1/g.dzF[zi1] + 2*w2/g.dzF[zi2] + 2*w3/g.dzF[zi3] + w4/g.dzF[zi4]) / 6 * ΔT
    x = periodic_domain(g.Nx, x + dx)
    y = periodic_domain(g.Ny, y + dy)
    z = max(1.1, min(g.Nz-0.1, z + dz))
    return x, y, z
end

"""
    agent_diffusionH(particle,grid,κh)
Using a random walk algorithm for horizontal diffusion
"""
function agent_diffusionH(particle,grid,κh)
    xi, yi = trunc(Int,particle.x[end]),trunc(Int,particle.y[end])
    x = particle.x[end] + rand(Uniform(-1.0,1.0)) * κh / (grid.xF[xi-1] - grid.xF[xi])
    y = particle.x[end] + rand(Uniform(-1.0,1.0)) * κh / (grid.yF[yi-1] - grid.yF[yi])
    x = periodic_domain(grid.Nx, x)
    y = periodic_domain(grid.Ny, y)
    return x, y
end

"""
    agent_diffusionV(particle,grid,κv)
Using a random walk algorithm for vertical diffusion
"""
function agent_diffusionV(particle,grid,κv)
    zi= trunc(Int,particle.z[end])
    z = particle.z[end] + rand(Uniform(-1.0,1.0)) * κv / (grid.zF[zi-1] - grid.zF[zi])
    z = max(1.1, min(grid.Nz-0.1, z))
    return z
end

end
