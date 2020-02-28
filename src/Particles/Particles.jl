module Particles
using Distributions
using Interpolations

mutable struct Particle
    x :: Number
    y :: Number
    z :: Number
    t :: Number
end

function Particle(grid)
    x = rand(Uniform(grid.xF[1],grid.xF[end]))
    y = rand(Uniform(grid.yF[1],grid.yF[end]))
    z = rand(Uniform(grid.zF[1],grid.zF[end]))
    t = 0.0
    return Particle(x,y,z,t)
end

"""
    generate_vel_itp(model)
Use Interpolations.jl to generate interpolation objects
"""
function generate_vel_itp(model)
    grid = model.grid
    u = model.velocities.u.data.parent
    v = model.velocities.v.data.parent
    w = model.velocities.w.data.parent
    # deal with halo points
    xF = collect(grid.xF)
    pushfirst!(xF,xF[1]-(xF[2]-xF[1]))
    yF = collect(grid.yF)
    pushfirst!(yF,yF[1]-(yF[2]-yF[1]))
    zF = collect(grid.zF)
    pushfirst!(zF,zF[1]-(zF[2]-zF[1]))
    pushfirst!(zF,zF[1]-(zF[2]-zF[1]))

    xC = collect(grid.xC)
    pushfirst!(xC,xF[2]-(xC[1]-xF[2]))
    push!(xC,xF[end]+(xF[end]-xC[end]))
    yC = collect(grid.yC)
    pushfirst!(yC,yF[2]-(yC[1]-yF[2]))
    push!(yC,yF[end]+(yF[end]-yC[end]))
    zC = collect(grid.zC)
    pushfirst!(zC,zF[2]-(zC[1]-zF[2]))
    push!(zC,zF[end]+(zF[end]-zC[end]))

    u_itp = interpolate((xF,yC,zC),u,Gridded(Linear()))
    v_itp = interpolate((xC,yF,zC),v,Gridded(Linear()))
    w_itp = interpolate((xC,yC,zF),w,Gridded(Linear()))
    return (u_itp, v_itp, w_itp)
end

"""
    get_vels(x, y, z, vel_itp)
Read and interpolate velocity at '(x,y,z)'
'(x,y,z)' is actual location of an individual
'vel_itp' is interpolation objects
"""
function get_vels(x, y, z, vel_itp)
    u_itp = vel_itp[1]; v_itp = vel_itp[2]; w_itp = vel_itp[3];
    u = u_itp(x, y, z); v = v_itp(x, y, z); w = w_itp(x, y, z);
    return u, v, w
end

"""
    periodic_domain(xF, x)
'xF' is the array containing the faces of grids
'x' is the coordinate of an individual
"""
function periodic_domain(xF, x)
    xF₀ = xF[1]; xFₜ= xF[end]
    if xF₀ < x < xFₜ
        return x
    elseif x ≥ xFₜ
        return x - xFₜ
    elseif x ≤ xF₀
        return x + xFₜ
    end
end

"""
    particle_advection(particle,vel_itp,g,ΔT)
Update grid indices of all the individuals according to velocity fields of each time step
Periodic domain is used
'vel_itp' is the tuple contains interpolations of u, v, w velocites of current time step
'g' is the grid information and 'ΔT' is time step
"""
function particle_advection(particle,vel_itp,g,ΔT)
    uvel, vvel, wvel = get_vels(particle.x, particle.y, particle.z, vel_itp)
    particle.x = particle.x + uvel*ΔT
    particle.y = particle.y + vvel*ΔT
    particle.z = max(g.zF[1],min(g.zF[end],particle.z + wvel*ΔT))
    # periodic domain
    particle.x = periodic_domain(g.xF, particle.x)
    particle.y = periodic_domain(g.yF, particle.y)
    return nothing
end
"""
    particle_advectionRK4(particle, vel_itps, g, ΔT::Int64)
'vel_itps' is an array of tuples containing interpolations of u, v, w velocites of current time step
"""
function particle_advectionRK4(particle, vel_itps, g, ΔT)
    u1,v1,w1 = get_vels(particle.x, particle.y, particle.z, vel_itps[1]) # velocites at t
    gx1 = periodic_domain(g.xF, particle.x + u1*0.5*ΔT)
    gy1 = periodic_domain(g.yF, particle.y + v1*0.5*ΔT)
    gz1 = particle.z + w1*0.5*ΔT
    gz1 = max(g.zF[1],min(g.zF[end],gz1))
    u2,v2,w2 = get_vels(gx1, gy1, gz1, vel_itps[2]) # velocites at t+0.5ΔT
    gx2 = periodic_domain(g.xF, particle.x + u2*0.5*ΔT)
    gy2 = periodic_domain(g.yF, particle.y + v2*0.5*ΔT)
    gz2 = particle.z + w2*0.5*ΔT
    gz2 = max(g.zF[1],min(g.zF[end],gz2))
    u3,v3,w3 = get_vels(gx2, gy2, gz2, vel_itps[2]) # velocites at t+0.5ΔT
    gx3 = periodic_domain(g.xF, particle.x + u3*0.5*ΔT)
    gy3 = periodic_domain(g.yF, particle.y + v3*0.5*ΔT)
    gz3 = particle.z + w3*0.5*ΔT
    gz3 = max(g.zF[1],min(g.zF[end],gz3))
    u4,v4,w4 = get_vels(gx3, gy3, gz3, vel_itps[3]) # velocites at t+ΔT
    dx = (u1 + 2*u2 + 2*u3 + u4) / 6 * ΔT
    dy = (v1 + 2*v2 + 2*v3 + v4) / 6 * ΔT
    dz = (w1 + 2*w2 + 2*w3 + w4) / 6 * ΔT
    particle.x = periodic_domain(g.xF, particle.x + dx)
    particle.y = periodic_domain(g.yF, particle.y + dy)
    particle.z = particle.z + dz
    particle.z = max(g.zF[1], min(g.zF[end], particle.z))
    return nothing
end

"""
    agent_diffusionX(particle,g,κh)
Using a random walk algorithm for horizontal diffusion
"""
function agent_diffusionX(particle,g,κh)
    particle.x += rand(Uniform(-1.0,1.0)) * κh
    particle.x = periodic_domain(g.xF, particle.x)
    return nothing
end

"""
    agent_diffusionX(particle,g,κh)
Using a random walk algorithm for horizontal diffusion
"""
function agent_diffusionY(particle,g,κh)
    particle.y += rand(Uniform(-1.0,1.0)) * κh
    particle.y = periodic_domain(g.yF, particle.y)
    return nothing
end
"""
    agent_diffusionZ(particle,g,κv)
Using a random walk algorithm for vertical diffusion
"""
function agent_diffusionZ(particle,g,κv)
    particle.z += rand(Uniform(-1.0,1.0)) * κv
    particle.z = max(g.zF[1], min(g.zF[end], particle.z))
    return nothing
end

end
