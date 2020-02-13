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

end
