using Random, Printf, Plots
using Oceananigans, Oceananigans.Utils

### Oceananigans Setup with heat induced convection ###
Nz = 32       # Number of grid points in x, y, z
Δz = 1.0      # Grid spacing in x, y, z (meters)
Qᵀ = 5e-5     # Temperature flux at surface
∂T∂z = 0.005    # Initial vertical temperature gradient
evaporation = 1e-7     # Mass-specific evaporation rate [m s⁻¹]
f = 1e-4     # Coriolis parameter
α = 2e-4     # Thermal expansion coefficient
β = 8e-4     # Haline contraction coefficient

grid = RegularCartesianGrid(size=(Nz, Nz, Nz), length=(Δz*Nz, Δz*Nz, Δz*Nz))
T_bcs = TracerBoundaryConditions(grid, top = BoundaryCondition(Flux, Qᵀ),
                                 bottom = BoundaryCondition(Gradient, ∂T∂z))

model = IncompressibleModel(
         architecture = CPU(),
                 grid = RegularCartesianGrid(size=(Nz, Nz, Nz), length=(Δz*Nz, Δz*Nz, Δz*Nz)),
             coriolis = FPlane(f=f),
             buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=α, β=β)),
              closure = AnisotropicMinimumDissipation(),
  boundary_conditions = (T=T_bcs,),
           parameters = (evaporation = evaporation,)
)

## Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise

## Temperature initial condition: a stable density tradient with random noise superposed.
T₀(x, y, z) = 20 + ∂T∂z * z + ∂T∂z * model.grid.Lz * 1e-6 * Ξ(z)

set!(model, T=T₀)
wizard = TimeStepWizard(cfl=0.2, Δt=1.0, max_change=1.1, max_Δt=5.0)
simulation = Simulation(model, Δt=wizard, stop_iteration=0, progress_frequency=1)

## generate particles
P = []
for i in 1:100
    particle = Particles.Particle(model.grid)
    push!(P, particle)
end

## run the model
for i in 1:200
    simulation.stop_iteration += 1
    run!(simulation)
    vel_itp = Particles.generate_vel_itp(model)
    for j in 1:100
        Particles.particle_advection(P[j],vel_itp,model.grid,1.0)
    end
end
