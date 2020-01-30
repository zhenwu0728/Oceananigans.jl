using Oceananigans, Oceananigans.OutputWriters, Printf

#####
##### Model set-up
#####

Nx, Ny, Nz = 128, 128, 20              # No. of grid points in x, y, and z, respectively.
Lx, Ly = 1000000, 1000000    # Length of the domain in x, y, and z, respectively (m).
const Lz = 4000
tf = 24*3600*360*3                     # Length of the simulation (s)


model = Model(     grid=RegularCartesianGrid(size=(Nx, Ny, Nz), x=(-500000, 500000), y=(-500000, 500000), z=(0, Lz)),
                   boundary_conditions = ShoeBoxSolutionBCs(),
                   coriolis=BetaPlane(latitude=0),
#                   forcing = forcing,
                   closure=ConstantAnisotropicDiffusivity(νh=1e2, κh=1e2, νv=1e-5, κv=1e-5))

# Add a cube-shaped warm temperature anomaly that takes up the middle 50%
# of the domain volume.
const radius = Ly/10

T(x,y,z) = 0.01*exp(-(x^2+y^2)/radius^2)*sin(π*z/Lz)
#T(x,y,z) = 0.01*exp(-(x^2+y^2)/radius^2)*(z > 0.96Lz)
set!(model, T=T)
#i1, i2 = round(Int, Nx/4), round(Int, 3Nx/4)
#j1, j2 = round(Int, Ny/4), round(Int, 3Ny/4)
#model.tracers.T.data[i1:i2, j1:j2, :] .+= 0.01

#####
##### Set up output
#####

# write_grid(model)

outputs =      Dict("u" => model.velocities.u,
                    "v" => model.velocities.v,
                    "w" => model.velocities.w,
                    "T" => model.tracers.T,
                    "S" => model.tracers.S)

outputattrib = Dict("u" => ["longname" => "Velocity in the x-direction", "units" => "m/s"],
                    "v" => ["longname" => "Velocity in the y-direction", "units" => "m/s"],
                    "w" => ["longname" => "Velocity in the z-direction", "units" => "m/s"],
                    "T" => ["longname" => "Temperature", "units" => "K"],
                    "S" => ["longname" => "Salinity", "units" => "g/kg"])

globalattrib = Dict("name" => "Equatorial waves expt 1")

# The following writer saves a yz plane at xC[5] for all fields that
# have xC as their dimension and at xF[6] for all fields that have xF
# as their dimension. Ranges also can be specified (e.g. xC=2:10)
#subsetwriter = NetCDFOutputWriter(model, outputs;
#                                  interval=10, filename="dump_subset.nc",
#                                  output_attributes=outputattrib,
#                                  global_attributes=globalattrib,
#                                  xC=5, xF=6)
#push!(model.output_writers, subsetwriter)

# The following writer saves a data from the entire domain
globalwriter = NetCDFOutputWriter(model, outputs, interval=24*3600*10,
                                  filename="equatorial_waves_global.nc")
push!(model.output_writers, globalwriter)

#####
##### Run the simulation
#####

function terse_message(model, walltime, Δt)
    cfl = Δt / Oceananigans.cell_advection_timescale(model)
    return @sprintf("i: %d, t: %.4f y, Δt: %.1f s, cfl: %.3f, wall time: %s\n",
                    model.clock.iteration, model.clock.time/24/3600/360, Δt, cfl, prettytime(walltime))
end

# A wizard for managing the simulation time-step.
wizard = TimeStepWizard(cfl=0.3, Δt=300.0, max_change=1.1, max_Δt=7200.0)

# Run the model
while model.clock.time < tf
    update_Δt!(wizard, model)
    walltime = @elapsed time_step!(model, 50, wizard.Δt)
    @printf "%s" terse_message(model, walltime, wizard.Δt)
end

# Close the NetCDFOutputWriter
# close(subsetwriter)
close(globalwriter)
