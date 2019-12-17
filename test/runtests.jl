using
    Oceananigans,
    Oceananigans.Operators,
    Oceananigans.Grids,
    Oceananigans.Solvers,
    Oceananigans.Diagnostics,
    Oceananigans.OutputWriters,
    Oceananigans.TurbulenceClosures,
    Oceananigans.AbstractOperations,
    Test,
    Random,
    JLD2,
    Printf,
    Statistics,
    OffsetArrays,
    FFTW

@hascuda begin
    import CUDAdrv
    using CUDAnative, CuArrays
end

using Statistics: mean
using LinearAlgebra: norm
using GPUifyLoops: @launch, @loop
using NCDatasets: Dataset

using Oceananigans: architecture, device, launch_config, datatuples, with_tracers,
                    Face, Cell, interiorparent, location, TracerFields, fill_halo_regions!

import Oceananigans: interior, datatuple

using Oceananigans.Solvers: PoissonSolver, PPN, PNN, solve_poisson_3d!

using Oceananigans.Diagnostics: run_diagnostic, velocity_div!

using Oceananigans.TimeSteppers: _compute_w_from_continuity!

using Oceananigans.AbstractOperations: Computation, compute!

# On CI servers select the GPU with the most available memory or with the
# highest capability if testing needs to be thorough).
# Source credit: https://github.com/JuliaGPU/CuArrays.jl/pull/526
@hascuda begin
    gpu_candidates = [(dev=dev, cap=CUDAdrv.capability(dev),
                       mem=CUDAdrv.CuContext(ctx->CUDAdrv.available_memory(), dev)) for dev in CUDAdrv.devices()]

    thorough = parse(Bool, get(ENV, "CI_THOROUGH", "false"))
    if thorough
        sort!(gpu_candidates, by=x->(x.cap, x.mem))
    else
        sort!(gpu_candidates, by=x->x.mem)
    end

    pick = last(gpu_candidates)
    device!(pick.dev)
end

datatuple(A) = NamedTuple{propertynames(A)}(Array(data(a)) for a in A)

function get_output_tuple(output, iter, tuplename)
    file = jldopen(output.filepath, "r")
    output_tuple = file["timeseries/$tuplename/$iter"]
    close(file)
    return output_tuple
end

float_types = (Float32, Float64)

archs = (CPU(),)
@hascuda archs = (CPU(), GPU())

closures = (
    :ConstantIsotropicDiffusivity,
    :ConstantAnisotropicDiffusivity,
    :AnisotropicBiharmonicDiffusivity,
    :TwoDimensionalLeith,
    :SmagorinskyLilly,
    :BlasiusSmagorinsky,
    :RozemaAnisotropicMinimumDissipation,
    :VerstappenAnisotropicMinimumDissipation
)

@testset "Oceananigans" begin
    include("test_grids.jl")
    include("test_fields.jl")
    include("test_halo_regions.jl")
    include("test_operators.jl")
    include("test_solvers.jl")
    include("test_coriolis.jl")
    include("test_surface_waves.jl")
    include("test_buoyancy.jl")
    include("test_models.jl")
    include("test_time_stepping.jl")
    include("test_boundary_conditions.jl")
    include("test_forcings.jl")
    include("test_turbulence_closures.jl")
    include("test_dynamics.jl")
    include("test_diagnostics.jl")
    include("test_output_writers.jl")
    include("test_regression.jl")
    include("test_examples.jl")
    include("test_abstract_operations.jl")
    include("test_verification.jl")
end
