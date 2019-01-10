var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#Oceananigans.jl-Documentation-1",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "#Oceananigans.RegularCartesianGrid",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.RegularCartesianGrid",
    "category": "type",
    "text": "RegularCartesianGrid{T::AbstractFloat} <: Grid\n\nA Cartesian grid with regularly spaces cells and faces so that Δx, Δy, and Δz are constants. Fields are stored using floating-point values of type T.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.RegularCartesianGrid-Tuple{Any,Any}",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.RegularCartesianGrid",
    "category": "method",
    "text": "RegularCartesianGrid(N, L; dim=3, FloatType=Float64)\n\nCreate a regular Cartesian grid with size N = (N_x N_y N_z) and domain size L = (L_x L_y L_z) where fields are stored using floating-point values of type T.\n\nExamples\n\njulia> g = RegularCartesianGrid((16, 16, 8), (2π, 2π, 2π))\n\n\n\n\n\n"
},

{
    "location": "#Grids-1",
    "page": "Oceananigans.jl Documentation",
    "title": "Grids",
    "category": "section",
    "text": "RegularCartesianGrid\nRegularCartesianGrid(N, L; dim=3, FloatType=Float64)"
},

{
    "location": "#Oceananigans.CellField",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.CellField",
    "category": "type",
    "text": "CellField{T<:AbstractArray} <: Field\n\nA cell-centered field defined on a grid G whose values are stored as floating-point values of type T.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.FaceFieldX",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.FaceFieldX",
    "category": "type",
    "text": "FaceFieldX{T<:AbstractArray} <: FaceField\n\nAn x-face-centered field defined on a grid G whose values are stored as floating-point values of type T.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.FaceFieldY",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.FaceFieldY",
    "category": "type",
    "text": "FaceFieldY{T<:AbstractArray} <: FaceField\n\nA y-face-centered field defined on a grid G whose values are stored as floating-point values of type T.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.FaceFieldZ",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.FaceFieldZ",
    "category": "type",
    "text": "FaceFieldZ{T<:AbstractArray} <: FaceField\n\nA z-face-centered field defined on a grid G whose values are stored as floating-point values of type T.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.CellField-Tuple{Grid}",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.CellField",
    "category": "method",
    "text": "CellField(grid::Grid)\n\nConstruct a CellField whose values are defined at the center of a cell.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.FaceFieldX-Tuple{Grid}",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.FaceFieldX",
    "category": "method",
    "text": "FaceFieldX(grid::Grid)\n\nA Field whose values are defined on the x-face of a cell.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.FaceFieldY-Tuple{Grid}",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.FaceFieldY",
    "category": "method",
    "text": "FaceFieldY(grid::Grid)\n\nA Field whose values are defined on the y-face of a cell.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.FaceFieldZ-Tuple{Grid}",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.FaceFieldZ",
    "category": "method",
    "text": "FaceFieldZ(grid::Grid)\n\nA Field whose values are defined on the z-face of a cell.\n\n\n\n\n\n"
},

{
    "location": "#Fields-1",
    "page": "Oceananigans.jl Documentation",
    "title": "Fields",
    "category": "section",
    "text": "CellField\nFaceFieldX\nFaceFieldY\nFaceFieldZ\nCellField(grid::Grid)\nFaceFieldX(grid::Grid)\nFaceFieldY(grid::Grid)\nFaceFieldZ(grid::Grid)"
},

{
    "location": "#Operators-1",
    "page": "Oceananigans.jl Documentation",
    "title": "Operators",
    "category": "section",
    "text": ""
},

{
    "location": "#Oceananigans.Operators.δx!-Tuple{RegularCartesianGrid,CellField,FaceField}",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.Operators.δx!",
    "category": "method",
    "text": "δx!(g::RegularCartesianGrid, f::CellField, δxf::FaceField)\n\nCompute the difference delta_x(f) = f_E - f_W between the eastern and western cells of a cell-centered field f and store it in a face-centered field δxf, assuming both fields are defined on a regular Cartesian grid g with periodic boundary condition in the x-direction.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.Operators.δx!-Tuple{RegularCartesianGrid,FaceField,CellField}",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.Operators.δx!",
    "category": "method",
    "text": "δx!(g::RegularCartesianGrid, f::FaceField, δxf::CellField)\n\nCompute the difference delta_x(f) = f_E - f_W between the eastern and western faces of a face-centered field f and store it in a cell-centered field δxf, assuming both fields are defined on a regular Cartesian grid g with periodic boundary conditions in the x-direction.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.Operators.δy!-Tuple{RegularCartesianGrid,CellField,FaceField}",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.Operators.δy!",
    "category": "method",
    "text": "δy!(g::RegularCartesianGrid, f::CellField, δyf::FaceField)\n\nCompute the difference delta_y(f) = f_N - f_S between the northern and southern cells of a cell-centered field f and store it in a face-centered field δyf, assuming both fields are defined on a regular Cartesian grid g with periodic boundary condition in the y-direction.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.Operators.δy!-Tuple{RegularCartesianGrid,FaceField,CellField}",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.Operators.δy!",
    "category": "method",
    "text": "δy!(g::RegularCartesianGrid, f::FaceField, δyf::CellField)\n\nCompute the difference delta_y(f) = f_N - f_S between the northern and southern faces of a face-centered field f and store it in a cell-centered field δyf, assuming both fields are defined on a regular Cartesian grid g with periodic boundary condition in the y-direction.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.Operators.δz!-Tuple{RegularCartesianGrid,CellField,FaceField}",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.Operators.δz!",
    "category": "method",
    "text": "δz!(g::RegularCartesianGrid, f::CellField, δzf::FaceField)\n\nCompute the difference delta_z(f) = f_T - f_B between the top and bottom cells of a cell-centered field f and store it in a face-centered field δzf, assuming both fields are defined on a regular Cartesian grid g with Neumann boundary condition in the z-direction.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.Operators.δz!-Tuple{RegularCartesianGrid,FaceField,CellField}",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.Operators.δz!",
    "category": "method",
    "text": "δz!(g::RegularCartesianGrid, f::FaceField, δzf::CellField)\n\nCompute the difference delta_z(f) = f_T - f_B between the top and bottom faces of a face-centered field f and store it in a cell-centered field δzf, assuming both fields are defined on a regular Cartesian grid g with Neumann boundary condition in the z-direction.\n\n\n\n\n\n"
},

{
    "location": "#Difference-operators-1",
    "page": "Oceananigans.jl Documentation",
    "title": "Difference operators",
    "category": "section",
    "text": "δx!(g::RegularCartesianGrid, f::CellField, δxf::FaceField)\nδx!(g::RegularCartesianGrid, f::FaceField, δxf::CellField)\nδy!(g::RegularCartesianGrid, f::CellField, δyf::FaceField)\nδy!(g::RegularCartesianGrid, f::FaceField, δyf::CellField)\nδz!(g::RegularCartesianGrid, f::CellField, δzf::FaceField)\nδz!(g::RegularCartesianGrid, f::FaceField, δzf::CellField)"
},

{
    "location": "#Oceananigans.Operators.avgx!-Tuple{RegularCartesianGrid,CellField,FaceField}",
    "page": "Oceananigans.jl Documentation",
    "title": "Oceananigans.Operators.avgx!",
    "category": "method",
    "text": "avgx(g::RegularCartesianGrid, f::CellField, favgx::FaceField)\n\nCompute the average overlinef^x = fracf_E + f_W2 between the eastern and western cells of a cell-centered field f and store it in a g face-centered field favgx, assuming both fields are defined on a regular Cartesian grid g with periodic boundary conditions in the x-direction.\n\n\n\n\n\n"
},

{
    "location": "#Averaging-operators-1",
    "page": "Oceananigans.jl Documentation",
    "title": "Averaging operators",
    "category": "section",
    "text": "avgx!(g::RegularCartesianGrid, f::CellField, favgx::FaceField)"
},

{
    "location": "#Divergence-operators-1",
    "page": "Oceananigans.jl Documentation",
    "title": "Divergence operators",
    "category": "section",
    "text": "Building on top of the differencing operators we can define operators that compute the divergencenablacdotpmathbff = frac1V left delta_x left( A_x f_x right)\n+ delta_yleft( A_y f_y right) + delta_zleft( A_z f_z right)right<!– @docs div!(g::RegularCartesianGrid, fx::FaceFieldX, fy::FaceFieldY, fz::FaceFieldZ, δfx::CellField, δfy::CellField, δfz::CellField, div::CellField) –>"
},

]}
