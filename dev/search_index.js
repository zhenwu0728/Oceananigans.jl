var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Oceananigans.jl-Documentation-1",
    "page": "Home",
    "title": "Oceananigans.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "#Oceananigans.RegularCartesianGrid",
    "page": "Home",
    "title": "Oceananigans.RegularCartesianGrid",
    "category": "type",
    "text": "RegularCartesianGrid\n\nA Cartesian grid with regularly spaces cells and faces so that Δx, Δy, and Δz are constants. Fields are stored using floating-point values of type T.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.RegularCartesianGrid-Tuple{ModelMetadata,Any,Any}",
    "page": "Home",
    "title": "Oceananigans.RegularCartesianGrid",
    "category": "method",
    "text": "RegularCartesianGrid(metadata::ModelMetadata, N, L)\n\nCreate a regular Cartesian grid with size N = (N_x N_y N_z) and domain size L = (L_x L_y L_z) where fields are stored using floating-point values of type T.\n\nExamples\n\njulia> g = RegularCartesianGrid((16, 16, 8), (2π, 2π, 2π))\n\n\n\n\n\n"
},

{
    "location": "#Grids-1",
    "page": "Home",
    "title": "Grids",
    "category": "section",
    "text": "RegularCartesianGrid\nRegularCartesianGrid(metadata::ModelMetadata, N, L)"
},

{
    "location": "#Oceananigans.CellField",
    "page": "Home",
    "title": "Oceananigans.CellField",
    "category": "type",
    "text": "CellField{T,G<:Grid{T}} <: Field{G}\n\nA cell-centered field defined on a grid G whose values are stored as floating-point values of type T.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.FaceFieldX",
    "page": "Home",
    "title": "Oceananigans.FaceFieldX",
    "category": "type",
    "text": "FaceFieldX{T,G<:Grid{T}} <: FaceField{G}\n\nAn x-face-centered field defined on a grid G whose values are stored as floating-point values of type T.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.FaceFieldY",
    "page": "Home",
    "title": "Oceananigans.FaceFieldY",
    "category": "type",
    "text": "FaceFieldY{T,G<:Grid{T}} <: FaceField{G}\n\nA y-face-centered field defined on a grid G whose values are stored as floating-point values of type T.\n\n\n\n\n\n"
},

{
    "location": "#Oceananigans.FaceFieldZ",
    "page": "Home",
    "title": "Oceananigans.FaceFieldZ",
    "category": "type",
    "text": "FaceFieldZ{T,G<:Grid{T}} <: FaceField{G}\n\nA z-face-centered field defined on a grid G whose values are stored as floating-point values of type T.\n\n\n\n\n\n"
},

{
    "location": "#Fields-1",
    "page": "Home",
    "title": "Fields",
    "category": "section",
    "text": "CellField\nFaceFieldX\nFaceFieldY\nFaceFieldZ"
},

{
    "location": "#Operators-1",
    "page": "Home",
    "title": "Operators",
    "category": "section",
    "text": ""
},

{
    "location": "#Difference-operators-1",
    "page": "Home",
    "title": "Difference operators",
    "category": "section",
    "text": "δx!(g::RegularCartesianGrid, f::CellField, δxf::FaceField)\nδx!(g::RegularCartesianGrid, f::FaceField, δxf::CellField)\nδy!(g::RegularCartesianGrid, f::CellField, δyf::FaceField)\nδy!(g::RegularCartesianGrid, f::FaceField, δyf::CellField)\nδz!(g::RegularCartesianGrid, f::CellField, δzf::FaceField)\nδz!(g::RegularCartesianGrid, f::FaceField, δzf::CellField)"
},

{
    "location": "#Averaging-operators-1",
    "page": "Home",
    "title": "Averaging operators",
    "category": "section",
    "text": "avgx!(g::RegularCartesianGrid, f::CellField, favgx::FaceField)"
},

{
    "location": "#Divergence-operators-1",
    "page": "Home",
    "title": "Divergence operators",
    "category": "section",
    "text": "Building on top of the differencing operators we can define operators that compute the divergencenablacdotpmathbff = frac1V left delta_x left( A_x f_x right)\n+ delta_yleft( A_y f_y right) + delta_zleft( A_z f_z right)right<!– @docs div!(g::RegularCartesianGrid, fx::FaceFieldX, fy::FaceFieldY, fz::FaceFieldZ, δfx::CellField, δfy::CellField, δfz::CellField, div::CellField) –>"
},

]}
