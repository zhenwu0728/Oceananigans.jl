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
    "location": "#Grids-1",
    "page": "Home",
    "title": "Grids",
    "category": "section",
    "text": "RegularCartesianGrid\nRegularCartesianGrid(N, L; dim=3, FloatType=Float64)"
},

{
    "location": "#Fields-1",
    "page": "Home",
    "title": "Fields",
    "category": "section",
    "text": "CellField\nFaceFieldX\nFaceFieldY\nFaceFieldZ\nCellField(grid::Grid)\nFaceFieldX(grid::Grid)\nFaceFieldY(grid::Grid)\nFaceFieldZ(grid::Grid)"
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
