"""
    TomoImages

Images that know where they are in space: grid geometries with an origin, a spacing, a direction and
a named coordinate space; images on such grids with named non-spatial axes; transforms between
spaces; and resampling from one grid onto another, on CPUs and GPUs.
"""
module TomoImages

using LinearAlgebra, StaticArrays
using Adapt
using KernelAbstractions
import Interpolations

const KA = KernelAbstractions

export AbstractSpace, DefaultSpace, LPS, RAS
export AbstractGeometry, GridGeometry, spacing, origin, direction, space, affine, index_to_world,
       world_to_index, fov, fov_center, axiscoords
export AbstractSpatialTransform, AffineTransform, source, target, transform_between, transform_geometry
export GridAxis, NumericAxis, CategoricalAxis
export AbstractSpatialImage, SpatialImage, geometry, meta, selectvalue
export AbstractResampler, NearestResampler, LinearResampler, BSplineResampler, resample

public spacename, worlddims, subgrid, compose, checkspaces
public dims, spatialdims, unit, indexof
public resample!, default_resampler, supports, resample_eltype, IndexMap, index_map

include("spaces/spaces.jl")
include("geometry/geometry.jl")
include("transforms/transforms.jl")
include("images/images.jl")
include("resampling/resampling.jl")

end
