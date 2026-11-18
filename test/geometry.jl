# Grid geometries. Spacings, origins and directions are dyadic (sums of powers of two) or 90-degree
# rotations, so positions are exact and compared with ==.

@testitem "Spaces" tags = [:spaces] begin
  struct Session <: AbstractSpace
    id::Int
  end
  # a space is a value: equal when it is the same coordinate system
  @test LPS() == LPS()
  @test LPS() != RAS()
  @test Session(1) == Session(1)
  @test Session(1) != Session(2)
  @test spacename(LPS()) == "LPS"
  @test spacename(Session(3)) == "Session"
  @test isbitstype(Session)
end

@testitem "GridGeometry: positions" tags = [:geometry] begin
  using StaticArrays
  g = GridGeometry((4, 3, 2); spacing=(0.5, 2.0, 1.0), origin=(1.0, -2.0, 0.25), space=LPS())
  # the origin is the center of the first voxel
  @test index_to_world(g, (1, 1, 1)) == SVector(1.0, -2.0, 0.25)
  # x = origin + direction * (spacing .* (i .- 1)), with the identity direction
  @test index_to_world(g, (3, 2, 2)) == SVector(2.0, 0.0, 1.25)
  # world_to_index inverts index_to_world, also for fractional indices
  for i in ((1, 1, 1), (4, 3, 2), (2.5, 1.5, 1.75))
    @test world_to_index(g, index_to_world(g, i)) == SVector(i)
  end
  # the affine matrix in homogeneous coordinates gives the same positions
  @test affine(g) * SVector(3, 2, 2, 1) == SVector(2.0, 0.0, 1.25, 1.0)
  # fov spans from edge to edge, the center lies midway between the outermost voxel centers
  @test fov(g) == SVector(2.0, 6.0, 2.0)
  @test fov_center(g) == (index_to_world(g, (1, 1, 1)) + index_to_world(g, (4, 3, 2))) / 2
  @test axiscoords(g, 2) == [-2.0, 0.0, 2.0]
  @test isbits(g)
end

@testitem "GridGeometry: direction" tags = [:geometry] begin
  using StaticArrays
  # index axis 1 along +y, index axis 2 along -x: a rotation by 90 degrees
  D = [0 -1; 1 0]
  g = GridGeometry((3, 2); spacing=(1.0, 2.0), origin=(0.0, 0.0), direction=D)
  @test index_to_world(g, (2, 1)) == SVector(0.0, 1.0)
  @test index_to_world(g, (1, 2)) == SVector(-2.0, 0.0)
  @test world_to_index(g, (-2.0, 1.0)) == SVector(2.0, 2.0)
  # axiscoords is only defined for axes that run along the axes of the space
  @test_throws ArgumentError axiscoords(g, 1)
  # sheared, scaled and non-unit directions are refused, as are non-positive spacings
  @test_throws ArgumentError GridGeometry((2, 2); direction=[1 1; 0 1])
  @test_throws ArgumentError GridGeometry((2, 2); direction=[2 0; 0 2])
  @test_throws ArgumentError GridGeometry((2, 2); spacing=(-1.0, 1.0))
  @test_throws ArgumentError GridGeometry((2, 2); spacing=(0.0, 1.0))
  # a direction rounded to Float32 is accepted: its columns are orthonormal up to a few eps(Float32)
  c, s = Float32(cos(0.3)), Float32(sin(0.3))
  @test GridGeometry((2, 2); direction=[c -s; s c]) isa GridGeometry
end

@testitem "GridGeometry: subgrid" tags = [:geometry] begin
  using StaticArrays
  g = GridGeometry((8, 6, 4); spacing=(0.5, 1.0, 2.0), origin=(1.0, 2.0, 3.0), space=RAS())
  # every selected voxel keeps its position: voxel k of the subgrid is voxel first + step (k - 1)
  for inds in ((2:3:8, 1:6, 4:-1:1), (8:-2:2, 3:3, 1:4))
    s = subgrid(g, inds)
    for k in CartesianIndices(size(s))
      full = ntuple(d -> inds[d][k[d]], 3)
      @test index_to_world(s, Tuple(k)) == index_to_world(g, full)
    end
    @test space(s) == RAS()
  end
  # a negative step flips the direction column instead of the spacing
  s = subgrid(g, (8:-2:2, 1:6, 1:4))
  @test spacing(s) == SVector(1.0, 1.0, 2.0)
  @test direction(s)[:, 1] == SVector(-1.0, 0.0, 0.0)
  # dropping a dimension gives a slice that stays in the 3D space
  slice = subgrid(g, (1:8, 3, 1:4))
  @test size(slice) == (8, 4)
  @test TomoImages.worlddims(slice) == 3
  for k in CartesianIndices(size(slice))
    @test index_to_world(slice, Tuple(k)) == index_to_world(g, (k[1], 3, k[2]))
  end
  @test_throws BoundsError subgrid(g, (0:2, 1:6, 1:4))
end

@testitem "GridGeometry: samegrid" tags = [:geometry] begin
  g = GridGeometry((4, 4); spacing=(1.0, 1.0), origin=(0.0, 0.0), space=LPS())
  @test samegrid(g, GridGeometry((4, 4); origin=(0.0, 0.0), space=LPS()))
  # within a thousandth of a voxel by default
  @test samegrid(g, GridGeometry((4, 4); origin=(0.0005, 0.0), space=LPS()))
  @test !samegrid(g, GridGeometry((4, 4); origin=(0.01, 0.0), space=LPS()))
  # a different space is a different grid, even with the same numbers
  @test !samegrid(g, GridGeometry((4, 4); origin=(0.0, 0.0), space=RAS()))
  @test !samegrid(g, GridGeometry((4, 5); space=LPS()))
  # the same voxels described from the other corner: flipped direction, moved origin
  flipped = GridGeometry((4, 4); origin=(3.0, 0.0), direction=[-1 0; 0 1], space=LPS())
  @test !samegrid(g, flipped)
  @test samegrid(subgrid(g, (4:-1:1, 1:4)), flipped)
end
