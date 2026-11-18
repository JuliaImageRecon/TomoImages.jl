@testitem "AffineTransform" tags = [:transforms] begin
  using StaticArrays
  struct Device <: AbstractSpace end
  # 90 degrees about z and a shift: exact in floating point
  t = AffineTransform([0 -1 0; 1 0 0; 0 0 1], [0, 0, 120]; source=LPS(), target=Device())
  @test t((1.0, 0.0, 0.0)) == SVector(0.0, 1.0, 120.0)
  @test source(t) == LPS() && target(t) == Device()
  @test isbits(t)
  # the inverse swaps the spaces and undoes the map
  ti = inv(t)
  @test source(ti) == Device() && target(ti) == LPS()
  @test ti(t((1.0, 2.0, 3.0))) == SVector(1.0, 2.0, 3.0)
  # compose applies the first transform, then the second
  c = compose(ti, t)
  @test source(c) == LPS() && target(c) == LPS()
  @test c((1.0, 2.0, 3.0)) == SVector(1.0, 2.0, 3.0)
  # composing transforms whose spaces do not chain is an error
  @test_throws ArgumentError compose(t, t)
end

@testitem "LPS and RAS" tags = [:transforms] begin
  using StaticArrays
  t = transform_between(LPS(), RAS())
  @test t((1.0, 2.0, 3.0)) == SVector(-1.0, -2.0, 3.0)
  @test transform_between(RAS(), LPS())(t((1.0, 2.0, 3.0))) == SVector(1.0, 2.0, 3.0)
  # a geometry written in RAS describes the same voxels: every voxel center maps to its RAS position
  g = GridGeometry((4, 3, 2); spacing=(0.5, 1.0, 2.0), origin=(1.0, 2.0, 3.0), space=LPS())
  h = transform_geometry(g, t)
  @test space(h) == RAS()
  for i in CartesianIndices(size(g))
    @test index_to_world(h, Tuple(i)) == t(index_to_world(g, Tuple(i)))
  end
  @test_throws ArgumentError transform_geometry(h, t)
end
