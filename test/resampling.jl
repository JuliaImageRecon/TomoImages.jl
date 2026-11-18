@testitem "Resampling: same grid and whole-voxel shifts" tags = [:resampling] begin
  g = GridGeometry((5, 4, 3); spacing=(1.0, 2.0, 0.5), origin=(0.0, 0.0, 0.0), space=LPS())
  data = reshape(Float32.(1:(2*5*4*3*3)), 2, 5, 4, 3, 3)
  img = SpatialImage(data, g; spatialdims=(2, 3, 4), extra=(a=CategoricalAxis(1:2), b=NumericAxis(1:3)))
  # on its own grid every output voxel is an input voxel center: the data, exactly
  for r in (NearestResampler(), LinearResampler(), BSplineResampler(0), BSplineResampler(1))
    out = resample(r, img, g)
    @test parent(out) == data
    @test geometry(out) === g && extra(out) === extra(img)
  end
  # a grid shifted by two voxels along x reads the input two voxels further; beyond the input: fill
  shifted = GridGeometry((5, 4, 3); spacing=(1.0, 2.0, 0.5), origin=(2.0, 0.0, 0.0), space=LPS())
  for r in (NearestResampler(), LinearResampler())
    out = parent(resample(r, img, shifted; fill=-1))
    @test out[:, 1:3, :, :, :] == data[:, 3:5, :, :, :]
    # the input extends half a voxel beyond its last center, so x = 4.5 still reads voxel 5
    @test all(==(-1), out[:, 4:5, :, :, :])
  end
end

@testitem "Resampling: linear images" tags = [:resampling] begin
  # f(i, j) = 3i - 2j + 7 at voxel (i, j); halfway between voxel centers the weights are 1/2 and
  # every sum is exact, so linear interpolation reproduces f exactly
  g = GridGeometry((6, 5); spacing=(1.0, 1.0), origin=(0.0, 0.0))
  f(i, j) = 3i - 2j + 7.0
  img = SpatialImage([f(i, j) for i in 1:6, j in 1:5], g)
  half = GridGeometry((10, 8); spacing=(0.5, 0.5), origin=(0.5, 0.5))
  out = parent(resample(LinearResampler(), img, half))
  for j in CartesianIndices(size(half))
    x = TomoImages.index_to_world(half, Tuple(j))
    @test out[j] == f(x[1] + 1, x[2] + 1)
  end
  # B-splines of degree 2 and 3 with linear extrapolation at the boundary also reproduce linear
  # images; their prefilter solves a diagonally dominant system with condition number below 3, so
  # the error stays within a few eps of the largest value (bound: 64 eps * maximum)
  for degree in (2, 3)
    s = parent(resample(BSplineResampler(degree), img, half))
    @test maximum(abs, s .- out) <= 64 * eps(Float64) * maximum(abs, out)
  end
end

@testitem "Resampling: transforms" tags = [:resampling] begin
  struct Device <: AbstractSpace end
  # the input on a device grid; the output on an LPS grid rotated by 90 degrees about z
  gi = GridGeometry((4, 3, 2); spacing=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0), space=Device())
  data = reshape(Float64.(1:24), 4, 3, 2)
  img = SpatialImage(data, gi)
  # t maps LPS positions (the output's space) into device positions (the input's space)
  t = AffineTransform([0 1 0; -1 0 0; 0 0 1], [0, 0, 0]; source=LPS(), target=Device())
  go = GridGeometry((3, 4, 2); spacing=(1.0, 1.0, 1.0), origin=(-2.0, 0.0, 0.0), space=LPS())
  out = parent(resample(img, go, t))
  # output voxel (i, j, k) sits at LPS (i - 3, j - 1, k - 1), which t maps to device (j - 1, 3 - i, k - 1)
  @test out == [data[j, 4-i, k] for i in 1:3, j in 1:4, k in 1:2]
  # spaces are checked: no transform across spaces, and the transform in the right direction
  @test_throws ArgumentError resample(img, go)
  @test_throws ArgumentError resample(img, go, inv(t))
end

@testitem "Resampling: label maps" tags = [:resampling] begin
  g = GridGeometry((4, 4); spacing=(1.0, 1.0), origin=(0.0, 0.0))
  labels = SpatialImage(UInt8[i > 2 ? 2 : 1 for i in 1:4, j in 1:4], g)
  # integer images default to nearest neighbor and keep their element type: no averaged labels
  @test default_resampler(labels) == NearestResampler()
  fine = GridGeometry((8, 8); spacing=(0.5, 0.5), origin=(-0.25, -0.25))
  out = resample(labels, fine)
  @test eltype(out) == UInt8
  @test Set(parent(out)) == Set(UInt8[1, 2])
  @test default_resampler(SpatialImage(rand(4, 4), g)) == LinearResampler()
  # B-splines need the data on the CPU in an Array
  @test supports(BSplineResampler(3), labels)
  @test !supports(BSplineResampler(3), SpatialImage(view(rand(4, 4), :, :), g))
end
