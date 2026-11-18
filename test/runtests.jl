# Test items (TestItemRunner): run all with `Pkg.test()`, a selection by tags with e.g.
# `Pkg.test(test_args=["geometry"])`, single items from VS Code.
#
# Tags: :spaces, :geometry, :transforms, :images, :resampling, :aqua

using TestItemRunner

const TAGS = Symbol.(ARGS)

@run_package_tests filter = ti -> isempty(TAGS) || any(in(ti.tags), TAGS)
