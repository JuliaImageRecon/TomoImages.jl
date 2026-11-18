using Documenter, Literate, TomoImages

# Examples are Literate scripts in literate/<section>/, rendered to src/generated/<section>/
# and run when the documentation is built.
const LITERATE = joinpath(@__DIR__, "literate")
const GENERATED = joinpath(@__DIR__, "src", "generated")
for dir in readdir(LITERATE), file in filter(endswith(".jl"), readdir(joinpath(LITERATE, dir)))
  Literate.markdown(joinpath(LITERATE, dir, file), joinpath(GENERATED, dir); documenter=true)
end

DocMeta.setdocmeta!(TomoImages, :DocTestSetup, :(using TomoImages); recursive=true)

makedocs(;
  modules=[TomoImages],
  sitename="TomoImages.jl",
  pagesonly=true,
  checkdocs=:exports,
  authors="Niklas Hackelberg and contributors",
  remotes=nothing,
  format=Documenter.HTML(;
    prettyurls=get(ENV, "CI", "false") == "true",
    edit_link="main",
    collapselevel=1,
  ),
  # Diátaxis: tutorials to learn, how-to guides for tasks, explanation for understanding,
  # reference for looking things up
  pages=[
    "Home" => "index.md",
    "Getting Started" => "generated/examples/getting_started.md",
    "Tutorials" => [
      "Overlay Images from Two Devices" => "generated/tutorials/overlay.md",
    ],
    "How to" => [
      "Define a Space" => "generated/howto/define_space.md",
      "Write a Resampler" => "generated/howto/write_resampler.md",
    ],
    "Explanation" => [
      "Grids, Spaces and Transforms" => "generated/explanation/conventions.md",
    ],
    "Reference" => [
      "Spaces and Grids" => "api/grids.md",
      "Transforms" => "api/transforms.md",
      "Images" => "api/images.md",
      "Resampling" => "api/resampling.md",
    ],
  ],
)
