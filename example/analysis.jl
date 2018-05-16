using DataFrames
using Gadfly

if length(ARGS) == 0
  println("CSV filename required.")
end

csvfilename = ARGS[1]
svgbase = join(split(csvfilename, ".")[1:end-1], ".")

df = readtable(csvfilename, makefactors=true)

df = by(df, [:generation, :simulationType]) do d
  DataFrame(
    meanFitness=mean(d[:meanFitness]),
    medianFitness=median(d[:medianFitness]),
    maxFitness=maximum(d[:maxFitness])
  )
end

# Plot development

function drawPlot(variable)
  p = plot(df, x="generation", y=variable, color="simulationType", Geom.line)
  draw(SVG("$(svgbase)_$(variable).svg", 12inch, 12inch), p)
end

drawPlot("meanFitness")
drawPlot("medianFitness")
drawPlot("maxFitness")

