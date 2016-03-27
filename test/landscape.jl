srand(1)

# TODO: Random links gets set up properly on instantiation.

context("landscape.jl") do
  context("NKLandscape") do
    l = NK.NKLandscape(4, 1)
    @fact length(l.links) --> l.n
    @fact length(l.contribs) --> l.n

    context("when near=true") do
      l = NK.NKLandscape(4, 2; near=true)
      @fact l.links[1] --> 0b1011
      @fact l.links[2] --> 0b0111
      @fact l.links[3] --> 0b1110
      @fact l.links[4] --> 0b1101
    end
  end

  context("NKqLandscape") do
    l = NK.NKqLandscape(4, 1, 2)
    @fact length(l.links) --> l.n
    @fact length(l.contribs) --> l.n

    context("when near=true") do
      l = NK.NKqLandscape(4, 2, 2; near=true)
      @fact l.links[1] --> 0b1011
      @fact l.links[2] --> 0b0111
      @fact l.links[3] --> 0b1110
      @fact l.links[4] --> 0b1101
    end
  end

  context("NKpLandscape") do
    l = NK.NKpLandscape(4, 1, 0.90)
    @fact length(l.links) --> l.n
    @fact length(l.contribs) --> l.n

    context("when near=true") do
      l = NK.NKpLandscape(4, 2, 0.90; near=true)
      @fact l.links[1] --> 0b1011
      @fact l.links[2] --> 0b0111
      @fact l.links[3] --> 0b1110
      @fact l.links[4] --> 0b1101
    end
  end
end
