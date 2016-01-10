import Base.Random: rand

export moransel, moransel!

function moransel(p::Population, ls::Landscape)
end

function moransel!(p::Population, ls::Landscape)
  fs = popfits(p, ls)
  fmax = maximum(fs)

  n = popsize(p)
  indkill = rand(1:n)
  indbirth = 0
  while indbirth == 0
    i = rand(1:n)
    w = fs[i] / fmax
    if rand() < w
      indbirth = i
    end
  end

  p[:,indkill] = p[:,indbirth]
end
