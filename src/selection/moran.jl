import Base.Random: rand

export moransel, moransel!

function moransel(p::Population, ls::Landscape, iters::Int64)
  np = copy(p)
  moransel!(np, ls)
  return np
end

function moransel!(p::Population, ls::Landscape, iters::Int64)
  fs = popfits(p, ls)
  fmax = maximum(fs)
  n = popsize(p)

  for _ = 1:iters
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
end
