import Base.Random: rand

export bwmutate, bwmutate!

function bwmutate(p::Population, ls::Landscape, mutprob::Float64)
  np = copy(p)
  bwmutate!(np, ls, mutprob)
  return np
end

function bwmutate!(p::Population, ls::Landscape, mutprob::Float64)
  for i = 1:popsize(p)
    for j = 1:ls.n
      # i = population member
      # j = locus to be mutated
      if rand() < mutprob
        # TODO: Should we exclude the current value?
        p[j,i] = rand(1:ls.a)
      end
    end
  end
end
