using NKLandscapes

const num_bins = 20
cond_mutsize_array = zeros(Int64,num_bins,num_bins)
abs_mutsize_array = zeros(Int64, num_bins)

function fitdiff(ls::Landscape)
  rg = rand(Genotype,ls)
  frg = fitness(rg,ls)
  frg,fitness(random_neighbor(rg,ls),ls)-frg
end

function bin(g) 
  convert(Int64,floor(div(num_bins,2)*g+div(num_bins,2))) + 1
end

function abin(g) 
  convert(Int64,floor(2*num_bins*g)) + 1
end


function abs_mutsize_increment!(mutsize_array::Array{Int64,1},ls::Landscape,reps::Int64)
  for i = 1:reps
    g, mutg = fitdiff(ls)
    mutsize_array[abin(abs(mutg))] += 1
  end
end

function abs_mutsize_tbl(ls_list::Array{NKLandscapes.NKLandscape,1},reps::Int64)
  abs_mutsize_array = zeros(Int64, num_bins+4,length(ls_list))
  j = 1
  for ls in ls_list
    for i = 1:reps
      g, mutg = fitdiff(ls)
      abs_mutsize_array[abin(abs(mutg)),j] += 1
    end
    j += 1
  end
  return abs_mutsize_array
end

function conditional_mutsize_increment!(mutsize_array::Array{Int64,2},ls::Landscape,reps::Int64)
  for i = 1:reps
    g, mutg = fitdiff(ls)
    mutsize_array[convert(Int64,floor(num_bins*g)),bin(mutg)] += 1
  end
end


