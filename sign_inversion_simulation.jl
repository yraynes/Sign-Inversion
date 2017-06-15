using Distributions

type Lineage
  fitness::Float64
  size::Int64
  state::Vector{Float64}
  mut_rate::Float64
end

function average_fitness(population::Dict{Array{Float64,1}, Lineage})
#calculate average fitness of the population. also outputs population size
  popN = Int64
  popN = 0
  popW = Float64
  popW = 0
  for line in values(population)
    popN+=line.size
    popW+=line.size * line.fitness
  end
  popW = popW/popN
  return popW, popN
end

function pop_size(population::Dict{Array{Float64,1}, Lineage})
  popN = Int64
  popN = 0
  for line in values(population)
    popN+=line.size
  end
  return popN
end

function assay_mutation_rate(population::Dict{Array{Float64,1}, Lineage})
#calculates mutator frequency
  popN = Int64
  popN = 0
  mutN = Int64
  mutN = 0
  for line in values(population)
    popN+=line.size
    if line.mut_rate>1 mutN+=line.size end
  end
  mutF = Float64
  mutF = mutN/popN
  return mutF
end

function fitness_calc(state::Array{Float64,1})
#calculates fitness from the lineage state = sum of all mutational effects - the mutator strength (100) - oriT marker (-100)
  w = Float64
  w = 1 + sum(state) - state[1]
  return w
end

function avail_sites(state::Array{Float64,1})
  sites = Int64[]
  for i = 1:length(state)
    if state[i] == 0 push!(sites,i) end
  end
  return sites
end

function mutate_population(population::Dict{Array{Float64,1}, Lineage}, Ub::Float64, sb::Float64, Ud::Float64, sd::Float64, mut_multiplier::Float64)
  new_population = Dict{Array{Float64,1}, Lineage}()
  #initialize new population dictionary
  for (key, line) in population
    mutations = rand(Poisson(line.mut_rate * (Ud + Ub) * line.size))
    if mutations > line.size mutations = line.size end
    bmutations = rand(Binomial(mutations, Ub/(Ub+Ud)))
    dmutations = mutations - bmutations

    open_sites =   avail_sites(line.state)
    #generates a list of indexes of all loci that are available for mutation
    mutation_positions = open_sites[rand(1:end,(dmutations+bmutations))] #generates positions of new mutations by randomly picking indexes from open_sites

      #for each new mutation, copies the state, adds mutation and generates the key.
      #check if it is alrady in the new_population. if it is, add 1 individual to it, if it's not - make another lineage with
      #the key and size 1
    c=1
    for i = 1:bmutations
      #assign beneficial mutations
      newstate = copy(line.state)#
      newstate[mutation_positions[c]] = sb
      c+=1
      if newstate in keys(new_population) new_population[newstate].size+=1
      else
        newfit::Float64 = 1 + sum(newstate[2:end])
        new_population[newstate] = Lineage(newfit, 1, newstate, newstate[1])
      end
    end

    for i = 1:dmutations
      #assign deleterious mutations
      newstate = copy(line.state)#
      newstate[mutation_positions[c]] = -sd
      c+=1
      if newstate in keys(new_population) new_population[newstate].size+=1
      else
        newfit::Float64 = 1 + sum(newstate[2:end])
        new_population[newstate] = Lineage(newfit, 1, newstate, newstate[1])
      end
    end

    new_size = line.size - bmutations - dmutations
    if  new_size > 0
      if key in keys(new_population) new_population[key].size+=new_size
      else new_population[key] = Lineage(line.fitness, new_size, line.state, line.mut_rate)
      end
    end
  end
  return new_population
end

#import Distributions.isprobvec
#isprobvec(p::Vector{Float64}) = true

function wright_fisher_reproduction(population::Dict{Array{Float64,1}, Lineage}, N0::Int64, Nnew::Int64, popw::Float64)
  proby_list = Float64[] #array of probabilities
  lineage_list = Lineage[] #initialize lineage liste
  for (key,line) in population
    push!(proby_list, max(line.size/N0 * line.fitness/popw, 0.0))
    #representation of a lineage in the next generations depends on its frequency and relative fitness
    push!(lineage_list, line)
  end
  new_counts_list = rand(Multinomial(Nnew, proby_list))
  new_population = Dict{Array{Float64,1}, Lineage}()
  for i = 1:length(new_counts_list)
    if new_counts_list[i] > 0
      new_population[lineage_list[i].state] = Lineage(lineage_list[i].fitness, new_counts_list[i], lineage_list[i].state, lineage_list[i].mut_rate)
    end
  end
  return new_population
end

function bottleneck(population::Dict{Array{Float64,1}, Lineage}, N0::Int64, Nnew::Int64)
  #subjects populations to random bottlenecks
  proby_list = Float64[] #array of probabilities
  lineage_list = Lineage[]
  for (key,line) in population
    push!(proby_list, line.size/N0)
    push!(lineage_list, line)
  end
  new_counts_list = rand(Multinomial(Nnew, proby_list))
  new_population = Dict{Array{Float64,1}, Lineage}()
  for i = 1:length(new_counts_list)
    if new_counts_list[i] > 0
      new_population[lineage_list[i].state] = Lineage(lineage_list[i].fitness, new_counts_list[i], lineage_list[i].state, lineage_list[i].mut_rate)
    end
  end
  return new_population
end


function simulate(pop_Ni)
  const init_mut_N = 1
  const sb = 0.1
  const sd = 0.1
  const mutator_strength = 100.0
  const Ub = 0.000001
  const Ud = 0.0001

  population = Dict{Array{Float64,1}, Lineage}()
  wt_state = zeros(Float64,100)
  wt_state[1] = 1.0
  population[wt_state] = Lineage(1,pop_Ni-init_mut_N,wt_state,1.0)
  #key = all non zero values
  mut_state = copy(wt_state)
  mut_state[1] = mutator_strength
  population[mut_state] = Lineage(1,init_mut_N,mut_state,mutator_strength)
  mut_f = assay_mutation_rate(population)
  popw, popN = average_fitness(population)
  generations = Int64
  generations = 0
  #initialize generations counter

  while 0.0<mut_f<1.0
    popw, pop_N = average_fitness(population)
    population = wright_fisher_reproduction(population, pop_N, pop_Ni, popw)
    population = mutate_population(population, Ub, sb, Ud, sd, mutator_strength)
    mut_f = assay_mutation_rate(population)
    popw, pop_N = average_fitness(population)
    generations+=1
  end
  return mut_f ==1.0, generations, popw
end

for n in [6, 8, 10, 12, 16, 20, 26, 32, 42, 52, 66, 84, 108, 136, 174, 220, 280, 356, 452, 574, 728, 924, 1172, 1486, 1886, 2394, 3040, 3856, 4894, 6210, 7880, 10000]
  time_to_mut= Float64[]
  #time_to_mut list records all times of successful mutator hitchhiking
  for run = 1:1000
    output = simulate(n)
    if output[1]
      push!(time_to_mut, output[2])
    else
      push!(time_to_mut, 0)
    end
    end
  outfile = open(string("time_to_mut.csv"), "a")
  write(outfile, join(time_to_mut, ","), "\n")
  close(outfile)
end
end
