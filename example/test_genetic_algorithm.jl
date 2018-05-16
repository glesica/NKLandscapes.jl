# Test various options of function genetic algorithm.

genetic_algorithm(ls,6,10,statistics_funct=simple_statistics_funct,selection_funct=tournsel!,tournsel_k=3)
genetic_algorithm(ls,6,10,statistics_funct=simple_statistics_funct,selection_funct=moransel!,moran_sel_iters=1.0)
genetic_algorithm(ls,6,10,statistics_funct=simple_statistics_funct,mutation_funct=bsmutate!)
genetic_algorithm(ls,6,10,statistics_funct=simple_statistics_funct,mut_prob=0.0)
genetic_algorithm(ls,6,100,statistics_funct=simple_statistics_funct,mut_prob=0.0,termination_funct=test_converged)
genetic_algorithm(ls,6,10,statistics_funct=simple_statistics_funct,init_genotype=rand(Genotype,ls))


