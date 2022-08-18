library(brkga)

# Leitura das instâncias
tspdl_data <-read_rds("inst/extdata/TSPDL_DATA.rds")
# Imprimer informações 
verbose = 2
# Número de réplicas por instância
Replicas = 1
# tibble para armazenar resultados
resultado <- tibble()

Instancia_initial = 1
Intanacia_final =  5

for (k in Instancia_initial:Intanacia_final){#
  distance <- tspdl_data[k,]$distance %>% unlist
  demand <- tspdl_data[k,]$demand %>% unlist
  draft <- tspdl_data[k,]$draft %>% unlist
  N <- tspdl_data[k,]$N
  instance_name <- tspdl_data$name[k]
  instance_n <- tspdl_data$n[k]
  instance_m <- tspdl_data$m[k]
  for (i in 1:Replicas){
    cat(sprintf("Instance %s, n = %s, m = %s. Replication %d of %d.             \n", instance_name, instance_n, instance_m, i, Replicas))
    seed <- 0 # as.integer(Sys.time())
    #
    rst <- brkga::calado_brkga(DistanceVector = distance,
                               DemandVector = demand,
                               DraftVector = draft,
                               N = N,
                               LS_INTVL = 1,
                               GEN_INTVL = 1,
                               MAX_TIME = 10, 
                               p = 250, 
                               pe = .2, 
                               pm = .2,
                               rhoe = .75,
                               THREADS = 8, 
                               K = 8, 
                               MAX_GENS = 1000,
                               RESET_AFTER = 20, 
                               verbose = verbose,
                               rngSeed = seed)
    
    resultado <- bind_rows(
      resultado,
      tibble(
        Replica = i, 
        Seed = seed,
        Tour = list(rst$Tour),
        MyTour = list(rst$MyTour),
        BRKGA = rst$BKFitness,
        N_gen = rst$`Generations Number`,
        N_bst = rst$`Best Generation`,
        Instancia = instance_name,
        n = instance_n,
        m = instance_m
      )
    )
  }  
  rm(i, rst, seed, distance, demand, draft, instance_name, instance_m, instance_n)
  write_rds(resultado, path = "inst/extdata/cadlado_brkga_resultado.rds")
}

cat("\nThis is the end. The Doors.")

# resultado %>% 
#   group_by(Instancia) %>% 
#   summarise(N = n(), 
#             MeanLSErr = mean(LSEr), 
#             MinLSErr = min(LSEr), 
#             MaxLSErr = max(LSEr))
# 
# 
# resultado %>% 
#   gather(Error, Value, -Order:-Target) %>% 
#   mutate(Error = if_else(Error == "BKEr", "BRKGA", "Local Search"))  %>% 
#   ggplot() + geom_boxplot(aes(x = Error, y = Value, fill = Error)) + facet_wrap(. ~ Instancia, ncol = 1) + coord_flip()
