library(tidyverse)
library(brkga)

lista <- read_rds(system.file("extdata", package = "brkga",  "mdplib.rds")) %>% 
  mutate(Name = paste0("conv_", Instance), 
         Type = str_sub(Name, 6,8), 
         Number = as.integer(str_remove(str_sub(Name, 12, 13), "_")),
         SubType = str_sub(Name, 10,10), 
         n = as.integer(str_remove(str_extract(Name,pattern = "n\\d+"), pattern = "n")),
         m = as.integer(str_remove(str_extract(Name,pattern = "m\\d+"), pattern = "m"))) %>% 
  filter(Type == "MDG")


N <- 10 # Número de rodadas


resultado <- tibble()

for (k in 1:nrow(lista)){#
  dist_matrix <- read_rds(system.file("extdata", package = "brkga", paste0(
                                      lista$Type[k], ".", lista$Number[k], ".", 
                                      lista$SubType[k], ".n", 
                                      lista$n[k], "m", 
                                      lista$m[k],".rds")                                    
                                      ))
    m <- lista$m[k]
    bestValue <- lista$`Best Value`[k]
    InstName <- lista$Name[k]
    for (i in 1:N){
      cat(sprintf("Instance %d of %d. Replication %d of %d.             \r", k, nrow(lista), i, N))
      seed <- as.integer(Sys.time())
      rst <- brkga::mdp_brkga(DistanceMatrix = dist_matrix,
                              m = m, 
                              LS_INTVL = 1,
                              GEN_INTVL = 1,
                              MAX_TIME = 7200, 
                              p = 150, 
                              pe = .2, 
                              pm = .2,
                              rhoe = .75,
                              THREADS = 4, 
                              K=4, 
                              MAX_GENS = 500,
                              RESET_AFTER = 20, 
                              verbose = FALSE,
                              rngSeed = seed)
      
      resultado <- bind_rows(
        resultado,
        tibble(
          Order = i, 
          Seed = seed,
          Tour = list(rst$Tour),
          BRKGA = rst$BKFitness,
          LS = rst$LSFitness,
          N_gen = rst$`Generations Number`,
          N_imp = rst$`Improvement Number`,
          N_bst = rst$`Best Generation`,
          Instancia = InstName,
          Duration = rst$Duration,
          Target = bestValue,
          LSEr = round(100*(LS - bestValue)/bestValue, digits = 2),
          BKEr = round(100*(BRKGA - bestValue)/bestValue, digits = 2)
        )
      )
      # try(if (i > 1) clear <- dev.off(dev.list()["RStudioGD"]), silent = TRUE)
      # print(
      #   resultado %>% filter(Instancia == InstName) %>% 
      #     gather(Error, Value, -Order:-Target) %>% 
      #     mutate(Error = if_else(Error == "BKEr", "BRKGA", "BRKGA+LS"))  %>% 
      #     ggplot() + 
      #     ggtitle(paste("Instance: ", InstName)) +
      #     geom_point(aes(x = as.factor(Order), y = Value, color = Error)) + 
      #     labs(x = "Execution", y = "Error", color = "Method") +
      #     #scale_x_discrete() +
      #     #geom_hline(yintercept = median(resultado$LSEr), color = "#00BFC4", alpha  = .25, size = 1) +
      #     #geom_hline(yintercept = median(resultado$BKEr), color = "#F8766D", alpha = .25, size = 1) + 
      #     ylim(min(min(resultado$BKEr),min(resultado$LSEr)) - .1, max(max(resultado$BKEr),max(resultado$LSEr)) + .1)
      # )
    }
    rm(i, rst, seed, bestValue, dist_matrix, InstName, m)
    write_rds(resultado, path = "~/data/brkga.rds")
  }


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
