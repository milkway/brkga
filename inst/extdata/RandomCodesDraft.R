library(tidyverse)

tabela <- tibble(path = list.files(path = 'inst/extdata/tspdl/', full.names = TRUE)) %>% 
  mutate(filename = str_sub(path, str_locate(path, '//')[[1]] + 2, end = -1),
         raw = sapply(path, read_file),
         name = str_trim(str_extract(raw, pattern = '(?<=\\!name:)(.)*(?=\\r\\n!type)')),
         n = as.integer(str_trim(str_extract(filename, pattern = '(?<=_)(.)*(?=_)'))),
         m = as.integer(str_trim(str_extract(filename, pattern = '(?<=_\\d{1,2}_)(.)*(?=\\.dat)'))),
         type = str_trim(str_extract(raw, pattern = '(?<=\\!type:)(.)*(?=\\r\\n!comment)')),
         comment = str_trim(str_extract(raw, pattern = '(?<=\\!comment:)(.)*(?=\\r\\nN)')),
         N = as.integer(str_trim(str_extract(raw, pattern = '(?<=N:)(.)*(?=\\r\\n!capacity)'))),
         capacity = as.numeric(str_trim(str_extract(raw, pattern = '(?<=\\!capacity)(.)*(?=\\r\\n\\!edgeWeightType)'))),
         edgeWeightType = str_trim(str_extract(raw, pattern = '(?<=\\!edgeWeightType:)(.)*(?=\\r\\n!edgeWeightFormat)')),
         edgeWeightFormat = str_trim(str_extract(raw, pattern = '(?<=\\!edgeWeightFormat:)(.)*(?=\\r\\n!edgeDataFormat)')),
         edgeDataFormat = str_trim(str_extract(raw, pattern = '(?<=\\!edgeDataFormat:)(.)*(?=\\r\\n!nodeCoordType)')),
         nodeCoordType = str_trim(str_extract(raw, pattern = '(?<=\\!nodeCoordType:)(.)*(?=\\r\\n!displayDataType)')),
         displayDataType = str_trim(str_extract(raw, pattern = '(?<=\\!displayDataType:)(.)*(?=\\r\\nNodes)')),
         nodes = lapply(str_split(str_trim(str_extract(raw, pattern = '(?<=Nodes: \\[\\r\\n)(.)*(?=\\r\\n)')), pattern = "\\s"), as.numeric),
         distance = lapply(str_split(str_remove_all(str_trim(str_extract(raw, pattern = '(?<=Distance\\:\\[\\r\\n)(.|\\n|\\r)*(?=\\s\\r\\n.\\r\\nPosX:)')), pattern = '\\n|\\r'), pattern = '\\s', simplify = FALSE), as.numeric),
         posX = lapply(str_split(str_trim(str_extract(raw, pattern = '(?<=PosX: \\[\\r\\n)(.)*(?=\\r\\n)')), pattern = '\\s', simplify = FALSE), as.numeric),
         posY = lapply(str_split(str_trim(str_extract(raw, pattern = '(?<=PosY: \\[\\r\\n)(.)*(?=\\r\\n)')), pattern = '\\s', simplify = FALSE), as.numeric),
         demand = lapply(str_split(str_trim(str_extract(raw, pattern = '(?<=Demand: \\[\\r\\n)(.)*(?=\\r\\n)')), pattern = '\\s', simplify = FALSE), as.numeric),
         draft = lapply(str_split(str_trim(str_extract(raw, pattern = '(?<=Draft: \\[\\r\\n)(.)*(?=\\r\\n)')), pattern = '\\s', simplify = FALSE), as.numeric)
         ) %>% arrange(name, n, m)
  

tabela[1,]$nodes
tabela[1,]$distance
tabela[1,]$demand
tabela[1,]$draft

teste <- tabela %>% 
  mutate(draft2 = sapply(draft, function(x) list(draft = as.integer(unlist(x))))) 

teste$draft2[1]

tabela %>% select(-path, -raw, -comment, -filename)



write_rds(tabela, file = "inst/extdata/TSPDL_RAW.rds")
write_rds(tabela %>% select(-path, -raw, -comment, -filename), file = "inst/extdata/TSPDL_DATA.rds")


library(tidyverse)
library(brkga)

tspdl_data <-read_rds("inst/extdata/TSPDL_DATA.rds")
##tspdl_raw <-read_rds("inst/extdata/TSPDL_RAW.rds")

#tabela
#apply(matrix(tabela$path, ncol = 1), 1L, function (x) length(read_file(x)))
#list(draft = as.integer(unlist(tspdl_data$draft[1])))

index <- 1
distance <- tspdl_data[index,]$distance %>% unlist
demand <- tspdl_data[index,]$demand %>% unlist
draft <- tspdl_data[index,]$draft %>% unlist
N <- tspdl_data[index,]$N

rst <- calado_brkga(DistanceVector = distance,
             DemandVector = demand,
             DraftVector = draft,
             N = N, 
             rngSeed = 0)
rst
table(rst$Tour)



sol_check <- function(rst, N, draft, demand, distance){
  if (sum((table(rst$Tour) == 1) == FALSE) > 0){
    cat("\nErro: Tour passa mais de uma vez por algum porto.")
    return(list(Contition = 'Erro!'))
  } 
  if (length(rst$Tour) != N) {
    cat("\nErro: Tour não passa por todos o nós, falta ", N - length(rst$Tour), "portos.")
    return(list(Contition = 'Erro!'))
  }
    
  if (rst$Tour[1] != 0) {
    cat("\nErro: Tour não inicia no porto 0")
    return(list(Contition = 'Erro!'))
  }
  
  mat_distance <- matrix(distance, nrow = N, ncol = N)
  
  dados <- tibble(N = 0:(N-1), Demand = demand, Draft = draft)
  
  table_tour <- tibble(Order = 1:length(rst$Tour), 
         Tour = rst$Tour) %>% 
    mutate(Next = lead(Tour, default = 0),
           #Distance = mat_distance[Tour, Next]
    ) %>% 
    group_by(Order, Tour, Next) %>% 
    summarise(Distance = mat_distance[Tour + 1,Next + 1], .groups = "drop") %>% 
    left_join(dados, by = c("Tour" = "N")) %>% 
    mutate(Demand_Cum = cumsum(Demand), Load = sum(dados$Demand) - Demand_Cum,
           Condition = if_else(Load <= Draft, "Success", "Fail")) 
  
  if (sum(table_tour$Condition == 'Fail') > 0) {
    cat("\nErro: Calado não atendido em alguma porto.")
    return(list(Contition = 'Erro!', table = table_tour))
  }
  
  if (sum(table_tour$Distance) != rst$BKFitness) {
    cat("\nErro: Fitness do tour não bate com resultado.")
    return(list(Contition = 'Erro!', table = table_tour))
  }
  
  return(list(Contition = 'Solução viável!', table = table_tour))
}

sol_check(rst = rst, N = N, draft = draft, demand = demand, distance = distance)



sol[1] <- 0
sol[6] <- 24

