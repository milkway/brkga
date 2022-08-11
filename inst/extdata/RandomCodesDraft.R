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
         capacity = str_trim(str_extract(raw, pattern = '(?<=\\!capacity)(.)*(?=\\r\\n\\!edgeWeightType)')),
         edgeWeightType = str_trim(str_extract(raw, pattern = '(?<=\\!edgeWeightType:)(.)*(?=\\r\\n!edgeWeightFormat)')),
         edgeWeightFormat = str_trim(str_extract(raw, pattern = '(?<=\\!edgeWeightFormat:)(.)*(?=\\r\\n!edgeDataFormat)')),
         edgeDataFormat = str_trim(str_extract(raw, pattern = '(?<=\\!edgeDataFormat:)(.)*(?=\\r\\n!nodeCoordType)')),
         nodeCoordType = str_trim(str_extract(raw, pattern = '(?<=\\!nodeCoordType:)(.)*(?=\\r\\n!displayDataType)')),
         displayDataType = str_trim(str_extract(raw, pattern = '(?<=\\!displayDataType:)(.)*(?=\\r\\nNodes)')),
         nodes = str_split(str_trim(str_extract(raw, pattern = '(?<=Nodes: \\[\\r\\n)(.)*(?=\\r\\n)')), pattern = '\\s', simplify = FALSE),
         distance = str_split(str_remove_all(str_trim(str_extract(raw, pattern = '(?<=Distance\\:\\[\\r\\n)(.|\\n|\\r)*(?=\\s\\r\\n.\\r\\nPosX:)')), pattern = '\\n|\\r'), pattern = '\\s', simplify = FALSE),
         posX = str_split(str_trim(str_extract(raw, pattern = '(?<=PosX: \\[\\r\\n)(.)*(?=\\r\\n)')), pattern = '\\s', simplify = FALSE),
         posY = str_split(str_trim(str_extract(raw, pattern = '(?<=PosY: \\[\\r\\n)(.)*(?=\\r\\n)')), pattern = '\\s', simplify = FALSE),
         demand = str_split(str_trim(str_extract(raw, pattern = '(?<=Demand: \\[\\r\\n)(.)*(?=\\r\\n)')), pattern = '\\s', simplify = FALSE),
         draft = str_split(str_trim(str_extract(raw, pattern = '(?<=Draft: \\[\\r\\n)(.)*(?=\\r\\n)')), pattern = '\\s', simplify = FALSE)
  ) %>% arrange(name, n, m)


teste <- tabela %>% 
  mutate(draft2 = sapply(draft, function(x) list(draft = as.integer(unlist(x))))) 

teste$draft2[1]

tabela %>% select(-path, -raw, -comment, -filename)



write_rds(tabela, file = "inst/extdata/TSPDL_RAW.rds")
write_rds(tabela %>% select(-path, -raw, -comment, -filename), file = "inst/extdata/TSPDL_DATA.rds")

tspdl_data <-read_rds("inst/extdata/TSPDL_DATA.rds")

tabela
apply(matrix(tabela$path, ncol = 1), 1L, function (x) length(read_file(x)))

list(draft = as.integer(unlist(tspdl_data$draft[1])))

