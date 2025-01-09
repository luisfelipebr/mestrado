
# open libraries
library(tidyverse)
library(sf)
library(klaR)
library(diceR)

# create support functions
kmodes_elbow <- function(data, modes, seeds) {
  withindiff <- list()
  for (seed in seeds) {
    set.seed(seed)
    model <- kmodes(data = data,
                    modes = modes)
    withindiff[[seed]] <- model$withindiff
  }
  sum <- lapply(FUN = sum, withindiff)
  output <- list()
  output[[1]] <- mean(unlist(sum))
  output[[2]] <- sd(unlist(sum))
  print(Sys.time())
  return(output)
}

kmodes_ensemble <- function(data, modes, seed) {
  set.seed(seed)
  sample <- sample(seq_len(nrow(data)),size = floor(0.8*nrow(data)))
  data2 <- data[sample,]
  set.seed(seed)
  model <- kmodes(data = as.matrix(data2),
                  modes = modes)
  df <- data.frame(nrow = sample,
                   cluster = model[[1]])
  data <- data %>% 
    mutate(nrow = seq(1:nrow(data))) %>%
    left_join(df)
  return(data$cluster)
}

# import data
polos_nuis <- readRDS("input/POLOS_NUIS.rds") %>%
  filter(V5a %in% c(1, 2))

# apply pre-processing steps
nuis <- polos_nuis %>%
  st_drop_geometry() %>%
  mutate(
    tipo_do_nui = case_when(V5a == 1 ~ "1 Favela",
                                V5a == 2 ~ "2 Loteamento irregular"),
    tempo_de_estabelecimento = case_when(V6 %in% c(1,2,3) ~ "1 Até 10 anos",
                                   V6 == 4 ~ "2 Acima de 10 anos",
                                   V6 == 5 ~ "3 NA"),
    dinamica_imobiliaria = case_when(V7 %in% c(1,2) ~ "1 Surgimento",
                            V7 == 3 ~ "2 Estável",
                            V7 %in% c(4,5) ~ "3 Diminuição",
                            V7 == 6 ~ "4 NA"),
    contiguidade_urbana = case_when(V8 == 3 ~ "1 Área urbana consolidada",
                                 V8 %in% c(2,1) ~ "2 Periferia ou remoto"),
    zeis = case_when(V9 %in% c(1,3) ~ "1 Sim",
                             V9 == 2 ~ "2 Não",
                             V9 == 4 ~ "3 NA"),
    areas_protegidas = case_when(V10a %in% c(2, 3, 4, 5, 6, 7, 8, 9) ~ "1 Sim",
                                V10a == 1 ~ "2 Não"),
    app = case_when(V11a == 1 ~ "1 Sim",
                                   V11a == 0 ~ "2 Não"),
    situacao_de_risco = case_when(V12a == 3 ~ "1 Sim",
                               V12a == 2 ~ "2 Não",
                               V12a == 1 ~ "3 NA"),
    suscetibilidade_a_risco = case_when(V13a %in% c(3,4) ~ "1 Sim",
                                    V13a == 2 ~ "2 Não",
                                    V13a == 1 ~ "3 NA"),
    tracado = case_when(V14a %in% c(1,2) ~ "1 Presente",
                               V14a %in% c(3,4) ~ "2 Ausente",
                               V14a == 5 ~ "3 NA"),
    definicao_dos_lotes = case_when(V15a %in% c(1,2) ~ "1 Bem-definidos",
                               V15a %in% c(3,4) ~ "2 Indefinidos",
                               V15a == 5 ~ "3 NA"),
    ocupacao_dos_lotes = case_when(V15a %in% c(1, 3) ~ "1 Baixa",
                                 V15a %in% c(2, 4) ~ "2 Alta",
                                 V15a == 5 ~ "3 NA"),
    condicao_das_construcoes = case_when(V16a == 1 ~ "1 Adequada",
                                    V16a == 2 ~ "2 Mista",
                                    V16a == 3 ~ "3 Precária",
                                    V16a == 4 ~ "4 NA"),
    urbanizacao_e_infraestrutura = case_when(V17 == 1 ~ "1 Adequada",
                                         V17 == 2 ~ "2 Parcial",
                                         V17 == 3 ~ "3 Precária",
                                         V17 == 4 ~ "4 NA")
  ) %>%
  rename(id = V0) %>%
  dplyr::select(
    id,
    tipo_do_nui,
    contiguidade_urbana,
    tempo_de_estabelecimento,
    dinamica_imobiliaria,
    zeis,
    areas_protegidas,
    app,
    situacao_de_risco,
    suscetibilidade_a_risco,
    tracado,
    definicao_dos_lotes,
    ocupacao_dos_lotes,
    condicao_das_construcoes,
    urbanizacao_e_infraestrutura
  ) %>%
  mutate(across(2:15, ~as.factor(.)))

# elbow plot
model1 <- lapply(FUN = kmodes_elbow,
                 X = c(2:10),
                 data = as.matrix(nuis)[,-c(1, 2, 5)],
                 seeds = seq(1:100)) %>%
  lapply(FUN = unlist) %>%
  bind_cols() %>%
  t() %>%
  as_tibble() %>%
  setNames(c("mean", "sd")) %>% 
  mutate(model = c(2:10))

model1

ggplot(model1, aes(x=model, y=mean)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  geom_vline(xintercept = 5, linetype="dashed", size = 0.25) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw() +
  scale_x_continuous(breaks = seq(2, 10, by = 1)) +
  xlab("k") +
  ylab("distância") +
  theme(panel.grid.minor.x = element_blank())

# final model
models <- lapply(FUN = kmodes_ensemble,
                  data = nuis[,-c(1, 2, 5)],
                  modes = 5,
                  X = seq(1, 5000, 1)) %>%
  bind_cols() %>%
  as.matrix()

final1 <- models %>% 
  majority_voting(is.relabelled = FALSE)

table(final1)

final1 <- case_when(final1 == 5 ~ "C4",
                    final1 == 2 ~ "C3",
                    final1 == 4 ~ "C1",
                    final1 == 1 ~ "C5",
                    final1 == 3 ~ "C2")

output1 <- polos_nuis %>%
  left_join(nuis, by = c("V0" = "id")) %>%
  mutate(clusters = final1)

results_table <- tibble(
  variable = c(
    names(nuis)[[3]],
    names(nuis)[[3]],
    names(nuis)[[4]],
    names(nuis)[[4]],
    names(nuis)[[4]],
    names(nuis)[[6]],
    names(nuis)[[6]],
    names(nuis)[[6]],
    names(nuis)[[7]],
    names(nuis)[[7]],
    names(nuis)[[8]],
    names(nuis)[[8]],
    names(nuis)[[9]],
    names(nuis)[[9]],
    names(nuis)[[9]],
    names(nuis)[[10]],
    names(nuis)[[10]],
    names(nuis)[[10]],
    names(nuis)[[11]],
    names(nuis)[[11]],
    names(nuis)[[11]],
    names(nuis)[[12]],
    names(nuis)[[12]],
    names(nuis)[[12]],
    names(nuis)[[13]],
    names(nuis)[[13]],
    names(nuis)[[13]],
    names(nuis)[[14]],
    names(nuis)[[14]],
    names(nuis)[[14]],
    names(nuis)[[14]],
    names(nuis)[[15]],
    names(nuis)[[15]],
    names(nuis)[[15]],
    names(nuis)[[15]]
  ),
  categories = c(
    names(prop.table(table(nuis[,3]))),
    names(prop.table(table(nuis[,4]))),
    names(prop.table(table(nuis[,6]))),
    names(prop.table(table(nuis[,7]))),
    names(prop.table(table(nuis[,8]))),
    names(prop.table(table(nuis[,9]))),
    names(prop.table(table(nuis[,10]))),
    names(prop.table(table(nuis[,11]))),
    names(prop.table(table(nuis[,12]))),
    names(prop.table(table(nuis[,13]))),
    names(prop.table(table(nuis[,14]))),
    names(prop.table(table(nuis[,15])))
  ),
  total = c(
    round(prop.table(table(nuis[,3]))[[1]],2),
    round(prop.table(table(nuis[,3]))[[2]],2),
    round(prop.table(table(nuis[,4]))[[1]],2),
    round(prop.table(table(nuis[,4]))[[2]],2),
    round(prop.table(table(nuis[,4]))[[3]],2),
    round(prop.table(table(nuis[,6]))[[1]],2),
    round(prop.table(table(nuis[,6]))[[2]],2),
    round(prop.table(table(nuis[,6]))[[3]],2),
    round(prop.table(table(nuis[,7]))[[1]],2),
    round(prop.table(table(nuis[,7]))[[2]],2),
    round(prop.table(table(nuis[,8]))[[1]],2),
    round(prop.table(table(nuis[,8]))[[2]],2),
    round(prop.table(table(nuis[,9]))[[1]],2),
    round(prop.table(table(nuis[,9]))[[2]],2),
    round(prop.table(table(nuis[,9]))[[3]],2),
    round(prop.table(table(nuis[,10]))[[1]],2),
    round(prop.table(table(nuis[,10]))[[2]],2),
    round(prop.table(table(nuis[,10]))[[3]],2),
    round(prop.table(table(nuis[,11]))[[1]],2),
    round(prop.table(table(nuis[,11]))[[2]],2),
    round(prop.table(table(nuis[,11]))[[3]],2),
    round(prop.table(table(nuis[,12]))[[1]],2),
    round(prop.table(table(nuis[,12]))[[2]],2),
    round(prop.table(table(nuis[,12]))[[3]],2),
    round(prop.table(table(nuis[,13]))[[1]],2),
    round(prop.table(table(nuis[,13]))[[2]],2),
    round(prop.table(table(nuis[,13]))[[3]],2),
    round(prop.table(table(nuis[,14]))[[1]],2),
    round(prop.table(table(nuis[,14]))[[2]],2),
    round(prop.table(table(nuis[,14]))[[3]],2),
    round(prop.table(table(nuis[,14]))[[4]],2),
    round(prop.table(table(nuis[,15]))[[1]],2),
    round(prop.table(table(nuis[,15]))[[2]],2),
    round(prop.table(table(nuis[,15]))[[3]],2),
    round(prop.table(table(nuis[,15]))[[4]],2)
  ),
  C1 = c(
    round(prop.table(table(nuis$contiguidade_urbana,final1),margin=2)[1,1],2),
    round(prop.table(table(nuis$contiguidade_urbana,final1),margin=2)[2,1],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[1,1],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[2,1],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[3,1],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[1,1],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[2,1],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[3,1],2),
    round(prop.table(table(nuis$areas_protegidas,final1),margin=2)[1,1],2),
    round(prop.table(table(nuis$areas_protegidas,final1),margin=2)[2,1],2),
    round(prop.table(table(nuis$app,final1),margin=2)[1,1],2),
    round(prop.table(table(nuis$app,final1),margin=2)[2,1],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[1,1],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[2,1],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[3,1],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[1,1],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[2,1],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[3,1],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[1,1],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[2,1],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[3,1],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[1,1],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[2,1],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[3,1],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[1,1],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[2,1],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[3,1],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[1,1],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[2,1],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[3,1],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[4,1],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[1,1],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[2,1],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[3,1],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[4,1],2)
  ),
  C2 = c(
    round(prop.table(table(nuis$contiguidade_urbana,final1),margin=2)[1,2],2),
    round(prop.table(table(nuis$contiguidade_urbana,final1),margin=2)[2,2],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[1,2],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[2,2],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[3,2],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[1,2],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[2,2],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[3,2],2),
    round(prop.table(table(nuis$areas_protegidas,final1),margin=2)[1,2],2),
    round(prop.table(table(nuis$areas_protegidas,final1),margin=2)[2,2],2),
    round(prop.table(table(nuis$app,final1),margin=2)[1,2],2),
    round(prop.table(table(nuis$app,final1),margin=2)[2,2],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[1,2],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[2,2],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[3,2],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[1,2],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[2,2],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[3,2],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[1,2],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[2,2],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[3,2],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[1,2],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[2,2],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[3,2],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[1,2],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[2,2],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[3,2],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[1,2],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[2,2],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[3,2],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[4,2],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[1,2],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[2,2],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[3,2],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[4,2],2)
  ),
  C3 = c(
    round(prop.table(table(nuis$contiguidade_urbana,final1),margin=2)[1,3],2),
    round(prop.table(table(nuis$contiguidade_urbana,final1),margin=2)[2,3],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[1,3],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[2,3],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[3,3],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[1,3],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[2,3],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[3,3],2),
    round(prop.table(table(nuis$areas_protegidas,final1),margin=2)[1,3],2),
    round(prop.table(table(nuis$areas_protegidas,final1),margin=2)[2,3],2),
    round(prop.table(table(nuis$app,final1),margin=2)[1,3],2),
    round(prop.table(table(nuis$app,final1),margin=2)[2,3],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[1,3],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[2,3],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[3,3],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[1,3],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[2,3],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[3,3],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[1,3],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[2,3],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[3,3],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[1,3],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[2,3],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[3,3],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[1,3],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[2,3],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[3,3],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[1,3],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[2,3],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[3,3],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[4,3],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[1,3],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[2,3],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[3,3],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[4,3],2)
  ),
  C4 = c(
    round(prop.table(table(nuis$contiguidade_urbana,final1),margin=2)[1,4],2),
    round(prop.table(table(nuis$contiguidade_urbana,final1),margin=2)[2,4],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[1,4],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[2,4],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[3,4],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[1,4],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[2,4],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[3,4],2),
    round(prop.table(table(nuis$areas_protegidas,final1),margin=2)[1,4],2),
    round(prop.table(table(nuis$areas_protegidas,final1),margin=2)[2,4],2),
    round(prop.table(table(nuis$app,final1),margin=2)[1,4],2),
    round(prop.table(table(nuis$app,final1),margin=2)[2,4],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[1,4],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[2,4],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[3,4],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[1,4],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[2,4],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[3,4],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[1,4],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[2,4],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[3,4],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[1,4],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[2,4],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[3,4],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[1,4],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[2,4],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[3,4],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[1,4],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[2,4],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[3,4],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[4,4],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[1,4],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[2,4],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[3,4],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[4,4],2)
    ),
  C5 = c(
    round(prop.table(table(nuis$contiguidade_urbana,final1),margin=2)[1,5],2),
    round(prop.table(table(nuis$contiguidade_urbana,final1),margin=2)[2,5],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[1,5],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[2,5],2),
    round(prop.table(table(nuis$tempo_de_estabelecimento,final1),margin=2)[3,5],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[1,5],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[2,5],2),
    round(prop.table(table(nuis$zeis,final1),margin=2)[3,5],2),
    round(prop.table(table(nuis$areas_protegidas,final1),margin=2)[1,5],2),
    round(prop.table(table(nuis$areas_protegidas,final1),margin=2)[2,5],2),
    round(prop.table(table(nuis$app,final1),margin=2)[1,5],2),
    round(prop.table(table(nuis$app,final1),margin=2)[2,5],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[1,5],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[2,5],2),
    round(prop.table(table(nuis$situacao_de_risco,final1),margin=2)[3,5],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[1,5],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[2,5],2),
    round(prop.table(table(nuis$suscetibilidade_a_risco,final1),margin=2)[3,5],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[1,5],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[2,5],2),
    round(prop.table(table(nuis$tracado,final1),margin=2)[3,5],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[1,5],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[2,5],2),
    round(prop.table(table(nuis$definicao_dos_lotes,final1),margin=2)[3,5],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[1,5],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[2,5],2),
    round(prop.table(table(nuis$ocupacao_dos_lotes,final1),margin=2)[3,5],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[1,5],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[2,5],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[3,5],2),
    round(prop.table(table(nuis$condicao_das_construcoes,final1),margin=2)[4,5],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[1,5],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[2,5],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[3,5],2),
    round(prop.table(table(nuis$urbanizacao_e_infraestrutura,final1),margin=2)[4,5],2)
  ))

table(output1$clusters)
round(prop.table(table(output1$clusters)),3)
table(output1$V1, output1$clusters)
round(prop.table(table(output1$V1, output1$clusters),margin=2),2)
table(polos_nuis$V2b, final1)
round(prop.table(table(polos_nuis$V2b, final1), margin = 2),2)
table(output1$V5a, output1$clusters)
round(prop.table(table(output1$V5a, output1$clusters),margin=2),3)

# export data and tables
write_sf(output1, "output/nuis_clusters.gpkg")
write_excel_csv2(results_table, "resultados/tabela_final.csv")
write.csv2(table(polos_nuis$V1, final1), "output/table30.csv")
write.csv2(round(prop.table(table(polos_nuis$V1, final1), margin = 1),2), "output/table3.csv")
write.csv2(table(polos_nuis$V2b, final1), "output/table40.csv")
write.csv2(round(prop.table(table(polos_nuis$V2b, final1), margin = 2),2), "output/table4.csv")
write.csv2(table(nuis$tipo_do_nui, final1), "output/table50.csv")
write.csv2(round(prop.table(table(nuis$tipo_do_nui, final1), margin = 2),2), "output/table5.csv")

