# ==============================================================================
# 02_eda.R
#
# Análise Exploratória dos Dados (EDA) do dataset airway.
#
# Objetivos:
# - Verificar distribuição das contagens por amostra
# - Filtrar genes com baixa expressão
# - Transformar dados para análise (vst)
# - Realizar PCA para investigar agrupamento das amostras
# - Construir heatmap dos genes mais variáveis
#
# Projeto: airway-eda
# Autor: Emmanuel Humberto
# Data: maio/2026
# ==============================================================================

# Pacotes ----------------------------------------------------------------------

library(airway)               # dataset
library(SummarizedExperiment) # classe de dados
library(ggplot2)              # visualização (vem com tidyverse)
library(tidyverse)            # manipulação de dados

# Carrega o dataset ------------------------------------------------------------

data(airway)

# ==============================================================================
# Bloco 1: Distribuição das contagens por amostra
# ==============================================================================

# Extrai a matraiz de contagens brutas do objeto airway
contagens <- assay(airway)

# Confere dimensão
dim(contagens)

# Resumo estatístico das primeiras 4 amostras (pra não poluir o console)
summary(contagens[, 1:4])

# Transformação log2(x + 1) - estabiliza a distribuição
contagens_log <- log2(contagens + 1)

# Resumo das mesmas 4 amostras após transformação
summary(contagens_log[, 1:4])

# Boxpot das contagens log por amostra. Aqui precisa de formato "long"
# Vamos converter usando tidyverse
contagens_long <- contagens_log |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  pivot_longer(
    cols = -gene_id,
    names_to = "amostra",
    values_to = "log2_contagem"
  )

# Verifica formato
head(contagens_long)

# Boxplot
ggplot(contagens_long, aes(x = amostra, y = log2_contagem)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  labs(
    title = "Distribuição de contagens (log2) por amostra",
    subtitle = "Dataset airway - 8 amostras", 
    x = "Amostra",
    y = "log2(contagem + 1)"
  ) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ==============================================================================
# Bloco 2: Filtragem de genes com baixa expressão
# ==============================================================================

# --- Antes da filtragem: quantos genes têm zero em todas as amostras? ---------

# Soma das contagens por gene (linha) em todas as 8 amostras
totais_por_gene <- rowSums(contagens)

# Quantos genes têm soma TOTAL = 0 (não expressos em nenhuma amostra)
sum(totais_por_gene == 0)

# Quantos têm soma > 0?
sum(totais_por_gene > 0)

# Percentual de genes não expressos 
mean(totais_por_gene == 0) * 100

# --- Aplica o filtro -----------------------------------------------------------

# Critério: pelo menos 4 amostras com contagem >= 10
genes_a_manter <- rowSums(contagens >= 10) >= 4

# Quantos genes passam no filtro?
sum(genes_a_manter)

# Quantos foram removidos?
sum(!genes_a_manter)

# Percentual de genes mantidos
mean (genes_a_manter) * 100

# Aplica o filtro: cria objeto airway filtrado
airway_filtrado <- airway[genes_a_manter, ]

# Confere dimensões antes/depois
dim(airway)
dim(airway_filtrado)
