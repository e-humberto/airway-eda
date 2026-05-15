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


# ==============================================================================
# Bloco 3: Transformação dos dados (VST)
# ==============================================================================

# Pacote DESeq2 para a transformação VSI
library(DESeq2)

# --- Cria objeto DESeqDataSet a partir do airway filtrado ---------------------

# Precisamos especificar o "design", ou seja, quais variáveis o modelo vai considerar.
# Por enquanto, ~1 (sem covariáveis). Vamos só transformar, não testar.
dds <- DESeqDataSet(airway_filtrado, design = ~ 1)

# Inspeciona
dds

# --- Aplica VST ---------------------------------------------------------------

# blind = TRUE: transformação independente do design experimental
# (apropriada para EDA, onde não queremos viés do "grupo conhecido")
vsd <- vst(dds, blind = TRUE)

# Inspeciona
vsd


# --- Compara distribuições: antes vs depois -----------------------------------

# Extrai matriz transformada
contagens_vst <- assay(vsd)

# Resumo estatístico das primeiras 4 amostras (compare com bloco 1)
summary(contagens_vst[, 1:4])


# --- Visualização: log2 vs VST ------------------------------------------------

# Já temos contagens_log do Bloco 1.
# Vamos fazer dois boxplots lado a lado para comparar

# Prepara dados em formato long - LOG2 (filtrado para os mesmos genes do VST)
log_long <- contagens_log[genes_a_manter, ] |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  pivot_longer(
    cols = -gene_id,
    names_to = "amostra",
    values_to = "valor"
  ) |>
  mutate(transformacao = "log2(count + 1)")

# Prepara dados em formato long - VST
vst_long <- contagens_vst |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  pivot_longer(
    cols = -gene_id,
    names_to = "amostra",
    values_to = "valor"
  ) |>
  mutate(transformacao = "VST")

# Combina os dois
comparacao <- bind_rows(log_long, vst_long)

# Boxplot comparativo dos dois painéis
ggplot(comparacao, aes(x = amostra, y = valor)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ transformacao, scales = "free_y") +
  labs(
    title = "Comparação: log2 simples vs VST",
    subtitle = "Após filtragem - 16.139 genes",
    x = "Amostra",
    y = "Valor transformado"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ==============================================================================
# Bloco 4: PCA (Principal Component Analysis)
# ==============================================================================

# --- PCA usando a função do DESeq2 --------------------------------------------

# DESeq2 tem função plotPCA() conveniente.
# Argumentos importantes:
#   ingroup: variável(is) do colData para colotrir os pontos
#   ntop: quantos genes mais variáveis usar (padrão = 500)
#   returnData: TRUE retorna a tabela (para customizar plot)

# PCA colorindo por tratamento (dex)
plotPCA(vsd, intgroup = "dex")

# PCA colorindo por linhagem celular (cell)
plotPCA(vsd, intgroup = "cell")

# PCa combinando ambos (tratamento + linhagem)
plotPCA(vsd, intgroup = c("dex", "cell"))


# --- PCA customizado: extrai os dados e plota com ggplot2 ---------------------

# Extrai os dados do PCA - permite customização total
pca_data <- plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE)

# Inspeciona estrutura
head(pca_data)

# Extrai percentual de variância explicada
percent_var <- round(100 * attr(pca_data, "percentVar"))
percent_var


# --- Plot customizado: tratamento na cor, doador na forma ---------------------

ggplot(pca_data, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(
    title = "PCA - Dataset airway",
    subtitle = "Top 500 mais variáveis (VST)",
    x = paste0("PC1: ", percent_var[1], "% variância"),
    y = paste0("PC2: ", percent_var[2], "% variância"),
    color = "Tratamento",
    shape = "Linhagem"
  ) +
  scale_color_manual(values = c("untrt" = "steelblue", "trt" = "tomato")) +
  theme_minimal() +
  theme(legend.position = "right")