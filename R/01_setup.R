# ==============================================================================
# 01_setup.R
#
# Carrega pacotes e o dataset airway. Inspeção inicial.
#
# Projeto: airway-eda
# Autor: Emmanuel Humberto
# Data: maio/2026
# ==============================================================================

# Pacotes ----------------------------------------------------------------------

library(airway)               # dataset de RNA-seq
library(SummarizedExperiment) # classe de dados (dependência do airway)

# Carrega o dataset ------------------------------------------------------------

data(airway)

# Inspeção inicial -------------------------------------------------------------

# Visão geral do objeto
airway

# Dimensões: número de genes (linhas) e amostras (colunas)
dim(airway)

# Metadados das amostras: tratamento, linhagem celular, etc.
colData(airway)

# Metadados dos genes: primeiras 6 linhas
head(rowData(airway))

# Matriz de contagens: primeiras 6 linhas, todas as 8 colunas
head(assay(airway))

# Informações da sessão (reprodutibilidade) ------------------------------------
sessionInfo()