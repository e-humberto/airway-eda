# ==============================================================================
# 03_export_figures.R
#
# Exporta as figuras finais da EDA para a pasta figures/.
# Requer que 02_eda.R tenha sido executado antes (objetos na memória).
#
# Projeto: airway-eda
# Autor: Emmanuel Humberto
# Data: maio/2026
# ==============================================================================


# Pacotes ----------------------------------------------------------------------

library(ggplot2)
library(pheatmap)



# --- Figura 1: Boxplot de distribuição (log2 vs VST) --------------------------

# Recria o plot do Bloco 3 e salva como PNG
plot_distribuicao <- ggplot(comparacao, aes(x = amostra, y = valor)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ transformacao, scales = "free_y") +
  labs(
    title = "Comparação de transformações: log2 vs VST",
    subtitle = "16.139 genes após filtragem - 8 amostras",
    x = "Amostra",
    y = "Valor transformado"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text (angle = 45, hjust = 1))

ggsave(
  filename = "figures/01_distribuicao_log2_vs_vst.png",
  plot = plot_distribuicao,
  width = 9,
  height = 5,
  dpi = 300,
  bg = "white"
)

# --- Figura 2: PCA customizado ------------------------------------------------

plot_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(
    title = "PCA - Dataset airway",
    subtitle = "Top 500 genes mais variáveis (VST)",
    x = paste0("PC1: ", percent_var[1], "% variância"),
    y = paste0("PC2: ", percent_var[2], "% variância"),
    color = "Tratamento",
    shape = "Linhagem"
  ) +
  scale_color_manual(values = c("untrt" = "steelblue", "trt" = "tomato")) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(
  filename = "figures/02_pca.png",
  plot = plot_pca,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


# --- Figura 3: Heatmap --------------------------------------------------------

# pheatmap não funciona com ggsave, pois tem mecanismo próprio de exportação
png(
  filename = "figures/03_heatmap_top50.png",
  width = 8,
  height = 8,
  units = "in",
  res = 300,
  bg = "white"
)

pheatmap(
  mat_top50_z,
  annotation_col = anotacao_amostras,
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Heatmap: 50 genes mais variáveis (z-score)",
  fontsize = 9
)

dev.off()
