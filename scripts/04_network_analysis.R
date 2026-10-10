# Regularized network estimation, bridge measures, community structure and NCT.


library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(qgraph)
library(bootnet)
library(networktools)
library(NetworkComparisonTest)
library(officer)
library(flextable)
library(igraph)
library(gridExtra)

OUTPUT_DIR <- Sys.getenv("CEBP_OUTPUT_DIR", unset = file.path(getwd(), "output"))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

df <- readRDS(file.path(OUTPUT_DIR, "cohort_clean_N805.rds"))
cat("Loaded N =", nrow(df), "\n")

set.seed(2026)
k_vars <- c("K_Smoke", "K_Obesity", "K_PA", "K_Alcohol", "K_HotFood",
            "K_Infection", "K_PlantFiber", "K_FriedMeat", "K_Sugar", "K_Breast")

p_vars <- c("smoke_score", "BMI_score", "pa_score", "alcohol_score",
            "hotfood_score", "hepb_score", "hpv_score", "fruvege_score",
            "wholflour_score", "dessert_score",
            "fried_score", "fastfood_score", "redmeat_score")

all_vars  <- c(k_vars, p_vars)
n_k       <- length(k_vars)
n_p       <- length(p_vars)
n_nodes   <- length(all_vars)

vars_T0 <- paste0(all_vars, "_T0")
vars_T1 <- paste0(all_vars, "_T1")

df_net <- df %>% select(all_of(c(vars_T0, vars_T1))) %>% na.omit()
cat("Network Analysis N =", nrow(df_net), "\n")

data_T0 <- df_net %>% select(all_of(vars_T0)) %>% setNames(all_vars)
data_T1 <- df_net %>% select(all_of(vars_T1)) %>% setNames(all_vars)

data_T1$hepb_score <- pmax(data_T0$hepb_score, data_T1$hepb_score)
data_T1$hpv_score  <- pmax(data_T0$hpv_score,  data_T1$hpv_score)

node_groups <- list(
  Knowledge = 1:n_k,
  Behavior  = (n_k + 1):n_nodes
)
comm_vec <- c(rep("Knowledge", n_k), rep("Behavior", n_p))

short_labels <- c(
  "K-Smoke",  "K-Obese",  "K-PA",    "K-Alc",   "K-Hot",
  "K-Infect", "K-Fiber",  "K-Meat",  "K-Sugar", "K-Breast",
  "B-Smoke",  "B-BMI",    "B-PA",    "B-Alc",   "B-Hot",
  "B-HepB",   "B-HPV",    "B-Veg",   "B-Grain", "B-Sweet",
  "B-Fried",  "B-Fast",   "B-RedM"
)

full_labels <- c(
  "K: Smoking",       "K: Obesity",      "K: Physical Act.", "K: Alcohol",
  "K: Hot Food",      "K: Infection",    "K: Plant Fiber",   "K: Fried Meat",
  "K: Added Sugar",   "K: Breastfeeding",
  "B: Smoking",       "B: BMI",          "B: Physical Act.", "B: Alcohol",
  "B: Hot Food",      "B: Hep B vacc.",  "B: HPV vacc.",     "B: Veg & Fruit",
  "B: Dietary Fiber",   "B: Dessert",
  "B: Fried Food",    "B: Fast Food",    "B: Red Meat"
)

COL_K   <- "#FDBF6F"
COL_P   <- "#CAB2D6"
COL_POS <- "darkgreen"
COL_NEG <- "darkred"
COL_GREY <- "grey85"

group_colors <- c(COL_K, COL_P)

cat("\n--- Correlation Matrix Diagnostic (T0) ---\n")

cor_mat_T0 <- cor_auto(data_T0)
cor_range_T0 <- range(cor_mat_T0[upper.tri(cor_mat_T0)])
cat("T0 correlation range:", round(cor_range_T0, 4), "\n")

suspect_T0 <- which(abs(cor_mat_T0) > 0.95 & abs(cor_mat_T0) < 1, arr.ind = TRUE)
cat("T0 high-correlation pairs (|r| > 0.95):", nrow(suspect_T0), "\n")
if (nrow(suspect_T0) > 0) {
  suspect_T0_named <- data.frame(
    Var1 = rownames(cor_mat_T0)[suspect_T0[, "row"]],
    Var2 = colnames(cor_mat_T0)[suspect_T0[, "col"]],
    r    = cor_mat_T0[suspect_T0]
  )
  print(suspect_T0_named)
}

cat("\n--- Correlation Matrix Diagnostic (T1) ---\n")

cor_mat_T1 <- cor_auto(data_T1)
cor_range_T1 <- range(cor_mat_T1[upper.tri(cor_mat_T1)])
cat("T1 correlation range:", round(cor_range_T1, 4), "\n")

suspect_T1 <- which(abs(cor_mat_T1) > 0.95 & abs(cor_mat_T1) < 1, arr.ind = TRUE)
cat("T1 high-correlation pairs (|r| > 0.95):", nrow(suspect_T1), "\n")
if (nrow(suspect_T1) > 0) {
  suspect_T1_named <- data.frame(
    Var1 = rownames(cor_mat_T1)[suspect_T1[, "row"]],
    Var2 = colnames(cor_mat_T1)[suspect_T1[, "col"]],
    r    = cor_mat_T1[suspect_T1]
  )
  print(suspect_T1_named)
}

cor_diagnostic_summary <- data.frame(
  Wave = c("T0", "T1"),
  Cor_Range_Min = c(cor_range_T0[1], cor_range_T1[1]),
  Cor_Range_Max = c(cor_range_T0[2], cor_range_T1[2]),
  N_Suspect_Edges = c(nrow(suspect_T0), nrow(suspect_T1))
)
print(cor_diagnostic_summary)


cat("\n--- Estimating T0 and T1 Networks ---\n")

net_T0 <- estimateNetwork(data_T0, default = "EBICglasso",
                          corMethod = "cor_auto", tuning = 0.5)
net_T1 <- estimateNetwork(data_T1, default = "EBICglasso",
                          corMethod = "cor_auto", tuning = 0.5)

get_edge_list <- function(net_obj, node_names) {
  W <- net_obj$graph
  W[lower.tri(W, diag = TRUE)] <- NA
  edges <- as.data.frame(as.table(W)) %>%
    filter(!is.na(Freq), Freq != 0) %>%
    rename(Node1 = Var1, Node2 = Var2, Weight = Freq)
  edges
}

edges_T0 <- get_edge_list(net_T0, all_vars)
edges_T1 <- get_edge_list(net_T1, all_vars)

edges_T0 <- get_edge_list(net_T0, all_vars) %>% rename(Weight_T0 = Weight)
edges_T1 <- get_edge_list(net_T1, all_vars) %>% rename(Weight_T1 = Weight)

edges_T0 <- edges_T0 %>% filter(abs(Weight_T0) >= 0.05)
edges_T1 <- edges_T1 %>% filter(abs(Weight_T1) >= 0.05)

edges_all <- full_join(edges_T0, edges_T1, by = c("Node1", "Node2")) %>%
  mutate(Weight_T0 = round(Weight_T0, 3),
         Weight_T1 = round(Weight_T1, 3))

label_lookup <- setNames(short_labels, all_vars)
edges_all <- edges_all %>%
  mutate(
    Node1 = label_lookup[as.character(Node1)],
    Node2 = label_lookup[as.character(Node2)],
    EdgeType = case_when(
      grepl("^K-", Node1) & grepl("^K-", Node2) ~ "K-K",
      grepl("^B-", Node1) & grepl("^B-", Node2) ~ "B-B",
      TRUE ~ "K-B"
    )
  ) %>%
  arrange(EdgeType, desc(pmax(abs(Weight_T0), abs(Weight_T1), na.rm = TRUE)))

ft_S10 <- flextable(edges_all) %>%
  set_header_labels(
    Node1 = "Node 1", Node2 = "Node 2",
    Weight_T0 = "Edge Weight (T0)", Weight_T1 = "Edge Weight (T1)",
    EdgeType = "Edge Type"
  ) %>%
  set_caption("Complete edge weight list for the estimated networks (T0 vs. T1)") %>%
  bg(part = "header", bg = "#1a4a6b") %>%
  color(part = "header", color = "white") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = c("Node1", "Node2"), align = "left", part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  add_footer_lines(
    "Note. Edge weights are partial correlation coefficients from the EBICglasso-estimated Gaussian graphical model (γ = 0.5, cor_auto). Only edges with |weight| >= 0.05 at either wave are listed, consistent with the visualization threshold (minimum = 0.05) used in the network visualization; a blank cell indicates the edge did not meet this threshold at that wave. Edge Type denotes whether the edge connects two Knowledge nodes (K-K), two Behavior nodes (B-B), or one of each (K-B)."
  ) %>%
  fontsize(size = 7, part = "footer") %>%
  italic(part = "footer") %>%
  autofit() %>%
  set_table_properties(layout = "autofit", opts_word = list(split = FALSE))
ft_S10

save_as_docx(
  ft_S10,
  path = file.path(OUTPUT_DIR, "Network_Edge_Weights.docx"),
  pr_section = officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(top = 0.5, bottom = 0.5, left = 0.5, right = 0.5)
  )
)
cat("Saved: Network_Edge_Weights.docx saved\n")

saveRDS(edges_all, file.path(OUTPUT_DIR, "Network_Edge_Weights.rds"))
cat("Saved: Network_Edge_Weights.rds saved\n")

cat("\n--- Plotting networks (combined PNG) ---\n")

set.seed(2026)
avg_layout <- averageLayout(net_T0, net_T1, layout = "spring", layout.par = list(niter = 500))

net_params <- list(
  layout       = avg_layout,
  groups       = node_groups,
  color        = group_colors,
  labels       = short_labels,
  label.color  = "black",
  label.cex    = 1.5,
  label.scale  = FALSE,
  vsize        = 10.5,
  border.color = "white",
  border.width = 2.5,
  minimum      = 0.05,
  cut          = 0,
  esize        = 8,
  posCol       = COL_POS,
  negCol       = COL_NEG,
  fade         = FALSE,
  legend       = FALSE
)

png(file.path(OUTPUT_DIR, "Network_Baseline_Followup.png"),
    width = 16, height = 8, units = "in", res = 300, type = "cairo")

par(cex = 1.3, mar = c(2, 2, 3, 2), oma = c(0, 0, 2, 0))

layout(matrix(c(1, 2, 3), nrow = 1), widths = c(4.5, 4.5, 2))

do.call(qgraph, c(list(input = net_T0$graph, title = "Baseline (T0)"), net_params))
do.call(qgraph, c(list(input = net_T1$graph, title = "Follow-up (T1)"), net_params))

par(mar = c(2, 1, 2, 1))
plot.new()

legend("center",
       legend = c("Knowledge nodes", "Behavior nodes", "Positive association", "Negative association"),
       pch    = c(21, 21, NA, NA),
       pt.bg  = c(COL_K, COL_P, NA, NA),
       col    = c("black", "black", COL_POS, COL_NEG),
       lwd    = c(NA, NA, 3, 3),
       lty    = c(NA, NA, 1, 1),
       pt.cex = 2.5,
       cex    = 1.3,
       bty    = "n",
       y.intersp = 1.5)

dev.off()
cat("Saved: Network_Combined_updated.png saved\n")


cat("\n--- Calculating Bridge Expected Influence ---\n")

bridge_T0 <- bridge(net_T0$graph, communities = node_groups, useCommunities = "all")
bridge_T1 <- bridge(net_T1$graph, communities = node_groups, useCommunities = "all")

ei_T0 <- as.numeric(bridge_T0$`Bridge Expected Influence (1-step)`)
ei_T1 <- as.numeric(bridge_T1$`Bridge Expected Influence (1-step)`)
bs_T0 <- as.numeric(bridge_T0$`Bridge Strength`)
bs_T1 <- as.numeric(bridge_T1$`Bridge Strength`)

bei_df <- data.frame(
  node      = full_labels,
  community = comm_vec,
  T0        = ei_T0,
  T1        = ei_T1,
  stringsAsFactors = FALSE
)

bei_df <- bei_df %>%
  arrange(T0) %>%
  mutate(node = factor(node, levels = node))

bei_long <- bei_df %>%
  pivot_longer(cols = c(T0, T1), names_to = "wave", values_to = "BEI") %>%
  mutate(
    wave = factor(wave, levels = c("T0", "T1"),
                  labels = c("Baseline (T0)", "Follow-up (T1)")),
    fill_col = ifelse(community == "Knowledge", COL_K, COL_P)
  )

p_bei <- ggplot(bei_long, aes(x = BEI, y = node, color = community)) +
  geom_segment(
    aes(x = 0, xend = BEI, y = node, yend = node),
    linewidth = 0.7, color = COL_GREY
  ) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.45,
             linetype = "solid") +
  geom_point(size = 3.2) +
  facet_wrap(~wave, nrow = 1) +
  scale_color_manual(
    values = c(Knowledge = COL_K, Behavior = COL_P),
    name   = NULL
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  labs(
    x = "Bridge Expected Influence (1-step)",
    y = NULL
  ) +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    panel.grid.major.x  = element_line(color = "grey93", linewidth = 0.4),
    panel.grid.major.y  = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.spacing       = unit(1.2, "lines"),

    strip.text          = element_text(face = "bold", size = 12,
                                       margin = margin(b = 6)),
    strip.background    = element_rect(fill = "grey97", color = NA),

    axis.text.y         = element_text(size = 9.5, color = "grey20",
                                       hjust = 1),
    axis.text.x         = element_text(size = 9),
    axis.title.x        = element_text(size = 10.5, margin = margin(t = 8)),
    axis.ticks.y        = element_blank(),

    legend.position     = "bottom",
    legend.key.size     = unit(0.9, "lines"),
    legend.text         = element_text(size = 10),
    legend.margin       = margin(t = 2),

    plot.margin         = margin(10, 15, 6, 10)
  )

ggsave(file.path(OUTPUT_DIR, "Bridge_EI_Combined.png"),
       p_bei, width = 11, height = 7, dpi = 180, bg = "white")
cat("Saved: Bridge_EI_Combined.png saved\n")


ei_T0 <- as.numeric(bridge_T0$`Bridge Expected Influence (1-step)`)
ei_T1 <- as.numeric(bridge_T1$`Bridge Expected Influence (1-step)`)

bridge_stability <- data.frame(
  Node        = short_labels,
  Community   = c(rep("Knowledge", n_k), rep("Behavior", n_p)),
  BS_T0       = round(bs_T0, 3),
  BEI_T0      = round(ei_T0, 3),
  Rank_T0     = rank(-ei_T0, ties.method = "min"),
  BS_T1       = round(bs_T1, 3),
  BEI_T1      = round(ei_T1, 3),
  Rank_T1     = rank(-ei_T1, ties.method = "min"),
  Rank_Change = rank(-ei_T1, ties.method = "min") - rank(-ei_T0, ties.method = "min")
) %>% arrange(Rank_T0)

print(bridge_stability)

cat("\n--- Bridge EI Stability (CS-coefficient) ---\n")

boot_bridge_T0 <- bootnet(net_T0,
                          nBoots     = 1000,
                          type       = "case",
                          statistics = "bridgeStrength",
                          communities = node_groups,
                          caseMin    = 0.05,
                          caseMax    = 0.95)
cs_bridge_T0 <- corStability(boot_bridge_T0, statistics = "bridgeStrength")
cat("Bridge Strength CS (T0):", cs_bridge_T0, "\n")

boot_bridge_T1 <- bootnet(net_T1,
                          nBoots     = 1000,
                          type       = "case",
                          statistics = "bridgeStrength",
                          communities = node_groups,
                          caseMin    = 0.05,
                          caseMax    = 0.95)
cs_bridge_T1 <- corStability(boot_bridge_T1, statistics = "bridgeStrength")
cat("Bridge Strength CS (T1):", cs_bridge_T1, "\n")


library(patchwork)
p1 <- plot(boot_bridge_T0, statistics = "bridgeStrength") +
  ggtitle("Bridge Strength Stability — Baseline (T0)") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12, hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93")
  )

p2 <- plot(boot_bridge_T1, statistics = "bridgeStrength") +
  ggtitle("Bridge Strength Stability — Follow-up (T1)") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12, hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93")
  )

p_bridge_combined <- p1 / p2 +
  plot_annotation(
    caption = "Dashed line: r = 0.70 (correlation threshold used to calculate CS). CS > 0.50 = acceptable; > 0.75 = excellent.",
    theme   = theme(plot.caption = element_text(size = 9, color = "grey50"))
  )

ggsave(file.path(OUTPUT_DIR, "Bridge_Stability_Combined.png"),
       p_bridge_combined, width = 9, height = 9, dpi = 180, bg = "white")
cat("Saved: Bridge_Stability_Combined.png saved\n")

cat("\n--- Bridge Strength CI per node (nonparametric bootstrap) ---\n")

boot_bridge_T0_np <- bootnet(net_T0,
                             nBoots = 1000,
                             type = "nonparametric",
                             statistics = "bridgeStrength",
                             communities = node_groups)

boot_bridge_T1_np <- bootnet(net_T1,
                             nBoots = 1000,
                             type = "nonparametric",
                             statistics = "bridgeStrength",
                             communities = node_groups)

p3 <- plot(boot_bridge_T0_np, statistics = "bridgeStrength",
           order = "sample", labels = TRUE) + ggtitle("Bridge Strength CI - T0")
p4 <- plot(boot_bridge_T1_np, statistics = "bridgeStrength",
           order = "sample", labels = TRUE) + ggtitle("Bridge Strength CI - T1")

png(file.path(OUTPUT_DIR, "Bridge_CI_Combined.png"),
    width = 2400, height = 1200, res = 180)
grid.arrange(p3, p4, ncol = 2)
dev.off()
cat("Saved: Bridge_CI_Combined.png saved\n")

summary_T0 <- summary(boot_bridge_T0_np, statistics = "bridgeStrength")
summary_T1 <- summary(boot_bridge_T1_np, statistics = "bridgeStrength")

print(summary_T0)
print(summary_T1)

supp_T0 <- summary_T0 %>%
  ungroup() %>%
  select(Node = node1, Sample = sample, Mean = mean,
         CIlower, CIupper, Prop0 = prop0) %>%
  mutate(
    Wave = "T0",
    Sig_NonZero = ifelse(CIlower > 0, "Yes", "No"),
    across(c(Sample, Mean, CIlower, CIupper), ~round(.x, 3)),
    Prop0 = round(Prop0, 3)
  )

supp_T1 <- summary_T1 %>%
  ungroup() %>%
  select(Node = node1, Sample = sample, Mean = mean,
         CIlower, CIupper, Prop0 = prop0) %>%
  mutate(
    Wave = "T1",
    Sig_NonZero = ifelse(CIlower > 0, "Yes", "No"),
    across(c(Sample, Mean, CIlower, CIupper), ~round(.x, 3)),
    Prop0 = round(Prop0, 3)
  )

bridge_CI_wide <- supp_T0 %>%
  select(Node, Sample_T0 = Sample, CIlower_T0 = CIlower,
         CIupper_T0 = CIupper, Prop0_T0 = Prop0, Sig_T0 = Sig_NonZero) %>%
  left_join(
    supp_T1 %>%
      select(Node, Sample_T1 = Sample, CIlower_T1 = CIlower,
             CIupper_T1 = CIupper, Prop0_T1 = Prop0, Sig_T1 = Sig_NonZero),
    by = "Node"
  ) %>%
  arrange(desc(Sample_T0))


cat("\n--- Global Expected Influence ---\n")

calc_EI <- function(adj_matrix) {
  diag(adj_matrix) <- 0

  ei_values <- rowSums(adj_matrix, na.rm = TRUE)

  return(ei_values)
}

global_ei_T0 <- calc_EI(net_T0$graph)
global_ei_T1 <- calc_EI(net_T1$graph)

stopifnot(
  exists("bridge_stability"),
  all(c("Node", "BS_T0", "BEI_T0", "BS_T1", "BEI_T1") %in% names(bridge_stability))
)

bridge_stability <- bridge_stability %>%
  mutate(
    GlobalEI_T0 = round(global_ei_T0[match(Node, short_labels)], 3),
    GlobalEI_T1 = round(global_ei_T1[match(Node, short_labels)], 3),
    Type = case_when(
      BS_T0 > median(BS_T0) & abs(GlobalEI_T0) > median(abs(GlobalEI_T0)) ~ "Core connector",
      BS_T0 > median(BS_T0) & abs(GlobalEI_T0) <= median(abs(GlobalEI_T0)) ~ "Bridge-only node",
      TRUE ~ "Other node"
    )
  )

print(bridge_stability)


node_map <- tibble::tribble(
  ~ShortLabel, ~FullName,
  "B-HepB",   "hepb_score",
  "B-PA",     "pa_score",
  "K-Smoke",  "K_Smoke",
  "B-HPV",    "hpv_score",
  "K-Fiber",  "K_PlantFiber",
  "K-Hot",    "K_HotFood",
  "B-Grain",  "wholflour_score",
  "K-PA",     "K_PA",
  "K-Sugar",  "K_Sugar",
  "K-Meat",   "K_FriedMeat",
  "K-Infect", "K_Infection",
  "B-Smoke",  "smoke_score",
  "K-Breast", "K_Breast",
  "B-Veg",    "fruvege_score",
  "K-Alc",    "K_Alcohol",
  "B-BMI",    "BMI_score",
  "B-Alc",    "alcohol_score",
  "B-Hot",    "hotfood_score",
  "B-Fried",  "fried_score",
  "K-Obese",  "K_Obesity",
  "B-Fast",   "fastfood_score",
  "B-Sweet",  "dessert_score",
  "B-RedM",   "redmeat_score"
)

tableS11_data <- bridge_stability %>%
  left_join(node_map, by = c("Node" = "ShortLabel")) %>%
  left_join(
    bridge_CI_wide %>% select(Node, CIlower_T0, CIupper_T0, CIlower_T1, CIupper_T1),
    by = c("FullName" = "Node")
  ) %>%
  arrange(desc(BS_T0)) %>%
  mutate(
    BS_T0  = round(BS_T0, 3), BS_T1  = round(BS_T1, 3),
    BEI_T0 = round(BEI_T0, 3), BEI_T1 = round(BEI_T1, 3),
    GlobalEI_T0 = round(GlobalEI_T0, 3),
    CI_T0 = sprintf("[%.3f, %.3f]", CIlower_T0, CIupper_T0),
    CI_T1 = sprintf("[%.3f, %.3f]", CIlower_T1, CIupper_T1)
  ) %>%
  select(Node, Community, BS_T0, CI_T0, BS_T1, CI_T1, BEI_T0, BEI_T1, GlobalEI_T0)

ft <- flextable(tableS11_data) %>%
  set_header_labels(
    Node = "Node", Community = "Community",
    BS_T0 = "Bridge Strength (T0)", CI_T0 = "95% CI (T0)",
    BS_T1 = "Bridge Strength (T1)", CI_T1 = "95% CI (T1)",
    BEI_T0 = "Bridge EI, 1-step (T0)", BEI_T1 = "Bridge EI, 1-step (T1)",
    GlobalEI_T0 = "Global EI (T0)"
  ) %>%
  set_caption("Node-level bridge strength, bridge expected influence, and global expected influence (T0 vs. T1)") %>%
  bg(part = "header", bg = "#1a4a6b") %>%
  color(part = "header", color = "white") %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "Node", align = "left", part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  add_footer_lines(
    "Note. Bridge Strength is the unsigned sum of absolute cross-community edge weights; 95% CIs derived from nonparametric bootstrap (1000 resamples). Bridge Expected Influence (1-step) is the signed sum of cross-community edge weights, reported as point estimate only (no bootstrap CI computed for this statistic). Correlation-stability (CS) coefficients from case-dropping bootstrap: Bridge Strength CS = 0.25 at both waves; Global Strength CS = 0.75 at both waves."
  ) %>%
  fontsize(size = 7, part = "footer") %>%
  italic(part = "footer") %>%
  autofit()
ft

save_as_docx(ft, path = file.path(OUTPUT_DIR, "Network_Bridge_Centrality.docx"))
cat("Saved: Network_Bridge_Centrality.docx saved\n")

saveRDS(tableS11_data, file.path(OUTPUT_DIR, "Network_Bridge_Centrality.rds"))
cat("Saved: Network_Bridge_Centrality.rds saved\n")

cat("\n--- Running NCT (1000 permutations) ---\n")
cat("This may take a few minutes...\n")

set.seed(2026)
nct_res <- NCT(data_T0, data_T1,
               it          = 1000,
               binary.data = FALSE,
               paired      = TRUE,
               test.edges  = TRUE,
               edges       = "all",
               gamma       = 0.5)

cat("\n[NCT Results]\n")
cat("Network Structure Invariance (M):  p =", round(nct_res$nwinv.pval, 3), "\n")
cat("Global Strength Invariance (S):    p =", round(nct_res$glstrinv.pval, 3), "\n")

if (!is.null(nct_res$einv.pvals)) {
  sig_edges <- nct_res$einv.pvals %>% filter(`p-value` < 0.05)

  if (nrow(sig_edges) > 0) {
    cat("\nSignificant edge differences:\n")
    print(sig_edges)
  } else {
    cat("No individual edges differ significantly after permutation.\n")
  }
}

network_results <- list(
  nct_res          = nct_res,
  bridge_stability = bridge_stability,
  cs_bridge_T0     = cs_bridge_T0,
  cs_bridge_T1     = cs_bridge_T1
)
saveRDS(network_results, file.path(OUTPUT_DIR, "04_NCT_results_N805.rds"))
cat("Saved: NCT results saved\n")

adjusted_p <- p.adjust(nct_res$einv.pvals$`p-value`, method = "holm")

nct_res$einv.pvals$p_holm <- adjusted_p

sig_after_correction <- nct_res$einv.pvals %>%
  filter(p_holm < 0.05)

cat(sprintf("\nNominally significant edges: %d\n", sum(nct_res$einv.pvals$`p-value` < 0.05)))
cat(sprintf("Edges significant after Holm correction: %d\n", nrow(sig_after_correction)))

if (nrow(sig_after_correction) > 0) {
  print(sig_after_correction[, c("Var1", "Var2", "p-value", "p_holm")])
} else {
  cat("No edges remained significant after Holm correction.\n")
}

write.csv(sig_edges, file.path(OUTPUT_DIR, "NCT_Edge_Nominal_Sig_Before_Holm.csv"), row.names = FALSE)

library(igraph)

diet_unhealthy_5 <- c("hotfood_score", "dessert_score", "fried_score",
                      "fastfood_score", "redmeat_score")
diet_core_4       <- c("hotfood_score", "dessert_score", "fried_score",
                       "fastfood_score")

detect_communities <- function(net_obj, node_names) {
  W <- net_obj$graph[node_names, node_names]
  g_signed <- graph_from_adjacency_matrix(W, mode = "undirected",
                                          weighted = TRUE, diag = FALSE)
  g_abs <- g_signed
  E(g_abs)$weight <- abs(E(g_abs)$weight)
  wt <- cluster_walktrap(g_abs, weights = E(g_abs)$weight, steps = 4)
  set.seed(2026)
  sg <- cluster_spinglass(g_signed, weights = E(g_signed)$weight,
                          implementation = "neg")
  list(walktrap = wt, spinglass = sg, graph = g_signed)
}

comm_T0 <- detect_communities(net_T0, p_vars)
comm_T1 <- detect_communities(net_T1, p_vars)

check_cluster <- function(comm_result, target_nodes, node_names) {
  idx <- match(target_nodes, node_names)
  list(
    walktrap_same  = length(unique(membership(comm_result$walktrap)[idx])) == 1,
    walktrap_mem   = membership(comm_result$walktrap)[idx],
    spinglass_same = length(unique(membership(comm_result$spinglass)[idx])) == 1,
    spinglass_mem  = membership(comm_result$spinglass)[idx]
  )
}

cluster_check_T0 <- check_cluster(comm_T0, diet_unhealthy_5, p_vars)
cluster_check_T1 <- check_cluster(comm_T1, diet_unhealthy_5, p_vars)

edge_density_signed <- function(W_sub) {
  diag(W_sub) <- NA
  mean(abs(W_sub), na.rm = TRUE)
}

perm_cluster_density <- function(net_obj, target_nodes, node_names, n_perm = 1000) {
  W_full <- net_obj$graph[node_names, node_names]
  obs <- edge_density_signed(W_full[target_nodes, target_nodes])
  k <- length(target_nodes)
  set.seed(2026)
  perm_dist <- replicate(n_perm, {
    rand_nodes <- sample(node_names, k)
    edge_density_signed(W_full[rand_nodes, rand_nodes])
  })
  list(observed = obs, perm_dist = perm_dist, p_value = mean(perm_dist >= obs))
}

perm_T0_full5 <- perm_cluster_density(net_T0, diet_unhealthy_5, p_vars)
perm_T1_full5 <- perm_cluster_density(net_T1, diet_unhealthy_5, p_vars)

perm_T0_core4 <- perm_cluster_density(net_T0, diet_core_4, p_vars)
perm_T1_core4 <- perm_cluster_density(net_T1, diet_core_4, p_vars)

group_density <- function(net_obj, group_vec) {
  W <- net_obj$graph; diag(W) <- NA
  c(KK = mean(abs(W[group_vec == "Knowledge", group_vec == "Knowledge"]), na.rm = TRUE),
    BB = mean(abs(W[group_vec == "Behavior",  group_vec == "Behavior"]),  na.rm = TRUE),
    KB = mean(abs(W[group_vec == "Knowledge", group_vec == "Behavior"]),  na.rm = TRUE))
}

perm_group_density <- function(net_obj, group_vec, n_perm = 1000) {
  obs <- group_density(net_obj, group_vec)
  set.seed(2026)
  perm_mat <- replicate(n_perm, group_density(net_obj, sample(group_vec)))
  raw_p <- c(KK_vs_BB = mean(abs(perm_mat["KK",] - perm_mat["BB",]) >= abs(obs["KK"] - obs["BB"])),
             KK_vs_KB = mean(abs(perm_mat["KK",] - perm_mat["KB",]) >= abs(obs["KK"] - obs["KB"])),
             BB_vs_KB = mean(abs(perm_mat["BB",] - perm_mat["KB",]) >= abs(obs["BB"] - obs["KB"])))
  list(observed = obs, perm_mat = perm_mat,
       p_values = raw_p, p_holm = p.adjust(raw_p, method = "holm"))
}

group_T0 <- perm_group_density(net_T0, comm_vec)
group_T1 <- perm_group_density(net_T1, comm_vec)


calc_add_one_p <- function(group_result, wave_label) {

  required_types <- c("KK", "BB", "KB")

  obs <- group_result$observed

  if (!all(required_types %in% names(obs))) {
    stop("group_result$observed must contain KK, BB, and KB.")
  }

  obs <- obs[required_types]

  perm_mat <- group_result$perm_mat

  if (
    !is.null(rownames(perm_mat)) &&
    all(required_types %in% rownames(perm_mat))
  ) {
    perm_values <- perm_mat[required_types, , drop = FALSE]

  } else if (
    !is.null(colnames(perm_mat)) &&
    all(required_types %in% colnames(perm_mat))
  ) {
    perm_values <- t(
      perm_mat[, required_types, drop = FALSE]
    )

  } else {
    stop(
      paste0(
        "Cannot identify KK, BB, and KB in perm_mat. ",
        "Please run rownames(group_result$perm_mat) and ",
        "colnames(group_result$perm_mat) to inspect its structure."
      )
    )
  }

  n_perm <- ncol(perm_values)

  observed_diff <- c(
    KK_vs_BB = abs(
      unname(obs["KK"]) - unname(obs["BB"])
    ),
    KK_vs_KB = abs(
      unname(obs["KK"]) - unname(obs["KB"])
    ),
    BB_vs_KB = abs(
      unname(obs["BB"]) - unname(obs["KB"])
    )
  )

  permuted_diff <- rbind(
    KK_vs_BB = abs(
      perm_values["KK", ] - perm_values["BB", ]
    ),
    KK_vs_KB = abs(
      perm_values["KK", ] - perm_values["KB", ]
    ),
    BB_vs_KB = abs(
      perm_values["BB", ] - perm_values["KB", ]
    )
  )

  extreme_count <- vapply(
    seq_along(observed_diff),
    function(i) {
      sum(
        permuted_diff[i, ] >= observed_diff[i],
        na.rm = TRUE
      )
    },
    integer(1)
  )

  names(extreme_count) <- names(observed_diff)

  raw_p_add1 <- (extreme_count + 1) / (n_perm + 1)

  adjusted_p_add1 <- p.adjust(
    raw_p_add1,
    method = "holm"
  )

  list(
    observed_means = data.frame(
      Wave = wave_label,
      Edge_type = required_types,
      Mean_absolute_edge_weight = as.numeric(obs),
      row.names = NULL
    ),

    tests = data.frame(
      Wave = wave_label,
      Comparison = names(observed_diff),
      Observed_difference = as.numeric(observed_diff),
      Extreme_count = as.integer(extreme_count),
      Total_permutations = n_perm,
      Raw_p_add1 = as.numeric(raw_p_add1),
      Holm_adjusted_p_add1 = as.numeric(adjusted_p_add1),
      row.names = NULL
    )
  )
}

edge_test_T0 <- calc_add_one_p(
  group_T0,
  "Baseline (T0)"
)

edge_test_T1 <- calc_add_one_p(
  group_T1,
  "Follow-up (T1)"
)

cat("\n============================================================\n")
cat("Observed mean absolute edge weights\n")
cat("============================================================\n")

print(
  rbind(
    edge_test_T0$observed_means,
    edge_test_T1$observed_means
  ),
  digits = 6,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("Add-one corrected pairwise permutation tests\n")
cat("============================================================\n")

print(
  rbind(
    edge_test_T0$tests,
    edge_test_T1$tests
  ),
  digits = 6,
  row.names = FALSE
)


community_results <- list(
  comm_T0 = comm_T0, comm_T1 = comm_T1,
  cluster_check_T0 = cluster_check_T0, cluster_check_T1 = cluster_check_T1,
  perm_T0_full5 = perm_T0_full5, perm_T1_full5 = perm_T1_full5,
  perm_T0_core4 = perm_T0_core4, perm_T1_core4 = perm_T1_core4,
  group_T0 = group_T0, group_T1 = group_T1
)
saveRDS(community_results, file.path(OUTPUT_DIR, "04b_community_results_N805.rds"))


cat("\n--- Case-dropping Bootstrap (1000 boots per network) ---\n")
cat("This takes significant computational time...\n")

set.seed(2026)

boot_T0 <- bootnet(net_T0, nBoots = 1000, type = "case", nCores = as.integer(Sys.getenv("CEBP_N_CORES", "12")))
boot_T1 <- bootnet(net_T1, nBoots = 1000, type = "case", nCores = as.integer(Sys.getenv("CEBP_N_CORES", "12")))

cs_T0 <- corStability(boot_T0)
cs_T1 <- corStability(boot_T1)
cat("\nCS-coefficient  T0:", round(cs_T0["strength"], 3),
    "  T1:", round(cs_T1["strength"], 3), "\n")

make_stability_plot <- function(boot_obj, wave_label, cs_val, col) {

  p <- plot(boot_obj, statistics = "strength") +
    geom_hline(yintercept = 0.7, linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    annotate("text", x = Inf, y = 0.72,
             label = "r = 0.70", hjust = 1.1, size = 3.2, color = "grey40") +
    annotate("text", x = 0.5, y = 0.05,
             label = paste0("CS = ", round(cs_val, 3)),
             hjust = 0.5, size = 3.5, fontface = "italic", color = col) +
    ggtitle(paste("Case-dropping Stability —", wave_label)) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12, hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey93")
    )

  return(p)
}

p_stab_T0 <- make_stability_plot(boot_T0, "Baseline (T0)", cs_T0["strength"], COL_K)
p_stab_T1 <- make_stability_plot(boot_T1, "Follow-up (T1)", cs_T1["strength"], COL_P)

p_stab_T0

p_stab_combined <- p_stab_T0 / p_stab_T1 +
  plot_annotation(
    caption = "Dashed line: r = 0.70 (correlation threshold used to calculate CS). CS > 0.50 = acceptable; > 0.75 = excellent.",
    theme   = theme(plot.caption = element_text(size = 9, color = "grey50"))
  )

ggsave(file.path(OUTPUT_DIR, "Network_Stability_Combined.png"),
       p_stab_combined, width = 9, height = 9, dpi = 180, bg = "white")
cat("Saved: Network_Stability_Combined.png saved\n")


cat("\n=== Summary ===\n")
cat(sprintf("T0  edges: %d  |  density: %.3f\n",
            sum(net_T0$graph[upper.tri(net_T0$graph)] != 0),
            sum(net_T0$graph[upper.tri(net_T0$graph)] != 0) / (n_nodes*(n_nodes-1)/2)))
cat(sprintf("T1  edges: %d  |  density: %.3f\n",
            sum(net_T1$graph[upper.tri(net_T1$graph)] != 0),
            sum(net_T1$graph[upper.tri(net_T1$graph)] != 0) / (n_nodes*(n_nodes-1)/2)))
cat(sprintf("NCT  M-stat p = %.3f  |  S-stat p = %.3f\n",
            nct_res$nwinv.pval, nct_res$glstrinv.pval))
cat("M-statistic (network invariance):", nct_res$nwinv.real, "\n")
cat("S-statistic (global strength invariance):", nct_res$glstrinv.real, "\n")
cat("Global strength T0 vs T1:", nct_res$glstrinv.sep, "\n")
cat(sprintf("CS   T0 = %.3f  |  T1 = %.3f\n",
            cs_T0["strength"], cs_T1["strength"]))

cat("\nOutput files (all PNG — RMD friendly):\n")
cat("  Network_Baseline_Followup.png\n")
cat("  Bridge_EI_Combined.png\n")
cat("  Network_Stability_Combined.png\n")
cat("  04_NCT_results_N805.rds\n")
cat("=== 04_network.R completed ===\n")
