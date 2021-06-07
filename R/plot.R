library(ggstatsplot)
library(tidyverse)
library(ggsignif)

# setwd("~/Desktop/iglR")

# global
IFOD2_control_global <- read_csv("matrices_excl/gretna/NetworkEfficiency/Group1/Eg.txt", col_names=F) %>%
  mutate(group = "IFOD2 c")

IFOD2_patho_global <- read_csv("matrices_excl/gretna/NetworkEfficiency/Group2/Eg.txt", col_names=F) %>%
  mutate(group = "IFOD2 p")

sd_control_global <- read_csv("matrices_sd/gretna/NetworkEfficiency/Group1/Eg.txt", col_names=F) %>%
  mutate(group = "SD_Stream c")

sd_patho_global <- read_csv("matrices_sd/gretna/NetworkEfficiency/Group2/Eg.txt", col_names=F) %>%
  mutate(group = "SD_Stream p")

all_global <- bind_rows(IFOD2_control_global, IFOD2_patho_global, sd_control_global, sd_patho_global)

all_global <- all_global %>%
  mutate(value = X1)

# local
IFOD2_control_local <- read_csv("matrices_excl/gretna/NetworkEfficiency/Group1/Eloc.txt", col_names=F) %>%
  mutate(group = "IFOD2 c")

IFOD2_patho_local <- read_csv("matrices_excl/gretna/NetworkEfficiency/Group2/Eloc.txt", col_names=F) %>%
  mutate(group = "IFOD2 p")

sd_control_local <- read_csv("matrices_sd/gretna/NetworkEfficiency/Group1/Eloc.txt", col_names=F) %>%
  mutate(group = "SD_Stream c")

sd_patho_local <- read_csv("matrices_sd/gretna/NetworkEfficiency/Group2/Eloc.txt", col_names=F) %>%
  mutate(group = "SD_Stream p")

all_local <- bind_rows(IFOD2_control_local, IFOD2_patho_local, sd_control_local, sd_patho_local)

all_local <- all_local %>%
  mutate(value = X1)

all_global %>%
  mutate(value = X1) %>%
  group_by(group) %>%
  summarize(m = mean(value))

all_local %>%
  mutate(value = X1) %>%
  group_by(group) %>%
  summarize(m = mean(value))

all_global

ggplot(all_global, aes(group, value)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(width=0.15, mapping=aes(color=group)) +
  geom_signif(
    comparisons = list(
      c("IFOD2 c", "IFOD2 p")
    ),
    y_position = 6600,
    test = t.test,
    map_signif_level = function(x) paste("p: ", signif(x, digits=2))
  ) +
  geom_signif(
    comparisons = list(
      c("SD_Stream c", "SD_Stream p")
    ),
    y_position = 6300,
    test = t.test,
    map_signif_level = function(x) paste("p: ", signif(x, digits=2))
  ) +
  geom_signif(
    comparisons = list(
      c("IFOD2 c", "SD_Stream c")
    ),
    y_position = 6900,
    test = t.test,
    map_signif_level = function(x) paste("p: ", signif(x, digits=2))
  ) +
  geom_signif(
    comparisons = list(
      c("IFOD2 p", "SD_Stream p")
    ),
    y_position = 7300,
    test = t.test,
    map_signif_level = function(x) paste("p: ", signif(x, digits=2))
  )
ggsave("global_eff.png")

ggplot(all_local, aes(group, value)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(width=0.15, mapping=aes(color=group)) +
  geom_signif(
    comparisons = list(
      c("IFOD2 c", "IFOD2 p")
    ),
    y_position = 6600,
    test = t.test,
    map_signif_level = function(x) paste("p: ", signif(x, digits=2))
  ) +
  geom_signif(
    comparisons = list(
      c("SD_Stream c", "SD_Stream p")
    ),
    y_position = 6300,
    test = t.test,
    map_signif_level = function(x) paste("p: ", signif(x, digits=2))
  ) +
  geom_signif(
    comparisons = list(
      c("IFOD2 c", "SD_Stream c")
    ),
    y_position = 6900,
    test = t.test,
    map_signif_level = function(x) paste("p: ", signif(x, digits=2))
  ) +
  geom_signif(
    comparisons = list(
      c("IFOD2 p", "SD_Stream p")
    ),
    y_position = 7300,
    test = t.test,
    map_signif_level = function(x) paste("p: ", signif(x, digits=2))
  )
ggsave("local_eff.png")


# ggwithinstats(
#   data = all_global,
#   x = group,
#   y = value,
#   type = "p",
#   comparison = list(
#     c("IFOD2 c","IFOD2 p")           
#   ),
#   pairwise = F,
#   plot.type = "box",
#   pairwise.display = "all",
#   p.adjust.method = "bonferroni",
#   title = "Global efficiency"
# )
# # ggsave("global_eff.png")
# 
# ggwithinstats(
#   data = all_local,
#   x = group,
#   y = value,
#   type = "p",
#   plot.type = "box",
#   pairwise.display = "all",
#   p.adjust.method = "bonferroni",
#   title = "Local efficiency"
# )
# # ggsave("local_eff.png")


