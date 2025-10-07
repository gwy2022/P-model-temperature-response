# ===================================================================
# load package
# ===================================================================
library(dplyr)
library(purrr)
library(tidyr)
library(smatr)
library(ggplot2)
library(ggnewscale)

# --- CHANGE: Loading the specified dataset ---
df <- read.csv("/Users/yongzhe/Documents/Reading/Topt_reading/ratio_vj.csv")
# ===================================================================
# manually set colour for each species 
# ===================================================================
species_colors <- c(
  "Cedrela odorata"    = "#b3004b",
  "Pinus nigra"         = "#f7a7ba",
  "Pinus pinaster"      = "#ffb300",
  "Pinus pinea"         = "#8c5c2c",
  "Pinus sylvestris"    = "#5c6643",
  "Quercus virginiana"  = "#227744",
  "Tamarindus indica"   = "#66bb66",
  "Ulmus americana"     = "#006400"
)

# ---- groupy by woody non woody pmodel and observation ----
df2_initial <- df %>%
  filter(Woody_Type %in% c("Woody", "Non-woody")) %>%
  mutate(
    DataType2 = case_when(
      tolower(DataType) == "obs"    ~ "Obs",
      tolower(DataType) == "pmodel" ~ "P model",
      TRUE                           ~ as.character(DataType)
    ),
    DataType2 = factor(DataType2, levels = c("Obs", "P model"))
  )

# --- set woody on left  ---
df2 <- df2_initial %>%
  mutate(
    Woody_Type = factor(Woody_Type, levels = c("Woody", "Non-woody"))
  )
df2$Species <- factor(df2$Species, levels = names(species_colors))


# ---- line color ----
line_cols <- c("Obs"="black", "P model"="red")
fill_cols <- c("Obs"="grey70", "P model"="red")

# ===================================================================
# model fit
# ===================================================================
fit_tbl <- df2 %>%
  group_by(DataType2, Woody_Type) %>%
  nest() %>%
  mutate(fit_data = map(data, ~{
    d <- .
    if (nrow(d) < 3 || n_distinct(d$T_growth) < 2) return(NULL)
    tryCatch({
      # --- CHANGE: The model formula now uses ratio_jv instead of topt ---
      fit <- sma(ratio_jv ~ T_growth, data = d, quiet = TRUE)
      
      coefs <- fit$coef[[1]]
      b0 <- coefs[1, "coef(SMA)"]; b1 <- coefs[2, "coef(SMA)"]
      r2 <- as.numeric(fit$r2); pv <- as.numeric(fit$p)
      sig <- pv < 0.05
      xr <- range(d$T_growth, na.rm = TRUE)
      xseq <- seq(xr[1], xr[2], length.out = 200)
      yseq <- b1 * xseq + b0
      lower_col <- grep("(?i)lower|lwr", colnames(coefs), value = TRUE)[1]
      upper_col <- grep("(?i)upper|upr", colnames(coefs), value = TRUE)[1]
      if (!is.na(lower_col) && !is.na(upper_col)) {
        b0_lo <- coefs[1, lower_col]; b0_hi <- coefs[1, upper_col]
        b1_lo <- coefs[2, lower_col]; b1_hi <- coefs[2, upper_col]
        y_low  <- b1_lo * xseq + b0_lo
        y_high <- b1_hi * xseq + b0_hi
      } else {
        y_low  <- rep(NA_real_, length(xseq)); y_high <- rep(NA_real_, length(xseq))
      }
      
      tibble(
        T_growth = xseq, y_fit = yseq, y_low = y_low, y_high = y_high,
        slope = b1, R2 = r2, pval = pv, signif = sig
      )
    }, error = function(e) NULL)
  })) %>%
  select(-data) %>%
  unnest(fit_data)



# ---- fitting for plots ----
fit_keep <- fit_tbl %>%
  filter(signif | (Woody_Type == "Non-woody" & DataType2 == "Obs"))

# ===================================================================
# rmse caculation regression fit
# ===================================================================
# observed
fit_params <- fit_tbl %>%
  group_by(DataType2, Woody_Type) %>%
  summarise(
    slope = first(slope), intercept = first(y_fit) - first(slope) * first(T_growth),
    .groups = "drop"
  )


rmse_own_fit_df <- df2 %>%
  left_join(fit_params, by = c("DataType2", "Woody_Type")) %>%
  filter(!is.na(slope)) %>%
  mutate(
    # --- CHANGE: Predicting ratio_jv and calculating residuals against ratio_jv ---
    ratio_jv_predicted = slope * T_growth + intercept,
    residual_sq = (ratio_jv - ratio_jv_predicted)^2
  ) %>%
  group_by(DataType2, Woody_Type) %>%
  summarise(rmse_own = sqrt(mean(residual_sq, na.rm = TRUE)), .groups = "drop")

# pmodel
pmodel_fits_df <- fit_tbl %>%
  filter(DataType2 == "P model") %>%
  group_by(Woody_Type) %>%
  summarise(
    pmodel_slope = first(slope),
    pmodel_intercept = first(y_fit) - first(slope) * first(T_growth),
    .groups = "drop"
  )

pmodel_ci_df_for_plot <- fit_tbl %>%
  filter(DataType2 == "P model") %>%
  mutate(DataType2 = "Obs")
# imposed
rmse_imposed_df <- df2 %>%
  filter(DataType2 == "Obs") %>%
  left_join(pmodel_fits_df, by = "Woody_Type") %>%
  mutate(
    # --- CHANGE: Predicting ratio_jv and calculating residuals against ratio_jv ---
    ratio_jv_predicted_by_pmodel = pmodel_slope * T_growth + pmodel_intercept,
    residual_sq = (ratio_jv - ratio_jv_predicted_by_pmodel)^2
  ) %>%
  group_by(DataType2, Woody_Type) %>%
  summarise(rmse_imposed = sqrt(mean(residual_sq, na.rm = TRUE)), .groups = "drop")

# ===================================================================
# creat label
# ===================================================================
n_points_df <- df2 %>% group_by(DataType2, Woody_Type) %>% summarise(n_points = n(), .groups = "drop")

label_df <- fit_keep %>%
  group_by(DataType2, Woody_Type) %>%
  summarise(
    slope = first(slope),
    R2 = first(R2),
    .groups = "drop"
  ) %>%
  left_join(n_points_df, by = c("DataType2", "Woody_Type")) %>%
  left_join(rmse_own_fit_df, by = c("DataType2", "Woody_Type")) %>%
  left_join(rmse_imposed_df, by = c("DataType2", "Woody_Type")) %>%
  mutate(
    label = case_when(
      DataType2 == "Obs" ~ sprintf(
        "Slope = %.2f\nR² = %.2f\nRMSE = %.2f\nRMSE_imp = %.2f",
        slope, R2, rmse_own, rmse_imposed
      ),
      DataType2 == "P model" ~ sprintf(
        "Slope = %.2f\nR² = %.2f\nRMSE_imp = %.2f",
        slope, R2, rmse_own
      )
    )
  )


# ===================================================================
# ---legend creation ---
# ===================================================================
legend_df_new <- tidyr::expand_grid(
  DataType2 = c("Obs", "P model"), Woody_Type = c("Woody", "Non-woody")
) %>%
  mutate(
    Woody_Type = factor(Woody_Type, levels = c("Woody", "Non-woody")),
    line_x_start = 15, line_x_end = 17, line_y = 3.8, # Example: Adjusted y value
    rect_x_start = 15, rect_x_end = 17, rect_y_center = 3.6, # Example: Adjusted y value
    rect_ymin = rect_y_center - 0.05, rect_ymax = rect_y_center + 0.05,
    text_x = 17.5, text_y_line = line_y, text_y_ci = rect_y_center
  )


# ===================================================================
# plot
# ===================================================================
df2 %>%
  mutate(Woody_Type = factor(trimws(Woody_Type), levels = c("Woody", "Non-woody"))) %>%
  # --- CHANGE: The y aesthetic is now ratio_jv ---
  ggplot(aes(x = T_growth, y = ratio_jv, color = Species)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = species_colors, guide = "none") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  # 1.  CI
  geom_ribbon(
    data = fit_keep %>% filter(!is.na(y_low), !is.na(y_high)),
    aes(x = T_growth, ymin = y_low, ymax = y_high, fill = DataType2, 
        group = interaction(DataType2, Woody_Type)),
    inherit.aes = FALSE, alpha = 0.18
  ) +
  # 2. Imposed CI
  geom_ribbon(
    data = pmodel_ci_df_for_plot,
    aes(x = T_growth, ymin = y_low, ymax = y_high, group = Woody_Type),
    fill = "lightblue", alpha = 0.3, inherit.aes = FALSE
  ) +
  # 3. 绘制原始的拟合线
  geom_line(
    data = fit_keep,
    aes(x = T_growth, y = y_fit, color = DataType2, linetype = signif,
        group = interaction(DataType2, Woody_Type)),
    inherit.aes = FALSE, linewidth = 1.2
  ) +
  # 4. Imposed Fit
  geom_abline(
    data = pmodel_fits_df %>% mutate(DataType2 = "Obs"),
    aes(slope = pmodel_slope, intercept = pmodel_intercept),
    color = "blue", linetype = "dotted", linewidth = 1.1
  ) +
  
  # --- CHANGE: In-plot legend commented out. Uncomment after adjusting y-coordinates in Step 3f ---
   # 5. legend creation
   geom_rect(
     data = legend_df_new %>% filter(Woody_Type == "Woody"),
     aes(xmin = rect_x_start, xmax = rect_x_end, ymin = rect_ymin, ymax = rect_ymax, fill = DataType2),
     inherit.aes = FALSE, alpha = 0.18
   ) +
   geom_segment(
     data = legend_df_new %>% filter(Woody_Type == "Woody"),
     aes(x = line_x_start, xend = line_x_end, y = line_y, yend = line_y, color = DataType2),
    inherit.aes = FALSE, linewidth = 1.2
  ) +
  geom_text(
    data = legend_df_new %>% filter(Woody_Type == "Woody"),
    aes(x = text_x, y = text_y_line, label = DataType2),
    inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 3
  )  +

scale_fill_manual(values = fill_cols, guide = "none") +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
  scale_color_manual(values = line_cols, guide = "none") +
  
  # --- CHANGE: Use robust, relative positioning for annotation text ---
  geom_text(
    data = label_df,
    aes(x = 35, y = Inf, label = label),
    hjust = 1, vjust = 1.1,  
    size = 3.5,
    lineheight = 0.9,
    inherit.aes = FALSE
  )+
  coord_cartesian(xlim = c(14, 36),ylim = c(0,3.5), expand = FALSE) +
  
  facet_grid(DataType2 ~ factor(Woody_Type, levels = c("Woody", "Non-woody"))) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black")
  ) +
  # --- CHANGE: Update the y-axis label ---
  labs(
    x = expression("T"["growth"] ~ "(°C)"),
    y = expression(V[cmax25]/J[max25] ~ Ratio) # A plausible label for ratio_jv
  )

# 只取 Woody 的数据
df2_woody <- df2 %>% filter(Woody_Type == "Woody")

# 重新获取 Woody 的物种颜色
species_woody <- unique(df2_woody$Species)
species_colors_woody <- species_colors[species_woody]

# 创建临时图，用于提取 legend
p_legend_woody <- ggplot(df2_woody, aes(x = T_growth, y = ratio_jv, color = Species)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = species_colors_woody, name = "Woody") +
  theme_void() +
  theme(legend.position = "right")

# 提取 legend
species_legend_woody <- cowplot::get_legend(p_legend_woody)

# 单独显示 legend
grid::grid.newpage()
grid::grid.draw(species_legend_woody)

