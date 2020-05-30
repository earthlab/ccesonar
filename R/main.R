library(tidyverse)
library(tidymodels)
library(sf)
library(here)
library(ranger)
library(parallel)

n_threads <- parallel::detectCores()


# Get and read the data ---------------------------------------------------

path_to_data <- here("analysis", "data", "sonar.csv")
if (!file.exists(path_to_data)) {
  download.file("https://ndownloader.figshare.com/files/22837592",
                destfile = path_to_data)
}

sonar_data <- read_csv(here::here("analysis", "data", "sonar.csv"))



# Split into spatial bins (K fold CV) and final year test set -------------

n_folds <- 10
data_to_split <- sonar_data %>%
  filter(year != max(year)) %>%
  mutate(lat_bin = as.character(cut(Lat,
                                    breaks = n_folds)))
d_split <- data_to_split %>%
  group_vfold_cv(group = "lat_bin")

get_id_left_out <- function(x) unique(assessment(x)$lat_bin)

heldout_bins <- map(d_split$splits, get_id_left_out) %>% unlist

crossval_pts <- data_to_split %>%
  distinct(Lon, Lat, lat_bin, year) %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326)

test_pts <- sonar_data %>%
  filter(year == 2015) %>%
  distinct(Lon, Lat, year) %>%
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326)


# Generate a plot of the study region -------------------------------------

basemap <- st_read("https://raw.githubusercontent.com/codeforamerica/click_that_hood/master/public/data/north-america.geojson") %>%
  summarize

ptsize <- .05
map_plot <- crossval_pts %>%
  ggplot() +
  geom_sf(data = basemap, color = NA, fill = "grey99") +
  geom_sf(aes(color = fct_rev(lat_bin), fill = fct_rev(lat_bin)), size = ptsize) +
  facet_wrap(~year, nrow = 1) +
  xlab("") +
  ylab("") +
  scale_color_discrete("CV fold") +
  scale_fill_discrete("CV fold") +
  theme(legend.position = "left") +
  geom_sf(data = test_pts, size = ptsize, alpha = .5) +
  geom_text(data = tibble(year = 2015),
            x = -Inf, y = -Inf, label = "2015 withheld as test data",
            inherit.aes = FALSE, vjust = 0, hjust = 0, size = 3) +
  coord_sf(crs = 4326, xlim = range(sonar_data$Lon),
           ylim = range(sonar_data$Lat)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
map_plot
ggsave(here("analysis", "figures", "map_plot.pdf"),
       map_plot, width = 9, height = 4)



# Define a function to do K fold cros validation ---------------
# syntax adapted from:
# https://www.benjaminsorensen.me/post/modeling-with-parsnip-and-tidymodels/

fit_kfold <- function(formula, folds) {
  rec <- folds$splits[[1]] %>%
    training() %>%
    recipe(formula)

  rf_mod <- rand_forest(
    mode = "regression",
    trees = 100
  ) %>%
    set_engine("ranger", num.threads = n_threads)

  folded <- folds %>%
    mutate(recipes = splits %>%
             # Prepper is a wrapper for `prep()` which handles `split` objects
             map(prepper, recipe = rec)) %>%
    mutate(test_data = splits %>% map(analysis)) %>%
    mutate(
      rf_fits =
        map2(
          recipes,
          test_data,
          ~ fit(
            rf_mod,
            formula(.x),
            data = bake(object = .x, new_data = .y)
          )
        )
    )


  predict_rf <- function(split, rec, model) {
    test <- bake(rec, assessment(split))
    predict(model, test) %>%
      bind_cols(tibble(obs = test$log_nasc))
  }

  predictions <-
    folded %>%
    mutate(
      pred =
        list(
          splits,
          recipes,
          rf_fits
        ) %>%
        pmap(predict_rf)
    )

  ## Evaluate
  eval <- predictions %>%
    mutate(
      metrics = pred %>% map(~ metrics(., truth = obs, estimate = .pred))
    ) %>%
    select(metrics) %>%
    mutate(fold = 1:n()) %>%
    unnest(metrics) %>%
    spread(.metric, .estimate) %>%
    mutate(lat_bin = heldout_bins[fold]) %>%
    arrange(lat_bin)

  pred_df <- predictions %>%
    mutate(id = parse_number(id),
           lat_bin = heldout_bins[id]) %>%
    select(pred, lat_bin) %>%
    unnest(cols = c(pred))

  rsq_df <- eval %>%
    rowwise %>%
    mutate(lat_bin = fct_rev(factor(paste("Lat:", lat_bin))),
           label = as.character(as.expression(substitute(italic(R)^2~"="~rsq,
                                                         list(rsq = format(rsq, digits = 2)))))) %>%
    ungroup

  return(list(eval = eval,
              pred = predictions,
              folded = folded,
              formula = formula,
              pred_df = pred_df,
              rsq_df = rsq_df))
}


# Fit models --------------------------------------------------------------

full_formula <- log_nasc ~ d7_chlor_a + d7_nflh + d7_poc +
  d7_sst + Day + distance_from_shore_m + Lat  +
  depth_m + wavelength + depth_bin

reduced_formula <- log_nasc ~ Day + distance_from_shore_m +
  Lat  + depth_m + wavelength + depth_bin


full_model <- fit_kfold(full_formula,
                        folds = d_split)
reduced_model <- fit_kfold(reduced_formula,
                           folds = d_split)


# Evaluate models ---------------------------------------------------------

rsquared <- function(y, yhat) {
  1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
}

full_rsq <- full_model$pred_df %>%
  rename(log_nasc = obs) %>%
  left_join(data_to_split %>%
              select(log_nasc, lat_bin, depth_bin, wavelength)) %>%
  distinct() %>%
  mutate(model = "Full")

reduced_rsq <- reduced_model$pred_df %>%
  rename(log_nasc = obs) %>%
  left_join(data_to_split %>%
              select(log_nasc, lat_bin, depth_bin, wavelength)) %>%
  mutate(model = "Reduced")

joined_rsq <- full_join(full_rsq, reduced_rsq)

rsq_split <- joined_rsq %>%
  group_by(depth_bin, wavelength, lat_bin, model) %>%
  summarize(rsq = rsquared(log_nasc, .pred)) %>%
  ungroup %>%
  mutate(model = factor(model, levels = c("Full", "Reduced")),
         rsq = ifelse(rsq < 0, 0, rsq))

rsq_split_plot <- rsq_split %>%
  mutate(depth_bin = factor(depth_bin, levels = c("Shallowwater",
                                                  "Midwater",
                                                  "Shelfwater",
                                                  "Deepwater")),
         wavelength_f = paste(wavelength, "kHz"),
         wavelength_f = fct_reorder(wavelength_f, wavelength)) %>%
  ggplot(aes(x = rsq, y = lat_bin, color = model)) +
  facet_grid(depth_bin~wavelength_f) +
  geom_point() +
  ylab("Latitude bin") +
  scale_color_manual("Model", values = c("black", "dodgerblue")) +
  theme(legend.position = "bottom") +
  xlab(expression(paste("Holdout ", R^2)))
rsq_split_plot
ggsave(here("analysis", "figures", "baseline-model-split.pdf"),
       plot = rsq_split_plot, width = 7, height = 7)

# print R-squared summary stats for each model
rsq_split %>%
  group_by(model) %>%
  summarize(mean_rsq = mean(rsq),
            min_rsq = min(rsq),
            max_rsq = max(rsq))




# Train a model with all but the 2015 data --------------------------------

training_data <- sonar_data %>%
  filter(year != max(year))
test_data <- sonar_data %>%
  filter(year == max(year))

full_test2015 <- ranger(full_formula,
                        num.threads = n_threads,
                        data = training_data)
full_pred2015 <- predict(full_test2015, data = test_data)

reduced_test2015 <- ranger(reduced_formula,
                        num.threads = n_threads,
                        data = training_data)
reduced_pred2015 <- predict(reduced_test2015, data = test_data)

test_data <- test_data %>%
  mutate(yhat_full = full_pred2015$predictions,
         yhat_reduced = reduced_pred2015$predictions,
         depth_bin = factor(depth_bin, levels = c("Shallowwater",
                                                  "Midwater",
                                                  "Shelfwater",
                                                  "Deepwater")),
         wavelength_f = paste(wavelength, "kHz"),
         wavelength_f = fct_reorder(wavelength_f, wavelength))

test_rsq <- test_data %>%
  group_by(depth_bin, wavelength) %>%
  summarize(full = rsquared(log_nasc, yhat_full),
            reduced = rsquared(log_nasc, yhat_reduced)) %>%
  ungroup %>%
  pivot_longer(cols = c("full", "reduced")) %>%
  mutate(Model = ifelse(name == "full", "Full", "Reduced"),
         Model = factor(Model, levels = c("Full", "Reduced")),
         value = ifelse(value < 0, 0, value))


# Evaluate performance on the 2015 test set ----------------------------
test_data %>%
  summarize(full = rsquared(log_nasc, yhat_full),
            reduced = rsquared(log_nasc, yhat_reduced)) %>%
  ungroup %>%
  pivot_longer(cols = c("full", "reduced")) %>%
  mutate(Model = ifelse(name == "full", "Full", "Reduced"),
         Model = factor(Model, levels = c("Full", "Reduced")),
         value = ifelse(value < 0, 0, value))
range(test_rsq$value)


test_cor <- test_data %>%
  group_by(depth_bin, wavelength_f) %>%
  summarize(full = cor(log_nasc, yhat_full),
            reduced = cor(log_nasc, yhat_reduced)) %>%
  ungroup %>%
  pivot_longer(cols = c("full", "reduced")) %>%
  mutate(Model = ifelse(name == "full", "Full", "Reduced"),
         Model = factor(Model, levels = c("Full", "Reduced")))
test_plot <- test_cor %>%
  ggplot(aes(x = value, y = fct_rev(depth_bin), color = Model)) +
  facet_wrap(~wavelength_f, nrow = 1) +
  geom_point() +
  ylab("Depth bin") +
  xlab(expression(paste("Holdout correlation"))) +
  scale_color_manual("Model", values = c("black", "dodgerblue")) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype = "dashed")
test_plot
ggsave(here("analysis", "figures", "test-plot.pdf"), plot = test_plot,
       width = 8, height = 2.5)


mean_df <- test_data %>%
  group_by(depth_bin, wavelength_f) %>%
  summarize(true_mean = mean(log_nasc),
            pred_mean = mean(yhat_full))

test_scatter <- test_data %>%
  ggplot(aes(x = log_nasc, y = yhat_full)) +
  stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon") +
  scale_fill_viridis_c() +
  facet_grid(depth_bin ~ wavelength_f) +
  geom_abline(linetype = "dashed") +
  geom_point(aes(x = true_mean, y = pred_mean),
            data = mean_df, color = "red") +
  theme(legend.position = "none") +
  xlab("True value") +
  ylab("Predicted value") +
  geom_text(aes(x = -Inf, y = Inf, label = paste0("r=", round(value, 2))),
            data = filter(test_cor, name == "full"),
            hjust = 0, vjust = 1, size = 3.5)
test_scatter
ggsave(here("analysis", "figures", "test-scatter.pdf"),
       plot = test_scatter, width = 7, height = 4.5)

