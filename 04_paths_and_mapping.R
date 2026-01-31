source("R/utils/helpers.R")

cfg <- read_config()
ensure_dir("data_intermediate")
ensure_dir("outputs/maps")

project_id <- cfg$project_id

df_in <- if (file.exists("data_intermediate/df_alltags_blocked.rds")) {
  readRDS("data_intermediate/df_alltags_blocked.rds")
} else {
  readRDS("data_intermediate/df_alltags_keep.rds")
}

df_paths <- get_path_by_run(df_in, project_id = project_id)
saveRDS(df_paths, "data_intermediate/df_paths.rds")

tag_keep <- cfg$tag_ids_keep
df_plot <- if (!is.null(tag_keep) && length(tag_keep) > 0) {
  dplyr::filter(df_paths, motusTagID %in% tag_keep)
} else df_paths

p1 <- ggplot(df_plot, aes(x = ts_h, y = recvLat)) +
  geom_point() + geom_path() +
  theme_bw() +
  labs(x = "Time (UTC)", y = "Receiver latitude", title = "Motus path (run-level)") +
  facet_wrap(~ motusTagID, scales = "free_x", ncol = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("outputs/maps/path_lat_time.png", p1, width = 11, height = 7, dpi = 300)

na_map <- map_data("world2") %>%
  dplyr::filter(region %in% c("Canada", "USA")) %>%
  dplyr::mutate(long = long - 360)

na_lakes <- map_data("lakes") %>% dplyr::mutate(long = long - 360)

xmin <- min(df_plot$recvLon, df_plot$tagDeployLon, na.rm = TRUE) - 2
xmax <- max(df_plot$recvLon, df_plot$tagDeployLon, na.rm = TRUE) + 2
ymin <- min(df_plot$recvLat, df_plot$tagDeployLat, na.rm = TRUE) - 2
ymax <- max(df_plot$recvLat, df_plot$tagDeployLat, na.rm = TRUE) + 2

p2 <- ggplot() +
  geom_polygon(data = na_map, aes(long, lat, group = group), colour = "grey60", fill = "grey98") +
  geom_polygon(data = na_lakes, aes(long, lat, group = group), colour = "grey70", fill = "white") +
  geom_path(data = df_plot, aes(recvLon, recvLat, group = as.factor(motusTagID), colour = as.factor(motusTagID))) +
  geom_point(data = df_plot, aes(tagDeployLon, tagDeployLat), colour = "black", shape = 4) +
  coord_map(projection = "mercator", xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  theme_bw() +
  labs(x = "", y = "", colour = "motusTagID", title = "Motus route (outline map)")

ggsave("outputs/maps/route_outline_map.png", p2, width = 11, height = 7, dpi = 300)

if (!is.null(cfg$google_api_key) && nzchar(cfg$google_api_key)) {
  ggmap::register_google(key = cfg$google_api_key)

  bbox <- c(left = xmin, right = xmax, bottom = ymin, top = ymax)
  gmap <- ggmap::get_stamenmap(bbox = bbox, maptype = "terrain-background", zoom = 6)

  p3 <- ggmap::ggmap(gmap) +
    geom_point(data = df_plot, aes(recvLon, recvLat), pch = 21, colour = "black", fill = "yellow") +
    geom_path(data = df_plot, aes(recvLon, recvLat, group = motusTagID, col = as.factor(motusTagID))) +
    theme_bw() +
    labs(colour = "motusTagID", title = "Motus route (ggmap)")

  ggsave("outputs/maps/route_ggmap.png", p3, width = 11, height = 7, dpi = 300)
}
