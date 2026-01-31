# Motus Data Analysis Pipeline (Project-based)

A reproducible, **data-free** R pipeline for downloading, filtering, summarizing, and mapping detections from a Motus project using the `motus` ecosystem.

This repository is designed for:
- updating a project `.motus` SQLite database
- removing dubious detections (Motus filter + project-specific blocking rules)
- building simplified movement paths from receiver runs
- creating maps and diagnostic plots
- computing transition distances/rates between receiver sites
- generating site and tag summaries

> **Data note:** This repo does **not** include Motus data. You will generate a local `.motus` file via `tagme()`.

---

## Project name suggestion

**MotusTrackR: Project-based Motus movement pipeline**

(Short, descriptive, and GitHub-friendly. You can rename the repo to whatever fits your paper/project.)

---

## Repository structure

- `scripts/`
  - `01_download_update_motus_db.R` — download/update `.motus` database for a project
  - `02_activity_filter_and_flatten.R` — apply `filterByActivity()`, flatten to a data.frame, keep “good” detections
  - `03_blocking_rules.R` — optional manual blocking rules (project-specific false detections)
  - `04_paths_and_mapping.R` — build paths (run-level summaries) and make maps/plots
  - `05_transitions_distance.R` — compute site-to-site transitions and movement distances
  - `06_site_tag_summaries.R` — site/tag summaries and helper plots
- `R/utils/` — reusable helper functions
- `config_example.yml` — copy to `config.yml` (ignored by git) and edit your paths/settings
- `outputs/` — figures and tables you generate (ignored by git)
- `data_intermediate/` — cached `.rds` objects (ignored by git)

---

## Quick start

### 1) Install packages

```r
install.packages(c(
  "tidyverse","lubridate","motus","motusData",
  "DBI","RSQLite","maps","rnaturalearth","rnaturalearthdata",
  "ggmap","yaml","plyr"
))
```

### 2) Create your config

Copy `config_example.yml` to `config.yml` and edit:
- `project_id`
- `motus_dir` (where the `.motus` database should be stored)
- `google_api_key` (optional; only if you use `ggmap`)
- `tag_ids_keep` (optional; tags to plot/compute transitions for)

### 3) Run scripts in order

```r
source("scripts/01_download_update_motus_db.R")
source("scripts/02_activity_filter_and_flatten.R")
source("scripts/03_blocking_rules.R")      # optional (if enabled in config.yml)
source("scripts/04_paths_and_mapping.R")
source("scripts/05_transitions_distance.R")
source("scripts/06_site_tag_summaries.R")
```

---

## License

MIT (see `LICENSE`).
