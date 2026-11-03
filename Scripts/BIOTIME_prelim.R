
# install.packages("BioTIMEr")
# citation(package = 'BioTIMEr')

library(BioTIMEr)
library(dplyr)
library(ggplot2)


# you can run the following commands to retrieve more information about the subsets.
?BTsubset_meta
?BTsubset_data

# Link to the location describing the columns (https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fgeb.70003&file=geb70003-sup-0001-FileS1.pdf)

BTsubset_meta <- BTsubset_meta
BTsubset_data <- BTsubset_data

# Studies are defined by a sampling methodology which is constant over time, and have a minimum of one site and two sampling events in different years, with at least one sampling location each


#####################################################
# 1. Gridding ----------------------------------------------------------------
####################################################

# Because the spatial extent varies among studies, an advisable step is to grid those studies that have large extents and multiple sampling locations into smaller equal-sized subsets.

# ***By default, the function separates studies into hexagonal cells of ~96 km2 (res=12), as this resolution was found to be the most appropriate when working on the whole BioTIME database.***



grid_samples <- gridding(
  meta = BTsubset_meta,
  btf = BTsubset_data,
  res = 12,
  resByData = FALSE
)

# Yet sometimes, particularly when working with a subset of BioTIME, you may want to explore alternative grid resolutions. The function has two options: you can either directly define a new resolution (see vignette(‘dggridR’)) or allow the function to find the best resolution based on the average study extent of your dataset. Let’s run over these alternative options with a few examples:

# A good rule of thumb is to never set res higher than the resolution you get when resByData = TRUE.

# check_res_2 <- grid_samples_14 |>
#   summarise(
#     n_cell = n_distinct(cell),
#     n_aID = n_distinct(assemblageID),
#     res = "res14",
#     .by = c(StudyMethod, STUDY_ID)
#   )
# 
# check_res_3 <- grid_samples_auto |>
#   summarise(
#     n_cell = n_distinct(cell),
#     n_aID = n_distinct(assemblageID),
#     res = "res15",
#     .by = c(StudyMethod, STUDY_ID)
#   )
# 
# checks <- rbind(check_res_1, check_res_2, check_res_3)
# 
# ggplot(data = checks) +
#   geom_bar(
#     mapping = aes(x = as.character(STUDY_ID), y = n_aID, fill = res),
#     stat = "identity",
#     position = "dodge"
#   ) +
#   scale_fill_discrete(type = c("#155f49", "#66c1d1", "#d9d956")) +
#   xlab("StudyID") +
#   ylab("Number of assemblages in a study") +
#   theme_bw()



# You can see that each sample is assigned a different combination of study ID and grid cell (based on its latitude and longitude) resulting in a unique identifier that determines our assemblage time series (assemblageID). This assemblageID can now be used as the study unit for analysis of biodiversity change based on similar spatial extent between them.

# Let’s take a look and see how many studies/cells/assemblages time series are there in this subset:

check_res_1 <- grid_samples %>% 
  summarise(
    n_cell = n_distinct(cell),
    n_aID = n_distinct(assemblageID),
    res = "res12",
    .by = c(STUDY_ID, StudyMethod)
  )

# How many samples were there in each study?
ggplot(data = check_res_1) +
  geom_bar(
    mapping = aes(x = as.character(STUDY_ID), y = n_aID, fill = res),
    stat = "identity"
  ) +
  scale_fill_discrete(type = c("#155f49")) +
  xlab("StudyID") +
  ylab("Number of assemblages in a study") +
  theme_bw() +
  theme(legend.position = "none")

# So for study 431 there are over 400 unique study units (unique spatial locations in different grid cells that fall under the one study)


#####################################################
# 2. Rarefaction ----------------------------------------------------------------
####################################################

# Now we will deal with another common issue when working with time series data: the fact that the number of samples may vary over time. When evaluating changes through time it’s important to be able to quantify change independently of the number of samples, i.e. we want to be able to assess whether e.g. alpha diversity has truly increased or decreased over time and that such signal is not a sampling artefact due to uneven sampling effort – i.e. resulting from having more or fewer samples in different years or periods of the time series being compared.

# Probably the most common and widely used approach to prevent temporal variation in sampling effort from affecting diversity estimate is to apply sample-based rarefaction. Here (add link) is a great introduction to it. BioTIMEr provides the resampling() function that implements sample-based rarefaction, already adapted to BioTIME data and designed to run smoothly over the previously gridded dataset as done above. The function identifies the minimum number of samples across all years in a time series and then uses this minimum to randomly resample each year down to that number.

# Translation - in one spatial grid cell for a given study, the number of samples (e.g. plots) for a given year might vary 

# QUESTION - must then lose a lot of data, so you'd want to possibly exclude certain low-sampling studies???
# ANSWER: they mention later that removing really low sampling events might be okay to do


# First, if you are not sure you need this step,
# you can always check how many samples there are in every year of the different time series:
check_samples <- grid_samples %>% 
  summarise(
    n_samples = n_distinct(SAMPLE_DESC),
    n_species = n_distinct(Species),
    .by = c(STUDY_ID, assemblageID, YEAR)
  )

head(check_samples)

# Let's apply resampling() then, using the data frame of the gridded data:
grid_rare <- resampling(
  x = grid_samples,
  measure = "ABUNDANCE",
  resamps = 1,
  conservative = FALSE
)


# The function returns a single long-form data frame containing the total currency of interest (i.e. summed across all samples) for each species in each year within each assemblage time series. You can specify what the currency field is: Abundance, Biomass or both (BioTIME studies can include records for either or both). Note that resampling() also tests and removes any NAs within the selected currency field(s) before implementing the sample-based rarefaction

# Remember that sample-based rarefaction involves randomly selecting a specified number of samples (i.e. sampling events) for each year, this number being equal to the number of samples in the year with the lowest sampling effort in each time series. If the sampling effort through time was equal or relatively stable, then you probably do not need to concern yourself any further. However, for many time series in BioTIME (or indeed other biodiversity data) sampling effort can vary substantially, and thus the sample-based rarefaction process can lead to different species diversity and composition simply by randomly selecting different samples. 

# SO AFTER THIS THE ABUNDANCE IS TOTALLED ACROSS SAMPLES

# One option to deal with this issue is to repeat this random sampling process multiple times, effectively creating multiple alternative datasets for the same time series (and which can be used in subsequent sensitivity analysis or to estimate the variability in diversity metrics that originates from the random sample selection process):

# Let's implement the sample-based rarefaction by resampling the dataset 10 times.
grid_rare_n10 <- resampling(
  x = grid_samples,
  measure = "ABUNDANCE",
  resamps = 10,
  conservative = FALSE
)

# Note that you may want to resample many more times (e.g. at least 30-100+ times, but up to 199
# if e.g. working with the whole BioTIME data), depending on how many iterations you want a
# subsequent bootstrap analysis to have.
# This may also take some computation time, so if you are working with a big subset of BioTIME
# is advisable to break it down in smaller subsets.

# Each resampling iteration will be identified as resamp = 1:n, in this case 1:10.
# Now we can check if there are differences across the first 3 of these iterations:

check_resamps <- grid_rare_n10[grid_rare_n10$resamp < 4, ] |>
  summarise(
    n_obs = n(),
    n_species = n_distinct(Species),
    n_year = n_distinct(YEAR),
    .by = c(STUDY_ID, assemblageID, resamp)
  )


# DATA ARE NOW READY FOR ANALYSIS -----------------------------------------



################################################
# PART 2 - Download entire database ---------------------------------------
################################################

# Downloaded on 2nd November 2026
# BioTIME citation (Dornelas et al., 2024):
# Maria Dornelas, Laura H. Antão, Amanda E. Bates, Viviana Brambilla, Jonathan M. Chase, Cher F. Y. Chow, Ada Fontrodona-Eslava, Anne E. Magurran, Inês S. Martins, Faye Moyes, Alban Sagouis, et al. BioTIME 2.0: Expanding and Improving a Database of Biodiversity Time Series. Global Ecol Biogeogr. 2025; 34(5):e70003. https://doi.org/10.1111/geb.70003

library(here)

# Only run this section once:

# dir.create(here("Data"), showWarnings = FALSE)
# 
# curl::curl_download("https://biotime.st-andrews.ac.uk/dl_request.php?dl=full_csv",
#                     here("Data", "biotime_full.zip"), quiet = FALSE)
# 
# unzip(here("Data", "biotime_full.zip"), exdir = here("Data"))

BT_data <- read.csv(here("Data/biotime_v2_rawdata_2025.csv"))
BT_meta <- read.csv(here("Data/biotime_v2_metadata_2025.csv"))

# DATA POINTS IS NUMBER OF YEARS COVERED BY THE STUDY
# NUMBER OF SAMPLES IS NUMBER OF UNIQUE SAMPLING EVENTS


#####################################################
# PART 2 - Clean the full BioTIME 2.0 database -------------------------------
#####################################################
# Pipeline (ordered so the biggest reductions happen first, keeping memory low):
#   1. Terrestrial realm                      (metadata)
#   2. Abundance or abundance + biomass only  (metadata, then record level)
#   3. Drop studies that can never reach 5 years (cheap pre-filter, study level)
#   4. Gridding (BioTIMEr::gridding)          -> assemblage time series
#   5. Look at studies / cells / assemblages per study
#   6. Keep assemblage time series with >= 5 sampled years (Antao et al. 2020)
#   7. Evaluation plots + summary tables
#
# Extra packages needed:  install.packages(c("maps"))   (world map outlines)

library(BioTIMEr)
library(dplyr)
library(ggplot2)
library(here)
library(maps)

dir.create(here("Data"),    showWarnings = FALSE)
dir.create(here("Plots/BioTIME_prelim"), showWarnings = FALSE)

min_years <- 5    # minimum number of distinct sampled years per time series
grid_res  <- 12   # ~96 km2 hexagons, as recommended for the whole database


#####################################################
# 0. Load data -----------------------------------------------------------------
#####################################################
# fread() is much faster than read.csv() on a file this size
BT_data <- data.table::fread(here("Data/biotime_v2_rawdata_2025.csv"),  data.table = FALSE)
BT_meta <- data.table::fread(here("Data/biotime_v2_metadata_2025.csv"), data.table = FALSE)

# Quick look at the categories we are about to filter on
table(BT_meta$REALM,  useNA = "ifany")
# Freshwater = 151 studies, Marine = 229 studies, Terrestrial = 328 studies
table(BT_meta$AB_BIO, useNA = "ifany")
# Abundance = 495 studies, AB = 143 studies, Biomass = 68 studies

# Running log of how much data each step removes
attrition <- list()
log_step <- function(step, studies, assemblages = NA_integer_, records = NA_integer_) {
  attrition[[length(attrition) + 1]] <<- data.frame(
    step = step, n_studies = studies, n_assemblages = assemblages, n_records = records
  )
}
log_step("All BioTIME 2.0", n_distinct(BT_meta$STUDY_ID), records = nrow(BT_data))


#####################################################
# 1-2. Terrestrial realm + abundance data -------------------------------------
#####################################################
# Done on the (small) metadata table first, so the big raw table is only
# subset once.
meta_sel <- BT_meta %>%
  filter(tolower(trimws(REALM)) == "terrestrial")
log_step("Terrestrial realm", nrow(meta_sel))

# AB_BIO: "A" = abundance only, "B" = biomass only, "AB" = both
meta_sel <- meta_sel %>%
  filter(trimws(AB_BIO) %in% c("A", "AB"))
log_step("Abundance (A) or abundance + biomass (AB) studies", nrow(meta_sel))

dat <- BT_data %>%
  filter(STUDY_ID %in% meta_sel$STUDY_ID)
log_step("  ...records in those studies", n_distinct(dat$STUDY_ID), records = nrow(dat))

# Free memory - the full table is no longer needed
rm(BT_data); invisible(gc())

# Record-level tidy-up. gridding() asserts that ABUNDANCE and BIOMASS are >= 0,
# so negative values (which are missing-value codes) must go.
dat$ABUNDANCE <- as.numeric(dat$ABUNDANCE)
dat$BIOMASS   <- as.numeric(dat$BIOMASS)
dat$BIOMASS[!is.na(dat$BIOMASS) & dat$BIOMASS < 0] <- NA   # biomass not used downstream

dat <- dat %>%
  filter(
    !is.na(ABUNDANCE), ABUNDANCE >= 0,     # a real abundance value (0 kept = recorded absence)
    !is.na(YEAR),
    !is.na(LATITUDE), !is.na(LONGITUDE)
  )
log_step("Records with valid abundance, year and coordinates",
         n_distinct(dat$STUDY_ID), records = nrow(dat))


#####################################################
# 3. Pre-filter: studies spanning < 5 sampled years ----------------------------
#####################################################
# An assemblage can only ever have a subset of its study's years, so any study
# with < 5 distinct years cannot contain a qualifying time series. Dropping
# them now makes gridding faster. The real 5-year filter is applied per
# assemblage in step 6, after gridding.
study_years <- dat %>%
  summarise(n_years = n_distinct(YEAR), .by = STUDY_ID)

dat <- dat %>%
  filter(STUDY_ID %in% study_years$STUDY_ID[study_years$n_years >= min_years])
log_step(paste0("Studies with >= ", min_years, " sampled years (pre-filter)"),
         n_distinct(dat$STUDY_ID), records = nrow(dat))


## Filtering out studies with < 5 years removes 100 studies out of 282

#####################################################
# 4. Gridding ---------------------------------------------------------------
#####################################################
# NOTE: gridding() classifies each study as single-location (SL) or
# multi-location (ML) using a threshold calculated from the metadata passed
# in (mean + SD of the area of single-point studies). I pass the FULL metadata
# (all realms) so that threshold matches the whole-database setting for which
# res = 12 was tuned. The inner join inside gridding() means only the studies
# present in `dat` come back out, so the filtering above is respected.
grid_samples <- gridding(
  meta      = BT_meta,
  btf       = dat,
  res       = grid_res,
  resByData = FALSE
)
rm(dat); invisible(gc())

log_step("Gridded (all assemblages)",
         n_distinct(grid_samples$STUDY_ID),
         assemblages = n_distinct(grid_samples$assemblageID),
         records = nrow(grid_samples))


#####################################################
# 5. Studies / cells / assemblages per study ----------------------------------
#####################################################
# (Same idea as check_res_1 in the tutorial, but on the full terrestrial data.)
# Within a study, one cell = one assemblage, so n_cell and n_aID are equal.
#  n_aID is the n_distinct(assemblageID)
check_res_full <- grid_samples %>%
  summarise(
    n_cell = n_distinct(cell),
    n_aID  = n_distinct(assemblageID),
    .by = c(STUDY_ID, StudyMethod, TAXA)
  ) %>%
  arrange(desc(n_aID))


cat("\nStudies:", n_distinct(grid_samples$STUDY_ID),
    "| Unique grid cells:", n_distinct(grid_samples$cell),
    "| Assemblage time series:", n_distinct(grid_samples$assemblageID), "\n")
print(table(check_res_full$StudyMethod)) # StudyMethod is: ML = multi-location, SL = single-location
print(head(check_res_full, 15))
print(summary(check_res_full$n_aID))

# Studies: 181 | Unique grid cells: 6843 | Assemblage time series: 7911 
# StudyMethod: ML = 44, SL = 137
# For number of distinct assemblages per study: median is one, max is 3075

# So unique grid cells are 6843, but there's more assemblage time series because there are assemblages that are from different studies but found in the same grid cell

#####################################################
# 6. Assemblage-level summaries + >= 5 year filter --------------------------
#####################################################
# Samples per assemblage per year (needed for sampling-effort diagnostics)
samples_per_year <- grid_samples %>%
  summarise(n_samples = n_distinct(SAMPLE_DESC),
            .by = c(STUDY_ID, assemblageID, YEAR))

# Year-level sampling effort and gaps for each assemblage
effort <- samples_per_year %>%
  arrange(assemblageID, YEAR) %>%
  summarise(
    min_samples    = min(n_samples),
    median_samples = median(n_samples),
    max_samples    = max(n_samples),
    total_samples  = sum(n_samples),
    # proportion of samples that would survive sample-based rarefaction, which
    # cuts every year down to the minimum number of samples in that series
    prop_retained  = min(n_samples) * n() / sum(n_samples),
    # longest run of missing years between two sampled years (0 = no gaps)
    max_gap        = if (n() > 1) max(diff(YEAR)) - 1 else 0,
    .by = assemblageID
  )

# Optional extra descriptors from the metadata (only those that exist)
meta_extra <- BT_meta %>%
  select(STUDY_ID, any_of(c("HABITAT", "BIOME_MAP", "ORGANISMS")))

assemblage_summary <- grid_samples %>%
  summarise(
    n_years    = n_distinct(YEAR),               # number of years SAMPLED
    start_year = min(YEAR),
    end_year   = max(YEAR),
    n_species  = n_distinct(Species),
    n_records  = n(),
    mean_lat   = mean(LATITUDE),
    mean_lon   = mean(LONGITUDE),
    .by = c(assemblageID, STUDY_ID, StudyMethod, TAXA, CLIMATE)
  ) %>%
  mutate(duration = end_year - start_year + 1) %>%  # number of years the series SPANS
  left_join(effort,     by = "assemblageID") %>%
  left_join(meta_extra, by = "STUDY_ID") %>%
  mutate(keep = n_years >= min_years)

cat("\nAssemblages with >=", min_years, "years:",
    sum(assemblage_summary$keep), "of", nrow(assemblage_summary), "\n")

# Assemblages with >= 5 years: 4262 of 7911 

# ---- apply the filter ----
ts5 <- assemblage_summary %>% filter(keep)

grid_samples_5yr <- grid_samples %>%
  filter(assemblageID %in% ts5$assemblageID)

log_step(paste0("Assemblages with >= ", min_years, " sampled years (FINAL)"),
         n_distinct(grid_samples_5yr$STUDY_ID),
         assemblages = n_distinct(grid_samples_5yr$assemblageID),
         records = nrow(grid_samples_5yr))

attrition_tbl <- bind_rows(attrition)
print(attrition_tbl)

# ---- save outputs ----
write.csv(grid_samples_5yr, here("Data", "BT_terrestrial_abund_gridded_5yr.csv"))
write.csv(assemblage_summary, here("Data", "BT_assemblage_summary.csv"), row.names = FALSE)
write.csv(attrition_tbl,      here("Data", "BT_cleaning_attrition.csv"), row.names = FALSE)


#####################################################
# 7. Evaluation plots ----------------------------------------------------------
#####################################################
theme_set(theme_bw(base_size = 11))

save_fig <- function(p, name, w = 9, h = 5) {
  ggsave(here("Plots/BioTIME_prelim", paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  print(p)
}

# Taxa palette used consistently across plots
taxa_levels <- ts5 %>% count(TAXA, sort = TRUE) %>% pull(TAXA)
ts5 <- ts5 %>% mutate(TAXA = factor(TAXA, levels = taxa_levels))
taxa_scale_fill   <- scale_fill_viridis_d(option = "turbo", name = "Taxa")
taxa_scale_colour <- scale_colour_viridis_d(option = "turbo", name = "Taxa")


# ---- 3. Assemblages per study (before vs after the 5-year filter) ------------
# Stacked bar for the 30 studies with the most assemblages: how many are
# retained vs dropped by the 5-year rule?
top_studies <- assemblage_summary %>%
  summarise(n_total = n(), n_kept = sum(keep), n_dropped = sum(!keep), .by = STUDY_ID) %>%
  arrange(desc(n_total)) %>%
  slice_head(n = 30) %>%
  mutate(STUDY_ID = factor(STUDY_ID, levels = rev(STUDY_ID)))

p_study <- bind_rows(
  top_studies %>% transmute(STUDY_ID, status = paste0("Kept (>= ", min_years, " yrs)"), n = n_kept),
  top_studies %>% transmute(STUDY_ID, status = paste0("Dropped (< ", min_years, " yrs)"), n = n_dropped)
) %>%
  ggplot(aes(x = n, y = STUDY_ID, fill = status)) +
  geom_col() +
  scale_fill_manual(values = c("#155f49", "#d9d956"), name = NULL) +
  labs(x = "Number of assemblage time series", y = "Study ID",
       title = "Assemblages per study (30 largest studies)") +
  theme(legend.position = "bottom")
save_fig(p_study, "03_assemblages_per_study_top30", h = 7)

p_study_hist <- ts5 %>%
  count(STUDY_ID, name = "n_aID") %>%
  ggplot(aes(n_aID)) +
  geom_histogram(bins = 30, fill = "#155f49", colour = "white") +
  scale_x_log10() +
  labs(x = "Assemblages per study (log scale, after filtering)", y = "Number of studies",
       title = "Most studies contribute only a few assemblages")
save_fig(p_study_hist, "03_assemblages_per_study_hist")


# ---- (a) Spatial location + density on a world map ---------------------------
world <- map_data("world")

base_map <- function() {
  list(
    geom_polygon(data = world, aes(long, lat, group = group),
                 fill = "grey92", colour = "grey70", linewidth = 0.1, inherit.aes = FALSE),
    coord_quickmap(xlim = c(-180, 180), ylim = c(-58, 85)),
    labs(x = NULL, y = NULL)
  )
}

# (a1) locations, coloured by taxa
p_map_pts <- ggplot() +
  base_map() +
  geom_point(data = ts5, aes(mean_lon, mean_lat, colour = TAXA),
             size = 0.9, alpha = 0.7) +
  taxa_scale_colour +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Location of terrestrial assemblage time series")
save_fig(p_map_pts, "a1_map_locations", w = 11, h = 6)

# (a2) density: number of ASSEMBLAGES per 5 x 5 degree bin
p_map_dens <- ggplot() +
  base_map() +
  stat_bin_2d(data = ts5, aes(mean_lon, mean_lat), binwidth = c(5, 5), alpha = 0.9) +
  scale_fill_viridis_c(option = "magma", direction = -1, trans = "log10",
                       name = "Assemblages\n(log10)") +
  labs(title = "Density of time series (assemblages per 5 x 5 degree bin)")
save_fig(p_map_dens, "a2_map_density_assemblages", w = 11, h = 6)

# (a3) density: number of independent STUDIES per bin. This removes the
# influence of a few studies that were split into hundreds of cells.
p_map_dens_study <- ggplot() +
  base_map() +
  stat_summary_2d(data = ts5, aes(mean_lon, mean_lat, z = STUDY_ID),
                  fun = function(z) length(unique(z)),
                  binwidth = c(5, 5), alpha = 0.9) +
  scale_fill_viridis_c(option = "viridis", direction = -1, name = "Studies",
                       breaks = function(x) unique(round(pretty(x)))) +
  labs(title = "Number of independent studies per 5 x 5 degree bin")
save_fig(p_map_dens_study, "a3_map_density_studies", w = 11, h = 6)


# ---- (b) Number of years of sampling per time series -------------------------
p_nyears <- ggplot(ts5, aes(n_years)) +
  geom_histogram(binwidth = 1, boundary = min_years - 0.5,
                 fill = "#155f49", colour = "white") +
  geom_vline(xintercept = median(ts5$n_years), linetype = "dashed") +
  annotate("text", x = median(ts5$n_years), y = Inf, vjust = 1.5, hjust = -0.1,
           label = paste("median =", median(ts5$n_years))) +
  labs(x = "Number of years sampled", y = "Number of time series",
       title = "Years of sampling per time series")
save_fig(p_nyears, "b_years_sampled")


# ---- (c) Temporal spread: years the time series covers -----------------------
p_duration <- ggplot(ts5, aes(duration)) +
  geom_histogram(binwidth = 2, fill = "#66c1d1", colour = "white") +
  geom_vline(xintercept = median(ts5$duration), linetype = "dashed") +
  labs(x = "Duration (last year - first year + 1)", y = "Number of time series",
       title = "Temporal extent of time series")
save_fig(p_duration, "c1_duration")

# Years sampled vs years spanned: points on the 1:1 line = sampled every year
p_dur_vs_n <- ggplot(ts5, aes(duration, n_years)) +
  geom_point(alpha = 0.25, size = 0.8, colour = "#155f49") +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  labs(x = "Duration (years spanned)", y = "Years sampled",
       title = "Sampling completeness (red line = sampled every year)")
save_fig(p_dur_vs_n, "c2_years_sampled_vs_duration")

# Which calendar years are covered, and by how many series / studies?
coverage_by_year <- samples_per_year %>%
  filter(assemblageID %in% ts5$assemblageID) %>%
  summarise(Assemblages = n_distinct(assemblageID),
            Studies     = n_distinct(STUDY_ID),
            .by = YEAR) %>%
  tidyr::pivot_longer(c(Assemblages, Studies), names_to = "metric", values_to = "n")

p_cov_year <- ggplot(coverage_by_year, aes(YEAR, n, colour = metric)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("#155f49", "#d98a1f"), name = NULL) +
  scale_y_log10() +
  labs(x = "Year", y = "Number sampled (log scale)",
       title = "Calendar-year coverage") +
  theme(legend.position = "bottom")
save_fig(p_cov_year, "c3_calendar_year_coverage")

# Start-year distribution
p_start <- ggplot(ts5, aes(start_year)) +
  geom_histogram(binwidth = 2, fill = "#66c1d1", colour = "white") +
  labs(x = "First year of time series", y = "Number of time series",
       title = "When do time series start?")
save_fig(p_start, "c4_start_year")


# ---- (d) Taxonomic spread ------------------------------------------------------
# Helper: studies vs assemblages per category. A category dominated by
# assemblages but not studies means a few big studies are doing the work.
count_plot <- function(df, var, title) {
  d <- df %>% mutate(category = coalesce(as.character(.data[[var]]), "Unknown"))
  counts <- bind_rows(
    d %>% summarise(n = n_distinct(STUDY_ID),      .by = category) %>% mutate(metric = "Studies"),
    d %>% summarise(n = n_distinct(assemblageID),  .by = category) %>% mutate(metric = "Assemblages")
  )
  ord <- counts %>% filter(metric == "Assemblages") %>% arrange(n) %>% pull(category)
  counts %>%
    mutate(category = factor(category, levels = ord),
           metric   = factor(metric, levels = c("Studies", "Assemblages"))) %>%
    ggplot(aes(n, category, fill = metric)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~metric, scales = "free_x") +
    scale_fill_manual(values = c("#d98a1f", "#155f49")) +
    labs(x = "Count", y = NULL, title = title)
}

p_taxa <- count_plot(ts5, "TAXA", "Taxonomic coverage")
save_fig(p_taxa, "d1_taxa_counts", h = 5)

p_taxa_years <- ggplot(ts5, aes(x = reorder(TAXA, n_years, FUN = median), y = n_years, fill = TAXA)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  scale_y_log10() +
  coord_flip() +
  taxa_scale_fill +
  labs(x = NULL, y = "Years sampled (log scale)",
       title = "Time-series length by taxonomic group")
save_fig(p_taxa_years, "d2_taxa_years_sampled", h = 5)

p_taxa_rich <- ggplot(ts5, aes(x = reorder(TAXA, n_species, FUN = median), y = n_species, fill = TAXA)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  scale_y_log10() +
  coord_flip() +
  taxa_scale_fill +
  labs(x = NULL, y = "Species recorded (log scale)",
       title = "Species richness (all years pooled) by taxonomic group")
save_fig(p_taxa_rich, "d3_taxa_species_richness", h = 5)


# ---- (e) Other coverage / quality checks -----------------------------------------
# e1. Latitudinal distribution
p_lat <- ggplot(ts5, aes(mean_lat, fill = TAXA)) +
  geom_histogram(binwidth = 5, boundary = 0, colour = "white", linewidth = 0.1) +
  coord_flip() +
  taxa_scale_fill +
  labs(x = "Latitude (5 degree bands)", y = "Number of time series",
       title = "Latitudinal coverage")
save_fig(p_lat, "e1_latitude", h = 6)

# e2. Climate and biome representation
p_climate <- count_plot(ts5, "CLIMATE", "Climate zones")
save_fig(p_climate, "e2_climate_counts", h = 4)

if ("BIOME_MAP" %in% names(ts5)) {
  p_biome <- count_plot(ts5, "BIOME_MAP", "Biomes (WWF)")
  save_fig(p_biome, "e3_biome_counts", w = 11, h = 6)
}

# e3. Taxa x climate matrix: where are the gaps?
p_tax_clim <- ts5 %>%
  mutate(CLIMATE = coalesce(as.character(CLIMATE), "Unknown")) %>%
  count(TAXA, CLIMATE) %>%
  ggplot(aes(CLIMATE, TAXA, fill = n)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n), size = 3) +
  scale_fill_viridis_c(direction = -1, trans = "log10", name = "Assemblages") +
  labs(x = NULL, y = NULL, title = "Taxa x climate coverage (number of assemblages)")
save_fig(p_tax_clim, "e4_taxa_by_climate", h = 5)

# e4. Dominance: do a handful of studies make up most of the dataset?
dominance <- ts5 %>%
  count(STUDY_ID, name = "n_aID") %>%
  arrange(desc(n_aID)) %>%
  mutate(rank = row_number(), cum_share = cumsum(n_aID) / sum(n_aID))

cat("\nShare of all assemblages from the 5 largest studies:",
    round(dominance$cum_share[5] * 100, 1), "%\n")

p_dom <- ggplot(dominance, aes(rank, cum_share)) +
  geom_line(colour = "#155f49", linewidth = 0.9) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Study rank (largest to smallest)", y = "Cumulative share of assemblages",
       title = "Dominance of large studies")
save_fig(p_dom, "e5_study_dominance")

# e5. Sampling effort and the cost of rarefaction
p_retained <- ggplot(ts5, aes(prop_retained)) +
  geom_histogram(bins = 40, fill = "#155f49", colour = "white") +
  labs(x = "Proportion of samples retained after sample-based rarefaction",
       y = "Number of time series",
       title = "How much data would rarefaction throw away?",
       subtitle = "Values near 1 = even sampling effort across years")
save_fig(p_retained, "e6_rarefaction_retention")

# e6. Temporal gaps
p_gap <- ggplot(ts5, aes(max_gap)) +
  geom_histogram(binwidth = 1, fill = "#66c1d1", colour = "white") +
  scale_y_sqrt() +   # sqrt (not log) axis: copes with empty bins
  labs(x = "Longest gap between sampled years (years missed)",
       y = "Number of time series (sqrt scale)",
       title = "Temporal gaps within time series")
save_fig(p_gap, "e7_max_gap")

# e7. Single-location vs multi-location studies
p_method <- ts5 %>%
  count(StudyMethod, TAXA) %>%
  ggplot(aes(TAXA, n, fill = StudyMethod)) +
  geom_col() +
  scale_fill_manual(values = c(SL = "#155f49", ML = "#d9d956"), name = "Study type") +
  coord_flip() +
  labs(x = NULL, y = "Number of assemblages",
       title = "Single- (SL) vs multi-location (ML) studies by taxa")
save_fig(p_method, "e8_SL_vs_ML", h = 5)


# ---- Quick numerical summary --------------------------------------------------
final_summary <- ts5 %>%
  summarise(
    n_assemblages     = n(),
    n_studies         = n_distinct(STUDY_ID),
    n_taxa            = n_distinct(TAXA),
    median_years      = median(n_years),
    median_duration   = median(duration),
    median_species    = median(n_species),
    median_retained   = median(prop_retained),
    pct_with_gaps     = mean(max_gap > 0) * 100
  )
print(t(final_summary))

# ---------------------------------------------------------------------------------
# NEXT STEP: rarefaction on the filtered, gridded data, e.g.
#   grid_rare <- resampling(x = grid_samples_5yr, measure = "ABUNDANCE",
#                           resamps = 1, conservative = FALSE)
# ---------------------------------------------------------------------------------






