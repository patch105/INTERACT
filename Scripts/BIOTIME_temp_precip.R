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

dir.create(here("Data/Raw"),    showWarnings = FALSE)
dir.create(here("Plots/BioTIME_prelim"), showWarnings = FALSE)

min_years <- 5    # minimum number of distinct sampled years per time series
grid_res  <- 12   # ~96 km2 hexagons, as recommended for the whole database

# Other potential filters
# Minimum time frame of the time series? Or time frame over which it's been measured?
# Freshwater?
# Biomass?

#####################################################
# 0. Load data -----------------------------------------------------------------
#####################################################
# fread() is much faster than read.csv() on a file this size
BT_data <- data.table::fread(here("Data/Raw/biotime_v2_rawdata_2025.csv"),  data.table = FALSE)
BT_meta <- data.table::fread(here("Data/Raw/biotime_v2_metadata_2025.csv"), data.table = FALSE)

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

# Remove records from first two filters
dat <- BT_data %>%
  filter(STUDY_ID %in% meta_sel$STUDY_ID)
log_step("  ...records in those studies", n_distinct(dat$STUDY_ID), records = nrow(dat))

# Free memory - the full table is no longer needed
rm(BT_data); gc()

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
# 2b. Species-level resolution only -------------------------------------------
#####################################################
# Drop records where the taxon is identified to a coarser resolution than
# species (e.g. genus, family). Records with a missing resolution are dropped
# too, because we can't confirm they are species-level.

cat("\nRecords by taxonomic resolution (before filtering):\n")
print(table(dat$resolution, useNA = "ifany"))

dat <- dat %>%
  filter(tolower(trimws(resolution)) == "species")
log_step("Records identified to species level",
         n_distinct(dat$STUDY_ID), records = nrow(dat))

# From 282 to 277 studies


#####################################################
# 2c. Other filters -----------------------------------------------------------
#####################################################
# Keep only records where no treatment was applied (GENERAL_TREAT is missing,
# "FALSE", blank or "Control"). GENERAL_TREAT is a study-level field, so if it
# isn't already in the record table it is joined on from the metadata.
if (!"GENERAL_TREAT" %in% names(dat)) {
  dat <- dat %>%
    left_join(BT_meta %>% select(STUDY_ID, GENERAL_TREAT), by = "STUDY_ID")
}

cat("\nRecords by GENERAL_TREAT (before filtering):\n")
print(table(dat$GENERAL_TREAT, useNA = "ifany"))

dat <- dat %>%
  filter(is.na(GENERAL_TREAT) | GENERAL_TREAT %in% c("FALSE", "", "Control"))
log_step("Records with no treatment applied",
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


## Filtering out studies with < 5 years removes ~100 studies out of 241

#####################################################
# 4. Gridding ---------------------------------------------------------------
#####################################################
# Species identifiers: gridding() drops ID_SPECIES and renames valid_name to
# `Species`. Both live in the raw record table, so capture the lookup here while
# `dat` still has them (restricted to the records that survived the filters
# above), and join it back onto the gridded data straight after gridding().
# One row per valid_name so the join can't duplicate records. If a valid_name is
# linked to several ID_SPECIES (e.g. synonyms), the smallest ID is kept and
# n_ID_SPECIES flags it (n_ID_SPECIES must be computed first).
sp_lookup <- dat %>%
  distinct(valid_name, ID_SPECIES) %>%
  summarise(n_ID_SPECIES = n_distinct(ID_SPECIES),
            ID_SPECIES   = min(ID_SPECIES),
            .by = valid_name)

cat("\nvalid_names linked to more than one ID_SPECIES:",
    sum(sp_lookup$n_ID_SPECIES > 1), "of", nrow(sp_lookup), "\n")

# valid_names linked to more than one ID_SPECIES: 521 valid names of 14,375

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

# Reconstruct the grid to save the centroids
dgg <- dggridR::dgconstruct(res = grid_res) # same call gridding() makes

cells <- grid_samples %>% distinct(STUDY_ID, cell, assemblageID)

# cell centres
ctr <- dggridR::dgSEQNUM_to_GEO(dgg, unique(cells$cell)) # Find the globally unique number identifying surface polygon
centres <- tibble(cell = unique(cells$cell),
                  cell_lon = ctr$lon_deg,
                  cell_lat = ctr$lat_deg)

# Cell outlines (sf polygons, with a `seqnum` column)
hexes <- dggridR::dgcellstogrid(dgg, unique(cells$cell))
hexes <- vect(hexes)

# Attach centres to the assemblage table so there's one centre per assemblage
assemblages <- cells %>% left_join(centres, by = "cell")

# store
writeVector(hexes, here("96km_hexagons_BioTIME.shp"))
saveRDS(assemblages, "assemblage_centres.rds")
st_write(hexes, "assemblage_hexes.gpkg")


rm(dat); gc()

# Add valid_name and ID_SPECIES back (Species is the renamed valid_name)
n_gridded <- nrow(grid_samples)
grid_samples <- grid_samples %>%
  mutate(valid_name = Species) %>%
  left_join(sp_lookup, by = "valid_name")
stopifnot(nrow(grid_samples) == n_gridded)   # join must not add/drop rows
cat("Gridded records with no ID_SPECIES match:", sum(is.na(grid_samples$ID_SPECIES)), "\n") # WAS ZERO

log_step("Gridded (all assemblages)",
         n_distinct(grid_samples$STUDY_ID),
         assemblages = n_distinct(grid_samples$assemblageID), # Count the number of unique assemblage+studyID points
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

# Studies: 152 | Unique grid cells: 6802 | Assemblage time series: 7848 
# StudyMethod: ML = 38, SL = 114
# For number of distinct assemblages per study: median is one, max is 3071

# So unique grid cells are 6802, but there's more assemblage time series because there are assemblages that are from different studies but found in the same grid cell

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

# Assemblages with >= 5 years - DON'T HAVE TO BE CONSECUTIVE: 4262 of 7911 

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
write.csv(grid_samples_5yr, here("Data/Intermediate", "BT_terrestrial_abund_gridded_5yr.csv"))
write.csv(assemblage_summary, here("Data/Intermediate", "BT_assemblage_summary.csv"), row.names = FALSE)
write.csv(attrition_tbl,      here("Data/Intermediate", "BT_cleaning_attrition.csv"), row.names = FALSE)



###########################################################
# SUBSETTING for temperature x precipitation study  -------------------------
############################################################

BT_terrestrial_abund_gridded_5yr <- read.csv(here("Data/Intermediate", "BT_terrestrial_abund_gridded_5yr.csv"))

length(unique(BT_terrestrial_abund_gridded_5yr$assemblageID))

count <- count(BT_terrestrial_abund_gridded_5yr, assemblageID, YEAR)

count <- count(BT_terrestrial_abund_gridded_5yr, assemblageID, YEAR, SAMPLE_DESC)

# Start off with 4208 unique assemblages, with 50,769 sampling occasions (assemblageID grid cell per year)

#### Filter to the year 2000
BT_microclim <- BT_terrestrial_abund_gridded_5yr %>% filter(YEAR >= 1941)

length(unique(BT_microclim$assemblageID))

count <- count(BT_microclim, assemblageID, YEAR)

count <- count(BT_microclim, assemblageID, YEAR, SAMPLE_DESC)

# Drops to 4161 assemblages, with 34,951 sampling occasions (assemblageID grid cell per year), and 287,364 samples


#### NOW, re-apply the >= 5 year filter :

BT_microclim_5yr <- BT_microclim %>% filter(n_distinct(YEAR) >= 5, .by = assemblageID)

length(unique(BT_microclim_5yr$assemblageID))

count <- count(BT_microclim_5yr, assemblageID, YEAR)

count <- count(BT_microclim_5yr, assemblageID, YEAR, SAMPLE_DESC)

# Drops to 3,600 assemblages, with 33,197 sampling occasions (assemblageID grid cell per year), 270,252 samples across these occasions


dgg <- dggridR::dgconstruct(res = res)

# Global grid of hexagonal cells derived from dgconstruct and assigning the individual samples to the cells across the grid based on its latitude and longitude.