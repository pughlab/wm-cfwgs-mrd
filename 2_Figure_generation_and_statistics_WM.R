# Figure_generation_and_statistics_WM.R -----------------------------------------------------------
# Reproducible analysis pipeline for cfWGS MRD detection in Waldenström
# Macroglobulinemia (WM).  Generates all tables & figures for:
#   Chow S., Abelman D. *et al.* (2025)
#   Author: D. Abelman | Date: May 5, 2025
# ------------------------------------------------------------------------------

# Libraries ---------------------------------------------------------------
### Load libraries
library(tidyverse)
library(ggrepel)
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(janitor)
library(RColorBrewer)
library(GGally)
library(glue)
library(knitr)






# Load MRDetect Data ------------------------------------------------------

## This starts after the main MRDetect .RDS file 
Merged_MRDetect_BRAWM_zscore <- readRDS("path/to/your/Merged_MRDetect_BRAWM_zscore.rds")


### Set dirs 
# Set output directory
outdir <- "path/to/output/directory/"



# Process MRDetect Data ---------------------------------------------------


# 1. Filter data for the Patient of interest
df_patient <- Merged_MRDetect_BRAWM_zscore %>%
  filter(
    Mut_source == "BM_cells",
    Filter_source == "STR_encode",
    plotting_type == "Matched_plasma",
    Sample_type_Bam == "P"
  )

df_patient$Mrd_by_WGS <- ifelse(df_patient$sites_rate_zscore_charm > 4.5, "Positive", "Negative")


# 2. Filter data for healthy controls (BAM starts with "TGL49")
df_healthy <- Merged_MRDetect_BRAWM_zscore %>%
  filter(
    Mut_source == "BM_cells",
    Filter_source == "STR_encode",
    grepl("^TGL49", BAM)
  ) %>%
  # Label these as "Healthy_controls" so they plot on the far right
  mutate(timepoint_info_Bam = "Healthy_controls")

# 3. Combine the two datasets
df_combined <- bind_rows(df_patient, df_healthy) %>%
  mutate(
    timepoint_info_plot = case_when(
      timepoint_info_Bam == "T0" ~ "Baseline",
      timepoint_info_Bam == "Healthy_controls" ~ "Healthy_controls",
      TRUE ~ paste0("Cycle ", timepoint_info_Bam)
    )
  ) %>%
  mutate(
    timepoint_info_plot = factor(
      timepoint_info_plot,
      levels = c("Baseline", "Cycle C7", "Cycle C12", "Cycle C18", "Healthy_controls")
    )
  )

## Export tables 
# 1) Filter & select
df_export <- df_combined %>%
  filter(plotting_type == "Matched_plasma") %>%
  select(
    # 1) patient identifiers & timepoints up front
    Patient,
    Patient_other_ID,
    
    # 2) technical/sample IDs
    Sample_ID_Bam,
    Sample_type_Bam,
    timepoint_info_Bam,
    timepoint_info_plot,
    
    # 3) the detection metrics
    sites_checked,
    reads_checked,
    sites_detected,
    reads_detected,
    total_reads,
    detection_rate,
    detection_rate_as_reads_detected_over_reads_checked,
    detection_rate_as_reads_detected_over_total_reads,
    sites_detection_rate,
    
    # 4) z‑score metrics and MRD call
    detection_rate_zscore_charm,
    sites_rate_zscore_charm,
    detection_rate_zscore_reads_checked_charm,
    detection_rate_zscore_total_reads_charm,
    Mrd_by_WGS
  )

# 2) Write it out
out_file <- file.path(outdir, "BRAWM_MRDetect_filtered_export_May2025.txt")
write.table(
  df_export,
  file      = out_file,
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

message("Exported ", nrow(df_export), " rows to:\n  ", out_file)

### Transform to current IDs
df_combined <- df_combined %>%
  mutate(
    Patient_number_plot = case_when(
      Patient_other_ID == "PATIENT_ID_1" ~ "Patient 1", ### Original IDs removed for confidentiality
      Patient_other_ID == "PATIENT_ID_2" ~ "Patient 2",
      Patient_other_ID == "PATIENT_ID_3" ~ "Patient 3",
      Patient_other_ID == "PATIENT_ID_4" ~ "Patient 4",
      Patient_other_ID == "PATIENT_ID_5" ~ "Patient 5",
      Patient_other_ID == "PATIENT_ID_6" ~ "Patient 6",
      Patient_other_ID == "PATIENT_ID_7" ~ "Patient 7",
      TRUE                               ~ NA_character_
    )
  )


# Plot

# Creating Figure 1A ------------------------------------------------------

ggplot(df_combined, aes(x = timepoint_info_plot, y = detection_rate_as_reads_detected_over_reads_checked)) +
  # Connect each patient's timepoints with a red line (exclude healthy controls)
  geom_line(
    data = filter(df_combined, timepoint_info_plot != "Healthy_controls"),
    aes(group = Patient),
    color = "red",
    size = 1
  ) +
  # Plot the Patient's points in red
  geom_point(
    data = filter(df_combined, timepoint_info_plot != "Healthy_controls"),
    color = "red",
    size = 3
  ) +
  # Plot healthy controls on the far right (jitter to avoid overlap)
  geom_jitter(
    data = filter(df_combined, timepoint_info_plot == "Healthy_controls"),
    width = 0.1,
    color = "black",
    size = 3
  ) +
  # Use a continuous y-axis that starts at 0 and multiply the labels by 100
  scale_y_continuous(
    limits = c(0, NA),
    labels = function(x) x * 100,
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Tumor-informed detection rate across timepoints and healthy controls",
    x = "Timepoint",
    y = "Tumor-informed cumulative VAF (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot to outdir
ggsave(
  filename = file.path(outdir, "Tumor_informed_detection_rate_BRAWM.png"),
  plot     = last_plot(),
  width    = 10,
  height   = 6,
  dpi      = 500
)




### Add boxplot for healthy control 
#### One only 
# 1) define a small helper
hc   <- df_combined %>% filter(timepoint_info_plot == "Healthy_controls")
pA <- ggplot(df_combined, aes(x = timepoint_info_plot, y = detection_rate_as_reads_detected_over_reads_checked)) +
  
  # patient red lines & points
  geom_line(
    data = filter(df_combined, timepoint_info_plot != "Healthy_controls"),
    aes(group = Patient),
    color = "red", size = 1
  ) +
  geom_point(
    data  = filter(df_combined, timepoint_info_plot != "Healthy_controls"),
    color = "red", size = 3
  ) +
  
  # boxplot of _all_ healthy controls at the far‐right x
  geom_boxplot(
    data            = hc,
    aes(x           = timepoint_info_plot,
        y           = detection_rate_as_reads_detected_over_reads_checked),
    width           = 0.3,
    fill            = "grey80",
    outlier.shape   = NA
  ) +
  # tiny points on top
  geom_jitter(
    data   = hc,
    aes(x    = timepoint_info_plot,
        y    = detection_rate_as_reads_detected_over_reads_checked),
    width  = 0.1,
    color  = "black",
    size   = 1.5
  ) +
  
  # axes, labels, theme (unchanged)
  scale_y_continuous(
    limits = c(0, NA),
    labels = function(x) x * 100,
    expand = expansion(mult = c(0, .05))
  ) +
  labs(
    title = "Tumor-informed detection rate across timepoints and healthy controls",
    x     = "Timepoint",
    y     = "Tumor-informed cumulative VAF (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(outdir, "figA_detection_rate_all_HC_box.png"),
  plot   = pA,
  width  = 10, height = 6, dpi = 500
)


### Alternatively, make one boxplot per pt 
pB_fixed <- ggplot(df_combined,
                   aes(x = timepoint_info_plot,
                       y = detection_rate_as_reads_detected_over_reads_checked)) +
  # patient trajectories (unchanged)
  geom_line(
    data = filter(df_combined, timepoint_info_plot != "Healthy_controls"),
    aes(group = Patient),
    color = "red", size = 1
  ) +
  geom_point(
    data  = filter(df_combined, timepoint_info_plot != "Healthy_controls"),
    color = "red", size = 3
  ) +
  
  # **per‐Patient boxplots at Healthy_controls:**  
  geom_boxplot(
    data          = hc,
    aes(fill      = Patient),                # ← map fill so dodge works
    width         = 0.4,
    outlier.shape = NA,
    position      = position_dodge(width = 0.8)
  ) +
  
  # **jitter on top, dodged by Patient via colour mapping**  
  geom_jitter(
    data          = hc,
    aes(colour    = Patient),                # ← map colour so dodge works
    position      = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.8
    ),
    size          = 2,
    alpha         = 0.8
  ) +
  
  # axes, labels, theme (unchanged)
  scale_y_continuous(
    limits = c(0, NA),
    labels = function(x) x * 100,
    expand = expansion(mult = c(0, .05))
  ) +
  labs(
    title = "Tumor-informed detection rate across timepoints and healthy controls",
    x     = "Timepoint",
    y     = "Tumor-informed cumulative VAF (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"   # hide the per‐patient legend if too busy
  )

# save
ggsave(
  file.path(outdir, "figB_detection_rate_perPatient_HC_box_fixed.png"),
  plot   = pB_fixed,
  width  = 10, height = 6, dpi = 500
)




## Now add labels
# 4. Prepare labels: select the last non-healthy timepoint for each patient
df_labels <- df_combined %>%
  filter(timepoint_info_plot != "Healthy_controls") %>%
  group_by(Patient) %>%
  filter(as.numeric(timepoint_info_plot) == max(as.numeric(timepoint_info_plot))) %>%
  ungroup()


# Compute the last timepoint for each patient (non-healthy)
p1 <- ggplot(df_combined, aes(x = timepoint_info_plot, y = detection_rate_as_reads_detected_over_reads_checked)) +
  # Connect each patient's timepoints with lines colored by MRD status (exclude healthy controls)
  geom_line(
    data = filter(df_combined, timepoint_info_plot != "Healthy_controls"),
    aes(group = Patient, color = Mrd_by_WGS),
    size = 1
  ) +
  # Plot the patient's points colored by MRD status
  geom_point(
    data = filter(df_combined, timepoint_info_plot != "Healthy_controls"),
    aes(color = Mrd_by_WGS),
    size = 3
  ) +
  # — replaced healthy‐control layer with boxplot + jitter —
  # 1) Jitter behind with reduced alpha
  geom_jitter(
    data   = filter(df_combined, timepoint_info_plot == "Healthy_controls"),
    aes(x    = timepoint_info_plot,
        y    = detection_rate_as_reads_detected_over_reads_checked),
    width  = 0.1,
    size   = 1,
    color  = "black",
    alpha  = 0.3         # make points semi‐transparent
  ) +
  
  # 2) Boxplot on top
  geom_boxplot(
    data          = filter(df_combined, timepoint_info_plot == "Healthy_controls"),
    aes(x         = timepoint_info_plot,
        y         = detection_rate_as_reads_detected_over_reads_checked),
    width         = 0.3,
    fill          = "grey80",
    color         = "black",
    outlier.shape = NA,
    alpha         = 0.7         # slight transparency so you still see points
  ) +
  # Add patient ID labels at the end of each line
  # replace geom_text(...) with:
  geom_text_repel(
    data            = df_labels,
    aes(
      x     = timepoint_info_plot,
      y     = detection_rate_as_reads_detected_over_reads_checked,
      label = Patient_number_plot
    ),
    nudge_x         = 0.37,          # push labels to the right
    direction       = "y",          # only repel vertically
    hjust           = 0,            # left‐align text
    segment.size    = 0.3,          # connector line thickness
    segment.color   = "grey50",     # connector line color
    box.padding     = 0.35,          # more space around each label
    point.padding   = 0.5,          # avoid overlapping the data point
    max.overlaps    = Inf,          # label everything (with repulsion)
    size            = 3
  ) +
  
  # Define manual colors: red for Positive, red-grey for Negative
  scale_color_manual(
    values = c("Positive" = "red", "Negative" = "gray50"),
    na.value = "black"
  ) +
  # Force x-axis to new names
  scale_x_discrete(
    limits = c("Baseline", "Cycle C7", "Cycle C12", "Cycle C18", "Healthy_controls"), 
    labels = c("Baseline", "C7",        "C12",        "M18",        "Healthy\ncontrols")
  ) +
  # Use a continuous y-axis that starts at 0 and multiply the labels by 100
  scale_y_continuous(
    limits = c(0, NA),
    labels = function(x) x * 100,
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Tumor-informed cfDNA Detection Rate Over Treatment Cycles by WGS",
    x = "Timepoint",
    y = "Tumor-informed cumulative VAF (%)",
    color = "MRD by WGS"
  ) +
  theme_classic() +
  theme(
    #   axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )


ggsave(filename = file.path(outdir, "Tumor_informed_detection_rate_BRAWM_updated7.png"), plot = p1, width = 10, height = 6, dpi = 1000)


## If want with right legend
test <- ggplot(df_combined, aes(x = timepoint_info_plot, y = detection_rate_as_reads_detected_over_reads_checked)) +
  # Connect each patient's timepoints with lines colored by MRD status (exclude healthy controls)
  geom_line(
    data = filter(df_combined, timepoint_info_plot != "Healthy_controls"),
    aes(group = Patient, color = Mrd_by_WGS),
    size = 1
  ) +
  # Plot the patient's points colored by MRD status
  geom_point(
    data = filter(df_combined, timepoint_info_plot != "Healthy_controls"),
    aes(color = Mrd_by_WGS),
    size = 3
  ) +
  # — replaced healthy‐control layer with boxplot + jitter —
  # 1) Jitter behind with reduced alpha
  geom_jitter(
    data   = filter(df_combined, timepoint_info_plot == "Healthy_controls"),
    aes(x    = timepoint_info_plot,
        y    = detection_rate_as_reads_detected_over_reads_checked),
    width  = 0.1,
    size   = 1,
    color  = "black",
    alpha  = 0.3         # make points semi‐transparent
  ) +
  
  # 2) Boxplot on top
  geom_boxplot(
    data          = filter(df_combined, timepoint_info_plot == "Healthy_controls"),
    aes(x         = timepoint_info_plot,
        y         = detection_rate_as_reads_detected_over_reads_checked),
    width         = 0.3,
    fill          = "grey80",
    color         = "black",
    outlier.shape = NA,
    alpha         = 0.7         # slight transparency so you still see points
  ) +
  # Add patient ID labels at the end of each line
  # replace geom_text(...) with:
  geom_text_repel(
    data            = df_labels,
    aes(
      x     = timepoint_info_plot,
      y     = detection_rate_as_reads_detected_over_reads_checked,
      label = Patient_number_plot
    ),
    nudge_x         = 0.37,          # push labels to the right
    direction       = "y",          # only repel vertically
    hjust           = 0,            # left‐align text
    segment.size    = 0.3,          # connector line thickness
    segment.color   = "grey50",     # connector line color
    box.padding     = 0.35,          # more space around each label
    point.padding   = 0.5,          # avoid overlapping the data point
    max.overlaps    = Inf,          # label everything (with repulsion)
    size            = 3
  ) +
  
  # Define manual colors: red for Positive, red-grey for Negative
  scale_color_manual(
    values = c("Positive" = "red", "Negative" = "gray50"),
    na.value = "black"
  ) +
  # Force x-axis to new names
  scale_x_discrete(
    limits = c("Baseline", "Cycle C7", "Cycle C12", "Cycle C18", "Healthy_controls"), 
    labels = c("Baseline", "C7",        "C12",        "M18",        "Healthy\ncontrols")
  ) +
  # Use a continuous y-axis that starts at 0 and multiply the labels by 100
  scale_y_continuous(
    limits = c(0, NA),
    labels = function(x) x * 100,
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Tumor-informed cfDNA Detection Rate Over Treatment Cycles by WGS",
    x = "Timepoint",
    y = "Tumor-informed cumulative VAF (%)",
    color = "MRD by WGS"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# 1) build the log‐scale plot from p1
p2 <- p1 +
  scale_y_log10(
    labels = function(x) paste0(round(x * 100, 1), "%"),
    expand = expansion(mult = c(0.01, 0.2))
  ) +
  # Force x-axis to new names
  scale_x_discrete(
    limits = c("Baseline", "Cycle C7", "Cycle C12", "Cycle C18", "Healthy_controls"), 
    labels = c("Baseline", "C7",        "C12",        "M18",        "Healthy\ncontrols")
  ) +
  annotation_logticks(sides = "l") +
  labs(y = "Tumor-informed cumulative VAF (%) [log₁₀]") +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
ggsave(filename = file.path(outdir, "Tumor_informed_detection_rate_right_ledend_updated.png"), plot = test, width = 10, height = 6, dpi = 1000)




### Now add the inlay of the log plot within Figure 1A

# 1) Build the no‑Baseline, log‑scale plot
p_noBase_log <- ggplot(df_combined, 
                       aes(x = timepoint_info_plot, 
                           y = detection_rate_as_reads_detected_over_reads_checked)) +
  
  # 1a) patient lines & points (but drop Baseline & Healthy_controls for these layers)
  geom_line(
    data = df_combined %>% filter(timepoint_info_plot %in% c("Cycle C7","Cycle C12","Cycle C18")),
    aes(group = Patient, color = Mrd_by_WGS),
    size = 1
  ) +
  geom_point(
    data = df_combined %>% filter(timepoint_info_plot %in% c("Cycle C7","Cycle C12","Cycle C18")),
    aes(color = Mrd_by_WGS),
    size = 3
  ) +
  # Force x-axis to new names
  scale_x_discrete(
    limits = c("Cycle C7", "Cycle C12", "Cycle C18", "Healthy_controls"), 
    labels = c("C7",        "C12",        "M18",        "Healthy\ncontrols")
  ) +
  
  # 1b) healthy‑control box + jitter (unchanged apart from dropping Baseline)
  geom_jitter(
    data   = df_combined %>% filter(timepoint_info_plot == "Healthy_controls"),
    aes(x    = timepoint_info_plot,
        y    = detection_rate_as_reads_detected_over_reads_checked),
    width  = 0.1,
    size   = 1,
    color  = "black",
    alpha  = 0.3
  ) +
  geom_boxplot(
    data          = df_combined %>% filter(timepoint_info_plot == "Healthy_controls"),
    aes(x         = timepoint_info_plot,
        y         = detection_rate_as_reads_detected_over_reads_checked),
    width         = 0.3,
    fill          = "grey80",
    color         = "black",
    outlier.shape = NA,
    alpha         = 0.7
  ) +
  
  # 1c) labels at last cycle (C18) only
  geom_text_repel(
    data            = df_labels,
    aes(
      x     = timepoint_info_plot,
      y     = detection_rate_as_reads_detected_over_reads_checked,
      label = Patient_number_plot
    ),
    nudge_x         = 0.3,
    direction       = "y",
    hjust           = 0,
    segment.size    = 0.3,
    segment.color   = "grey50",
    box.padding     = 0.4,
    point.padding   = 0.5,
    max.overlaps    = Inf,
    size            = 3
  ) +
  
  # 2) log‑scale and ticks
  scale_y_log10(
    labels = function(x) paste0(round(x * 100, 1), "%"),
    expand = expansion(mult = c(0.01, 0.2))
  ) +
  annotation_logticks(sides = "l") +
  
  # 3) colors, axes, labs
  scale_color_manual(
    values   = c("Positive" = "red", "Negative" = "gray50"),
    na.value = "black"
  ) +
  labs(
    title = "Tumor-informed detection rate across cycles (no Baseline) and healthy controls",
    x     = "Timepoint",
    y     = "Tumor-informed cumulative VAF (%) [log₁₀]",
    color = "MRD by WGS"
  ) +
  
  # 4) theme
  theme_classic() +
  theme(
    #    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# 2) Save it
ggsave(
  filename = file.path(outdir, "Tumor_informed_detection_rate_log_noBaseline.png"),
  plot     = p_noBase_log,
  width    = 8, 
  height   = 5, 
  dpi      = 500
)

# 3) export both
ggsave(
  filename = file.path(outdir, "Tumor_informed_detection_rate_log_withBaseline.png"),
  plot     = p2,
  width    = 8, height = 5, dpi = 500
)













# Making Figure 1B - VDJ Rearrangements ------------------------------------

##### Now load in the clinical data - the VAFS by NGS, clonoSEQ, and anything else she wants me to include
# 1) path to your file
xlsx_file <- "/path/WM-BRAWM-VDJclonalfract-time-updated.xlsx"

# 1) Read both sheets as all text
cf_wide <- read_excel(
  xlsx_file,
  sheet     = "cfDNA",
  col_types = "text"
) %>% rename(SampleID = 1)

bm_wide <- read_excel(
  xlsx_file,
  sheet     = "bone marrow",
  col_types = "text"
) %>% rename(SampleID = 1)

# 2) Pivot each to long
cf_long <- cf_wide %>%
  pivot_longer(
    cols      = -SampleID,
    names_to  = "Chain",
    values_to = "ClonalFraction"
  ) %>%
  mutate(Compartment = "cfDNA")

bm_long <- bm_wide %>%
  pivot_longer(
    cols      = -SampleID,
    names_to  = "Chain",
    values_to = "ClonalFraction"
  ) %>%
  mutate(Compartment = "BoneMarrow")

# 3) Bind the two long tables
dat_long <- bind_rows(cf_long, bm_long)

# 4) Clean up missing and convert to numeric
dat_clean <- dat_long %>%
  mutate(
    ClonalFraction = na_if(ClonalFraction, "-"),
    ClonalFraction = na_if(ClonalFraction, "?"),
    ClonalFraction = as.numeric(ClonalFraction)
  ) %>%
  # 5) Extract Patient & Timepoint, now using M12/M18
  mutate(
    Patient   = str_remove(SampleID, "-c7$|-m12$|-m18$"),
    Timepoint = case_when(
      str_detect(SampleID, "-c7$")  ~ "C7",
      str_detect(SampleID, "-m12$") ~ "C12",    # was "C12"
      str_detect(SampleID, "-m18$") ~ "M18",    # was "C18"
      TRUE                           ~ "Baseline"
    ),
    Timepoint = factor(
      Timepoint,
      levels = c("Baseline", "C7", "C12", "M18")
    )
  )

# replace any NA with 0
dat_clean <- dat_clean %>% mutate(
  ClonalFraction = tidyr::replace_na(ClonalFraction, 0)
)
dat_clean <- dat_clean %>% filter(!is.na(SampleID))

# 1) Identify each patient’s dominant clone at Baseline
dominant_clone <- dat_clean %>%
  filter(Timepoint == "Baseline") %>%
  group_by(Patient) %>%
  filter(!is.na(ClonalFraction)) %>%
  slice_max(order_by = ClonalFraction, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Patient, Chain)

### Edit this to the ones Signy wanted when multiple dominant clones present
dominant_clone <- dominant_clone %>%
  mutate(
    Chain = case_when(
      Patient == "PATIENT_ID_3"  ~ "IGL-1",
      Patient == "PATIENT_ID_1"  ~ "IGH-1",
      Patient == "PATIENT_ID_7"   ~ "IGL-1",
      TRUE                  ~ Chain
    )
  )

# 2) Join back to the full data so you only have those chains per patient
dat_to_plot <- dat_clean %>%
  inner_join(dominant_clone, by = c("Patient", "Chain"))

# 3) Quick check: one chain per patient
dat_to_plot %>% 
  distinct(Patient, Chain) %>% 
  print()

dat_to_plot <- dat_to_plot %>%
  mutate(
    Patient_plot = recode(
      Patient,
      "PATIENT_ID_1"  = "Patient 1",
      "PATIENT_ID_2"   = "Patient 2",
      "PATIENT_ID_3"  = "Patient 3",
      "PATIENT_ID_4"  = "Patient 4",
      "PATIENT_ID_5"   = "Patient 5",
      "PATIENT_ID_6"  = "Patient 6",
      "PATIENT_ID_7"   = "Patient 7"
    ),
    # now make it a factor in the exact order you want
    Patient_plot = factor(
      Patient_plot,
      levels = paste("Patient", 1:7)
    )
  )

# 4) Now plot dat_to_plot instead of dat_clean
ggplot(dat_to_plot, aes(x = Timepoint, y = ClonalFraction, color = Patient)) +
  geom_line(aes(group = Patient), size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Dominant V(D)J Clone Fraction Over Time",
    y     = "Clonal Fraction (%)",
    x     = "Timepoint"
  ) +
  theme_classic()


# 1) Time‐series plots for each compartment
p_cfDNA <- dat_to_plot %>%
  filter(Compartment == "cfDNA") %>%
  ggplot(aes(Timepoint, ClonalFraction, group = Patient_plot, color = Patient_plot)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(add = c(0.01, 0.05))
  ) +
  scale_color_brewer(
    palette = "Dark2",
    name    = "Patient ID"
  ) +
  labs(
    title = "Dominant cfDNA V(D)J Clone Fraction",
    x     = "Treatment Cycle",
    y     = "cfDNA Clonal Fraction (%)"
  ) +
  theme_classic() +
  theme(
    axis.line.x        = element_line(color = "black"),
    axis.ticks         = element_line(color = "black"),
    legend.position    = "right",
    plot.title         = element_text(face = "bold", hjust = 0.5)
  )


p_BM <- dat_to_plot %>%
  filter(Compartment == "BoneMarrow") %>%
  ggplot(aes(Timepoint, ClonalFraction,
             group = Patient_plot,
             color = Patient_plot)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(add = c(0.01, 0.05))
  ) +
  scale_color_brewer(
    palette = "Dark2",
    name    = "Patient ID"
  ) +
  labs(
    title = "Dominant BM V(D)J Clone Fraction",
    x     = "Treatment Cycle",
    y     = "Bone Marrow Clonal Fraction (%)"
  ) +
  theme_classic() +
  theme(
    axis.line.x     = element_line(color = "black"),
    axis.ticks      = element_line(color = "black"),
    legend.position = "none",
    plot.title      = element_text(face = "bold", hjust = 0.5)
  )

# 2) Build a wide table so each row has BM & cfDNA for the same Patient_plot/Chain/Timepoint
df_wide <- dat_to_plot %>%
  filter(Compartment %in% c("cfDNA", "BoneMarrow")) %>%
  select(Patient_plot, Chain, Timepoint, Compartment, ClonalFraction) %>%
  pivot_wider(
    names_from  = Compartment,
    values_from = ClonalFraction
  )

# 3) Correlation plot
p_corr <- ggplot(df_wide, aes(x = BoneMarrow, y = cfDNA, color = Patient_plot)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Correlation of Clonal Fraction: BM vs cfDNA",
    x     = "Bone Marrow Clonal Fraction (%)",
    y     = "cfDNA Clonal Fraction (%)"
  ) +
  theme_classic() +
  theme(legend.position = "right")

## Combine with Patchwork 
## Version A: all three plots in one row
layout_all3 <- p_cfDNA | p_BM | p_corr
ggsave(
  filename = file.path(outdir, "VDJ_clonalFrac_all3_sideBySide_updated.png"),
  plot     = layout_all3,
  width    = 15,   # adjust width so each panel is ~5 inches
  height   = 4,
  dpi      = 500
)

## Version B: only cfDNA and BM side‑by‑side
layout_2 <- p_BM | p_cfDNA
ggsave(
  filename = file.path(outdir, "VDJ_clonalFrac_cfDNA_BM_sideBySide_updated.png"),
  plot     = layout_2,
  width    = 10,   # two panels ~5 inches each
  height   = 4,
  dpi      = 500
)

ggsave(
  filename = file.path(outdir, "VDJ_clonalFrac_cfDNA_BM_sideBySide2_updated.png"),
  plot     = layout_2,
  width    = 10,   # two panels ~5 inches each
  height   = 3.5,
  dpi      = 500
)

ggsave(
  filename = file.path(outdir, "VDJ_clonalFrac_cfDNA_BM_sideBySide3_updated.png"),
  plot     = layout_2,
  width    = 10,   # two panels ~5 inches each
  height   = 3,
  dpi      = 500
)

## Export table 
# Export dat_clean
write.csv(dat_clean, file = file.path(outdir, "VDJ_dat_clean_updated.csv"), row.names = FALSE)
saveRDS(dat_clean, file = file.path(outdir, "VDJ_dat_clean_updated.rds"))

# Export dat_to_plot
write.csv(dat_to_plot, file = file.path(outdir, "VDJ_dat_to_plot_updated.csv"), row.names = FALSE)
saveRDS(dat_to_plot, file = file.path(outdir, "VDJ_dat_to_plot_updated.rds"))






# Making Figure 1C - Mutations --------------------------------------------


##### Now process the mutation data 
# ── 1. load the raw sheet (everything as text) ──────────────────────────
# ── libraries ───────────────────────────────────────────────────────────
xlsx_mut <- "/path/Cleaned-table-for-targetedseq-WM-Dory-edited.xlsx"

raw       <- read_excel(xlsx_mut, col_names = FALSE, col_types = "text")

# ─────────────────────────────────────────────────────────────────────────────
# 2) Extract and clean the two header rows
# ─────────────────────────────────────────────────────────────────────────────
pat <- as.character(raw[1, ])    # e.g. "", "Patient", "1", "", "2", ... 
tp  <- as.character(raw[2, ])    # e.g. "", "Timepoint", "Screen", "c7", ...

# fill forward the patient numbers across blanks (only for columns ≥3)
for(i in 3:length(pat)) {
  if(is.na(pat[i]) || pat[i] == "") pat[i] <- pat[i-1]
}

# build new column names: first two are fixed, the rest "PatientX_Timepoint"
new_names <- c(
  "SampleType",
  "Descriptor",
  paste0(
    "Patient", pat[-c(1,2)], "_",
    ifelse(
      str_detect(tp[-c(1,2)], regex("screen", TRUE)),
      "Baseline",
      toupper(tp[-c(1,2)])
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Drop the header rows and assign new_names
# ─────────────────────────────────────────────────────────────────────────────
df_wide <- raw[-c(1,2), ]
colnames(df_wide) <- new_names

# ─────────────────────────────────────────────────────────────────────────────
# 4) Pivot long & split out Patient / Timepoint
# ─────────────────────────────────────────────────────────────────────────────
long <- df_wide %>%
  pivot_longer(
    cols      = -c(SampleType, Descriptor),
    names_to  = "PatientTimepoint",
    values_to = "Value"
  ) %>%
  separate(
    PatientTimepoint,
    into = c("Patient", "Timepoint"),
    sep  = "_"
  ) %>%
  mutate(
    SampleType = if_else(str_detect(SampleType, "^Bone", TRUE),
                         "BoneMarrow", "cfDNA"),
    Patient    = factor(Patient, levels = paste0("Patient", 1:7)),
    # rename “M12” → “C12” 
    Timepoint = case_when(
      Timepoint == "M12" ~ "C12",
      TRUE               ~ as.character(Timepoint)
    ),
    Timepoint  = factor(Timepoint,
                        levels = c("Baseline","C7","C12","M18")),
    Descriptor = str_trim(Descriptor),
    Value      = str_trim(Value)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 5) Within each (SampleType, Patient, Timepoint, Gene), pick Mutation + VAF
# ─────────────────────────────────────────────────────────────────────────────
tidy <- long %>%
  # identify VAF rows vs mutation rows
  mutate(
    isVAF = toupper(Descriptor) == "VAF"
  ) %>%
  # rename Descriptor => Gene to keep the gene name for mutations
  rename(Gene = Descriptor) %>%
  group_by(SampleType, Patient, Timepoint, Gene) %>%
  summarise(
    Mutation = Value[!isVAF][1] %||% NA_character_,
    VAF      = as.numeric(Value[ isVAF ][1])  %||% 0,
    .groups = "drop"
  ) %>%
  select(Patient, Timepoint, SampleType, Gene, Mutation, VAF)

# 1) Split tidy into two tables: real mutations vs the VAF‐rows
mut_rows <- tidy %>%
  filter(!str_ends(Gene, "VAF")) %>%
  rename(Gene = Gene)   # Gene is “MYD88” or “CXCR4”

vaf_rows <- tidy %>%
  filter(str_ends(Gene, "VAF")) %>%
  mutate(
    Gene = str_remove(Gene, "VAF$"),            # drop the VAF suffix
    VAF_num = as.numeric(Mutation) %>% replace_na(0)  # parse the numeric string
  ) %>%
  select(-Mutation)  # we only need VAF_num

# 2) Left‐join them back together on Patient/Timepoint/SampleType/Gene
df_fixed <- mut_rows %>%
  left_join(
    vaf_rows,
    by = c("Patient", "Timepoint", "SampleType", "Gene")
  ) %>%
  # pick columns and rename
  select(
    Patient, Timepoint, SampleType, Gene,
    Mutation,
    VAF = VAF_num
  )

df_fixed <- df_fixed %>%
  mutate(
    # if the called mutation is NA, "NA", or "not identified", drop the VAF back to NA
    VAF = if_else(
      is.na(Mutation) |
        Mutation == "NA" |
        str_to_lower(Mutation) == "not identified",
      NA_real_,     # set VAF to NA in these cases
      VAF       # otherwise keep the numeric VAF
    )
  ) 

# custom 7-colour qualitative palette
pal <- brewer.pal(7, "Dark2")

# Bone Marrow panel
p_BM <- df_fixed %>%
  filter(SampleType == "BoneMarrow") %>%
  ggplot(aes(x = Timepoint, y = VAF, group = Patient, color = Patient)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Gene, ncol = 1, scales = "free_y") +
  scale_color_manual(values = pal, name = "Patient ID") +
  scale_y_continuous(
    labels = percent_format(1),
    expand = expansion(add = c(0.01, 0.05))
  ) +
  labs(
    title = "Bone Marrow Mutation VAF",
    x     = "Treatment Cycle",
    y     = "VAF (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line.x     = element_line(color = "black"),
    axis.ticks      = element_line(color = "black"),
    legend.position = "none",
    plot.title      = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "grey95", color = NA)
  )

# cfDNA panel
p_cf <- df_fixed %>%
  filter(SampleType == "cfDNA") %>%
  ggplot(aes(x = Timepoint, y = VAF, group = Patient, color = Patient)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Gene, ncol = 1, scales = "free_y") +
  scale_color_manual(values = pal, name = "Patient ID") +
  scale_y_continuous(
    labels = percent_format(1),
    expand = expansion(add = c(0.01, 0.05))
  ) +
  labs(
    title = "cfDNA Mutation VAF",
    x     = "Treatment Cycle",
    y     = "VAF (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line.x     = element_line(color = "black"),
    axis.ticks      = element_line(color = "black"),
    legend.position = "right",
    plot.title      = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "grey95", color = NA)
  )

# combine side by side
layout_2 <- p_BM | p_cf

# export
ggsave(
  filename = file.path(outdir, "Mutation_VAF_BM_vs_cfDNA.png"),
  plot     = layout_2,
  width    = 12,
  height   = 4,
  dpi      = 500
)


### Updating to have fixed scale within a specific sample type comparison 
# 1) pre‐compute per–sample‐type VAF maxima
max_bm <- df_fixed %>%
  filter(SampleType == "BoneMarrow") %>%
  pull(VAF) %>%
  max(na.rm = TRUE)


top_break <- ceiling(max_bm * 100) / 100  # round up to next percent, e.g. 0.12 for 12%

max_cf <- df_fixed %>%
  filter(SampleType == "cfDNA") %>%
  pull(VAF) %>%
  max(na.rm = TRUE)

# 2) redraw the two panels, fixing the y–limits
# redraw BM with a little bottom padding
p_BM <- df_fixed %>%
  filter(SampleType == "BoneMarrow") %>%
  ggplot(aes(Timepoint, VAF, group=Patient, color=Patient)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Gene, ncol=1) +
  scale_color_manual(values=pal) +
  scale_y_continuous(
    labels = percent_format(1),
    limits = c(0, max_bm*1.05),
    expand = expansion(add = c(max_bm*0.04, 0))  # 3% of the BM‐max below zero
  ) +
  labs(title="Bone Marrow Mutation VAF", x="Treatment Cycle", y="VAF (%)") +
  theme_classic(base_size=12) +
  theme(
    axis.line.x     = element_line(color="black"),
    axis.ticks      = element_line(color="black"),
    legend.position = "none",
    plot.title      = element_text(face="bold", hjust=0.5),
    strip.background=element_rect(fill="grey95", color=NA)
  )

# 3) redraw cfDNA with similar bottom padding
p_cf <- df_fixed %>%
  filter(SampleType == "cfDNA") %>%
  ggplot(aes(Timepoint, VAF, group=Patient, color=Patient)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Gene, ncol=1) +
  scale_color_manual(values=pal) +
  scale_y_continuous(
    labels = percent_format(1),
    breaks = seq(0, 1, by = 0.25),    # ticks at 0,25,50,75,100
    limits = c(0, 1),
    expand = expansion(add = c(max_cf*0.05, 0))  # 5% of the cfDNA‐max below zero
  ) +
  labs(title="cfDNA Mutation VAF", x="Treatment Cycle", y="VAF (%)") +
  theme_classic(base_size=12) +
  theme(
    axis.line.x     = element_line(color="black"),
    axis.ticks      = element_line(color="black"),
    legend.position = "right",
    plot.title      = element_text(face="bold", hjust=0.5),
    strip.background=element_rect(fill="grey95", color=NA)
  )


# 3) stitch & save
layout_2 <- p_BM | p_cf

ggsave(
  filename = file.path(outdir, "Mutation_VAF_BM_vs_cfDNA_fixed_scales.png"),
  plot     = layout_2,
  width    = 12,
  height   = 4,
  dpi      = 500
)





## Export 
# Export as CSV
write.csv(df_fixed, file = file.path(outdir, "WM_mutation_VAF_table.csv"), row.names = FALSE)

# Export as RDS
saveRDS(df_fixed, file = file.path(outdir, "WM_mutation_VAF_table.rds"))







# Making figure 1D: ClonoSEQ ----------------------------------------------


### Now load in the clonoSEQ
# 1) read  full ClonoSEQ assay5 summary
cls_raw <- read_csv("/path/brawm_blood_evaluable_summary_assay5.csv")

# 2) extract one baseline-per-million value per patient
baseline_summary <- cls_raw %>%
  select(CustomerSubjectId, Scn_PerMillionCount) %>%
  distinct(CustomerSubjectId, .keep_all = TRUE) %>%
  rename(
    Patient_other_ID = CustomerSubjectId,
    per_million      = Scn_PerMillionCount
  ) %>%
  # compute everything you need for a baseline row
  mutate(
    SampleName    = NA_character_,
    Timepoint_raw = "Baseline",
    mrd1e6        = NA_character_,
    Patient = case_when(
      Patient_other_ID == "PATIENT_ID_1" ~ "Patient 1",
      Patient_other_ID == "PATIENT_ID_2"  ~ "Patient 2",
      Patient_other_ID == "PATIENT_ID_3" ~ "Patient 3",
      Patient_other_ID == "PATIENT_ID_4" ~ "Patient 4",
      Patient_other_ID == "PATIENT_ID_5"  ~ "Patient 5",
      Patient_other_ID == "PATIENT_ID_6" ~ "Patient 6",
      Patient_other_ID == "PATIENT_ID_7"  ~ "Patient 7",
      TRUE                          ~ NA_character_
    ) %>% factor(levels = paste("Patient", 1:7)),
    Timepoint = factor("Baseline", levels = c("Baseline","C7","C12","M18","M24")),
    MRD_status = NA_character_,
    freq_pct   = per_million / 1e6 * 100
  ) %>%
  select(Patient_other_ID, SampleName, Timepoint_raw, per_million, mrd1e6,
         Patient, Timepoint, MRD_status, freq_pct)

# 3) now process your real data exactly as before, stopping just before the final filter:
cls <- cls_raw %>%
  rename(
    Patient_other_ID = CustomerSubjectId,
    SampleName       = SampleName,
    Timepoint_raw    = CustomerTimepoint,
    per_million      = PerMillionCount,
    mrd1e6           = Mrd1E6
  ) %>%
  mutate(
    per_million = as.numeric(per_million),
    Patient = case_when(
      Patient_other_ID == "PATIENT_ID_1" ~ "Patient 1",
      Patient_other_ID == "PATIENT_ID_2"  ~ "Patient 2",
      Patient_other_ID == "PATIENT_ID_3" ~ "Patient 3",
      Patient_other_ID == "PATIENT_ID_4" ~ "Patient 4",
      Patient_other_ID == "PATIENT_ID_5"  ~ "Patient 5",
      Patient_other_ID == "PATIENT_ID_6" ~ "Patient 6",
      Patient_other_ID == "PATIENT_ID_7"  ~ "Patient 7"
    ) %>% factor(levels = paste("Patient", 1:7)),
    Timepoint = case_when(
      str_detect(Timepoint_raw, regex("screen", TRUE))        ~ "Baseline",
      str_detect(Timepoint_raw, regex("cycle[_ ]*7", TRUE))   ~ "C7",
      str_detect(Timepoint_raw, regex("cycle[_ ]*12", TRUE))  ~ "C12",
      str_detect(Timepoint_raw, regex("month[_ ]*18", TRUE))  ~ "M18",
      str_detect(Timepoint_raw, regex("month[_ ]*24", TRUE))  ~ "M24",
      str_detect(Timepoint_raw, regex("18\\s*Months?", TRUE)) ~ "M18",
      str_detect(Timepoint_raw, regex("24\\s*Months?", TRUE)) ~ "M24",
      TRUE                                                    ~ Timepoint_raw
    ) %>% factor(levels = c("Baseline","C7","C12","M18","M24")),
    MRD_status = case_when(
      str_to_lower(mrd1e6) == "positive"       ~ "Positive",
      str_to_lower(mrd1e6) == "negative"       ~ "Negative",
      str_detect(str_to_lower(mrd1e6),"indet") ~ "Indeterminate",
      TRUE                                     ~ NA_character_
    ),
    freq_pct = per_million / 1e6 * 100
  ) %>%
  filter(!is.na(Patient)) %>%
  select(Patient_other_ID, SampleName, Timepoint_raw, per_million, mrd1e6,
         Patient, Timepoint, MRD_status, freq_pct)

# 4) bind the baseline rows on top
cls_final <- bind_rows(baseline_summary, cls) %>%
  arrange(Patient, Timepoint)

cls_final <- cls_final %>% filter(!is.na(Patient))

## Now get baseline to be filled 
cls <- cls_final %>%
  tidyr::replace_na(list(MRD_status = "Positive")) 

pal   <- brewer.pal(7, "Dark2")
shps  <- c(Positive=16, Negative=15, Indeterminate=17)

# --- 3) Linear‐scale plot (p1) ---
p1 <- ggplot(cls, aes(Timepoint, freq_pct, group=Patient, color=Patient)) +
  geom_line(size=1) +
  geom_point(aes(shape=MRD_status), size=3) +
  scale_color_manual(values=pal, name="Patient ID") +
  scale_shape_manual(values=shps, name="MRD status") +
  scale_y_continuous(
    labels = label_number(suffix="%"),
    expand = expansion(add = c(0.5, 2))      # add 0.5% below zero, 2% above max
  ) +
  labs(
    title = "ClonoSEQ Frequency Over Treatment",
    x     = "Treatment Cycle",
    y     = "ClonoSEQ freq. (%)"
  ) +
  theme_classic(base_size=12) +
  theme(
    # axis.text.x     = element_text(angle=45, hjust=1),
    legend.position = "none",
    plot.title      = element_text(face="bold", hjust=0.5)
  )

# --- 4) Log‐scale version (p2) ---
# 1) Prepare data: drop Baseline and map zeros → 1e‑6
cls2 <- cls %>%
  filter(Timepoint != "Baseline") %>%
  mutate(
    freq_pct_plot = if_else(freq_pct == 0, 1e-6, freq_pct)
  )

# 2) Define your color + shape palettes (reuse from before)
pal  <- brewer.pal(7, "Dark2")
shps <- c(Positive=16, Negative=15, Indeterminate=17)

# 3) Build the log plot from scratch
p2_new <- ggplot(cls2, 
                 aes(x = Timepoint, 
                     y = freq_pct_plot, 
                     group = Patient, 
                     color = Patient)) +
  geom_line(size = 1) +
  geom_point(aes(shape = MRD_status), size = 3) +
  
  # 4) Custom log scale: breaks and labels
  scale_y_log10(
    breaks = c(1e-6,1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1),
    labels = c("Not detected", "0.00001%", "0.0001%","0.001%", "0.01%", "0.1%", "1%", "10%"),
    expand = expansion(add = c(0.1, 0.2))
  ) +
  annotation_logticks(sides = "l") +
  
  # 5) Colors, shapes, and theme
  scale_color_manual(values = pal, name = "Patient ID") +
  scale_shape_manual(values = shps, name = "MRD status") +
  labs(
    title = "ClonoSEQ Frequency at MRD Timepoints [log]",
    x     = "Treatment Cycle",
    y     = "ClonoSEQ freq. (%) [log₁₀]"
  ) +
  theme_classic(base_size = 12) +
  theme(
    #    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.title      = element_text(face = "bold", hjust = 0.5)
  )

# Save both ---
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

ggsave(file.path(outdir, "ClonoSEQ_freq_linear_long.png"), p1, width=11, height=3.75, dpi=500)
ggsave(file.path(outdir, "ClonoSEQ_freq_linear_narrow.png"), p1, width=4.5, height=3.75, dpi=500)

ggsave(file.path(outdir, "ClonoSEQ_freq_log2.png"),    p2, width=10, height=6, dpi=500)
ggsave(file.path(outdir, "ClonoSEQ_freq_log2_new_narrow.png"),    p2_new, width=7, height=3.5, dpi=750)
ggsave(file.path(outdir, "ClonoSEQ_freq_log2_new_narrow2.png"),    p2_new, width=4.5, height=3.75, dpi=750)


# 1) Assemble p1 and p2_new side by side
clono_layout <- p1 + p2_new + 
  plot_layout(ncol = 2, widths = c(1, 1)) & 
  theme(
    plot.margin = margin(5, 5, 5, 5)
  )

# 2) Export at narrow aspect (two panels ~4.5" each)
ggsave(
  filename = file.path(outdir, "ClonoSEQ_freq_linear_and_log_narrow.png"),
  plot     = clono_layout,
  width    = 12,      # two 4.5" panels
  height   = 4,
  dpi      = 500
)





# Making Figure 2 - Scatterplot Matrix  -----------------------------------



#### Now load in the BM clonoseq data for use in Figure 2
# 1) read your full ClonoSEQ assay5 summary
cls_raw_bm <- read_csv("/path/brawm_bone_marrow_evaluable_summary_assay5.csv")

# 2) extract one baseline-per-million value per patient
baseline_summary <- cls_raw_bm %>%
  select(CustomerSubjectId, Scn_PerMillionCount) %>%
  distinct(CustomerSubjectId, .keep_all = TRUE) %>%
  rename(
    Patient_other_ID = CustomerSubjectId,
    per_million      = Scn_PerMillionCount
  ) %>%
  # compute everything you need for a baseline row
  mutate(
    SampleName    = NA_character_,
    Timepoint_raw = "Baseline",
    mrd1e6        = NA_character_,
    Patient = case_when(
      Patient_other_ID == "PATIENT_ID_1" ~ "Patient 1",
      Patient_other_ID == "PATIENT_ID_2"  ~ "Patient 2",
      Patient_other_ID == "PATIENT_ID_3" ~ "Patient 3",
      Patient_other_ID == "PATIENT_ID_4" ~ "Patient 4",
      Patient_other_ID == "PATIENT_ID_5"  ~ "Patient 5",
      Patient_other_ID == "PATIENT_ID_6" ~ "Patient 6",
      Patient_other_ID == "PATIENT_ID_7"  ~ "Patient 7",
      TRUE                          ~ NA_character_
    ) %>% factor(levels = paste("Patient", 1:7)),
    Timepoint = factor("Baseline", levels = c("Baseline","C7","C12","M18","M24")),
    MRD_status = NA_character_,
    freq_pct   = per_million / 1e6 * 100
  ) %>%
  select(Patient_other_ID, SampleName, Timepoint_raw, per_million, mrd1e6,
         Patient, Timepoint, MRD_status, freq_pct)

# 3) now process your real data exactly as before, stopping just before the final filter:
tmp <- cls_raw_bm %>%
  rename(
    Patient_other_ID = CustomerSubjectId,
    SampleName       = SampleName,
    Timepoint_raw    = CustomerTimepoint,
    per_million      = PerMillionCount,
    mrd1e6           = Mrd1E6
  ) %>%
  mutate(
    per_million = as.numeric(per_million),
    Patient = case_when(
      Patient_other_ID == "PATIENT_ID_1" ~ "Patient 1",
      Patient_other_ID == "PATIENT_ID_2"  ~ "Patient 2",
      Patient_other_ID == "PATIENT_ID_3" ~ "Patient 3",
      Patient_other_ID == "PATIENT_ID_4" ~ "Patient 4",
      Patient_other_ID == "PATIENT_ID_5"  ~ "Patient 5",
      Patient_other_ID == "PATIENT_ID_6" ~ "Patient 6",
      Patient_other_ID == "PATIENT_ID_7"  ~ "Patient 7"
    ) %>% factor(levels = paste("Patient", 1:7)),
    Timepoint = case_when(
      str_detect(Timepoint_raw, regex("screen", TRUE))        ~ "Baseline",
      str_detect(Timepoint_raw, regex("cycle[_ ]*7", TRUE))   ~ "C7",
      str_detect(Timepoint_raw, regex("cycle[_ ]*12", TRUE))  ~ "C12",
      str_detect(Timepoint_raw, regex("month[_ ]*18", TRUE))  ~ "M18",
      str_detect(Timepoint_raw, regex("month[_ ]*24", TRUE))  ~ "M24",
      str_detect(Timepoint_raw, regex("18\\s*Months?", TRUE)) ~ "M18",
      str_detect(Timepoint_raw, regex("24\\s*Months?", TRUE)) ~ "M24",
      TRUE                                                    ~ Timepoint_raw
    ) %>% factor(levels = c("Baseline","C7","C12","M18","M24")),
    MRD_status = case_when(
      str_to_lower(mrd1e6) == "positive"       ~ "Positive",
      str_to_lower(mrd1e6) == "negative"       ~ "Negative",
      str_detect(str_to_lower(mrd1e6),"indet") ~ "Indeterminate",
      TRUE                                     ~ NA_character_
    ),
    freq_pct = per_million / 1e6 * 100
  ) %>%
  filter(!is.na(Patient)) %>%
  select(Patient_other_ID, SampleName, Timepoint_raw, per_million, mrd1e6,
         Patient, Timepoint, MRD_status, freq_pct)

# 4) bind the baseline rows on top
cls_final_bm <- bind_rows(baseline_summary, tmp) %>%
  arrange(Patient, Timepoint)

cls_final_bm <- cls_final_bm %>% filter(!is.na(Patient))

cls_bm <- cls_final_bm %>%
  tidyr::replace_na(list(MRD_status = "Positive")) 









#### Now make correlation matrix with everyting together 
# 1a) MRDetect: rename timepoint to factor & Patient to match
md <- df_combined %>%
  # overwrite Patient, then drop the helper column
  mutate(
    Patient   = Patient_number_plot
  ) %>%
  select(-Patient_number_plot, -Timepoint) %>%
  
  # rename timepoint_info_plot to Timepoint
  rename(
    Timepoint = timepoint_info_plot
  ) %>%
  
  # drop healthy controls
  filter(Timepoint != "Healthy_controls") %>%
  
  # make Timepoint a clean factor
  mutate(
    Timepoint = factor(Timepoint,
                       levels = c("Baseline","Cycle C7","Cycle C12","Cycle C18"))
  ) %>%
  
  # pick exactly the columns you want
  select(
    Patient,
    Timepoint,
    
    sites_checked,
    reads_checked,
    sites_detected,
    reads_detected,
    total_reads,
    detection_rate,
    detection_rate_as_reads_detected_over_reads_checked,
    detection_rate_as_reads_detected_over_total_reads,
    sites_detection_rate,
    
    detection_rate_zscore_charm,
    sites_rate_zscore_charm,
    detection_rate_zscore_reads_checked_charm,
    detection_rate_zscore_total_reads_charm,
    
    Mrd_by_WGS
  )

# 1b) VDJ: pivot your df_wide so you get bone & plasma fractions
df_wide <- dat_to_plot %>%
  filter(Compartment %in% c("cfDNA", "BoneMarrow")) %>%
  select(Patient_plot, Chain, Timepoint, Compartment, ClonalFraction) %>%
  pivot_wider(
    names_from  = Compartment,
    values_from = ClonalFraction
  )

vdj <- df_wide %>%
  rename(
    Patient = Patient_plot,
    Timepoint = Timepoint
  ) %>%
  pivot_longer(
    cols=c(BoneMarrow, cfDNA),
    names_to="Compartment",
    values_to="VDJ_frac"
  ) %>%
  pivot_wider(
    names_from=Compartment,
    values_from=VDJ_frac,
    names_prefix="VDJ_"
  ) %>%
  select(Patient, Timepoint, VDJ_BoneMarrow, VDJ_cfDNA)

# 1c) Mutations: each row is one Gene; pivot it wide so MYD88/CXCR4 get their own columns
mut <- df_fixed %>%
  select(Patient,Timepoint,SampleType,Gene,VAF) %>%
  pivot_wider(
    names_from = c(SampleType,Gene),
    values_from = VAF,
    names_sep = "_",
    values_fill = NA_real_
  )
# this gives e.g. BoneMarrow_MYD88, BoneMarrow_CXCR4, cfDNA_MYD88, cfDNA_CXCR4

# 1d) ClonoSEQ: we already have freq_pct + MRD_status
cs <- cls %>%
  select(Patient,Timepoint,freq_pct,MRD_status) %>%
  rename(
    ClonoSEQ_pct  = freq_pct,
    ClonoSEQ_MRD  = MRD_status
  )

cs_bm <-  cls_bm %>%
  select(Patient,Timepoint,freq_pct,MRD_status) %>%
  rename(
    ClonoSEQ_pct  = freq_pct,
    ClonoSEQ_MRD  = MRD_status
  )

# 2) Full join them all together
clean_labels <- function(df) {
  df %>%
    mutate(
      #––– Fix Patient: ensure a space before the digit and factorize –––
      Patient = str_replace(as.character(Patient),
                            "^\\s*Patient\\s*(\\d+)\\s*$",
                            "Patient \\1"),
      Patient = factor(Patient, levels = paste("Patient", 1:7)),
      
      #––– Now clean Timepoint exactly the same way –––
      Timepoint = str_trim(Timepoint),
      Timepoint = str_remove(Timepoint,
                             regex("^Cycle\\s+", ignore_case = TRUE)),
      Timepoint = str_replace(Timepoint,
                              "^M(\\d+)$",
                              "C\\1"),
      Timepoint = if_else(str_to_lower(Timepoint) == "screening",
                          "Baseline",
                          Timepoint),
      Timepoint = str_to_title(Timepoint),
      Timepoint = factor(Timepoint,
                         levels = c("Baseline","C7","C12","C18","C24"))
    )
}

# Apply to each of your tables:
md_clean  <- clean_labels(md)
vdj_clean <- clean_labels(vdj)
mut_clean <- clean_labels(mut)
cs_clean  <- clean_labels(cs)
cs_clean_bm  <- clean_labels(cs_bm)

cs_clean <- cs_clean %>% 
  rename(
    ClonoSEQ_pct_cfDNA = ClonoSEQ_pct,
    ClonoSEQ_MRD_cfDNA = ClonoSEQ_MRD
  )

cs_clean_bm <- cs_clean_bm %>% 
  rename(
    ClonoSEQ_pct_BM = ClonoSEQ_pct,
    ClonoSEQ_MRD_BM = ClonoSEQ_MRD
  )

all_metrics <- list(
  md_clean,
  vdj_clean,
  mut_clean,
  cs_clean,
  cs_clean_bm
) %>% 
  reduce(full_join, by = c("Patient","Timepoint"))


all_metrics <- all_metrics %>%
  mutate(
    BM_cells_target_VAF   = pmax(BoneMarrow_MYD88,
                                 BoneMarrow_CXCR4,
                                 na.rm = TRUE),
    cfDNA_target_VAF      = pmax(cfDNA_MYD88,
                                 cfDNA_CXCR4,
                                 na.rm = TRUE)
  )


#### Now make figure 

# 1) Define palette
pal <- brewer.pal(7, "Dark2") # Same as earlier 

# Pick Metrics + Patient, renaming the detection rate
metrics_df <- all_metrics %>%
  # compute % detection rate
  mutate(
    cfWGS_pct = detection_rate_as_reads_detected_over_reads_checked * 100
  ) %>%
  # pick just the columns we need
  select(
    Patient,
    cfWGS_pct,
    VDJ_BoneMarrow,
    VDJ_cfDNA,
    ClonoSEQ_pct_BM,
    ClonoSEQ_pct_cfDNA,
    BM_cells_target_VAF,
    cfDNA_target_VAF
  ) %>%
  # rename them to a consistent short form
  rename(
    cfWGS        = cfWGS_pct,
    VDJ_BM       = VDJ_BoneMarrow,
    VDJ_cfDNA    = VDJ_cfDNA,
    ClonoSEQ_BM  = ClonoSEQ_pct_BM,
    ClonoSEQ_cfDNA = ClonoSEQ_pct_cfDNA,
    VAF_BM       = BM_cells_target_VAF,
    VAF_cfDNA    = cfDNA_target_VAF
  ) %>%
  # and make sure Patient has the right factor levels
  mutate(
    Patient = factor(Patient, levels = paste("Patient", 1:7))
  )


# metrics_df <- metrics_df %>% filter(!is.na(VDJ_BM)) # Remove ones missing many lab tests

# patient colours
patient_cols <- brewer.pal(7, "Dark2")
names(patient_cols) <- paste("Patient", 1:7)

# -- 1) Make safe spearman tests wrapper 
safe_spearman <- function(data, mapping, method="spearman",
                          size=3, colour="black", ...) {
  x_var <- as_label(mapping$x)
  y_var <- as_label(mapping$y)
  x <- as.numeric(data[[x_var]])
  y <- as.numeric(data[[y_var]])
  
  # only complete pairs
  valid_idx <- which(!is.na(x) & !is.na(y))
  n_pairs   <- length(valid_idx)
  
  # too few pairs
  if(n_pairs < 2) {
    return(
      ggally_text(
        label   = paste0("n=", n_pairs),
        mapping = aes(), colour = colour, size = size
      ) + theme_void()
    )
  }
  
  # zero variance in either vector
  if(var(x[valid_idx]) == 0 || var(y[valid_idx]) == 0) {
    return(
      ggally_text(
        label   = paste0("n=", n_pairs, "\nρ undefined"),
        mapping = aes(), colour = colour, size = size
      ) + theme_void()
    )
  }
  
  # otherwise compute Spearman
  ct   <- cor.test(x[valid_idx], y[valid_idx], method = method)
  rho  <- round(ct$estimate, 2)
  stars <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols   = c("***", "**", "*", "'", " ")
  )
  lbl  <- paste0("n=", n_pairs, "\nρ=", rho, stars)
  
  ggally_text(
    label   = lbl,
    mapping = aes(), colour = colour, size = size
  ) + theme_void()
}

# -- 2) A custom lower‐panel that draws points and a box around each cell --
lower_with_border <- function(data, mapping, ...) {
  ggally_points(data, mapping, ...) +
    theme_minimal() +
    theme(
      panel.border     = element_rect(color="black", fill=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# -- 3) A simple diagonal (all‐data) histogram --
diag_hist <- function(data, mapping, ...) {
  x_var <- as_label(mapping$x)
  ggplot(data, aes(x=.data[[x_var]])) +
    geom_histogram(bins=15, fill="grey80", color="black") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x  = element_blank()
    )
}

# -- 4) Assemble  pairwise plot --
corr_plot <- ggpairs(
  metrics_df,
  columns = 2:8,                        # numeric metrics
  mapping = aes(color = Patient),
  upper   = list(continuous = safe_spearman),
  lower   = list(continuous = wrap(lower_with_border,
                                   size=1.5, alpha=0.6)),  
  diag  = list(continuous = diag_hist),
  axisLabels = "show",
  progress   = FALSE
)  +
  scale_color_manual(values = patient_cols) +
  theme(
    strip.text      = element_text(size=8),
    legend.position = "none"
  ) 

# -- 5) Save it --
ggsave(
  filename = file.path(outdir, "AllMetrics_correlation_matrix_fixed_updated6.png"),
  plot     = corr_plot,
  width    = 10,
  height   = 10,
  dpi      = 500
)






# Performing Additional Analysis - Summarizing Mutations Detected by WGS  -------------------------


#### Summarise the mutation info - amount called at each step pre and post filtering
# 1) Read the Excel sheet
df_raw <- read_excel("/path/Mutation_counts.xlsx", sheet = 1, col_types = "text") %>%
  # ensure Mutation_Count is numeric
  mutate(Mutation_Count = as.integer(Mutation_Count))

# 2) Extract patient IDs from the filename and recode
df <- df_raw %>%
  mutate(
    # 1) extract either "WAL-1234-5" or "1234-5"
    raw_id = str_extract(VCF_File, "(?:WAL-)?[0-9]+-[0-9]+"),
    # 2) remove any leading "WAL-" so raw_id2 is always "1234-5"
    raw_id2 = str_remove(raw_id, "^WAL-"),
    # 3) rebuild a canonical patient_other_ID = "WAL-1234-5"
    patient_other_ID = paste0("WAL-", raw_id2),
    # 4) map to Patient 1–7
    Patient = case_when(
      patient_other_ID == "WAL-PATIENT_ID_1" ~ "Patient 1",
      patient_other_ID == "WAL-PATIENT_ID_2"  ~ "Patient 2",
      patient_other_ID == "WAL-PATIENT_ID_3" ~ "Patient 3",
      patient_other_ID == "WAL-PATIENT_ID_4" ~ "Patient 4",
      patient_other_ID == "WAL-PATIENT_ID_5"  ~ "Patient 5",
      patient_other_ID == "WAL-PATIENT_ID_6" ~ "Patient 6",
      patient_other_ID == "WAL-PATIENT_ID_7"  ~ "Patient 7",
      TRUE                             ~ NA_character_
    )
  ) %>%
  select(-raw_id, -raw_id2)


# 3) Compute per‑stage summary stats, now including min & max (and a human‑readable range)
# A) Compute per‐stage summary stats, including ranges + absolute differences
stats <- df %>%
  group_by(Stage) %>%
  summarize(
    mean_muts = mean(Mutation_Count),
    sd_muts   = sd(Mutation_Count),
    min_muts  = min(Mutation_Count),
    max_muts  = max(Mutation_Count),
    .groups = "drop"
  ) %>%
  arrange(factor(Stage, levels = c("10_VAF", "STR", "Post_RS"))) %>%
  mutate(
    range_muts = paste0(min_muts, "–", max_muts),
    prev_mean  = lag(mean_muts),
    diff_muts  = mean_muts - prev_mean
  )

# B) Glue it all together into one concise sentence
cat(glue(
  "After filtering for somatic variants with VAF ≥10% (\"10_VAF\"), the mean mutation count was ",
  "{round(stats$mean_muts[1],1)} (SD {round(stats$sd_muts[1],1)}, range {stats$range_muts[1]}). ",
  "Exclusion of known polymorphisms and simple repeats (\"STR\" stage) ",
  "reduced the mean by {round(stats$diff_muts[2],1)} to ",
  "{round(stats$mean_muts[2],1)} (SD {round(stats$sd_muts[2],1)}, range {stats$range_muts[2]}). ",
  "Finally, removal of dbSNP‐flagged variants (\"Post_RS\") ",
  "further decreased the mean by {round(stats$diff_muts[3],1)} to ",
  "{round(stats$mean_muts[3],1)} (SD {round(stats$sd_muts[3],1)}, range {stats$range_muts[3]})."
), "\n")


# C) Save the final per-patient, per-stage data and the summary stats
write_csv(
  df,
  file.path(outdir, "mutation_counts_per_patient_per_stage.csv")
)
write_csv(
  stats,
  file.path(outdir, "mutation_counts_stage_summary.csv")
)

message("Exported plots and tables to:\n  ", outdir)








# Performing Additional Statistical Comparisons for MS Text ---------------------------


### First get the number of mutations detected and range at baseline 
# 1. Filter to the 6 baseline samples
baseline_df <- df_combined %>%
  filter(timepoint_info_plot == "Baseline")

# 2. Sites detected: mean, median, sd and range
sites_stats <- baseline_df %>%
  summarise(
    mean_sites   = mean(sites_detected),
    median_sites = median(sites_detected),
    sd_sites     = sd(sites_detected),
    min_sites    = min(sites_detected),
    max_sites    = max(sites_detected)
  )

# 3. Detection‐rate (reads_detected / reads_checked): mean, median, sd and range, expressed as %
rate_stats <- baseline_df %>%
  summarise(
    mean_rate   = mean(detection_rate_as_reads_detected_over_reads_checked) * 100,
    median_rate = median(detection_rate_as_reads_detected_over_reads_checked) * 100,
    sd_rate     = sd(detection_rate_as_reads_detected_over_reads_checked)   * 100,
    min_rate    = min(detection_rate_as_reads_detected_over_reads_checked)  * 100,
    max_rate    = max(detection_rate_as_reads_detected_over_reads_checked)  * 100
  )

# 4. (Optional) Combine into a single table for easy reporting
summary_table <- bind_rows(
  sites_stats %>% mutate(metric = "sites_detected", across(everything(), ~ .)),
  rate_stats  %>% mutate(metric = "tumour_fraction_%", across(everything(), ~ .))
) %>%
  select(metric, everything())

# Print out
sites_stats
rate_stats
summary_table


### Now get comparison to MYD88 at baseline
# 1. Restrict to baseline
baseline <- all_metrics %>%
  filter(Timepoint == "Baseline")

# 2. Define the four comparisons you care about
corr_pairs <- tibble::tibble(
  comparison = c("MYD88 (BM)"      ,
                 "MYD88 (cfDNA)"   ,
                 "clonoSEQ (BM)"   ,
                 "clonoSEQ (PB)")  ,
  x           = "detection_rate_as_reads_detected_over_reads_checked",
  y           = c("BoneMarrow_MYD88",
                  "cfDNA_MYD88"     ,
                  "ClonoSEQ_pct_BM" ,
                  "ClonoSEQ_pct_cfDNA")
)

# 3. A helper that drops NAs and only runs cor.test if there are >=3 points
safe_spearman2 <- function(x, y) {
  df <- tibble::tibble(x = x, y = y) %>% filter(!is.na(x) & !is.na(y))
  if(nrow(df) < 3) {
    return(tibble::tibble(r = NA_real_, p = NA_real_, n = nrow(df)))
  }
  test <- cor.test(df$x, df$y, method = "spearman")
  tibble::tibble(
    r = unname(test$estimate),
    p = test$p.value,
    n = nrow(df)
  )
}

# 4. Run them all and tidy up
results <- corr_pairs %>%
  mutate(
    out = map2(x, y, ~ safe_spearman2(
      baseline[[.x]], baseline[[.y]]
    )),
    out = map(out, ~ .x %>%
                mutate(
                  r = round(r, 2),
                  p = if_else(p < 0.001, "<0.001", as.character(round(p, 3)))
                ))
  ) %>%
  select(comparison, out) %>%
  unnest(out)

results


### Now compare to the followup samples
c7_summary <- all_metrics %>%
  filter(Timepoint == "C7") %>%
  # compute % tumour fraction
  mutate(tf_pct = detection_rate_as_reads_detected_over_reads_checked * 100) %>%
  # split positive vs negative
  group_by(Mrd_by_WGS) %>%
  summarise(
    patients = paste(unique(Patient), collapse = ", "),
    sites    = paste(sites_detected, collapse = " & "),
    tf       = paste0(round(tf_pct, 2), collapse = "% & "),
    n        = n()
  ) %>%
  ungroup()

consistency <- all_metrics %>%
  filter(Timepoint == "C7", Mrd_by_WGS == "Negative") %>%
  transmute(
    Patient,
    panel_BM    = if_else(BoneMarrow_MYD88 > 0,   "Positive", "Negative"),
    panel_cfDNA = if_else(cfDNA_MYD88 > 0,        "Positive", "Negative"),
    clono_BM    = ClonoSEQ_MRD_BM,
    clono_PB    = ClonoSEQ_MRD_cfDNA,
    vdj_BM      = if_else(VDJ_BoneMarrow > 0,     "Positive", "Negative"),
    vdj_cfDNA   = if_else(VDJ_cfDNA > 0,          "Positive", "Negative")
  ) %>%
  summarise(
    n                   = n(),
    panel_concordance   = sum(panel_BM    == "Negative" & panel_cfDNA == "Negative"),
    clono_concordance   = sum(clono_BM    == "Negative" & clono_PB    == "Negative"),
    vdj_concordance     = sum(vdj_BM      == "Negative" & vdj_cfDNA   == "Negative")
  )






# Export the all metrics table to be included as supplemental info --------
### Export all metrics 
write_csv(all_metrics, "all_metrics.csv")
write_tsv(all_metrics, "all_metrics.tsv")
saveRDS(all_metrics, "all_metrics.rds")

# HTML table — for embedding in an R Markdown report
kable(all_metrics, format = "html", table.attr = "class='table table-striped'")


# ── Save session info ─────────────────────────────────────────────────────
writeLines(
  capture.output(sessionInfo()),
  file = "sessionInfo.txt"
)

# ── Done! ─────────────────────────────────────────────────────────────────
message("All metrics exported and session info saved. Script complete.")
