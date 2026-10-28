# paper_fig_WL.R
#
# Assembles the multi-panel manuscript figures (PNG and 300 dpi LZW TIFF) from the individual panels.
# Run after all figure scripts  |  Output: plots1/paper_fig/

getwd()
library(magick)
if (!requireNamespace("pdftools", quietly = TRUE)) install.packages("pdftools")
library(pdftools)

# REUSABLE FUNCTIONS (define once, reuse for every figure)
prep_panel_unlabeled <- function(path, top_pad = 120, border_size = "20x20", pdf_density = 300) {
  if (grepl("\\.pdf$", path, ignore.case = TRUE)) {
    img <- image_read_pdf(path, density = pdf_density)
  } else {
    img <- image_read(path)
  }
  img <- tryCatch(image_trim(img), error = function(e) img)
  img <- image_border(img, "white", border_size)
  
  info <- image_info(img)
  w <- info$width; h <- info$height
  
  # Reserve blank space on top for the label -- but don't write it yet
  canvas <- image_blank(width = w, height = h + top_pad, color = "white")
  image_composite(canvas, img, offset = paste0("+0+", top_pad))
}

# Add a height_weights parameter to give some rows extra room
combine_figure <- function(row_paths, row_labels, output_basename,
                           output_dir = "plots1/paper_fig",
                           target_height = 1320,
                           label_size = 85,
                           top_pad = 120,
                           quantize_png = TRUE,
                           outer_border = "40x40",
                           row_gap = 0,
                           height_weights = NULL) {   # <-- NEW: e.g. c(1, 1, 2) to make row 3 twice as tall
  
  stopifnot(length(row_paths) == length(row_labels))
  if (is.null(height_weights)) height_weights <- rep(1, length(row_paths))
  stopifnot(length(height_weights) == length(row_paths))
  
  prep_panel_unlabeled <- function(path, th, top_pad, border_size = "20x20", pdf_density = 300) {
    if (grepl("\\.pdf$", path, ignore.case = TRUE)) {
      img <- image_read_pdf(path, density = pdf_density)
    } else {
      img <- image_read(path)
    }
    img <- tryCatch(image_trim(img), error = function(e) img)
    img <- image_border(img, "white", border_size)
    img <- image_resize(img, paste0("x", th))
    info <- image_info(img)
    canvas <- image_blank(width = info$width, height = info$height + top_pad, color = "white")
    image_composite(canvas, img, offset = paste0("+0+", top_pad))
  }
  
  rows_raw <- list()
  for (i in seq_along(row_paths)) {
    row_target_height <- round(target_height * height_weights[i])  # scale this row's base height by its weight
    imgs <- lapply(row_paths[[i]], function(p) prep_panel_unlabeled(p, row_target_height, top_pad))
    
    if (row_gap > 0 && length(imgs) > 1) {
      spacer <- image_blank(width = row_gap, height = image_info(imgs[[1]])$height, color = "white")
      imgs_with_gaps <- list(imgs[[1]])
      for (j in 2:length(imgs)) imgs_with_gaps <- c(imgs_with_gaps, list(spacer), list(imgs[[j]]))
      imgs <- imgs_with_gaps
    }
    
    rows_raw[[i]] <- image_append(do.call(c, imgs))
  }
  
  widths <- sapply(rows_raw, function(r) image_info(r)$width)
  max_w <- max(widths)
  rows_scaled <- lapply(rows_raw, function(r) image_resize(r, paste0(max_w, "x")))
  
  panel_positions <- list()
  y_cursor <- 0
  for (i in seq_along(row_paths)) {
    scale_factor <- max_w / widths[i]
    row_target_height <- round(target_height * height_weights[i])
    imgs_unscaled <- lapply(row_paths[[i]], function(p) prep_panel_unlabeled(p, row_target_height, top_pad))
    panel_widths <- sapply(imgs_unscaled, function(im) image_info(im)$width)
    panel_widths_scaled <- panel_widths * scale_factor
    
    if (row_gap > 0 && length(panel_widths_scaled) > 1) {
      gap_scaled <- row_gap * scale_factor
      x_offsets <- cumsum(c(0, panel_widths_scaled[-length(panel_widths_scaled)] + gap_scaled))
    } else {
      x_offsets <- cumsum(c(0, panel_widths_scaled[-length(panel_widths_scaled)]))
    }
    
    panel_positions[[i]] <- data.frame(x = x_offsets, y = y_cursor, label = row_labels[[i]], stringsAsFactors = FALSE)
    y_cursor <- y_cursor + image_info(rows_scaled[[i]])$height
  }
  
  figure_final <- image_append(do.call(c, rows_scaled), stack = TRUE)
  
  all_positions <- do.call(rbind, panel_positions)
  for (k in seq_len(nrow(all_positions))) {
    figure_final <- image_annotate(
      figure_final, all_positions$label[k],
      size = label_size, weight = 700, color = "black",
      location = paste0("+", round(all_positions$x[k]) + 15, "+", round(all_positions$y[k]) + 20)
    )
  }
  
  figure_final <- image_border(figure_final, "white", outer_border)
  figure_final_flat <- image_flatten(image_background(figure_final, "white"))
  figure_final_flat <- image_background(figure_final_flat, "white", flatten = TRUE)
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  png_path  <- file.path(output_dir, paste0(output_basename, ".png"))
  tiff_path <- file.path(output_dir, paste0(output_basename, ".tiff"))
  png_out <- if (quantize_png) image_quantize(figure_final_flat, max = 256, colorspace = "sRGB") else figure_final_flat
  
  image_write(png_out, png_path, format = "png", density = "300x300")
  image_write(figure_final_flat, tiff_path, format = "tiff", density = "300x300", compression = "LZW")
  
  cat("Saved:", png_path, "\nSaved:", tiff_path, "\n")
  print(figure_final_flat)
  invisible(figure_final_flat)
}

# Figure 1: Derivation of the 9-gene stromal-associated prognostic signature
#   A - Volcano plot (stromal DEGs)
#   B - LASSO cross-validation curve
#   C - LASSO coefficient path
#   D - Final 9-gene signature bar plot
path_volcano    <- "plots1/Stromal-Associated Differentially Expressed Genes.png"
path_lasso_cv   <- "plots1/Cross Validation Curve lasso.png"
path_lasso_path <- "plots1/LASSO Coefficient Path.png"
path_barplot    <- "plots1/lasso gene Prognostic Signature.png"

combine_figure(
  row_paths  = list(c(path_volcano, path_lasso_cv), c(path_lasso_path, path_barplot)),
  row_labels = list(c("A", "B"), c("C", "D")),
  output_basename = "Figure1",
  output_dir = "plots1/paper_fig"
)

# Figure 2: Discrimination of the risk score across cohorts
path_risk_dist   <- "plots1/Risk Score Distr.png"       # Panel A
path_surv_status <- "plots1/survival status.png"        # Panel B
path_km_train    <- "plots1/Train KM.png"                # Panel C
path_km_test     <- "plots1/Test KM.png"                 # Panel D
path_km_full     <- "plots1/Full KM.png"                 # Panel E
path_auc_train   <- "plots1/Train AUROC.png"             # Panel F
path_auc_test    <- "plots1/Test AUROC.png"              # Panel G
path_auc_full    <- "plots1/Full AUROC.png"              # Panel H

combine_figure(
  row_paths  = list(
    c(path_risk_dist, path_surv_status),
    c(path_km_train, path_km_test, path_km_full),
    c(path_auc_train, path_auc_test, path_auc_full)
  ),
  row_labels = list(
    c("A", "B"),
    c("C", "D", "E"),
    c("F", "G", "H")
  ),
  output_basename = "Figure2",
  output_dir = "plots1/paper_fig"
)

# Figure 3: Multivariable Cox models
#   A - Age + stage model
#   B - Subtype-adjusted model
path_multicox <- "plots1/multicox subtype-adjusted.png"    # Panel B
path_multicox_age_stage <- "plots1/multicox age+stage.png"      # Panel A

combine_figure(
  row_paths  = list(
    c(path_multicox_age_stage, path_multicox)
  ),
  row_labels = list(
    c("A", "B")
  ),
  output_basename = "Figure3",
  output_dir = "plots1/paper_fig"
)

# FIGURE 4 -- Model calibration and clinical utility
#   A - Calibration plot, 1-year
#   B - Calibration plot, 3-year
#   C - Calibration plot, 5-year
#   D - Nomogram
#   E - Decision curve analysis (DCA)
path_calib_1yr <- "plots1/cal1.png"      # Panel A
path_calib_3yr <- "plots1/cal3.png"      # Panel B
path_calib_5yr <- "plots1/cal5.png"      # Panel C
path_nomogram <- "plots1/Nomogram1.png"              # Panel D
path_dca       <- "plots1/DCA.png"                  # Panel E

combine_figure(
  row_paths  = list(
    c(path_calib_1yr, path_calib_3yr, path_calib_5yr),
    c(path_nomogram, path_dca)
  ),
  row_labels = list(
    c("A", "B", "C"),
    c("D", "E")
  ),
  output_basename = "Figure4",
  output_dir = "plots1/paper_fig",
  row_gap = 60   # <-- adds a 60px white gap between D and E only
)

# FIGURE 5 -- GSVA Hallmark pathway + CIBERSORT immune infiltration
#   A - GSVA barplot (significant pathways)
#   B - GSVA heatmap
#   C - CIBERSORT cell-type differences by risk group
#   D - CIBERSORT heatmap
path_gsva_bar    <- "plots1/GSVA BARPLOT (ALL SIGNIFICANT PATHWAYS).png"
path_gsva_heat   <- "plots1/GSVA HEATMAP (ALL SIGNIFICANT PATHWAYS).png"
path_ciber_vio   <- "plots1/IMMUNE CELL DIFFERENCES ACROSS RISK GROUPS.png"
path_ciber_heat  <- "plots1/IMMUNE HEATMAP.png"

combine_figure(
  row_paths  = list(
    c(path_gsva_bar, path_gsva_heat),
    c(path_ciber_vio, path_ciber_heat)
  ),
  row_labels = list(c("A", "B"), c("C", "D")),
  output_basename = "Figure5",
  output_dir = "plots1/paper_fig",
  row_gap = 60   # <-- adds a 60px white gap between A/B and between C/D
)

# FIGURE 6: Single-cell resolution validation of the 9-gene signature
#   A - UMAP with 28 annotated cell type clusters
#   B - UMAP colored by combined 9-gene signature module score
#   C - Module score violin (risk-associated vs. protective) across 24 cell types
path_umap_clusters <- "plots1/scRNA/FINAL_Celltype_Annotation_UMAP.png"                    # Panel A
path_umap_9gene    <- "plots1/scRNA/FINAL_Cellular_Distribution_9Gene_Signature.png"       # Panel B
path_modscore      <- "plots1/scRNA/FINAL_Risk_vs_Protective_Violin.png"                   # Panel C

combine_figure(
  row_paths  = list(
    c(path_umap_clusters, path_umap_9gene),
    c(path_modscore)
  ),
  row_labels = list(c("A", "B"), c("C")),
  output_basename = "Figure6",
  output_dir = "plots1/paper_fig"
)

# FIGURE 7: Molecular subtype, TMB, and HRD confound resolution
#   A - Molecular subtype composition by risk group
#   B - Significance summary (before/after adjustment)
#   C - TMB by risk group, within subtype
#   D - HRD by risk group, within subtype
path_subtype     <- "plots1/tmb_hrd_figures_clean/A3_panelA_subtype.png"       # Panel A
path_significance <- "plots1/tmb_hrd_figures_clean/A3_panelD_significance.png" # Panel D
path_tmb         <- "plots1/tmb_hrd_figures_clean/A3_panelB_TMB.png"           # Panel B
path_hrd         <- "plots1/tmb_hrd_figures_clean/A3_panelC_HRD.png"           # Panel C

combine_figure(
  row_paths  = list(
    c(path_subtype, path_significance),
    c(path_tmb),
    c(path_hrd)
  ),
  row_labels = list(
    c("A", "B"),
    c("C"),
    c("D")
  ),
  output_basename = "Figure7",
  output_dir = "plots1/paper_fig"
)

# FIGURE 8: Immunotherapy-response evidence (TIDE, Responder Rate, IPS)
path_tide      <- "plots1/tide_figures_clean/TIDE.png"
path_responder <- "plots1/tide_figures_clean/Responder_Rate_Barplot.png"
path_ips       <- "plots1/tide_figures_clean/IPS_Boxplot.png"

combine_figure(
  row_paths  = list(
    c(path_tide, path_responder, path_ips)
  ),
  row_labels = list(
    c("A", "B", "C")
  ),
  output_basename = "Figure8",
  output_dir = "plots1/paper_fig"
)

# FIGURE 9: Exclusion mechanism decomposition
path_dysexc  <- "plots1/tide_figures_clean/DysfunctionExclusion.png"
path_cells   <- "plots1/tide_figures_clean/CellTypes.png"
path_mech    <- "plots1/tide_figures_clean/Mechanism.png"
path_checkpt <- "plots1/tide_figures_clean/Checkpoints.png"

combine_figure(
  row_paths  = list(
    c(path_dysexc, path_cells),
    c(path_mech, path_checkpt)
  ),
  row_labels = list(
    c("A", "B"),
    c("C", "D")
  ),
  output_basename = "Figure9",
  output_dir = "plots1/paper_fig"
)

# FIGURE 10: Drug sensitivity across risk groups
path_drug_forest <- "plots1/drug_figures_clean/drug_panelA_forest.png"
path_drug_boxplots <- "plots1/drug_figures_clean/drug_panelB_boxplots.png"

combine_figure(
  row_paths  = list(
    c(path_drug_forest),
    c(path_drug_boxplots)
  ),
  row_labels = list(
    c("A"),
    c("B")
  ),
  output_basename = "Figure10",
  output_dir = "plots1/paper_fig"
)

# FIGURE 11: Checkpoint and TME replication (CPTAC-UCEC)
path_cptac_checkpoints <- "plots1/cptac_figures/cptac_panelA_checkpoints.png"
path_cptac_tme         <- "plots1/cptac_figures/cptac_panelB_TME.png"

combine_figure(
  row_paths  = list(
    c(path_cptac_checkpoints),
    c(path_cptac_tme)
  ),
  row_labels = list(
    c("A"),
    c("B")
  ),
  output_basename = "Figure11",
  output_dir = "plots1/paper_fig"
)

# FIGURE 12: Immunotherapy-response and pathway replication (CPTAC-UCEC)
path_cptac_tide        <- "plots1/cptac_figures/CPTAC_TIDE.png"
path_cptac_responder   <- "plots1/cptac_figures/CPTAC_Responder_Rate.png"
path_cptac_gsva        <- "plots1/cptac_figures/CPTAC_GSVA_Replication.png"

combine_figure(
  row_paths  = list(
    c(path_cptac_tide, path_cptac_responder),
    c(path_cptac_gsva)
  ),
  row_labels = list(
    c("A", "B"),
    c("C")
  ),
  output_basename = "Figure12",
  output_dir = "plots1/paper_fig"
)

standardize_figure <- function(input_path, output_basename, output_dir = "plots1/paper_fig",
                               quantize_png = TRUE) {
  
  img <- image_read(input_path)
  
  img_flat <- image_flatten(image_background(img, "white"))
  img_flat <- image_background(img_flat, "white", flatten = TRUE)
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  png_path  <- file.path(output_dir, paste0(output_basename, ".png"))
  tiff_path <- file.path(output_dir, paste0(output_basename, ".tiff"))
  
  png_out <- if (quantize_png) image_quantize(img_flat, max = 256, colorspace = "sRGB") else img_flat
  
  image_write(png_out, png_path, format = "png", density = "300x300")
  image_write(img_flat, tiff_path, format = "tiff", density = "300x300", compression = "LZW")
  
  cat("Saved:", png_path, "\n")
  cat("Saved:", tiff_path, "\n")
  
  print(img_flat)
  invisible(img_flat)
}

# FIGURE 13 -- CRISPR essentiality (DepMap)
standardize_figure(
  input_path = "plots1/depmap_figures/DepMap_signature_gene_essentiality.png",
  output_basename = "Figure13",
  output_dir = "plots1/paper_fig" 
)

# SUPPLEMENTARY FIGURE S1: Risk-associated and protective program UMAPs
#   A - UMAP colored by risk-associated program score
#   B - UMAP colored by protective program score
path_umap_risk       <- "plots1/scRNA/FINAL_RiskAssociated_Program_FeaturePlot.png"    # Panel A
path_umap_protective <- "plots1/scRNA/FINAL_ProtectiveImmune_Program_FeaturePlot.png"  # Panel B

combine_figure(
  row_paths  = list(
    c(path_umap_risk, path_umap_protective)
  ),
  row_labels = list(
    c("A", "B")
  ),
  output_basename = "SupplementaryFigureS1",
  output_dir = "plots1/paper_fig"
)
