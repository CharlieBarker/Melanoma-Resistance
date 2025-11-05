# Install if necessary
library(magick)

# Specify your PDF files (vector of file paths)
pdf_files <- list.files(
  "~/Desktop/Melanoma_Resistance/paper/Supplementary_plots/evangelia_edited/",
  pattern = "\\.pdf$",
  full.names = TRUE
)
pdf_files<- c("/home/charlie/Desktop/Melanoma_Resistance/paper/Supplementary_plots/evangelia_edited//raw_data_figures_new.pdf",
              "/home/charlie/Desktop/Melanoma_Resistance/paper/Supplementary_plots/evangelia_edited//differential_analysis_supp_new.pdf",
              "/home/charlie/Desktop/Melanoma_Resistance/paper/Supplementary_plots/evangelia_edited//factor1_supp_new.pdf",
              "/home/charlie/Desktop/Melanoma_Resistance/paper/Supplementary_plots/evangelia_edited//combo_supplementary_new.pdf",
              "/home/charlie/Desktop/Melanoma_Resistance/paper/Supplementary_plots/evangelia_edited//Factor3_supplementary_new.pdf",
              "/home/charlie/Desktop/Melanoma_Resistance/paper/Supplementary_plots/evangelia_edited//combined_network_supplement_new.pdf",
              "/home/charlie/Desktop/Melanoma_Resistance/paper/Supplementary_plots/evangelia_edited//differential_analysis_supp_neusethis.pdf")
# Set output folders for JPEGs and PNGs
out_dir_jpg <- "/home/charlie/Desktop/Melanoma_Resistance/paper/Supplementary_plots/evangelia_edited/"
out_dir_png <- "/home/charlie/Desktop/Melanoma_Resistance/paper/Supplementary_plots/evangelia_edited/"
dir.create(out_dir_jpg, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_png, showWarnings = FALSE, recursive = TRUE)

# Loop through each PDF and convert to JPEG + PNG
for (pdf in pdf_files) {
  message("Converting: ", pdf)

  # Read PDF (can contain multiple pages)
  img <- image_read_pdf(pdf, density = 300)  # increase density for higher resolution

  # Loop over pages
  n_pages <- length(img)
  for (i in seq_len(n_pages)) {
    base_name <- tools::file_path_sans_ext(basename(pdf))

    # Define output paths
    out_jpg <- file.path(out_dir_jpg, paste0(base_name, "_page", i, ".jpg"))
    out_png <- file.path(out_dir_png, paste0(base_name, "_page", i, ".png"))

    # Convert and save
    image_write(image = img[i], path = out_jpg, format = "jpg", quality = 95)
    image_write(image = img[i], path = out_png, format = "png")
  }
}

