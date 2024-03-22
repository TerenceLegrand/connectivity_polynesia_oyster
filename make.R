#' connectivity_polynesia_oyster: A Research Compendium
#' 
#' @description 
#' A paragraph providing a full description of the project and describing each 
#' step of the workflow.
#' 
#' @author Terence Legrand \email{legrandterence@gmail.com}
#' 
#' @date 2024/02/28



## Install Dependencies (listed in DESCRIPTION) ----
devtools::install_deps(upgrade = "never")


## Load Project Addins (R Functions and Packages) ----
devtools::load_all(here::here())

## Add dependencies 
rcompendium::add_dependencies(".")


## Global Variables ----

# Biophysical model site 
site <- utils::read.csv(here::here("data","raw-data","Genetics","df_islands_processed.csv"),
                         sep = ",")

# Genetic site
pop.coordinates <- openxlsx::read.xlsx(here::here("data","raw-data","Genetics","pop_coordinates.xlsx"))


# Diffrentiation matrix
matrix_genetic <- utils::read.csv(here::here("data","raw-data","Genetics","fst_matrix.csv"),
                                  sep = " ")

# Nbr of generation
nbr.gen <- c(seq(1,10,1),seq(15,500,5))

## Run Project ----

# Compute a unique matrix
temporal_mean_matrix(site = site)

# Produce the three scenarios
filter_matrix(site = site)

# Transform label indexed matrices to numeric indexed matrices
index_matrix(site = site, archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas")
index_matrix(site = site, archipelago = "tuamotuWest_tuamotuEast_gambier_society")
index_matrix(site = site, archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas_austral")

# Correlation between connectivity outputs and genetic differentiation
correlation(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas",
            nbr.gen = nbr.gen, matrix.diff = matrix_genetic, pop.coordinates = pop.coordinates)

correlation(archipelago = "tuamotuWest_tuamotuEast_gambier_society",
            nbr.gen = nbr.gen, matrix.diff = matrix_genetic, pop.coordinates = pop.coordinates)

correlation(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas_austral",
            nbr.gen = nbr.gen, matrix.diff = matrix_genetic, pop.coordinates = pop.coordinates)

# Saturation of connectivity outputs
saturation(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas_austral",
           nbr.gen = nbr.gen, site_data = site)
saturation(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas",
           nbr.gen = nbr.gen, site_data = site)
saturation(archipelago = "tuamotuWest_tuamotuEast_gambier_society",
           nbr.gen = nbr.gen, site_data = site)


## ------------
# Report of results

# Network SI
source("function_plot.R")
panel_A <- network_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas")
panel_B <- network_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas_austral")
panel_C <- network_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society")
panel_D <- legend_network_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society")
legend <- cowplot::get_legend(panel_D)

layout <- "
    AB
    CD
    "

figure_netwok_SI <- panel_A + panel_B + panel_C + legend +
  patchwork::plot_annotation(tag_levels = 'A') + 
  patchwork::plot_layout(design = layout, guides = "collect") & theme(legend.position = 'bottom')

ggsave(filename = paste0("network_plot",".pdf"),
       plot = figure_netwok_SI,
       device = cairo_pdf(),
       path = here::here("outputs")
)

# Saturation SI

panel_A <- saturation_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas")
panel_B <- saturation_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas_austral")
panel_C <- saturation_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society")

layout <- "
    AB
    CD
    "

figure_saturation_SI <- panel_A + panel_B + panel_C + patchwork::guide_area() +
  patchwork::plot_annotation(tag_levels = 'A') +
  patchwork::plot_layout(design = layout, guides = "collect") & theme(legend.position = 'left')

ggsave(filename = paste0("saturation_plot",".pdf"),
       plot = figure_saturation_SI,
       device = cairo_pdf(),
       path = here::here("outputs")
)

# Find whereis the extra connections 

connectivity.study <- filial_vs_coalescent_connection(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas", nbr_gen = 500)



# Optimal AIC SI
panel_A <- aic_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas")
panel_B <- aic_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas_austral")
panel_C <- aic_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society")

layout <- "
    AB
    CD
    "

figure_AIC_SI <- panel_A + panel_B + panel_C + patchwork::guide_area() +
  patchwork::plot_annotation(tag_levels = 'A') +
  patchwork::plot_layout(design = layout, guides = "collect") & theme(legend.position = 'left')

ggsave(filename = paste0("AIC_plot",".pdf"),
       plot = figure_AIC_SI,
       device = cairo_pdf(),
       path = here::here("outputs")
)

# Centrality
plot_hm <- centrality_hm(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas")
plot_btw <- centrality_btw(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas")

layout <- "
    AB
    ##
    "

figure_centrality <- plot_hm + plot_btw +
  patchwork::plot_annotation(tag_levels = 'A') +
  patchwork::plot_layout(design = layout, guides = "collect") & theme(legend.position = 'right')

ggsave(filename = paste0("centrality",".pdf"),
       plot = figure_centrality,
       device = cairo_pdf(),
       path = here::here("outputs")
)

# Save centrality per sampled pop
centrality_population(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas",site = site, pop.coordinates = pop.coordinates)

# Summary
final_table_summary <- rbind(table_summary(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas"),
                             table_summary(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas_austral"),
                             table_summary(archipelago = "tuamotuWest_tuamotuEast_gambier_society"))

file_name <- here::here("outputs","summary_table.xlsx")
openxlsx::write.xlsx(final_table_summary, file = file_name, rowNames = FALSE, colNames = TRUE)



# Plot scatter
panel_A <- scatter_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas", connection = "coalescent")
panel_B <- scatter_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society_marquesas_austral", connection = "coalescent")
panel_C <- scatter_plot(archipelago = "tuamotuWest_tuamotuEast_gambier_society", connection = "coalescent")
panel_D <- deviation_island(archipelago = "tuamotuWest_tuamotuEast_gambier_society", connection = "coalescent")

layout <- "
    AB
    CD
    "

figure_scatter <- panel_A + panel_B + panel_C + panel_D +
  patchwork::plot_annotation(tag_levels = 'A') 

ggsave(filename = paste0("scatter_plot",".pdf"),
       plot = figure_scatter,
       device = cairo_pdf(),
       path = here::here("outputs")
)


Kruskal_Wallis <- deviation_island_test(archipelago = "tuamotuWest_tuamotuEast_gambier_society", connection = "coalescent")


# List all R scripts in a sequential order and using the following form:
# source(here::here("analyses", "script_X.R"))
