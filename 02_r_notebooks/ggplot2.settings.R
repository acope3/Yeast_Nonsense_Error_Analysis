library(ggpubr)
library(cowplot)
library(ggrepel)
library(broom)
library(patchwork)
library(gridExtra)  ## grid.arrange()
library(ggside)
library(viridis)
library(scales)

theme_set(theme_cowplot())
theme_update(
  #ensure white background (default is transparent)
  plot.background = element_rect(fill="white"),
  ## Center text as default
  text = element_text(hjust = 0.5),
  ## Legend text is left justified
  legend.text = element_text(hjust = 0),
  plot.caption = element_text(hjust = 0),
  plot.margin = margin(1, 1, 1, 1, "cm") ## Ensure axes labels are printed.
  ## ggsave3 will get rid of excess margins in pdfs.
)


switch(plotting.target,
       "paper" = {
       fig.dir = "Figures_Paper"},
       "talk" = {
         theme_update(
           text = element_text(
             face = "bold"),
           line = element_line(
             linewidth = 2)
         )
         fig.dir = "Figures_Talk"
       }
       )

#' Ensure that PDF files are properly cropped when saving. Uses cowplot::ggsave2() knitr::plot_crop
#' 
#' @param file name of file to save in
#' @param object graphic object to save
#' @param crop = TRUE option for turning of cropping
#' @param ... arguments passed to ggsave2()
#' 

ggsave3 <- function(file, plot, crop = TRUE, dpi = 600,...) {
  cowplot::ggsave2(file, plot, dpi = dpi, ...)
  if(stringr::str_ends(file, "(PDF|pdf)") & crop) knitr::plot_crop(file)
}
