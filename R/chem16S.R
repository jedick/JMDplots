# JMDplots/chem16S.R
# Make figure for chem16S paper
# 20230704 first version
# 20230709 moved to JMDplots

# Required packages:
#library(chem16S)
#library(phyloseq)
#library(ggplot2)
#library(patchwork)

chem16S_1 <- function(pdf = FALSE) {

  if(!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Pleast first install phyloseq from Bioconductor")
  }

  # This is needed for composing plots with arithmetic operators (/, -, |)
  if(!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Pleast first install patchwork from CRAN")
  }

  # Colors from scales::hue_pal()(3)
  red <- "#F8766D"
  green <- "#00BA38"
  blue <- "#619CFF"

  # Setup ggplot2 theme
  theme_set(theme_bw())

  ## Panel A: Code example

  # Define text
  input_text <- data.frame(label = '  library(chem16S)
  data(GlobalPatterns, package = "phyloseq")
  ps_metrics(GlobalPatterns, refdb = "RefSeq_206") |> head(n = 3)











  data(mouse.GTDB, package = "chem16S")
  ps_metrics(mouse.GTDB) |> head(n = 3)





')

  output_text <- data.frame(label = '>
>
>
[1] "map_taxa: using these manual mapping(s) to NCBI RefSeq:"
order_Rhizobiales --> order_Hyphomicrobiales (0.2%)
order_Clostridiales --> order_Eubacteriales (0.6%)
family_Ruminococcaceae --> family_Oscillospiraceae (3.1%)
[1] "map_taxa: can\'t map groups order_Stramenopiles (12.94%),
family_ACK-M1 (3.27%), 374 others (11.75%)"
[1] "map_taxa: mapping rate to RefSeq taxonomy is 71.9%"
            Zc        nO2       nH2O
CL3 -0.1439756 -0.6697709 -0.7475996
CC1 -0.1447548 -0.6701114 -0.7450737
SV1 -0.1403848 -0.6623117 -0.7428770
>
>
[1] "map_taxa: mapping rate to GTDB taxonomy is 100.0%"
               Zc        nO2       nH2O
F3D0   -0.1572185 -0.7051765 -0.7704916
F3D141 -0.1549137 -0.7020583 -0.7755448
F3D142 -0.1540525 -0.7008292 -0.7777001
')
                            
  # NULL assignments are needed to avoid "Undefined global functions or variables" in R CMD check
  label <- NULL
  # Start with blank plot then add text
  # https://stackoverflow.com/questions/12518387/can-i-create-an-empty-ggplot2-plot-in-r
  pA <- ggplot() + theme_void() +
    geom_text(data = input_text, aes(label = label), x = -0.2, y = 0.5, family = "Courier", hjust = 0, size = 2.7, colour = green) +
    geom_text(data = output_text, aes(label = label), x = -0.2, y = 0.5, family = "Courier", hjust = 0, size = 2.7) +
    coord_cartesian(clip = "off") +
    labs(title = "(A) Chemical metrics")

  ## Panel B: Mouse gut dataset
  mouse.GTDB <- NULL
  data(mouse.GTDB, package = "chem16S", envir = environment())
  pB1 <- plot_ps_metrics(mouse.GTDB, x = "Day", metrics = "Zc") +
    facet_wrap(~When, scales = "free_x") + geom_line(colour = red) +
    labs(title = "(B) Mouse gut dataset")
  pB2 <- plot_ps_metrics(mouse.GTDB, x = "Day", metrics = "nH2O", quiet = TRUE) +
    facet_wrap(~When, scales = "free_x") + geom_line(colour = blue)
  pB <- pB1 / pB2

  ## Panel C: GlobalPatterns dataset
  data(GlobalPatterns, package = "phyloseq", envir = environment())
  Human = phyloseq::get_variable(GlobalPatterns, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
  phyloseq::sample_data(GlobalPatterns)$Human <- factor(Human)
  SampleType <- NULL
  pC1 <- plot_ps_metrics2(GlobalPatterns, color = "SampleType", shape = "Human", refdb = "RefSeq_206") +
    geom_polygon(aes(fill = SampleType), alpha = 0.5) + geom_point(size = 3) +
    labs(title = "(C) GlobalPatterns dataset") +
    theme(legend.position = "none")
  pC2 <- plot_ps_metrics2(GlobalPatterns, c("O/C", "H/C"), color = "SampleType", shape = "Human", refdb = "RefSeq_206", quiet = TRUE) +
    geom_polygon(aes(fill = SampleType), alpha = 0.5) + geom_point(size = 3) +
    theme(legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 9)) +
    # Change legend title 20230716
    labs(color = "Sample type", fill = "Sample type") +
    # Change the point symbols in the color legend 20230707
    # https://stackoverflow.com/questions/13456765/ggplot2-custom-legend-shapes
    guides(colour = guide_legend(override.aes = list(shape = c(17, 19, 19, 17, 19, 19, 17, 19, 17))))
  pC <- pC1 + pC2 + patchwork::plot_layout(nrow = 1)

  ## Panel D: Zero-width spacer
  pD <- patchwork::plot_spacer()

  # Compose panels
  ## Not used: panels are aligned, but there is empty space above the guide (legend)
  #pall <- (pA | pB) / pC + patchwork::plot_layout(heights = c(3, 2))
  pall <- (pA | pB) / (pC - pD + patchwork::plot_layout(widths = c(1, 0))) + patchwork::plot_layout(heights = c(3, 2))

  # Make PDF or return ggplot2 object
  if(pdf) {
    pdf("Figure_1.pdf", width = 9, height = 7)
    print(pall)
    dev.off()
  } else pall

}
