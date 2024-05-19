# This is a script for plotting multidimensional scaling components to visualize population stratification

# load necessary packages
source(file = "../scripts/load_necessary_packages.R")

# load the MDS data
MDSdata <- read.table(
                      file="MDS_merge2.mds",
                      header=TRUE
                      )

# load the 1000 Genomes population codes
populationCodes <- 
  read.table(
    file="populationCodefile.txt",
    header=TRUE
    )

# join the MDS data with the population code data by individual ID (IID)
MDSdata <- 
  MDSdata %>%
  select(
    .data = .,
    -FID
    ) %>%
  left_join(
    x = .,
    y = populationCodes %>%
        select(
          .data = .,
          -FID
          ),
    by = "IID"
      ) %>%
   arrange(
     .data = .,
     populationCode
     )

# create a scatter plot of the first 2 MDS coordinates
# labeled by 1000 Genomes continental population code
MDSplot <-
  MDSdata %>% 
  ggplot(
    data = .,
    mapping = aes(
                  x = C1,
                  y = C2,
                  color = populationCode
                    )
    ) +
    geom_point(
               size = 1.2, 
               alpha = 1,
               colour = "black",
               pch = 21,
               aes(fill = populationCode)
               ) +
    theme_classic() +
    scale_fill_discrete(name = paste0("1000 genomes","\n","   population code")) +
    xlab(label = "MDS component 1") +
    ylab(label = "MDS component 2") +
    theme(
      legend.title = element_text(family = "Arial",hjust = 0.5,face = "bold",size = 11),
      legend.text = element_text(family = "Arial",hjust = 0.5,face = "bold",size = 11),
      axis.title = element_text(family = "Arial",face = "bold",size = 11),
      axis.text = element_text(family = "Arial",size = 11,face = "bold"),
      legend.position = c(0.8,0.3),
      legend.box.background = element_rect(colour = "black",size = 1.5)
    ) 

tiff(
  filename = "./MDS.tiff",
  compression = "none",
  res = 300,
  units = "in",
  width = 7,
  height = 5
  )

MDSplot

dev.off()
