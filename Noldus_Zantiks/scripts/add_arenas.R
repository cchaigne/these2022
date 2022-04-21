
# BE SURE TO EMPTY THE 'WITHOUT' FOLDER BEFORE RUNNING THIS SCRIPT#
# CHECK THAT THE N CORRESPOND TO THE ARENA COLUMN #

rm(list=ls())  # clean R environment

# 0. LIBRAIRIES ----
library(readxl)
library(writexl)

n <- 3  # define the column arena
m <- 96  # number of arenas


# 1. DEFINE FOLDERS ----
setwd("Y://Clair//for_Elise//add_arenas//without")  # input files
dossier_sortie <- "Y://Clair//for_Elise//add_arenas//with"  # output files

# 2. LOAD FILES ----
nom_fichiers <- list.files()
nom_fichiers
fichiers <- list()
for (i in 1:length(nom_fichiers))
{
  fichiers[[i]] <- read_excel(nom_fichiers[i], col_names = FALSE)
}
names(fichiers) <- sub('\\.xlsx$', '', nom_fichiers)


# 3. ADD ARENAS  ----
for (i in 1:length(nom_fichiers))
{
  j <- 1
  for(j in 1:m)
  {
    fichiers[[i]][, n][fichiers[[i]][, n]==j] <- paste("Arena", j, sep = " ")
  }
}


# 4. EXPORT FILES -----
nom_fichiers <- sub('\\.xlsx$', '', nom_fichiers)
for (i in 1:length(nom_fichiers))
{
  write_xlsx(fichiers[[i]], path = paste(dossier_sortie, "//",nom_fichiers[i], "-arena.xlsx", sep = ""), col_names = FALSE)
}


# 5. MANUALLY MODIFY THE WITHOUT FILES ----

  # open the files in the without folder and manually copy/paste the corresponding column arena from files within the "with" folder
  # use the modified files from the without folder for wakefish






