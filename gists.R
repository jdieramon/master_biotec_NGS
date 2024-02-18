# Install libraries
install.packages("BiocManager")
BiocManager::install("Rbowtie2")

# Load libraries
library(Rbowtie2)

# Manipulation of Directories 
dir.exists()
dir.create()

# List the files in a directory
list.files()
dir()

# File manipulation
file.exists()
file.create()
file.rename()
file.remove()

# Download File from the Internet
download.file()

# Remove object
rm()


# Object features
class()
length()


#Test Objects for Exact Equality
identical()

# Scatter-plot
plot()

# Histogram
hist()


# Save object : .csv file
write.csv()

# Load object 
read.csv()

# Timing code 
start_time <- Sys.time()
 # do something
end_time <- Sys.time()
end_time - start_time



# Create dataframe (table data) 
# ----------------------------------------------------------------------------- 
# Create sample vector
sampleID = paste0("SRR93364", 68:76)

# Create condition vector
conditions = rep(c("basal", "indust1", "indust2"), each = 3)

# Create description vector
description = rep(c("pH5_0.04CO2", "pH5_5CO2", "pH3_0.04CO2"), each = 3)

# Create data frame
saccha_metadata <- data.frame(sampleID, conditions, description)

# Assign the row names of the data frame
rownames(saccha_metadata) = paste(conditions, 1:3, sep = "_")


# Save the df into a csv file 
# ----------------------------------------------------------------------------- 
write.csv(saccha_metadata, file = "saccha_metadata.csv")

