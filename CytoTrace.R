# CytoTRACE Analysis (Trajectory Analysis)

print ("R is running!")

install.packages("BiocManager")  # If not installed already
BiocManager::install("sva")

devtools::install_local("C:\\Users\\annan\\OneDrive\\Documents\\BENg 204 CRC Project\\CytoTRACE_0.3.3.tar.gz")

# Calling CytoTRACE Library
library(CytoTRACE)
library(data.table)

# Setting Working Directory
setwd("C:\\Users\\annan\\OneDrive\\Documents\\BENg 204 CRC Project")

# Unzip file
unzip("cytotrace_input_normalized.csv.zip", exdir = "C:\\Users\\annan\\OneDrive\\Documents\\BENg 204 CRC Project")

# Accessing File
print("reading CSV")
ematrix <- fread("cytotrace_input_normalized.csv", data.table = FALSE)

rownames(ematrix) <- ematrix[,1]  # Set first column as row names
ematrix <- ematrix[,-1] # Remove the first column from data

# Conversion of Matrix
print("Matrix conversion in progress...")
ematrix <- as.matrix(ematrix)
mode(ematrix)

# Matrix QC
print("Running QC checks...")
print(paste("Matrix dimensions:", paste(dim(ematrix), collapse = " x ")))  # Print (genes, cells)
print(paste("Missing values:", sum(is.na(ematrix))))  # Count NA values

# Run CytoTRACE if the matrix is valid
if (ncol(ematrix) > 0 && nrow(ematrix) > 0) {
  print("Running CytoTRACE...")
  rna_ta <- CytoTRACE(ematrix)
  
  # Output results
  print("CytoTRACE Analysis Complete!")
  print(head(rna_ta$CytoTRACE))
  
  # Saving output
  write.csv(rna_ta$CytoTRACE, "CytoTRACE_scores.csv")
} else {
  print("Error: Expression matrix is empty or incorrectly loaded.")
}

# Plotting Trajectory Results
plotCytoTRACE(rna_ta)

