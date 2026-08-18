#### Create a folder/project and set it as your working directory####

#### Create project setup
library(SciDataReportR)
library(here)
CreateProjectFolders()

# Now, Manually Copy over your data to the Data/Raw folder#######


### Load Data and Use it to Create Codebook that you can manually edit ####
Data <- readxl::read_xlsx(here("Data", "Raw", "DataFileName.xlsx"))
CreateVariableTypesTemplate(Data, here("Data", "Raw", "VariableTypesUnedited.csv"))

# Now, Manually Open the VariableTypesUnedited csv file and edit it.
# You can add columns for filtering if needed
# Save it in Data/Clean as a csv or xlsx file depending on your preferences
# For best practice, append a date to it

### Ready to Analyze!!! #######

# Start from scratch, or with a template

use_EDATemplate()

