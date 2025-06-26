#-----------------------------------------
# Open Libraries
#-----------------------------------------
library(tidyverse)
library(openxlsx)
#-----------------------------------------

# Load the data
valdata_VPA_POS <- read.table("data/Apolar_validation_fake/p2024_051_POS_VPA.txt", header = TRUE, sep = "\t")
valdata_VPB_POS <- read.table("data/Apolar_validation_fake/p2024_051_POS_VPB.txt", header = TRUE, sep = "\t")
valdata_VPCD_POS <- read.table("data/Apolar_validation_fake/p2024_051_POS_VPC_VPD.txt", header = TRUE, sep = "\t")

# Load aliases
compound_alias <- read.xlsx("data/Apolar_validation_fake/20250625_compound_alias.xlsx") %>%
  rename(Component.Name = compound)

# Give aliases to the compounds
valdata_VPA_POS2 <- left_join(valdata_VPA_POS, compound_alias) %>%
  select(-Component.Name) %>%
  rename(Component.Name = compound_alias)

# export
write.csv(valdata_VPA_POS2, "data/Apolar_validation_fake/p2024_051_POS_VPA_fake.csv", row.names = FALSE)
  
# Give aliases to the compounds
valdata_VPB_POS2 <- left_join(valdata_VPB_POS, compound_alias) %>%
  select(-Component.Name) %>%
  rename(Component.Name = compound_alias)

# export
write.csv(valdata_VPB_POS2, "data/Apolar_validation_fake/p2024_051_POS_VPB_fake.csv", row.names = FALSE)

# Give aliases to the compounds
valdata_VPCD_POS2 <- left_join(valdata_VPCD_POS, compound_alias) %>%
  select(-Component.Name) %>%
  rename(Component.Name = compound_alias)

# export
write.csv(valdata_VPCD_POS2, "data/Apolar_validation_fake/p2024_051_POS_VPCD_fake.csv", row.names = FALSE)
#----------------------------
# Load valtargets
valtargets <- read.csv("data/Apolar_validation_fake/20250315_valtargets_XAB_positive.csv")
compound_alias2 <- compound_alias %>%
  rename(compound = Component.Name)
# use the aliases
valtargets_aliased <- left_join(valtargets, compound_alias2) %>%
  select(-compound, -PubChem_CID) %>%
  rename(compound = compound_alias) %>%
  select(compound, rt, mz, MAC_mix)

#export
write.csv(valtargets_aliased, "data/Apolar_validation_fake/20250315_valtargets_XAB_positive_fake.csv")
