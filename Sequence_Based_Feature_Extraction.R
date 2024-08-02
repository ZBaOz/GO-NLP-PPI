install.packages("protr")
library("protr")

install.packages("xlsx")
library(xlsx)

data <- read.xlsx(file = 'VeriSeti.xlsx', 1, header=TRUE)
pozProt =as.character(as.vector(data[,8]))

negProt =as.character(as.vector(data[,9]))


#EXTRACTING CONJOINT TRIAD BASED FEATURES
sekans=unlist(negProt[1],use.names=F)
CTriadOzellik=extractCTriad(sekans)
for (i in c(2:332))
{
  sekans=unlist(negProt[i],use.names=F)
  y=extractCTriad(sekans)
  CTriadOzellik = rbind(CTriadOzellik, y)
}
write.xlsx(CTriadOzellik, 'CTriadOzellikNeg.xlsx', sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)



#EXTRACTING AMINO ACID COMPOSITION BASED FEATURES

sekans=unlist(negProt[1],use.names=F)
AACOzellik=extractAAC(sekans)
for (i in c(1:332))
{
  sekans=unlist(negProt[i],use.names=F)
  y=extractAAC(sekans)
  AACOzellik = rbind(AACOzellik, y)
}
write.xlsx(AACOzellik, 'AACozellik.xlsx', sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)


#EXTRACTING PSEUDO AMINO ACID BASED FEATURES

sekans=unlist(pozProt[1],use.names=F)
PAACOzellik=extractPAAC(sekans, props = c("Hydrophobicity", "Hydrophilicity","SideChainMass"), lambda = 30, w = 0.05, customprops = NULL)
for (i in c(1:332))
{
  sekans=unlist(negProt[i],use.names=F)
  y=extractPAAC(sekans, props = c("Hydrophobicity", "Hydrophilicity","SideChainMass"), lambda = 30, w = 0.05, customprops = NULL)
  PAACOzellik = rbind(PAACOzellik, y)
}
write.xlsx(PAACOzellik, 'PAACozellik.xlsx', sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)


#EXTRACTING AMPHIPHILIC AMINO ACID COMPOSITION BASED FEATURES
sekans=unlist(pozProt[1],use.names=F)
APAACOzellik=extractAPAAC(sekans, props = c("Hydrophobicity", "Hydrophilicity"), lambda = 30, w = 0.05, customprops = NULL)
for (i in c(1:332))
{
  sekans=unlist(negProt[i],use.names=F)
  y=extractAPAAC(sekans, props = c("Hydrophobicity", "Hydrophilicity"), lambda = 30, w = 0.05, customprops = NULL)
  APAACOzellik = rbind(APAACOzellik, y)
}
write.xlsx(APAACOzellik, 'APAACozellik.xlsx', sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)


#EXTRACTING DIPEPTIDE COMPOSITION BASED FEATURES
sekans=unlist(negProt[1],use.names=F)
DCOzellik=extractDC(sekans)
for (i in c(2:332))
{
  sekans=unlist(negProt[i],use.names=F)
  y=extractDC(sekans)
  DCOzellik = rbind(DCOzellik, y)
}
write.xlsx(DCOzellik, 'DCozellikNeg.xlsx', sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)

