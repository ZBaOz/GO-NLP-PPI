install.packages("rlist")
library("rlist")

install.packages("stringr")
library("stringr")

install.packages("xlsx")
library(xlsx)

install.packages("writexl")
library("writexl")

library("protr")


#read positive and negative proteins from the file
data <- read.xlsx(file = 'VeriSeti.xlsx', 1, header=TRUE)
pozProt =as.character(as.vector(data[,8]))

negProt =as.character(as.vector(data[,9]))

proteinSayi=332

aminoasit <- c('A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V')


komb=vector()
indis=1
for (i in 1:20)
{
  for (j in 1:20)
  {
    tmp <- c(aminoasit[i],aminoasit[j])
    str=paste(tmp, collapse="")
    komb[indis]=str
    indis=indis+1
  }
}

proteinSayi=length(negProt)
  
CKSAAPNeg=matrix(nrow=proteinSayi,ncol=400)
colnames(CKSAAPNeg) = komb
CKSAAPNeg <- data.frame(CKSAAPNeg)

#k=0 When there is no residue between 2 amino acids, the possible combinations are identified
k=0
for (i in 1:proteinSayi)
{
  proteinID=negProt[i]
  protein<-getUniProt(proteinID)
  protUz=str_length(protein)
  kombSayisi=protUz-k-1;
  for (j in 1:400)
  {
    str=komb[j]
    eslesme=str_extract_all(protein,komb[j])
    eslesme=unlist(eslesme,use.names=F)
    CKSAAPNeg[i,j]=length(eslesme)/kombSayisi
  }
}
write_xlsx(CKSAAPNeg,"CKSAAPNeg.xlsx")


#k=1 When there is no residue between 2 amino acids, the possible combinations are identified
komb=vector()
indis=1
for (i in 1:20)
{
  for (j in 1:20)
  {
    tmp <- c(aminoasit[i],".",aminoasit[j])
    str=paste(tmp, collapse="")
    komb[indis]=str
    indis=indis+1
  }
}

k=1
for (i in 1:proteinSayi)
{
  proteinID=negProt[i]
  protein<-getUniProt(proteinID)
  protUz=str_length(protein)
  kombSayisi=protUz-k-1
  for (j in 1:400)
  {
    str=komb[j]
    eslesme=str_extract_all(protein,komb[j])
    eslesme=unlist(eslesme,use.names=F)
    CKSAAPNeg[i,j]=length(eslesme)/kombSayisi
  }
}

write_xlsx(CKSAAPNeg,"CKSAAPNeg.xlsx")


#k=2 When there is no residue between 2 amino acids, the possible combinations are identified
komb=vector()
indis=1
for (i in 1:20)
{
  for (j in 1:20)
  {
    tmp <- c(aminoasit[i],"..",aminoasit[j])
    str=paste(tmp, collapse="")
    komb[indis]=str
    indis=indis+1
  }
}

k=2
for (i in 1:proteinSayi)
{
  proteinID=negProt[i]
  protein<-getUniProt(proteinID)
  protUz=str_length(protein)
  kombSayisi=protUz-k-1
  for (j in 1:400)
  {
    str=komb[j]
    eslesme=str_extract_all(protein,komb[j])
    eslesme=unlist(eslesme,use.names=F)
    CKSAAPNeg[i,j]=length(eslesme)/kombSayisi
  }
}

write_xlsx(CKSAAPNeg,"CKSAAPNeg2.xlsx")


#k=3 When there is no residue between 2 amino acids, the possible combinations are identified
komb=vector()
indis=1
for (i in 1:20)
{
  for (j in 1:20)
  {
    tmp <- c(aminoasit[i],"...",aminoasit[j])
    str=paste(tmp, collapse="")
    komb[indis]=str
    indis=indis+1
  }
}

CKSAAPOzellik=matrix(nrow=proteinSayi,ncol=400)
colnames(CKSAAPOzellik) = komb
CKSAAPOzellik <- data.frame(CKSAAPOzellik)
k=3
for (i in 1:proteinSayi)
{
  protein<-negProt[i]
  protUz=str_length(protein)
  kombSayisi=protUz-k-1
  for (j in 1:400)
  {
    str=komb[j]
    eslesme=str_extract_all(protein,komb[j])
    eslesme=unlist(eslesme,use.names=F)
    CKSAAPOzellik[i,j]=length(eslesme)/kombSayisi
  }
}

write_xlsx(CKSAAPOzellik,"CKSAAPNeg.xlsx")




