install.packages("xlsx")
library(xlsx)
setwd("C:/Users/zeyne/Desktop/Covid-PPI Proje/Veriler")

#pozitif insan proteinleri dosyadan okunur
pozitifProteinler <- read.xlsx(file = 'Proteinler.xlsx', 1, header=TRUE)
pozitifProteinler =as.character(as.vector(pozitifProteinler[,2]))

#t�m aday proteinleri dosyadan okunur
adayProteinler <- read.xlsx(file = 'AdayProteinler.xlsx', 1, header=TRUE)
adayProteinler <- adayProteinler[,1]

#aday proteinler ile pozitif proteinlerin sekanslar� al�n�r.
install.packages("protr")
library("protr")
pozProtSeq<-getUniProt(pozitifProteinler)
adayProtSeq<-getUniProt(adayProteinler)

#pairwise alignment i�in k�t�phane tan�mlar�
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("pairwiseAlignment")
library(Biostrings)
data(BLOSUM62)

install.packages("writexl")
library("writexl")

#pozitif veri k�mesi ile aday proteinlerin identity de�erleri 'identity' matrisinde tutulur.
#matrisin sat�r isimleri pozitif protein id'ler, s�tun isimleri ise aday protein id ler olarak atand�
identity = matrix(NA,nrow=length(pozitifProteinler),ncol=length(adayProteinler))
rownames(identity) = pozitifProteinler
colnames(identity) = adayProteinler
identity <- data.frame(identity)
dim(identity)

#Pozitif veri k�mesindeki proteinler ile aday veri k�mesindeki proteinlerin amino asit sekanslar�
#uniprotdan �ekilir. For d�ng�s�nde, her bir pozitif proteinin her bir aday protein ile sequence identity
#skoru hesaplan�r. sonu�lar 'identity' matrisine al�n�r.
for (i in 332:length(pozitifProteinler))
{
  print(i)
  #s1 <- getUniProt(pozitifProteinler[183])
  s1=pozProtSeq[i]
  s1=unlist(s1,use.names=F)
  for (j in 1:length(adayProteinler))
  {
    #s2 <- getUniProt(adayProteinler[j])
    s2=adayProtSeq[j]
    s2=unlist(s2,use.names=F)
    palign1 <- pairwiseAlignment(s1, s2, substitutionMatrix = "BLOSUM62",type="global",gapOpening=8, gapExtension = 4)
    identity[i,j]=pid(palign1)
  }
  write_xlsx(identity,"C:/Users/zeyne/Desktop/identity.xlsx")
}

identity2 <- read.xlsx(file = 'identity.xlsx', 1, header=TRUE)

#pozitif etkile�imlerden, her bir sars proteini ile ka� tane insan proteinin etkile�ime girdi�ini tutan veri
#okunur. bu bir sars proteini ile ka� tane pozitif etkile�im varsa o kadar negatif etkile�im verisi olu�turmak i�in
library("readxl")
etkilesimSayi <- read.xlsx(file = 'Proteinler.xlsx', 7, header=TRUE)

#�rne�in sars cov 2 e ile etkile�imli 6 insan proteini var. Bu durumda 6 negatif protein bulunacak. identity
#matrisindeki ilk 6 sat�r, cov2 e ile etikle�imli proteinlerin aday proteinler ile benzerli�ini i�erir.
#o 6 sat�r protein matrisine al�nd�. Ama� cov2 e ile en d���k ortalama benzerli�e sahip 6 aday proteini belirlemek
#protein matrsinde her bir s�tunun ortalamas� hesaplan�r. min ilk 6 h�crenin s�tun isimleri, sars cov2 e ile
#etkile�imsiz olarak etiketlenir.
negatifProtein = matrix(NA,nrow=332,ncol=2)
indis=1
for (j in 1:length(etkilesimSayi))
{
  sayi=etkilesimSayi[1,j]
  protein=identity[1:sayi,]
  protein=protein[,2:5708]
  toplam=colSums(protein, na.rm = FALSE, dims = 1)
  ortalama=toplam/sayi
  # 
  # for(i in 1:(sayi*2))
  # {
  #   min=which.min(ortalama)
  #   ortalama[min]=10000
  # }
  
  for(i in 1:sayi)
  {
    min=which.min(ortalama)
    print(min)
    negatifProtein3[indis,1]=names(ortalama[min])
    negatifProtein3[indis,2]=ortalama[min]
    ortalama[min]=10000
    indis=indis+1
  }
  protein=NULL
}
#negatifProtein matrisinde ilk s�tun proteinler 2. s�tun ortalama benzerliklerdir.

negatifProtein <- data.frame(negatifProtein)
write_xlsx(negatifProtein,"C:/Users/zeyne/Desktop/negatifData.xlsx")



