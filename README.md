# GO-NLP-PPI
NLP and Gene Ontology-based PPI for SARS-CoV-2

->Download the positive proteins dataset from the related study. (D.E. Gordon, G.M. Jang, M. Bouhaddou, J. Xu, K. Obernier, K.M. White, M.J. O’Meara, V.V. Rezelj, J.Z. Guo, D.L. Swaney, A SARS-CoV-2 protein interaction map reveals targets for drug repurposing, Nature, 583(7816) (2020) 459-468.)

->Generate the negative dataset (Negative_Dataset_Generation).

->Extract GO terms both positive and negative datasets (GO_Term_Extraction). With this process, protein vectors are represented with GO terms.

->GO term vectors are converted to protein vectors using NLP techniques(NLP_Based_Vector_Generation.R). There are 4 datasets for each GO sub ontology (BP, CC, MF). 

->After feature selection, machine learnşng algorithms is applied
