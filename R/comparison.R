library(readr)
library(readxl)
library(tidyverse)
library(matrixStats)

pep2prot = read_delim("~/Projects/pep2prot/pep2prot/data/2019-019_HYE/2019-019 HYE neu_user designed 20190621-083010_protein_quantification_report.csv", ";", escape_double = FALSE, trim_ws = TRUE)
jorg_hf = read_excel("~/Projects/pep2prot/pep2prot/data/2019-019_HYE/2019-019 HYE neu_user designed 20190621-083010_quantification_report_Jorg.xlsx", sheet = "TOP3 quantification", skip = 1)

colnames(jorg_hf)

P2P = pep2prot[,c(1,8:13)]
JHF = jorg_hf[,c(5,10:15)]

A = c('A_1','A_2','A_3')
B = c('B_1','B_2','B_3')
AB = c(A,B)
colnames(P2P) = c('id', A, B)
colnames(JHF) = colnames(P2P)

row_meds = function(x, cols) rowMedians(as.matrix(x[,cols]))
P2P$A_med = row_meds(P2P,A)
P2P$B_med = row_meds(P2P,B)

JHF$A_med = row_meds(JHF,A)
JHF$B_med = row_meds(JHF,B)

JHF$R = log(JHF$A_med) - log(JHF$B_med)
P2P$R = log(P2P$A_med) - log(P2P$B_med)

par(mfrow=c(1,2))
plot(log(P2P$A_med), P2P$R)
plot(log(JHF$A_med), JHF$R)

P2P$algo = 'pep2prot'
JHF$algo = 'homology_filter'

D = rbind(P2P, JHF)

ggplot(D, aes(x=A_med, y=R)) + geom_point() + scale_x_log10() +
  facet_grid(.~algo) + theme_bw() + ylab('med(A)/med(B)') + xlab('med(A)')
