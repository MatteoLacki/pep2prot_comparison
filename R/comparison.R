library(readr)
library(readxl)
library(tidyverse)
library(matrixStats)
library(stringr)

pep2prot = read_delim("data/2019-019-hye-prot-pep2prot.csv", ";", escape_double = FALSE, trim_ws = TRUE)
jorg_hf = read_excel("data/2019-019-hye-prot-jorg.xlsx", sheet = "TOP3 quantification", skip = 1)
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
plot(log(P2P$A_med), P2P$R, pch='.')
plot(log(JHF$A_med), JHF$R, pch='.')

P2P$algo = 'pep2prot'
JHF$algo = 'homology_filter'
D = rbind(P2P, JHF)

ggplot(D, aes(x=A_med, y=R)) + 
  geom_hline(yintercept=c(-2,0,1), color='red') +
  geom_point(size=.1) + 
  scale_x_log10() +
  facet_grid(.~algo, scales='free_x') + 
  theme_bw() + 
  ylab('med(A)/med(B)') +
  xlab('med(A)')

W = data.frame(str_split_fixed(P2P$id, "_", 2))
colnames(W) = c('prot', 'origin')
P2P = cbind(P2P, W)

setdiff(P2P$prot, JHF$id)
setdiff(JHF$id, P2P$prot)
intersect(JHF$id, P2P$prot)

DD = D %>% filter(!str_detect(D$id, 'CONTA'))

tail(DD)
W = 
sum(W[,2] == '')
table(W[,2])
D$id
