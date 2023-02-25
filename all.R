source('data_generator.R')
source('method.R')
library("energy")
p=6 
data_method = c('MDA-1','MDA-2','LDA-1-Ind','LDA-1-Cor', 'LDA-2-Ind','LDA-2-Cor','SIM-1','SIM-2','LR-Ind','LR-Cor','MIM-1','MIM-2','MIM-3','LDA-Ind','LDA-Cor')   
estimate_method = c("mases","energy",'sir','save','dr')
basis = generate_basis(p)
for (i in data_method) {
  print(i)
  result =Simulation(basis,estimate_method,i,2)
}



library(ggplot2)
result =Simulation(basis,estimate_method,'MDA-1',1)
Data_MDA_1 =MDA_2(basis)
Y = Data_MDA_1[,1]
X = Data_MDA_1[,-1]
b1 =t(basis[,1:1])%*% t(Data_MDA_1[,-1])
b2 =t(basis[,2:2])%*% t(Data_MDA_1[,-1])
df = data.frame(t(rbind(Y,b1,b2)))

ggplot(data = df, aes(x=V2, y=V3,col = factor(Y))) +
  geom_point(size=2) +theme_classic() + 
  labs(y="beta2", x = "beta1",title='MDA-1')


