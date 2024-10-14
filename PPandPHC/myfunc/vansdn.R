#分母がnの分散
#x:データベクトル
van<-function(x){mean((x-mean(x))^2)} 

#分母がnの標準偏差
#x:データベクトル
sdn<-function(x){sqrt(van(x))} 

