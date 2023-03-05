library(fitdistrplus)
library(data.table)
library(stats)
library(pheatmap)

args = commandArgs(trailingOnly=TRUE)
input = args[1]
output = args[2]

getDiag=function(m,d){
    tmp=c()
    for(i in c(1:(ncol(m)-d))){
        j=i+d
        tmp=c(tmp,m[i,j])
    }
    if(length(tmp[tmp!=0])>1){
        fw <- fitdist(tmp[tmp!=0], "weibull")
        tmp2 = pweibull(tmp, shape=fw$estimate[1], scale = fw$estimate[2], lower.tail = FALSE, log.p = TRUE)
    }else{
        tmp2=rep(0,length(tmp))
    }
    #print(fw$estimate)
    return(-1*tmp2)
}

m=as.matrix(fread(input))
m1=unlist(as.list(m[m!=0]))
fw <- fitdist(m1, "weibull")
tmp2 = pweibull(unlist(as.list(m)), shape=fw$estimate[1], scale = fw$estimate[2], lower.tail = FALSE, log.p = TRUE)
m3=-1*matrix(tmp2,nrow=nrow(m),ncol=ncol(m))
m3[is.infinite(m3)] = 0
m3[is.na(m3)] = 0
m3[is.nan(m3)] = 0

m4=m3
for(d in c(0:(ncol(m3)-100))){
    p1=getDiag(m,d)
    c=0
    #print(d)
    for(i in c(1:(ncol(m4)-d))){
        c=c+1
        j=i+d
        m4[i,j]=p1[c]
        m4[j,i]=p1[c]
    }
}

for(d in c(0:(ncol(m3)-100))){
    c=0
    #print(d)
    for(i in c(1:(ncol(m4)-d))){
        c=c+1
        j=i+d
        m4[i,j]=max(m4[i,j],m3[i,j])
        m4[j,i]=max(m4[j,i],m3[j,i])
    }
}

write.table(m4, file = output, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)

