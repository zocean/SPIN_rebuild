library(data.table)
#library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
input = args[1]
method = args[2]
size = as.numeric(args[3])
output_matrix = args[4]

smoothMatrix=function(m,method="mean",size=1){
    m1=m
    for(i in c(1:nrow(m))){
        for(j in c(1:ncol(m))){
            x1=max(1,i-size)
            x2=min(nrow(m),i+size)
            y1=max(1,j-size)
            y2=min(ncol(m),j+size)
            t=as.numeric(unlist(m[x1:x2,y1:y2]))
            #print(c(x1,x2,y1,y2))
            if(method=="mean"){score=mean(t)}
            if(method=="median"){score=median(t)}
            m1[i,j]=score
        }
    }
    return(m1)
}


m=fread(input)
m1=smoothMatrix(m=as.matrix(m),method=method,size=size)
write.table(m1, file = output_matrix, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE)

