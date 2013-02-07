library("reshape");

data <- read.table("vector_bench_memcpy.dat", header=TRUE);


adata <- aggregate((mpi_pack / farc_pack) ~ size + blklen + stride + count + pack_count, data=data, FUN=median);
names(adata)[6]<-"performance"


plot(x=adata$blklen, y=adata$performance, log="xy");