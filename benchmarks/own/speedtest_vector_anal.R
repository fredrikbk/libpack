library("reshape");

postscript("speedtest_vector_anal.ps");

data <- read.table("kepler_vector.dat.bz2", header=TRUE);


adata <- aggregate((farc_pack / mpi_pack) ~ size + blklen + stride + count + pack_count, data=data, FUN=median);
names(adata)[6]<-"performance"


plot(x=adata$size, y=adata$performance, log="x", type="p",
     xlab="Datasize [B]", ylab="Packing time relative to MPI", xaxt="n");
axis(1, at=c(2e4, 5e4, 1e5, 2e5), labels=c("20K", "50K", "100K", "200K"))

dev.off();
