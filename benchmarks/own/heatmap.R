svg("out_always_memcpy.svg")

library("ggplot2");

data <- read.table("farc_vs_mpi_original.dat", header=TRUE);

p <- ggplot(data=data, aes(x=blklen, y=count)) + 
            geom_tile(aes(fill = farc_pack/mpi_pack), colour = "white") + 
            scale_fill_gradient(low = "green", high = "red", limits=c(0,2))

print(p);

dev.off();