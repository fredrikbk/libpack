svg("speedtest_vector_pack.svg")

library("ggplot2");

data <- read.table("speedtest_vector_pack.dat", header=TRUE);


p <- ggplot(data=data, aes(x=blklen, y=count)) + 
            geom_tile(aes(fill = farc_pack/mpi_pack)) + 
            geom_text(aes(label=sprintf("%.2f",farc_pack/mpi_pack)), color="gray20", size=2) +
            scale_fill_gradient(low = "green", high = "red", limits=c(0,2))

print(p);

dev.off();
