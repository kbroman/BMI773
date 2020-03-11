# gough x wsb F2 body weights

library(broman)
library(qtl2)

url <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/Gough/gough_pheno.csv"
file <- basename(url)
if(!file.exists(file)) download.file(url, file)
phe <- qtl2::read_csv_numer(file)
sm <- phe[,1:16+16]
de <- phe[,1:16+16*2]
phe <- phe[,1:16]

x <- seq_len(ncol(phe))
yli <- c(0, 1.05*max(phe,na.rm=TRUE))
xli <- c(min(x)-0.5, max(x)+0.5)

ind <- 417 # individual of interest

pdf("../Figs/gough_weights1.pdf", height=5.5, width=9.75)
par(mar=c(4.1, 4.1, 0.6, 0.6))
grayplot(x, phe[1,], ylim=yli, xlim=xli, xaxs="i", xat=c(1, 5, 10, 15),
         yaxs="i", type="n", xlab="Week", ylab="Body weight (g)")
points(x, phe[ind,], pch=21, bg="lightblue")
dev.off()


pdf("../Figs/gough_weights2.pdf", height=5.5, width=9.75)
par(mar=c(4.1, 4.1, 0.6, 0.6))
grayplot(x, phe[1,], ylim=yli, xlim=xli, xaxs="i", xat=c(1, 5, 10, 15),
         yaxs="i", type="n", xlab="Week", ylab="Body weight (g)")
points(x, phe[ind,], pch=21, bg="lightblue")
points(x[4], phe[ind,4], pch=21, bg=brocolors("web")["orange"], cex=1.3)
dev.off()


pdf("../Figs/gough_weights3.pdf", height=5.5, width=9.75)
par(mar=c(4.1, 4.1, 0.6, 0.6))
grayplot(x, phe[1,], ylim=yli,  xlim=xli, xaxs="i", xat=c(1, 5, 10, 15),
         yaxs="i", type="n", xlab="Week", ylab="Body weight (g)")
for(i in 1:nrow(phe)) lines(x, phe[i,], col="gray")
points(x, phe[ind,], pch=21, bg="lightblue")
points(x[4], phe[ind,4], pch=21, bg=brocolors("web")["orange"], cex=1.3)
dev.off()


# calculate differences and 2nd-differences

diffs <- function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    c(d1=max(abs(diff(x))),
      d2=max(abs(diff(diff(x)))))
}

d <- t(apply(phe, 1, diffs))

pdf("../Figs/gough_weights4.pdf", height=5.5, width=6.5)
par(mar=c(4.1, 4.1, 0.6, 0.6))
grayplot(d[,1], d[,2], xlab="Max absolute change",
         ylab="Max absolute 2nd-difference")
points(d[ind,1], d[ind,2], pch=21, bg="orange", cex=1.3)
dev.off()



pdf("../Figs/gough_weights5.pdf", height=5.5, width=9.75)
par(mar=c(4.1, 4.1, 0.6, 0.6))
grayplot(x, phe[1,], ylim=yli,  xlim=xli, xaxs="i", xat=c(1, 5, 10, 15),
         yaxs="i", type="n", xlab="Week", ylab="Body weight (g)")
points(x, phe[ind,], pch=21, bg="lightblue")
lines(x, sm[ind,])
dev.off()


resid <- apply( abs(phe-sm)/sm*100, 1, max, na.rm=TRUE)
pdf("../Figs/gough_weights6.pdf", height=5.5, width=9.75)
par(mar=c(4.1, 4.1, 0.6, 0.6))
grayplot(resid, xlab="Index", ylab="Absolute value of relative residual (%)")
points(ind, resid[ind], pch=21, bg="orange", cex=1.3)
dev.off()
