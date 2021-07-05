qvalue <-function(pval) {

lambda=seq(0,0.90,0.05);

pi0 <- rep(0,length(lambda))
       for(i in 1:length(lambda)) {
            pi0[i] <- mean(pval >= lambda[i])/(1-lambda[i]);
}
spi0 <- smooth.spline(lambda,pi0,df=3);
pi0 <- predict(spi0,x=max(lambda))$y;
pi0 <- min(pi0,1);
m <- length(pval);

u <- order(pval);
idx <- sort.list(pval);
fc <- factor(pval);
nl <- length(levels(fc));
bin <- as.integer(fc);
tbl <- tabulate(bin);
cs <- cumsum(tbl);
tbl <- rep(cs, tbl);
tbl[idx] <- tbl;

qvalue <- pi0*m*pval/tbl;
qvalue[u[m]] <- min(qvalue[u[m]],1);
for(i in (m-1):1) {
	qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1);
}
return (qvalue);

}

