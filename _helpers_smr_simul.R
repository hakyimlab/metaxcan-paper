###############################################################################
smr_stat_simul_ <- function(n=20000) {
  x <- rnorm(n)
  y <- rnorm(n)
  d <- data.frame(x, y)
  
  d$t <- (d$x**2 * d$y**2)/(d$x**2 + d$y**2)
  
  d$chi2_t <- pchisq(d$t, 1, lower.tail=FALSE)
  d$chi2_i <- rchisq(n, 1)
  d
}


# chisqzy = rchisq(10e4,df=1,ncp=0)
# chisqzx = rchisq(10e4,df=1,ncp=29)
# instruments = which(chisqzx > qchisq(5e-8, df = 1, lower.tail=F))
# chisqzy = chisqzy[instruments]
# chisqzx = chisqzx[instruments]
# smrstat = chisqzy * chisqzx / (chisqzy + chisqzx)
# ##added two lines to simulate a chi2 random variable and generated qqplot against smrstat.
# kk = rchisq(length(instruments),df=1,ncp=0)
# qqplot(smrstat,kk);abline(0,1)

smr_stat_simul_2_ <- function(n=1e5) {
  chisqzy <- rchisq(n, df=1,ncp=0)  #GWAS
  chisqzx <- rchisq(n, df=1,ncp=29) #eQTL
  
  d <- data.frame(x=chisqzx, y=chisqzy)
  d <- d %>% dplyr::filter(x > qchisq(5e-8, df=1, lower.tail = FALSE)) # genome wide significance
  d$t <- d$y * d$x / (d$y + d$x)
  d$chi2_i <- rchisq(nrow(d),df=1,ncp=0)
  d$chi2_t <- pchisq(d$t, 1, lower.tail=FALSE)
  d
}

smr_stat_simul <- function(n=20000) {
  d <- smr_stat_simul_(n)
  d$chi2_xy <- pchisq(d$x**2+d$y**2,2, lower.tail=FALSE)
  d
}

smr_stat_simul_2 <- function(n=1e5) {
  d <- smr_stat_simul_2_(n)
  d$chi2_xy <- pchisq(d$x+d$y, 2, lower.tail=FALSE)
  d
}

plot_smr_simul_qq_ <- function (d) {
  d_ <- data.frame(t_smr=sort(d$t), x_chi=sort(d$chi2_i))
  r_ <- max(d_$t_smr, d_$x_chi)
  ggplot2::ggplot(data=d_) + 
    ggplot2::geom_point(mapping=ggplot2::aes(y=t_smr, x=x_chi, size=1)) + 
    ggplot2::geom_abline(slope=1, intercept=0) + 
    ggplot2::scale_x_continuous(limits=c(0,r_)) +
    ggplot2::scale_y_continuous(limits=c(0,r_)) +
    ggplot2::ggtitle(expression(bold(paste("QQ-Plot for ",T[SMR], " and a sample ", chi^{2}, " distribution")))) +
    ggplot2::ylab(expression(bold(T[SMR]))) + 
    ggplot2::xlab(expression(bold(chi^{2}))) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=29),
                   axis.title = ggplot2::element_text(size=30, face="bold"),
                   axis.text = ggplot2::element_text(size=27, face="bold"),
                   legend.position="none")
}

plot_smr_simul_qq <- function (d, slope=NULL) {
  p <- plot_smr_simul_qq_(d) 
  if (!is.null(slope)) {
    p <- p + ggplot2::geom_abline(slope=slope, intercept=0, colour="darkgray")  
  }
  p
}

plot_smr_simul <- function (d) {
  r_ <- max(-log10(d$chi2_xy))
  ggplot2::ggplot(data=d) + 
    ggplot2::geom_smooth(mapping=ggplot2::aes(x=-log10(chi2_xy), y=-log10(chi2_t)), method="lm") + 
    ggplot2::geom_point(mapping=ggplot2::aes(x=-log10(chi2_xy), y=-log10(chi2_t)), alpha=0.1, size=1) + 
    ggplot2::geom_abline(slope=1, intercept=0) + 
    ggplot2::ggtitle(expression(bold(paste("Simulation of ",T[SMR], " for 20000 sample genes")))) +
    ggplot2::xlab(expression(bold(paste("-log10[p",chi[2](Z[eQTL]^{2}+Z[GWAS]^{2},"dof=2"),"]")))) + 
    ggplot2::ylab(expression(bold(paste("-log10[p",chi[2](T[SMR],"dof=1"),"]")))) +
    ggplot2::scale_x_continuous(limits=c(0,r_)) +
    ggplot2::scale_y_continuous(limits=c(0,r_)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=32),
                   axis.title = ggplot2::element_text(size=30, face="bold"),
                   axis.text = ggplot2::element_text(size=27, face="bold"))
}

###############################################################################
meta_smr_simul <- function(n=1e5, s=1000) {
  d <- data.frame()
  for (i in 1:s) {
    simul <- smr_stat_simul_(n)
    tr <- mean(simul$t)
    cr <- mean(simul$chi2_i)
    d <- rbind(d, data.frame(chi2r=cr, tr=tr))
  }
  d
}

meta_smr_simul_2 <- function(n=1e5, s=1000) {
  d <- data.frame()
  for (i in 1:s) {
    simul <- smr_stat_simul_(n)
    tr <- mean(simul$t)
    cr <- mean(simul$chi2_i)
    
    x <- sort(simul$chi2_i)
    y <- sort(simul$t)
    k <- lm(y~x-1)
    e <- coefficients(summary(k))[1]
    
    d <- rbind(d, data.frame(chi2r=cr, tr=tr, e=e))
  }
  d
}

meta_smr_simul_yang <- function(n=1e5, s=1000) {
  d <- data.frame()
  for (i in 1:s) {
    simul <- smr_stat_simul_2_(n)
    tr <- mean(simul$t)
    cr <- mean(simul$chi2_i)
    
    x <- sort(simul$chi2_i)
    y <- sort(simul$t)
    k <- lm(y~x-1)
    e <- coefficients(summary(k))[1]
    
    d <- rbind(d, data.frame(chi2r=cr, tr=tr, e=e))
  }
  d
}

plot_meta_smr_simul <- function (d, nn=1000) {
  nn <- paste(nn)
  ggplot2::ggplot(data=d) + 
    ggplot2::geom_smooth(mapping=ggplot2::aes(x=chi2r, y=tr)) + 
    ggplot2::geom_point(mapping=ggplot2::aes(x=chi2r, y=tr), alpha=0.1) + 
    ggplot2::ggtitle(paste("Mean(",expression(T[SMR]), ") for", nn,"simulations")) +
    ggplot2::ylab(expression(paste("Mean(",T[SMR], ")"))) +
    ggplot2::xlab(expression(paste("Mean(",X[chi^{2}], ")"))) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=30),
                   axis.title = ggplot2::element_text(size=25, face="bold"),
                   axis.text = ggplot2::element_text(size=20, face="bold"))
}

plot_meta_smr_simul_density <- function(d, title, column="tr", xlabel=expression(bold(paste("mean", T[SMR],"")))) {
  ggplot2::ggplot(data=d) +
    ggplot2::geom_density(ggplot2::aes_string(x=column), size=1) + 
    ggplot2::scale_x_continuous(limits=c(0,1), breaks= c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(xlabel) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=37),
                   axis.title = ggplot2::element_text(size=30, face="bold"),
                   axis.text = ggplot2::element_text(size=27, face="bold"))
}

plot_meta_smr_simul_histogram <- function(d, title, column="tr", xlabel=expression(bold(paste("mean", T[SMR],"")))) {
  ggplot2::ggplot(data=d) +
    ggplot2::geom_histogram(ggplot2::aes_string(x=column), size=1, fill="grey", color="black") + 
    ggplot2::ggtitle(title) +
    ggplot2::xlab(xlabel) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=37),
                   axis.title = ggplot2::element_text(size=30, face="bold"),
                   axis.text = ggplot2::element_text(size=27, face="bold"))
}
