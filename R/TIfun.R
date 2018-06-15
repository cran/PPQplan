#' Probability of Passing PPQ Test using Tolerance Interval
#'
#' The function for calculating the probability of passing critical quality attributes (CQA) PPQ test .
#'
#' @usage ti.pp(Llim, Ulim, mu, sigma, n, n.batch, alpha, coverprob, side)
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param n sample size (number of locations) per batch
#' @param n.batch number of batches for passing PPQ during validation
#' @param alpha significant level for constructing the tolerance interval
#' @param coverprob coverage probability for constructing the tolerance interval
#' @param side whether a 1-sided or 2-sided tolerance interval is required (determined by side = 1 or side = 2, respectively).
#' @return
#' A numeric value of the passing/acceptance probability
#' @seealso \code{rl.pp}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' ti.pp(sigma=0.5, mu=2.5, n=10, n.batch=1, Llim=1.5, Ulim=3.5, alpha=0.05)
#'
#' sapply(X=c(0.1,0.5, 1,2,3,4,5,10), FUN = ti.pp, mu=97, n=10, Llim=95, Ulim=105,
#' n.batch=1, alpha=0.05)
#' sapply(X=c(0.1,0.5, 1,2,3,4,5,10), FUN = ti.pp, mu=100, n=10, Llim=95, Ulim=105,
#' n.batch=1, alpha=0.05)
#' @importFrom tolerance K.factor
#' @import stats
#' @export
ti.pp <- function(Llim, Ulim, mu, sigma, n=10, n.batch=1, alpha=0.05, coverprob = 0.675, side = 2){
  k <- tolerance::K.factor(n = n, alpha = alpha, P = coverprob, side = side)
  Func <- function(V){
    (pnorm(q = Ulim-k*sqrt(V), mean = mu, sd = sigma/sqrt(n))-pnorm(q = Llim+k*sqrt(V), mean = mu, sd = sigma/sqrt(n)))*dchisq(x = (n-1)*V/sigma^2, df = n-1)
  }
  (min(1, integrate(Func, lower = 0, upper = ((Ulim-Llim)/(2*k))^2, rel.tol = 1e-10)$value*(n-1)/sigma^2))^n.batch
}



#' Operating Characteristic (OC) Curves for the PPQ Plan using Tolerance Interval.
#'
#' The function for  plotting the OC curve to show the PPQ plan based on the specification test, given lower and upper specification limits.
#'
#' @usage ti.occurve(attr.name, attr.unit, Llim, Ulim, mu, sigma, n, n.batch, alpha,
#' coverprob, side, add.reference, NV)
#' @param attr.name user-defined attribute name
#' @param attr.unit user-defined attribute unit
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param n sample size (number of locations) per batch
#' @param n.batch number of batches for passing PPQ during validation
#' @param alpha significant level for constructing the tolerance interval.
#' @param coverprob converage probability for constructing the tolerance interval
#' @param side whether a 1-sided or 2-sided tolerance interval is required (determined by side = 1 or side = 2, respectively).
#' @param add.reference logical; if \code{TRUE}, then add reference OC curves (Baseline and High Performance) in the plot.
#' @param NV nominal volume for the specification test.
#' @return
#' OC curves for specification test and PPQ plan.
#' @seealso \code{ti.pp} and \code{rl.pp}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#'
#'
#' ti.occurve(attr.name = "Sterile Concentration Assay", attr.unit="%",
#' mu=97, sigma=seq(0.1, 10, 0.1), Llim=95, Ulim=105, n=10, add.reference=TRUE)
#'
#' ti.occurve(attr.name = "Sterile Concentration Assay", attr.unit="%",
#' mu=100, sigma=seq(0.1, 10, 0.1), Llim=95, Ulim=105, n=10, add.reference=TRUE)
#'
#' ti.occurve(attr.name = "Extractable Volume", attr.unit = "% of NV=3mL",
#' Llim = 100, Ulim = Inf, mu=102.5, sigma=seq(0.2, 6 ,0.05), n=40,
#' alpha = 0.05, coverprob = 0.97, side=1, NV=3)
#'
#' ti.occurve(attr.name = "Extractable Volume", attr.unit = "% of NV=3mL",
#' Llim = 100, Ulim = Inf, mu=102.5, sigma=seq(0.2, 6 ,0.05), n=40,
#' alpha = 0.05, coverprob = 0.992, side=1, NV=3)
#'
#'
#' @export
ti.occurve <- function(attr.name="", attr.unit="", Llim, Ulim, mu, sigma, n=10, n.batch=1, alpha=0.05, coverprob=0.675, side=2, add.reference=FALSE, NV=10){
  rlpp  <- rl.pp(Llim=Llim, Ulim=Ulim, mu = mu, sigma = sigma, NV=NV)
  if (length(mu)==1 & length(sigma)>1){
    tipp <- sapply(X=sigma, FUN = ti.pp, mu=mu, n=n, Llim=Llim, Ulim=Ulim, n.batch=n.batch, alpha=alpha, coverprob=coverprob, side=side)
    plot(sigma, rlpp, xlab=paste0("Standard Deviation (",attr.unit, ")"),
         ylab="Probability of Passing", type="l", lty=1, col=1, lwd=2,
         ylim=c(0,1)) # plot release test
    lines(sigma, tipp, lty=1, col=2, lwd=2) # add PPQ Assessment line
    bl.sigma <- sigma[which.min(abs(rlpp-0.9))]  # search the std dev for baseline plan
    hp.sigma <- sigma[which.min(abs(rlpp-0.95))] # search the std dev for high performance plan
    abline(h=0.95, lty=2, col=16)
    abline(h=0.90, lty=2, col=16)
    abline(h=0.20, lty=2, col=16)
    abline(h=0.05, lty=2, col=16)
    abline(v=bl.sigma, lty=2, col=16)
    abline(v=hp.sigma, lty=2, col=16)
    points(bl.sigma,0.2, pch="*", col=4, cex=3) # mark baseline
    points(hp.sigma,0.05, pch="*", col=3, cex=3) # mark high performance

    if (add.reference == TRUE){
      bl.level <- optimize(f = function(x){abs(ti.pp(sigma = bl.sigma, mu = mu, n = n, Llim = Llim, Ulim = Ulim, n.batch = n.batch, alpha = x, coverprob=coverprob, side=side)-0.2)}, interval = c(0,0.5), tol = 1e-6)$minimum # optimize the alpha for baseline
      bl.tipp <- sapply(X=sigma, FUN = ti.pp, mu=mu, n=n, Llim=Llim, Ulim=Ulim, n.batch=n.batch, alpha=bl.level, coverprob=coverprob, side=side)
      lines(sigma, bl.tipp, lty=2, col=4, lwd=2) # add PPQ Assessment baseline
      hp.level <- optimize(f = function(x){abs(ti.pp(sigma = hp.sigma, mu = mu, n = n, Llim = Llim, Ulim = Ulim, n.batch = n.batch, alpha = x, coverprob=coverprob, side=side)-0.05)}, interval = c(0,0.5), tol = 1e-6)$minimum # optimize the alpha for baseline
      hp.tipp <- sapply(X=sigma, FUN = ti.pp, mu=mu, n=n, Llim=Llim, Ulim=Ulim, n.batch=n.batch, alpha=hp.level, coverprob=coverprob, side=side)
      lines(sigma, hp.tipp, lty=2, col=3, lwd=2) # add PPQ Assessment high performance line
      legend("topright", legend=c(paste0("Spec: ", Llim, "-", Ulim, attr.unit),
                                  paste0("TI ", (1-alpha)*100, "% Confidence", coverprob*100, "% Coverage"),
                                  paste0("TI ", round((1-bl.level)*100,2), "%/", coverprob*100, "% Baseline"),
                                  paste0("TI ", round((1-hp.level)*100,2), "%/", coverprob*100, "% High Perform")), bg="white",
             col=c(1,2,4,3), lty=c(1,1,2,2), cex=0.6)
    } else{
      legend("topright", legend=c(paste0("Spec: ", Llim, "-", Ulim, attr.unit),
                                  paste0("TI ", (1-alpha)*100, "% Confidence", coverprob*100, "% Coverage")), bg="white",
             col=c("1","2"), lty=c(1,1), cex=0.6)
    }
    title(paste0(attr.name, " OC curves for ", n.batch, " future result \n Assume process mean = ", mu, attr.unit) , cex=0.9)
  }

  else if (length(mu)!=1 & length(sigma)==1){
    tipp <- sapply(X=mu, FUN = ti.pp, sigma=sigma, n=n, Llim=Llim, Ulim=Ulim, n.batch=n.batch, alpha=alpha, side=side)
    plot(mu, rlpp, xlab=paste0("Mean Value (",attr.unit, ")"),
         ylab="Probability of Passing", type="l", lty=1, col=1, lwd=2,
         ylim=c(0,1)) # plot release test
    lines(mu, tipp, lty=1, col=2, lwd=2) # add PPQ Assessment line


    bl.mu <- mu[order(abs(rlpp-0.9))[1:2]]  # search the std dev for baseline plan
    hp.mu <- mu[order(abs(rlpp-0.95))[1:2]] # search the std dev for high performance plan
    abline(h=0.95, lty=2, col=16)
    abline(h=0.90, lty=2, col=16)
    abline(h=0.20, lty=2, col=16)
    abline(h=0.05, lty=2, col=16)
    abline(v=bl.mu, lty=2, col=16)
    abline(v=hp.mu, lty=2, col=16)
    points(bl.mu[1],0.2, pch="*", col=4, cex=3) # mark baseline
    points(bl.mu[2],0.2, pch="*", col=4, cex=3) # mark baseline
    points(hp.mu[1],0.05, pch="*", col=3, cex=3) # mark high performance
    points(hp.mu[2],0.05, pch="*", col=3, cex=3) # mark high performance

    if (add.reference == TRUE){
      bl.level <- optimize(f = function(x){abs(ti.pp(mu = bl.mu[1], sigma=sigma, n = n, Llim = Llim, Ulim = Ulim, n.batch = n.batch, alpha = x, coverprob=coverprob, side=side)-0.2)}, interval = c(0,0.5), tol = 1e-6)$minimum # optimize the alpha for baseline
      bl.tipp <- sapply(X=mu, FUN = ti.pp, sigma=sigma, n=n, Llim=Llim, Ulim=Ulim, n.batch=n.batch, alpha=bl.level, coverprob=coverprob, side=side)
      lines(mu, bl.tipp, lty=2, col=4, lwd=2) # add PPQ Assessment baseline
      hp.level <- optimize(f = function(x){abs(ti.pp(mu = hp.mu[1], sigma=sigma, n = n, Llim = Llim, Ulim = Ulim, n.batch = n.batch, alpha = x, coverprob=coverprob, side=side)-0.05)}, interval = c(0,0.5), tol = 1e-6)$minimum # optimize the alpha for baseline
      hp.tipp <- sapply(X=mu, FUN = ti.pp, sigma=sigma, n=n, Llim=Llim, Ulim=Ulim, n.batch=n.batch, alpha=hp.level, coverprob=coverprob, side=side)
      lines(mu, hp.tipp, lty=2, col=3, lwd=2) # add PPQ Assessment high performance line
      legend("center", legend=c(paste0("Spec: ", Llim, "-", Ulim, attr.unit),
                                paste0("TI ", (1-alpha)*100, "% Confidence", coverprob*100, "% Coverage"),
                                paste0("TI ", round((1-bl.level)*100,2), "%/", coverprob*100, "% Baseline"),
                                paste0("TI ", round((1-hp.level)*100,2), "%/", coverprob*100, "% High Perform")), bg="white",
             col=c(1,2,4,3), lty=c(1,1,2,2), cex=0.6)
    } else{
      legend("center", legend=c(paste0("Spec: ", Llim, "-", Ulim, attr.unit),
                                paste0("TI ", (1-alpha)*100, "% Confidence", coverprob*100, "% Coverage")), bg="white",
             col=c("1","2"), lty=c(1,1), cex=0.6)
    }
    title(paste0(attr.name, " OC curves for ", n.batch, " future result \n Assume process standard deviation = ", sigma, attr.unit) , cex=0.9)
  } else {stop("mu and sigma should be one single value and one vector!")}
}



#' Heatmap/Contour Plot for Assessing Power of the PPQ Plan using Tolerance Interval.
#'
#' The function for plotting the heatmap to evaluate the PPQ plan based on the specification test, given lower and upper specification limits.
#'
#' @usage ti.ctplot(attr.name, attr.unit, Llim, Ulim, mu, sigma, n, n.batch,
#' alpha, coverprob, side, test.point)
#' @param attr.name user-defined attribute name for PPQ assessment
#' @param attr.unit user-defined attribute unit
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param n sample size (number of locations) per batch
#' @param n.batch number of batches for passing PPQ during validation
#' @param alpha significant level for constructing the tolerance interval.
#' @param coverprob converage probability for constructing the tolerance interval
#' @param side whether a 1-sided or 2-sided tolerance interval is required (determined by side = 1 or side = 2, respectively).
#' @param test.point (optional) actual process data points for testing whether the processes pass PPQ
#' @return
#' Heatmap (or Countour Plot) for PPQ Assessment.
#' @seealso \code{ti.pp} and \code{ti.occurve}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#'
#' mu <- seq(95,105,0.1)
#' sigma <- seq(0.1,2.5,0.1)
#' ti.ctplot(attr.name = "Sterile Concentration Assay", attr.unit = "%",
#' mu = mu, sigma = sigma, Llim=95, Ulim=105)
#'
#' ti.ctplot(attr.name = "Extractable Volume", attr.unit = "% of NV=1mL",
#' Llim = 100, Ulim = Inf, mu=seq(100, 110, 0.5), sigma=seq(0.2, 15 ,0.5), n=40,
#' alpha = 0.05, coverprob = 0.675, side=1)
#'
#' @export

ti.ctplot <- function(attr.name="", attr.unit="", Llim, Ulim, mu, sigma, n=10, n.batch=1, alpha=0.05, coverprob=0.675, side=2, test.point=c()){

  para <- expand.grid(mu,sigma)
  ct.df <- data.frame(para, sapply(1:nrow(para), function(i) ti.pp(mu = para[i,1], sigma=para[i,2], n = n, n.batch = n.batch, Llim = Llim, Ulim = Ulim, alpha = alpha, coverprob=coverprob, side=side)))
  colnames(ct.df) <- c("Mean", "Std Dev", "Passing Probability")

  ct.mat<-matrix(ct.df$`Passing Probability`, nrow=length(mu), ncol=length(sigma), byrow=FALSE)

  filled.contour(mu,sigma,ct.mat,levels=c(0, 0.8, 0.9, 0.95, 0.99,1), # ylim=c(0.05,0.27),
                 col=c("red", "orange", "yellow", "light green", "green"), ylab=paste0("Standard Deviation",  " (", attr.unit, ")"),
                 xlab=paste0("Mean for ", attr.name, " (", attr.unit, ")"),
                 main=paste0("Probability of Passing ", round(100*(1-alpha),2), "%/", coverprob*100, "% TI Criteria \n", attr.name, " Spec: ", Llim, "-", Ulim, attr.unit),
                 plot.axes={axis(1); axis(2); points(x=test.point[,1], y=test.point[,2], pch=8)},
                 key.axes=axis(4, at=c(0.80, 0.90, 0.95, 0.99)))
}
