#' Heatmap/Contour Plot for Assessing Power of the CQA PPQ Plan Using General Multiplier.
#'
#' The function for plotting the heatmap to evaluate the PPQ plan based on the specification test, given lower and upper specification limits.
#'
#' @usage PPQ_ctplot(attr.name, attr.unit, Llim, Ulim, mu, sigma, n, n.batch, k, test.point)
#' @param attr.name (optional) user-defined attribute name for PPQ assessment
#' @param attr.unit (optional) user-defined attribute unit
#' @param Llim lower specification limit
#' @param Ulim upper specification limit
#' @param mu hypothetical mean of the attribute
#' @param sigma hypothetical standard deviation of the attribute
#' @param n sample size (number of locations) per batch
#' @param n.batch number of batches for passing PPQ during validation
#' @param k general multiplier for constructing the specific interval
#' @param test.point (optional) actual process data points for testing whether the processes pass PPQ
#' @return
#' Heatmap (or Contour Plot) for PPQ Assessment.
#' @seealso \code{PPQ_pp} and \code{PPQ_occurve}.
#' @references
#' Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
#' Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
#' \emph{Springer}.
#' @author Yalin Zhu
#' @examples
#' \dontrun{
#' mu <- seq(1.6,3.4,0.05)
#' sigma <- seq(0.05,0.8,0.01)
#' PPQ_ctplot(attr.name = "Total Protein", attr.unit = "mg/mL", Llim=1.5, Ulim=3.5,
#' mu = mu, sigma = sigma, k=2.373)
#'
#' ## Example verifying simulation resutls in the textbook page 249
#' mu <- seq(95, 105, 0.1)
#' sigma <- seq(0.2, 5, 0.1)
#' PPQ_ctplot(attr.name = "Composite Assay", attr.unit = "%LC", Llim=95, Ulim=105,
#' mu = mu, sigma = sigma, k=2.373)
#' mu <- seq(90, 110, 0.5)
#' PPQ_ctplot(attr.name = "Composite Assay", attr.unit = "%LC", Llim=90, Ulim=110,
#' mu = mu, sigma = sigma, k=2.373)
#'
#' mu <- seq(95,105,0.1)
#' sigma <- seq(0.1,2.5,0.1)
#' PPQ_ctplot(attr.name = "Sterile Concentration Assay", attr.unit = "%", Llim=95, Ulim=105,
#' mu = mu, sigma = sigma, k=2.373)
#' test <- data.frame(mean=c(97,98.3,102.5), sd=c(0.55, 1.5, 1.2))
#' PPQ_ctplot(attr.name = "Sterile Concentration Assay", attr.unit = "%", Llim=95, Ulim=105,
#' mu = mu, sigma = sigma, k=2.373, test.point=test)
#' }
#' @import graphics
#' @export

PPQ_ctplot <- function(attr.name="", attr.unit="", Llim, Ulim, mu, sigma, n=10, n.batch=1, k, test.point=c()){

  para <- expand.grid(mu,sigma)
  ct.df <- data.frame(para, sapply(1:nrow(para), function(i) PPQ_pp(mu = para[i,1], sigma=para[i,2], n = n, n.batch = n.batch, Llim = Llim, Ulim = Ulim, k = k)))
  colnames(ct.df) <- c("Mean", "Std Dev", "Passing Probability")
  ct.mat<-matrix(ct.df$`Passing Probability`, nrow=length(mu), ncol=length(sigma), byrow=FALSE)

  filled.contour(mu,sigma,ct.mat,levels=c(0, 0.8, 0.9, 0.95, 0.99,1), # ylim=c(0.05,0.27),
                 col=c("red", "orange", "yellow", "light green", "green"), ylab=paste0("Standard Deviation",  " (", attr.unit, ")"),
                 xlab=paste0("Mean for ", attr.name, " (", attr.unit, ")"),
                 main=paste0("Heatmap for ", attr.name, "\nLSL = ", Llim, attr.unit, ", USL = ", Ulim, attr.unit,  ", k = ", k),
                 plot.axes={axis(1); axis(2); points(x=test.point[,1], y=test.point[,2], pch=8)},
                 key.axes=axis(4, at=c(0.80, 0.90, 0.95, 0.99)))
}

