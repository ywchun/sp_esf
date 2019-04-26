library(RColorBrewer)
library(classInt)

mapping.seq <- function(polys, x, n.class, main="") {  
  pal.red <- brewer.pal(n.class, "Reds")
  q.n <- classIntervals(x, n.class, style="quantile") 
  cols.red <- findColours(q.n, pal.red)
  plot(polys, col=cols.red)
  brks <- round(q.n$brks,3)
  leg <- paste(brks[-(n.class+1)], brks[-1], sep=" - ")
  legend("bottomright", fill=pal.red, legend=leg, bty="n")
  if (!missing(main)) title(main)
}  

mapping.bipolar <- function(polys, x, n.class, main="", legend.loc="bottomright"){
  if ((n.class %% 2) == 1){
    cat(paste("The number of classes is set to an even number." , n.class+1,
                  " classes will be used!!"), "\n",sep="")
      n.class <- n.class + 1
  }
  ql <- quantile(x[x < mean(x)], seq(0,1,2/n.class))
  qr <- quantile(x[x >= mean(x)], seq(0,1,2/n.class))
  qn <- c(ql[-length(ql)], mean(x), qr[-1])

  ncl <- length(qn)-1
  n.clp <- length(qn) 
  pal <- rev(brewer.pal(ncl,"RdYlGn"))

  brks <- round(qn,3)
  cols <- pal[findInterval(x, qn, rightmost.closed=T)]
  legs <- paste(brks[-(n.class+1)], brks[-1], sep=" - ")

  plot(polys, col=cols, main=main, axes=FALSE)
  legend(legend.loc, legend=legs, fill=pal, bty="n",ncol=1)
}

stepwise <- function(full, initial, sle=0.1, sls=0.1001, verbose=T) {
      cut.string <- function(string) {
        if (length(string) > 1L)
            string[-1L] <- paste("\n", string[-1L], sep = "")
        string
      }
      msef <- sum(full$residuals^2)/full$df
      n <- length(full$residuals)
      allvars <- attr(full$terms, "predvars")
      current <- initial
      while (TRUE) {
        temp <- summary(current)
        rnames <- names(current$coefficient)
        if (verbose) print(current$coefficients)
        p <- current$rank
        mse <- sum(current$residuals^2)/current$df
        if (verbose) {
          cat("S = ", temp$sigma, " R-sq = ", temp$r.squared, "\n\n")
        }
		
		if (p > initial$rank) {  
          d <- drop1(current, test="F")
          pmax <- max(d[-c(1:initial$rank),6])
          if (pmax > sls) {
            var.del <- rownames(d)[d[,6] == pmax];  # name of variable to delete
            if (length(var.del) > 1) var.del <- var.del[2]			
            current <- update(current, paste("~ . - ", var.del), evaluate=FALSE)
            current <- eval.parent(current)
			if (verbose) {
              cat("\nStep:  \n", cut.string(deparse(as.vector(formula(current)))), "\n\n", sep = "")
              utils::flush.console()
			}
            next
          }
        }

        a <- tryCatch(add1(current, scope=full, test="F"), error=function(e) NULL);
        if (is.null(a)) break;

        pmin <- min(a[-1,6])
        if (pmin < sle) {
          var <- rownames(a)[a[,6] == pmin];
          if (length(var) > 1) {
            var <- var[2]
          }
          current <- update(current, paste("~ . + ", var), evaluate=FALSE)
          current <- eval.parent(current)
          if (verbose) {
            cat("\nStep:  \n", cut.string(deparse(as.vector(formula(current)))), "\n\n", sep = "")
            utils::flush.console()
          }
          next
        }
        break
      }
      return(current)
}


