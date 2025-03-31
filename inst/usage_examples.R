library(devtools)
install_github("C2S2-HF/InverseL0TF", repos = NULL, type = "source")

library(L0TFinv)

### Figure 1 and Figure 2---------------------------

## generate data
tau = c(0.1, 0.25, 0.3, 0.4, 0.7, 0.85, 0.95)
h = c(-1, 5, 3, 0, -1, 2, 0, -1)
BlocksData <- SimuBlocksInv(n = 500, sigma = 1, seed = 25, tau = tau ,h = h)
print(BlocksData$setA)
## fit the model: up to 20 change points, using the BIC penalty
res <- L0TFinv.opt(y=BlocksData$y, kmax=20, q=0, penalty="bic")
print(res$Aopt)

df <- data.frame(num = 1:500, input = BlocksData$y,
                 y0 = BlocksData$y0, fit = res$y.all[, res$kopt])

matplot(df$num, df[, c("y0", "fit", "input")],
        type = c("l","l","p"),
        lty = c("solid", "dashed", "blank"),
        col = c("black", "red", "grey"),
        lwd = 2,
        pch = c(NA, NA, 1),
        cex = c(1, 1, 0.1),
        xlab = "Sample Size", ylab = "Function value",
        main = "")

legend("topright",
       legend = c("True trend", "Fitted trend", "Observation"),
       col = c("black", "red", "grey"),
       lty = c("solid", "dashed", "blank"),
       pch = c(NA, NA, 1),
       lwd = 2,
       bty = "n",
       cex = 0.8)




