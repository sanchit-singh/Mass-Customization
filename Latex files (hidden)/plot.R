pdf("comparsion_all_4.pdf",width=6, height=4)
fun1 <- function(x) 175.0+0*x
fun2 <- function(x) 166+10*x
fun3 <- function(x) 163+8*x
fun4 <- function(x) 174.78+0.44*x
x <- seq(0, 5, 0.001)
matplot(x,fun1(x),type="l",lwd='1',col="black",ylim=c(140,240),xlab="h",ylab="Cost")
curve(fun2(x), add = TRUE, col = "red")
curve(fun3(x), add = TRUE, col = "green")
curve(fun4(x), add = TRUE, col = "blue")
legend(0,240, c("MTO", "MTS", "ATO", "ATO-CS-hybrid"), col = c("black", "red", "green", "blue"), lwd = 1, pch = c(NA, NA, NA, NA), bg = "gray90")
dev.off()

