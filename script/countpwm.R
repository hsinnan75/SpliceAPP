suppressMessages({
library(Biostrings)

args <- commandArgs(T)

 human_PWM <- scan("data/selected_pwm.txt", what = "", quiet = TRUE)
 RBP.name <- as.character(sub(">","",human_PWM[grep(">",human_PWM)]))
 Sep.length <- as.numeric(human_PWM[grep(">",human_PWM)+1])
 grp.index <- unlist(lapply(seq_along(RBP.name), FUN = function(x){return(rep(RBP.name[x], Sep.length[x]*4+2))}))
 RBP.list <- split(human_PWM, grp.index)[RBP.name]
 RBP.matrix <- sapply(RBP.list, function(x){
   temp <- matrix(as.numeric(x[-c(1,2)]), nrow = 4, byrow = F)
   row.names(temp) <- c("A","C","G","T")
   return(temp)
 })

#RBP.matrix <- readRDS(RBP_matrix.rds)

pwm.res <- lapply(RBP.matrix, function(x) countPWM(pwm = x, subject = args[1], min.score = "95%"))
pwm.res <- do.call(rbind, pwm.res)
pwm.res <- pwm.res[pwm.res[, 1] > 0, ]
write.table(pwm.res, "pwm_motif.tsv", sep = "\t", col.names = F, quote = F)

})
