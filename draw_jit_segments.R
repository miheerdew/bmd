draw_jit_segments <- function (Df, jitscale = 1, 
                               colors = rep("black", nrow(Df)), ...) {
  
  reorder <- sample(nrow(Df))
  Df <- Df[reorder, ]
  colors <- colors[reorder]
  
  for (i in 1:nrow(Df)) {
    
    x1 <- Df$x1[i]
    
    newcp <- c(x1/2, 0)
    while (newcp[1] <= x1 || newcp[1] >= x2) {
      x1 <- Df$x1[i]
      y1 <- Df$y1[i]
      x2 <- Df$x2[i]
      y2 <- Df$y2[i]
      switch <- FALSE
      if (y2 <= y1) {
        switch <- TRUE
        y0 <- y1
        y1 <- y2
        y2 <- y0
      }
      
      x <- rnorm(1)
      jit <- x * jitscale
      #cat(jit, "\n")
      hypot <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
      cpx <- hypot / 2
      cpy <- jit
      #cat(cpy, "\n")
      hypotv <- c(x2 - x1, y2 - y1)
      theta <- acos(hypotv[1] / hypot)
      rmat <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol = 2)
      newcp <- rmat %*% (c(cpx, cpy)) + c(x1, y1)
      #cat(newcp, "\n")
      
      if (switch) {
        y0 <- y1
        y1 <- y2
        y2 <- y0
        newcp[1] <- x1 + x2 - newcp[1]
      }
    }
    
    segments(x0 = c(x1, newcp[1]),
             y0 = c(y1, newcp[2]),
             x1 = c(newcp[1], x2),
             y1 = c(newcp[2], y2), col = colors[i], ...)
    #cat("----\n")
  }
  
}