calc_pcareduceres <- function(obj, reduce) {
  if (!reduce) {
    return(t(obj))
  }
  df <- data.frame(x = 1:20)
  df$sdev <- prcomp(t(obj), scale = TRUE)$sdev[1:20]
  optpoint <- which.min(
    vapply(
      2:10,
      function(i) {
        data$x2 <- pmax(0, data$x - i)
        sum(lm(sdev ~ x + x2, df)$residuals^2)
      },
      0
    )
  )
  pcadim <- optpoint + 1
  tmpdata <- t(apply(obj, 1, scale))
  colnames(tmpdata) <- colnames(obj)
  tmppc <- prcomp(t(tmpdata), scale = TRUE)
  return(t(tmpdata) %*% tmppc$rotation[, 1:pcadim])
}
