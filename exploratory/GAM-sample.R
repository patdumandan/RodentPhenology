tmpf <- tempfile()
download.file("https://gist.github.com/gavinsimpson/b52f6d375f57d539818b/raw/2978362d97ee5cc9e7696d2f36f94762554eefdf/load-process-cet-monthly.R",
              tmpf)
source(tmpf)
ls()

ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
m2 <- gamm(Temperature ~ s(nMonth, bs = "cc", k = 12) + s(Time, k = 20),
           data = cet, correlation = corARMA(form = ~ 1|Year, p = 2),
           control = ctrl)

want <- seq(1, nrow(cet), length.out = 200)
pdat <- with(cet,
             data.frame(Time = Time[want], Date = Date[want],
                        nMonth = nMonth[want]))
p2 <- predict(m2$gam, newdata = pdat, type = "terms", se.fit = TRUE)
pdat <- transform(pdat, p2 = p2$fit[,2], se2 = p2$se.fit[,2])

df.res <- df.residual(m2$gam)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat <- transform(pdat,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))

Term <- "Time"
m2.d <- Deriv(m2)
m2.dci <- confint(m2.d, term = Term)
m2.dsig <- signifD(pdat$p2, d = m2.d[[Term]]$deriv,
                     +                    m2.dci[[Term]]$upper, m2.dci[[Term]]$lower)

ylim <- with(pdat, range(upper, lower, p2))
ylab <- expression(Temperature ~ (degree*C * ":" ~ centred))
 
 plot(p2 ~ Date, data = pdat, type = "n", ylab = ylab, ylim = ylim)
lines(p2 ~ Date, data = pdat)
lines(upper ~ Date, data = pdat, lty = "dashed")
lines(lower ~ Date, data = pdat, lty = "dashed")
lines(unlist(m2.dsig$incr) ~ Date, data = pdat, col = "blue", lwd = 3)
lines(unlist(m2.dsig$decr) ~ Date, data = pdat, col = "red", lwd = 3)