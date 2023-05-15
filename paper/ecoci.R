## ----"knitr-setup", echo = FALSE----------------------------------------------
## knitr options
library(knitr)
opts_chunk$set(fig.height = 4,
               echo = FALSE,
               warning = FALSE,
               message = FALSE,
               cache = FALSE,
               eval = TRUE)

## should sessionInfo be printed at the end?
Reproducibility <- TRUE

## size of the code font
codeFont <- "normalsize"


## ----"main-setup", include = FALSE--------------------------------------------
## Packages
library(ciCalibrate) # support intervals
library(lamW) # Lambert W function
library(colorspace) # color palettes
library(ggplot2) # plotting

## Function to nicely format Bayes factors
.formatBF_ <- function(BF, digits = "default") {
    ## check inputs
    stopifnot(
        length(BF) == 1,
        is.numeric(BF),
        (is.finite(BF) && 0 < BF) || is.na(BF),

        length(digits) == 1,
        (is.character(digits) && digits == "default") ||
        (is.numeric(digits) && 0 <= digits)
    )
    ## return NA if input NA/NaN
    if (is.na(BF) || is.nan(BF))
        result <- NA
    else {
        ## format BF
        if (digits == "default") {
            if (BF < 1/1000)
                result <- "< 1/1000"
            if ((BF >= 1/1000) & (BF <= 1/10))
                result <- paste0("1/", as.character(round(1/BF)))
            if ((BF > 1/10) & (BF < 1))
                result <- paste0("1/", as.character(round(1/BF, digits = 1)))
            if ((BF < 10) & (BF >= 1))
                result <- as.character(round(BF, digits = 1))
            if ((BF >= 10) & (BF <= 1000))
                result <- as.character(round(BF))
            if (BF > 1000)
                result <- "> 1000"
        } else {
            if (BF < 1)
                result <- paste0("1/", as.character(round(1/BF, digits = digits)))
            else
                result <- as.character(round(BF, digits = digits))
        }
        ## when 1/1 return 1
        if (result == "1/1") result <- "1"
    }
    return(result)
}
formatBF <- Vectorize(FUN = .formatBF_)


## ----"figure-BFfun-RECOVERY", fig.height = 3.5--------------------------------
## data from RECOVERY trial abstract
## https://www.nejm.org/doi/10.1056/NEJMoa2021436
HR <- 0.83
ciHR <- c(0.75, 0.93)
logHR <- log(HR)
selogHR <- diff(log(ciHR))/(2*qnorm(p = 0.975))

## "power to detect a clinically relevant proportional reduction of 20% (an
## absolute difference of 4 percentage points)"
m <- log(0.8)
v <- 4 # unit information variance for a logHR
nevents <- 4/selogHR^2 # implicit event count

## compute k SI for different k
k <- c(1/10, 1, 10)
siList <- lapply(X = k, FUN = function(k) {
    si <- ciCalibrate(estimate = logHR, se = selogHR, siLevel = k,
                      method = "SI-normal", priorMean = m, priorSD = sqrt(v))
})
siDF <- do.call("rbind", lapply(X = siList, FUN = function(x) {
    data.frame(k = x$siLevel, lower = x$si[1], upper = x$si[2],
               est = x$estimate)
}))

## compute BF function
logHRseq <- seq(log(0.65), log(1.05), length.out = 500)
bfFun <- siList[[1]]$bfFun
plotDF <- data.frame(logHR = logHRseq, HR = exp(logHRseq),
                     BF = bfFun(x = logHRseq))

## plot BFfun and k SI
bfBks <- c(1, 3, 10, 30, 100, 300, 1000)
bfBks <- c(1, 10, 100, 1000)
bfBks <- c(1/bfBks, bfBks)
bfMax <- bfFun(x = logHR)
plot <- ggplot(data = plotDF, aes(x = logHR, y = BF)) +
    geom_ribbon(aes(ymin = 0, ymax = BF), alpha = 0.05, fill = 1) +
    geom_line(size = 0.3) +
    annotate(geom = "segment", x = rep(-0.42, 2), xend = rep(-0.42, 2),
             y = 1.1^(c(-1, 1)), yend = 5^(c(-1, 1)),
             arrow = arrow(length = unit(0.07, "inches")),
             lineend = "round", linejoin = "round", alpha = 0.6) +
    annotate(geom = "text", x = rep(-0.4075, 2), y = 2^(c(-1, 1)),
             label = c("italic('H')[1]", "italic('H')[0]"), alpha = 0.6,
             parse = TRUE, size = 4) +
    geom_errorbarh(data = siDF, aes(xmin = lower, xmax = upper, y = k),
                   height = 0.15, size = 0.3) +
    annotate(geom = "point", x = logHR, y = bfMax, size = 0.7) +
        annotate(geom = "text", x = logHR, y = bfMax*1.5, size = 3,
                 label = "'MLE' ~ hat(theta)", parse = TRUE) +
    geom_text(data = siDF,
              aes(x = est, y = k*1.4,
                  label = paste0("italic(k) == ", formatBF(BF = k))),
              parse = TRUE,
              size =  4) +
    scale_y_log10(breaks = bfBks, labels = formatBF(BF = bfBks)) +
    coord_cartesian(ylim = c(1/100, 100), xlim = c(-0.41, 0.035)) +
    labs(x = bquote("Log hazard ratio" ~ theta[scriptstyle("0")]),
         y = bquote("BF"["01"] * "(data; " * theta[scriptstyle("0")] * ")")) +
    theme_bw()
plot +
    theme(panel.grid.minor = element_blank())


## ----"existence-nonlocal-SI"--------------------------------------------------
## W_0 function has to be larger than 1/2 for support interval to exist
rootFun <- function(x) lamW::lambertW0(x = x) - 0.5
lArg <- uniroot(f = rootFun, interval = c(0, exp(1)))$root


## ----"compare-priors-growth", fig.height = 5----------------------------------
## names of priors, levels, and colors for plot
priorlevels <- c("Normal (wrong mean)", "Normal (correct mean)", "Local normal",
                 "Nonlocal normal moment")
priorcols <- c("#000000", "#009E73", "#E69F00", "#56B4E9")
names(priorcols) <- priorlevels
priorlty <- c("dotted", "dashed", "solid", "dotdash")
names(priorlty) <- priorlevels

## function to compute width of support interval
Width. <- function(n, priorvar, unitvar, priormean = NULL, t, method, k) {
    width <- diff(ciCalibrate(estimate = t, se = sqrt(unitvar/n), siLevel = k,
                              priorMean = priormean, priorSD = sqrt(priorvar),
                              method = method)$si)
}
Width <- Vectorize(FUN = Width., vectorize.args = c("n", "priorvar"))

## find prior parameters such that the same width for a certain sample size
n <- 2 # sample size where width should match
k <- 1 # support level
unitvar <- 4 # unit variance
priorvarNorm <- unitvar # unit-information prior
t <- 0 # observed effect estimate
priormean <- t + sqrt(unitvar) # prior mean off by one unit standard deviation
priormean2 <- t # prior mean for correct normal prior
priorvarNorm2 <- unitvar ## uniroot(f = function(priorvarNorm2) {
##     Width(n = n, k = k, priorvar = priorvarNorm, unitvar = unitvar,
##           priormean = priormean, t = t, method = "SI-normal") -
##         Width(n = n, k = k, priorvar = priorvarNorm2, unitvar = unitvar,
##           priormean = priormean2, t = t, method = "SI-normal")
## }, interval = c(0.001, 30))$root
priorvarLocal <- unitvar ## uniroot(f = function(priorvarLocal) {
##     Width(n = n, k = k, priorvar = priorvarNorm, unitvar = unitvar,
##           priormean = priormean, t = t, method = "SI-normal") -
##         Width(n = n, k = k, priorvar = priorvarLocal, unitvar = unitvar,
##               t = t, method = "SI-normal-local")
## }, interval = c(0.1, 10))$root
priorvarNonlocal <- unitvar ## uniroot(f = function(priorvarNonlocal) {
##     Width(n = n, k = k, priorvar = priorvarNorm, unitvar = unitvar,
##           priormean = priormean, t = t, method = "SI-normal") -
##         Width(n = n, k = k, priorvar = priorvarNonlocal, unitvar = unitvar,
##               t = t, method = "SI-normal-nonlocal")
## }, interval = c(0.1, 5))$root

## compute SI width
nseq <- exp(seq(log(n), log(3000), length.out = 1000))
widthNormWrong <- Width(n = nseq, k = k, priorvar = priorvarNorm,
                        unitvar = unitvar, priormean = priormean, t = t,
                        method = "SI-normal")
widthNormCorrect <- Width(n = nseq, k = k, priorvar = priorvarNorm2,
                          unitvar = unitvar, priormean = priormean2, t = t,
                          method = "SI-normal")
widthLocal <- Width(n = nseq, k = k, priorvar = priorvarLocal,
                    unitvar = unitvar, t = t, method = "SI-normal-local")
widthNonlocal <- Width(n = nseq, k = k, priorvar = priorvarNonlocal,
                       unitvar = unitvar, t = t, method = "SI-normal-nonlocal")
plotDF <- rbind(data.frame(n = nseq, k = k, width = widthNormWrong,
                           type = priorlevels[1]),
                data.frame(n = nseq, k = k, width = widthNormCorrect,
                           type = priorlevels[2]),
                data.frame(n = nseq, k = k, width = widthLocal,
                           type = priorlevels[3]),
                data.frame(n = nseq, k = k, width = widthNonlocal,
                           type = priorlevels[4]))
plotDF$type <- factor(x = plotDF$type, levels = priorlevels)
plotA <- ggplot(data = plotDF, aes(x = n, y = width, color = type, linetype = type)) +
    geom_line(alpha = 0.7, size = 0.7) +
    scale_x_log10(breaks = c(2, 10, 100, 1000)) +
    expand_limits(y = 0) +
    labs(color = "", linetype = "", x = bquote("Sample size" ~ italic(n)),
         y = bquote("Width of" ~ italic(k) == .(k) ~ "support interval")) +
    scale_color_manual(values = priorcols) +
    scale_linetype_manual(values = priorlty) +
    guides(linetype = guide_legend(override.aes = list(linewidth = 0.5)),
           color = guide_legend(override.aes = list(linewidth = 0.5))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

## for a given sample size compute the highest support level for which an SI is
## non-empty (i.e., the BF function evaluated at the observed effect estimate)
kMax. <- function(n, priorvar, unitvar, priormean, t, method) {
    k <- 1
    ciCalibrate(estimate = t, se = sqrt(unitvar/n), siLevel = k,
                priorMean = priormean, priorSD = sqrt(priorvar),
                method = method)$bfFun(t)
}
kMax <- Vectorize(FUN = kMax., vectorize.args = c("n", "priorvar"))

## choose prior parameters so that the same highest support level for n
n <- 2
priorvarNorm2 <- unitvar ## uniroot(f = function(priorvarNorm2) {
##     kMax(n = n, priorvar = priorvarNorm, unitvar = unitvar,
##          priormean = priormean, t = t, method = "SI-normal") -
##         kMax(n = n, priorvar = priorvarNorm2, unitvar = unitvar,
##              priormean = priormean2, t = t, method = "SI-normal")
## }, interval = c(0.001, 20))$root
priorvarLocal <- unitvar ## uniroot(f = function(priorvarLocal) {
##     kMax(n = n, priorvar = priorvarNorm, unitvar = unitvar,
##          priormean = priormean, t = t, method = "SI-normal") -
##         kMax(n = n, priorvar = priorvarLocal, unitvar = unitvar,
##              priormean = priormean, t = t, method = "SI-normal-local")
## }, interval = c(0.1, 10))$root
priorvarNonlocal <- unitvar ## uniroot(f = function(priorvarNonlocal) {
##     kMax(n = n, priorvar = priorvarNorm, unitvar = unitvar,
##          priormean = priormean, t = t, method = "SI-normal") -
##         kMax(n = n, priorvar = priorvarNonlocal, unitvar = unitvar,
##              priormean = priormean, t = t, method = "SI-normal-nonlocal")
## }, interval = c(0.1, 5))$root

## compute highest support level
kMaxnorm1 <- kMax(n = nseq, priorvar = priorvarNorm, unitvar = unitvar,
                  priormean = priormean, t = t, method = "SI-normal")
kMaxnorm2 <- kMax(n = nseq, priorvar = priorvarNorm2, unitvar = unitvar,
                  priormean = priormean2, t = t, method = "SI-normal")
kMaxloc <- kMax(n = nseq, priorvar = priorvarLocal, unitvar = unitvar,
                priormean = priormean2, t = t, method = "SI-normal-local")
kMaxnonloc <- kMax(n = nseq, priorvar = priorvarNonlocal, unitvar = unitvar,
                   priormean = priormean2, t = t, method = "SI-normal-nonlocal")
plotDF2 <- rbind(data.frame(n = nseq, kMax = kMaxnorm1, type = priorlevels[1]),
                 data.frame(n = nseq, kMax = kMaxnorm2, type = priorlevels[2]),
                 data.frame(n = nseq, kMax = kMaxloc, type = priorlevels[3]),
                 data.frame(n = nseq, kMax = kMaxnonloc, type = priorlevels[4]))
plotDF2$type <- factor(x = plotDF2$type, levels = priorlevels)
plotB <- ggplot(data = plotDF2, aes(x = n, y = kMax, color = type, linetype = type)) +
    geom_line(alpha = 0.7, size = 0.7) +
    scale_x_log10(breaks = c(2, 10, 100, 1000)) +
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
    expand_limits(y = 1) +
    labs(color = "Prior", linetype = "Prior", x = bquote("Sample size" ~ italic(n)),
         y = bquote("Highest support level" ~ italic(k))) +
    scale_color_manual(values = priorcols) +
    scale_linetype_manual(values = priorlty) +
    guides(linetype = guide_legend(override.aes = list(linewidth = 0.25)),
           color = guide_legend(override.aes = list(linewidth = 0.25))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
ggpubr::ggarrange(plotA, plotB, ncol = 1, common.legend = TRUE, align = "hv")


## ----"mapping-conf-minsupport"------------------------------------------------
## functions to compute minimum support level k from confidence level such that
## the same standard error multiplier
## -----------------------------------------------------------------------------
## class of all priors
conf2kall <- function(conf) {
    exp(-0.5*qnorm(p = (1 + conf)*0.5)^2)
}
## class of local normal priors
conf2knormallocal <- function(conf) {
    z <- qnorm(p = (1 + conf)*0.5)
    abs(z)*exp(-0.5*(z^2 - 1))
}
## -eplogp
conf2keplogp <- function(conf) {
    k <- -exp(1)*log(1 - conf)*(1 - conf)
    return(k)
}
conf2k <- list(conf2kall, conf2keplogp, conf2knormallocal)
conf2kName <- c("italic(k) ~ 'minimum support (all priors)'",
                "italic(k) ~ 'minimum support ('*-italic('e') *italic('p') * 'log' *italic('p')* ')'",
                "italic(k) ~ 'minimum support (local normal priors)'")

## which k minimum SI correspond to 95% CI?
k95 <- sapply(X = conf2k, FUN = function(f) f(0.95))

## functions to compute confidence level from  minimum support level such that
## the same standard error multiplier
## -----------------------------------------------------------------------------
## class of all priors
k2confall <- function(k) {
    2*pnorm(q = sqrt(-2*log(k))) - 1
}
## class of local normal priors
k2confnormallocal <- function(k) {
    2*pnorm(q = sqrt(-lambertWm1(x = -k^2/exp(1)))) - 1
}
## -eplogp
k2confeplogp <- function(k) {
    2*pnorm(q = qnorm(p = 1 - exp(lambertWm1(-k/exp(1)))/2)) - 1
}
k2conf <- list(k2confall, k2confeplogp, k2confnormallocal)

## which confidence level corresponds to k = 1/10 minimum support level?
conf10 <- sapply(X = k2conf, FUN = function(f) f(1/10))


## ----"figure-mapping-conf-minsupport", fig.height = 3.3-----------------------
## computing mapping for a grid of confidence levels
confseq <- seq(0.8, 0.99999, length.out = 2000)
mseq <- qnorm(p = 0.5*(1 + confseq))
plotDF <- do.call("rbind", lapply(X = seq_along(conf2k), FUN = function(i) {
    f <- conf2k[[i]]
    k <- f(conf = confseq)
    data.frame(type = conf2kName[i], conflevel = confseq, m = mseq, k = k)
}))
plotDF$type <- factor(x = plotDF$type, levels = conf2kName)

## plotting conflevel vs minimum support level
bfBks <- c(1/300, 1/100, 1/30, 1/10, 1/3, 1)
confbks <- seq(50, 100, 0.5)
ggplot(data = plotDF, aes(x = conflevel*100, y = k, color = type, linetype = type)) +
    geom_line(alpha = 0.8, size = 0.7) +
    ## scale_x_continuous(sec.axis = sec_axis(~ conf2z(./100), name = "SE multiplier",
    ##                                        breaks = c(1.7, 1.8, 2, 2.3, 3))) +
    scale_x_continuous(breaks = confbks) +
    scale_y_log10(breaks = bfBks, labels = formatBF(BF = bfBks),
                  expand = c(0, 0)) +
    expand_limits(y = 1.5) +
    coord_cartesian(xlim = c(95, 99.9), ylim = c(1/250, 1.3)) +
    labs(y = bquote("Minimum support level" ~ italic(k)),
         x = bquote("Confidence level" ~ (1 - alpha) * "100%"),
         color = "", linetype = "") +
    scale_color_brewer(palette = "Dark2", labels = scales::parse_format()) +
    scale_linetype_discrete(labels = scales::parse_format()) +
    theme_bw() +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          legend.text = element_text(size = 7.5))


## ----"eliciation-spread-nonlocal"---------------------------------------------
## elicitation of a suitable value for the spread parameter of the nonlocal
## moment prior as in Pramanik and Johnson (2022)
## -----------------------------------------------------------------------------

## density function of normal moment prior centered around m with spread s
dnormMoment <- function(x, m, s) {
    dnorm(x = x, mean = m, sd = s) * (x - m)^2/s^2
}
## cumulative distribution function of normal moment prior centered around m
## with spread s
pnormMoment_ <- function(q, m, s, lower.tail = TRUE) {
    p <- integrate(f = dnormMoment, lower = -Inf, upper = q, m = m, s = s)$value
    if (lower.tail == FALSE) p <- 1 - p
    return(p)
}
pnormMoment <- Vectorize(FUN = pnormMoment_)

## determine spread s such that 90% probability within [m - log2, m + log2]
pnonlocal <- 0.9
rootFun <- function(s) {
    p <- pnormMoment(q = log(2), m = 0, s = s) -
        pnormMoment(q = -log(2), m = 0, s = s)
    return(p - pnonlocal)
}
snonlocal <- uniroot(f = rootFun, lower = 0.1, upper = 10)$root
vnonlocal <- snonlocal^2

## ## look at resulting prior
## xseq <- seq(-1.5, 1.5, 0.001)
## plot(xseq, dnormMoment(x = xseq, m = 0, s = sNonlocal), type = "l",
##      xlab = bquote("Log hazard ratio" ~ theta), ylab = "Density")


## ----"compute-SI-RECOVERY"----------------------------------------------------
## compute different types of support intervals for different support levels k
## -----------------------------------------------------------------------------

## set up grid of methods and support levels
k <- c(1/10, 1/3, 1, 3, 10)
methods <- c("SI-normal","SI-normal-local","SI-normal-nonlocal", "mSI-all",
             "mSI-normal-local", "mSI-eplogp")
methodName <- c("italic(k) ~ 'SI (normal prior)'",
                "italic(k) ~ 'SI (local normal prior)'",
                "italic(k) ~ 'SI (nonlocal normal moment prior)'",
                "italic(k) ~ 'minSI (all priors)'",
                "italic(k) ~ 'minSI (local normal priors)'",
                "italic(k) ~ 'minSI ('*-italic('e') *italic('p') * 'log' *italic('p')* ')'")
names(methods) <- methodName
ciName <- "'95% CI'"
typeLevels <- rev(c(ciName, methodName))
applyGrid <- expand.grid(k = k, methodi = seq(1, length(methods)))
applyGrid$method <- methods[applyGrid$methodi]
applyGrid$name <- methodName[applyGrid$methodi]

## compute SIs
m <- log(0.8) # mean of normal prior
v <- 4 # unit information variance normal and local normal prior
vnonlocal <- 0.28^2 # spread parameter of nonlocal normal prior
siDF <- do.call("rbind", lapply(X = seq(1, nrow(applyGrid)), FUN = function(i) {
    k <- applyGrid$k[i]
    method <- applyGrid$method[i]
    name <- applyGrid$name[i]
    ## HACK prior parameters are ignored for methods that do not depend on them
    ## so we can set priorMean = m for all methods since only "SI-normal"
    ## depends on it. Is there a cleaner solution?
    if (method == "SI-normal-nonlocal") {
        priorSD <- sqrt(vnonlocal)
    } else {
        priorSD <- sqrt(v)
    }
    si <- ciCalibrate(estimate = logHR, se = selogHR, method = method,
                      siLevel = k, priorMean = m, priorSD = priorSD)
    data.frame(k = k, method = method, type = name, lower = si$si[1],
               upper = si$si[2], est = logHR)
}))

## combine SI data frame with CI
ciDF <- data.frame(k = NA, method = "CI", type = ciName,
                   lower = log(ciHR)[1], upper = log(ciHR[2]), est = logHR)
plotDF <- rbind(siDF, ciDF)
plotDF$type <- ordered(x = plotDF$type, levels = typeLevels)

## select some SI to discuss in the text
si10norm <- subset(x = plotDF, k == 10 & method == "SI-normal")
si10normlocal <- subset(x = plotDF, k == 10 & method == "SI-normal-local")
si10normnonlocal <- subset(x = plotDF, k == 10 & method == "SI-normal-nonlocal")


## ----"figure-SIcomparison-RECOVERY", fig.height = 3.25------------------------
## define diverging color palette
plotDF$kfac <- factor(plotDF$k, levels = k, labels = formatBF(BF = k))
cols <- divergingx_hcl(n = length(k), palette = "Zissou 1", rev = TRUE)
names(cols) <- levels(plotDF$kfac)

## plot CI and SIs
plotA <- ggplot() +
    geom_point(data = plotDF, aes(x = type, y = est), size = 0, alpha = 0) +
    geom_errorbar(data = ciDF, aes(x = type, ymin = lower, ymax = upper),
                  col = 1, alpha = 1, width = 0.35, size = 0.4)
## HACK draw SIs in correct order so that narrower SI on top of longer ones
for (i in seq_along(k)) {
    ki <- k[i]
    dat <- subset(plotDF, k == ki)
    plotA <- plotA +
        geom_errorbar(data = dat, key_glyph = "rect",  alpha = 1,
                      size = 0.4 + i*0.02,
                      aes(x = type, ymin = lower, ymax = upper, color = kfac),
                      width = 0.35)
}
plotA +
    geom_point(data = subset(plotDF, type == ciName),
               aes(x = type, y = est), size = 0.9, color = 1, alpha = 1) +
               #shape = "square") +
    labs(x = "", y = bquote("Log hazard ratio" ~ theta)) +
    scale_x_discrete(labels = function(l) parse(text = l)) +
    scale_color_manual(breaks = names(cols), values = cols, name = bquote("Support level" ~ italic(k))) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(hjust = 0, color = "black"))


## ----"sample-size-based-on-support"-------------------------------------------
## determine sample size such that k SI based on JAB has width w
k <- 10
w <- 0.2
lambda <- 2
expterm1 <- lamW::lambertW0(x = -k^2*w^2/4/lambda^2)
expterm2 <- lamW::lambertWm1(x = -k^2*w^2/4/lambda^2)
n <- k^2*exp(-c(expterm1, expterm2))

## ## check that both solutions correct
## estCheck <- 0
## seCheck <- lambda/sqrt(n)
## par(mfrow = c(1, 2))
## plot(ciCalibrate(siLevel = k, estimate = estCheck, se = seCheck[1],
##                  priorMean = estCheck, priorSD = lambda))
## plot(ciCalibrate(siLevel = k, estimate = estCheck, se = seCheck[2],
##                  priorMean = estCheck, priorSD = lambda))


## ----"code-snippet-ciCalibrate", echo = TRUE, fig.height = 4, size = codeFont----
## 95% CI from RECOVERY trial
logHRci <- c(-0.29, -0.07)
## compute a support interval with level k = 10
library("ciCalibrate") # install with install.packages("ciCalibrate")
si10 <- ciCalibrate(ci = logHRci, ciLevel = 0.95, siLevel = 10,
                    method = "SI-normal", priorMean = 0, priorSD = 2)
si10
## plot Bayes factor function with support interval
plot(si10)



## ----"sessionInfo2", echo = Reproducibility, results = Reproducibility, size = codeFont----
cat(paste(Sys.time(), Sys.timezone(), "\n"))
sessionInfo()

