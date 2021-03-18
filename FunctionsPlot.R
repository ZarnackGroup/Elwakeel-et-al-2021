plot.correlationsTak1l <- function(corrDf) {
    p1 = ggplot(corrDf, aes(x = corpval)) +
        geom_histogram(bins = 50,
                       color = "lightgrey",
                       fill = "darkgrey") +
        theme_bw() +
        xlab("P-value") +
        ylab("Count") +
        labs(tag = "A")

    p2 = ggplot(corrDf, aes(x = corcoef, fill = sigAll)) +
        geom_histogram(bins = 100) +
        scale_fill_manual(values = c("darkgrey", "#7E6148FF", "#B09C85FF")) +
        theme_bw() +
        xlab("Correlation coefficient") +
        ylab("Count") +
        labs(tag = "B")

    p3 = ggplot(corrDf, aes(
        x = log2(baseMean),
        y = corcoef,
        color = sigAll
    )) +
        geom_hline(yintercept = 0) +
        geom_point_rast(alpha = .3, size = 3) +
        scale_color_manual(values = c("darkgrey", "#7E6148FF", "#B09C85FF")) +
        theme_bw() +
        xlab("Mean gene expression [log2]") +
        ylab("Correlation coefficient") +
        labs(tag = "C")

    p4 = ggplot(corrDf, aes(
        x = (corcoef),
        y = -log10(padj),
        color = sigAll
    )) +
        geom_vline(xintercept = -0.3) +
        geom_vline(xintercept = 0.3) +
        geom_hline(yintercept = -log10(.05)) +
        geom_point_rast(
            alpha = 0.3,
            size = 3,
            position = position_jitter(width = 0.05, height = 0.1)
        ) +
        scale_color_manual(values = c("darkgrey", "#7E6148FF", "#B09C85FF")) +
        theme_bw() +
        xlab("Correlation coefficient") +
        ylab("Adjusted p-value [-log10]") +
        labs(tag = "D")
    ret = list(p1, p2, p3, p4)
    return(ret)
}

get_density <- function(x, y, n) {
    dens <- MASS::kde2d(x = x, y = y, n = n)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}

pairsBsScatterRast <-
    function(data,
             mapping,
             N = 100,
             diag = FALSE,
             reg = FALSE,
             ellipse = FALSE,
             ...) {
        get_density <- function(x, y, n) {
            dens <- MASS::kde2d(x = x, y = y, n = n)
            ix <- findInterval(x, dens$x)
            iy <- findInterval(y, dens$y)
            ii <- cbind(ix, iy)
            return(dens$z[ii])
        }

        X <- GGally::eval_data_col(data, mapping$x)
        Y <- GGally::eval_data_col(data, mapping$y)

        data$density <- get_density(x = X, y = Y, n = N)

        p <- ggplot(data, mapping) +
            geom_point_rast(aes(colour = density),
                            size = 4,
                            raster.dpi = 50) +
            scale_color_viridis(option = "cividis") +
            theme_bw() +
            theme(axis.text.x = element_text(
                angle = 90,
                hjust = 1,
                vjust = 0.5
            ))

        return(p)
    }

pairsDensity <- function(data, mapping, ...) {
    require(tidyr)
    # clean input
    data = data %>% drop_na()
    # compute quantiles to display
    x <- GGally::eval_data_col(data, mapping$x)
    x = as.numeric(x)
    dt = data.frame(x = c(1:length(x)), y = x)
    dens <- density(dt$y)
    df <- data.frame(x = dens$x, y = dens$y)
    probs <- c(0, 0.5, 1)
    quantiles <- quantile(dt$y, prob = probs)
    quant <- factor(findInterval(df$x, quantiles))
    df$quant <- factor(findInterval(df$x, quantiles))

    # make plot
    p = ggplot(df, aes(x, y)) +
        geom_line() + geom_ribbon(aes(ymin = 0, ymax = y), fill = "#B09C85FF") +
        scale_x_continuous(breaks = quantiles) +
        scale_fill_brewer(guide = "none") +
        theme_bw() +
        theme(axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
        ))
    return(p)
}


pairsCorColor <-
    function(data,
             mapping,
             color = I("black"),
             sizeRange = c(1, 5),
             ...) {
        # set P value cutoff
        pvalCo = 0.05

        # get the x and y data to use the other code
        x <- GGally::eval_data_col(data, mapping$x)
        y <- GGally::eval_data_col(data, mapping$y)

        # calc pairwise pearson correlation
        ct <- cor.test(x, y, method = "pearson")

        # extract P value
        pval = ct$p.value

        # exctract correlation coefficient
        r <- unname(ct$estimate)
        rt <- format(r, digits = 1)[1]
        tt <- as.character(rt)

        # plot the cor value
        p <- ggally_text(
            label = tt,
            mapping = aes(),
            xP = 0.5,
            yP = 0.5,
            size = 4,
            color = color,
            ...
        ) +
            theme_void() +
            theme(
                panel.background = element_rect(fill = "white"),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank()
            )

        corColors1 = sequential("#E64B35FF", percentage = 12, plot = F)
        corColors2 = sequential("#3C5488FF", percentage = 12, plot = F)

        if (r >= .8) {
            corCol = corColors2[9]
        } else if (r >= .7) {
            corCol = corColors2[8]
        } else if (r >= .6) {
            corCol = corColors2[7]
        } else if (r >= .5) {
            corCol = corColors2[6]
        } else if (r >= .4) {
            corCol = corColors2[5]
        } else if (r >= .3) {
            corCol = corColors2[4]
        } else if (r >= .2) {
            corCol = corColors2[3]
        } else if (r >= .1) {
            corCol = corColors2[2]
        }
        else if (r <= -.8) {
            corCol = corColors1[9]
        } else if (r <= -.7) {
            corCol = corColors1[8]
        } else if (r <= -.6) {
            corCol = corColors1[7]
        } else if (r <= -.5) {
            corCol = corColors1[6]
        } else if (r <= -.4) {
            corCol = corColors1[5]
        } else if (r <= -.3) {
            corCol = corColors1[4]
        } else if (r <= -.2) {
            corCol = corColors1[3]
        } else if (r <= -.1) {
            corCol = corColors1[2]
        }
        else {
            corCol = "white"
        }

        # set background color to grey for not significant correlations
        if (pval >= pvalCo) {
            corCol = "#A9A9A9"
        }

        p <- p + theme(panel.background = element_rect(fill = corCol,
                                                       color = corCol))
        p
    }
