# Definition for function used for generating the correlation matrices
generateMatrixPlot <- function(plotFormula, diagLabels, ageType) {
    png(
        filename = paste0(myEnv$baseDirectory, "/matrixplot", ageType, ".png"),
        width = 2500, height = 2500, res = 300, type = "cairo"
    )
    pairs(
        plotFormula,
        data = myEnv$pdataSVs,
        upper.panel = twolines,
        lower.panel = panel.cor,
        diag.panel = function(x, ...) mydiag.panel(x, labels = diagLabels, ...),
        labels = rep("", length(diagLabels)),
        gap = 0.3,
        main = ""
    )
    dev.off()
}

# Definition for function used for generating the grouped bar chart
createGroupedBarChart <- function(data, x, y, fill, title) {
    # All measures are every column except the x column (preserve incoming order)
    measureVars <- setdiff(colnames(data), x)

    # Long format with tidyverse; stable names via all_of()
    long <- tidyr::pivot_longer(
        data,
        cols = tidyselect::all_of(measureVars),
        names_to = fill,
        values_to = y
    )

    # Force fill order to follow original column order
    long[[fill]] <- factor(long[[fill]], levels = measureVars)

    # Palette sized to the number of measures
    nFills <- length(measureVars)
    colorPalette <- scales::hue_pal()(nFills)
    names(colorPalette) <- measureVars

    p <- ggplot(
        long,
        aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])
    ) +
        geom_bar(
            stat = "identity",
            width = 0.6,
            position = ggplot2::position_dodge(width = 0.6)
        ) +
        labs(x = x, y = "Age", title = title) +
        scale_fill_manual(values = colorPalette, drop = FALSE) +
        theme_minimal() +
        theme(
            plot.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 65, hjust = 1)
        )

    ggsave(
        filename = file.path(myEnv$baseDirectory, "GroupedAge.png"),
        plot = p, width = 35, height = 10, units = "in", dpi = 300
    )

    # Optional: scatterplots of clocks vs age (wide data is fine)
    if ("Age" %in% colnames(myEnv$pdataSVs) || "age" %in% colnames(data)) {
        clocksOnly <- setdiff(measureVars, c("age", x))
        for (nm in clocksOnly) {
            p2 <- ggpubr::ggscatter(
                data, x = "age", y = nm,
                add = "reg.line",
                add.params = list(color = "blue", fill = "lightgray")
            ) +
                ggpubr::stat_cor(method = "pearson")

            ggsave(
                filename = file.path(myEnv$baseDirectory, paste0("plot_", nm, ".png")),
                plot = p2, width = 1000, height = 1000, units = "px"
            )
        }
    }
}

# Definition for function used in matrix generation
panel.cor <- function(x, y, digits = 2, prefix = "",
                      cex.min = 0.05, cex.max = 2, shrink = 1, ...) {
    usr <- par("usr"); on.exit(par(usr = usr))
    par(usr = c(0, 1, 0, 1))

    r <- abs(stats::cor(x, y, use = "complete.obs")); if (is.na(r)) r <- 0
    txt <- paste0(prefix, formatC(r, format = "f", digits = digits))

    # Measure at cex = 1 so size math is stable
    w <- strwidth(txt, cex = 1); h <- strheight(txt, cex = 1)
    cex_fit <- 0.9 / max(w, h)          # fits inside unit box
    cex <- max(cex.min, min(cex_fit * shrink, cex.max))

    text(0.5, 0.5, txt, cex = cex)
}


# Definition for function used in matrix generation
mydiag.panel <- function(x, labels,
                         pad = 0.9,
                         cex.min = 0.2,
                         cex.max = 3,
                         box.col = "#CC7178",
                         txt.col = "black",
                         font = 2,
                         ...) {
    usr0 <- par("usr"); on.exit(par(usr = usr0))
    par(usr = c(0, 1, 0, 1))

    rect(0, 0, 1, 1, col = box.col, border = NA)

    idx <- par("mfg")[2]
    lab <- labels[idx]

    # measure at cex = 1, then scale to fit
    w <- strwidth(lab,  cex = 1, font = font)
    h <- strheight(lab, cex = 1, font = font)
    cex.fit <- pad / max(w, h)           # uniform scale to keep aspect
    cex <- max(cex.min, min(cex.fit, cex.max))

    text(0.5, 0.5, lab, cex = cex, col = txt.col, font = font)
}

# Definition for function used in matrix generation
twolines <- function(x, y) {
    points(x, y, pch = 20)
    abline(lm(y ~ x), col = "#CC7178")
}

# Making data frame for plotting
preparePlotDf <- function(results, pdataSVs) {
    clockCols <- pickPlotClocks(results)
    df <- data.frame(sample = colnames(myEnv$bVals),
                     age = pdataSVs$Age,
                     results[, clockCols, drop = FALSE],
                     check.names = FALSE)
    df
}
