setwd('~/Dropbox/msdbPaper')
#### Figure 1 ####
{
    rm(list=ls())
    h2 <- 0.5
    library('wesanderson')
    source('scripts/figureFuncs.R')
    my.cols <- c('forestgreen', 'firebrick4','dodgerblue4')
    my.cols <- c(wes_palette('Rushmore1')[c(3,4)],wes_palette('FantasticFox1')[5])[c(1,3,2)]
    png(
        'figures/paperFiguresForRealForRealThisTime/SchematicFigure1.png',
        height = 9 ,
        width = 22,
        units = 'cm',
        res = 500
    )
    nf <- layout(matrix(c(1, 2), nrow = 1, byrow = TRUE))
    
    ## 1A
    ## mid zoom liability distribution
    op1 <- par(mar = c(2, 2, 1, 1.3))
    {
        makePhenLi(
            xlim = c(-3, 4),
            prev = 0.01,
            shade = 'disease',
            w.thr = FALSE,
            cols = my.cols
        )
        makeGenLi(
            h2 = h2,
            xlim = c(-3, 4),
            plot.x.axis = FALSE,
            plot.y.axis = FALSE,
            add = TRUE,
            y.modifier = dnorm(0) / dnorm(0, sd = sqrt(h2)),
            my.color = my.cols[1]
        )
    }
    
    ## lines drawn on the plot
    {
        abline(v = 0 , lty = 2 , lwd = 2, col = my.cols[1])
        lines(
            x = rep(qnorm(0.99), 2),
            y = c(0, dnorm(0) * 0.81) ,
            lty = 1 ,
            lwd = 2 ,
            col = my.cols[2]
        )
        arrows(
            x0 = 0,
            x1 = -h2 * 2,
            y0 = dnorm(0) * 0.35 ,
            y1 = dnorm(0) * 0.35 ,
            length = 0.05
        )
        arrows(
            x0 = 0 ,
            x1 = -1.53,
            y0 = dnorm(0) * 0.3 ,
            y1 = dnorm(0) * 0.3 ,
            length = 0.05
        )
    }
    
    ## text drawn on the plot
    {
        text(
            x = -0.45 ,
            y = dnorm(0) * 0.4 ,
            labels = expression(sqrt(V[A])) ,
            cex = 0.86
        )
        text(
            x = -0.6 ,
            y = dnorm(0) * 0.25 ,
            labels = expression(sqrt(V[A] + V[E])) ,
            cex = 0.86
        )
        text(
            x = -1.4 ,
            y = dnorm(0) * 1.1 ,
            labels = expression(paste('mean liability, ', bar(G), sep = '')),
            cex = 1
        )
        text(
            x = 3.3 ,
            y = dnorm(0) * 0.78 ,
            labels = 'threshold, T',
            cex = 1 ,
            col = my.cols[2]
        )
        text(
            x = 3.3 ,
            y = dnorm(0) * 0.27 ,
            labels = 'individuals' ,
            cex = 0.7 ,
            col = my.cols[2]
        )
        text(
            x = 3.4 ,
            y = dnorm(0) * 0.22 ,
            labels = 'that develop' ,
            cex = 0.7 ,
            col = my.cols[2]
        )
        text(
            x = 3.37 ,
            y = dnorm(0) * 0.18 ,
            labels = 'the disease' ,
            cex = 0.7 ,
            col = my.cols[2]
        )
        arrows(
            x0 = 2.9 ,
            y0 = dnorm(0) * 0.15 ,
            x1 = 2.77 ,
            y = dnorm(0)  * 0.04 ,
            col = my.cols[2] ,
            length = 0.1
        )
    }
    
    ## axis labels and  legend
    {
        mtext(
            side = 1 ,
            text = 'Liability',
            line = 1,
            cex = 1.2
        )
        mtext(
            side = 2 ,
            text = 'Density',
            line = 1,
            cex = 1.2
        )
        # mtext(
        #     side = 3 ,
        #     text = 'A' ,
        #     line = -0.4,
        #     at = -3.5,
        #     cex = 1.5
        # )
        legend(
            'topright' ,
            lty = 1 ,
            lwd = 1.5 ,
            col = c(my.cols[3], my.cols[1]) ,
            legend = c('total liability' , 'genetic liability') ,
            bty = 'n'
        )
    }
    
    par(op1)
    
    {## 1B
        
        this.col = 'darkorchid4'#my.cols[2]
        y.max.factor = 2
        op4 <- par(mar = c(2, 2, 1, 1))
        this.ymax <- 1.128378
        makeEmptyFrame(
            xlim = c(-3, 6.2) ,
            y.max.factor = 2,
            plot.x.axis = FALSE
        )
        lines(
            x = rep(qnorm(1 - 0.01),2),
            y = c(dnorm(qnorm(0.99),0,sqrt(h2)),this.ymax*1.14/1.2),
            lty = 1 ,
            lwd = 2 ,
            col = my.cols[2],
            )
        axis(side = 1 , at = seq(-3,6,by=1) , labels = rep('',6-(-3)+1))
        this.ymax <- dnorm(0, sd = sqrt(h2)) * y.max.factor
        abline(
            h = (this.ymax / 1.2) / (1 - 0.01 * 0.5) ,
            lty = 2 ,
            col = this.col
        )
        lines(
            x = c(-4,-0.43) ,
            y = rep((0.5 * this.ymax / 1.2) / (1 - 0.01 * 0.5) , 2 ) ,
            lty = 2 ,
            col = this.col
        )
        lines(
            x = c(0.44,7) ,
            y = rep((0.5 * this.ymax / 1.2) / (1 - 0.01 * 0.5) , 2 ) ,
            lty = 2 ,
            col = this.col
        )
        text(
            x = 3.58 ,
            y = this.ymax*1.11/1.2 ,
            labels = 'threshold, T',
            cex = 1 ,
            col = my.cols[2]
        )
        text(
            x = -1 ,
            y = this.ymax*1.06/1.2 ,
            labels = expression(paste('1/', (1-bar(R)*C),sep = '') ),
            cex = 1 ,
            col = this.col
        )
        text(
            x = 5.17 ,
            y = this.ymax*0.56/1.2 ,
            labels = expression(paste('(1-C)/', (1-bar(R)*C),sep = '') ),
            cex = 1 ,
            col = this.col
        )
        {
            
            this.ymax <-
                makeGenLi(
                    xlim = c(-3, 4),
                    h2 = h2,
                    plot.y.axis = FALSE ,
                    plot.x.axis = FALSE ,
                    y.max.factor = 2,
                    add = TRUE ,
                    my.color=my.cols[1]
                )

            addRiskCurve(
                h2 = h2,
                prev = 0.01,
                xlim = c(-3, 6) ,
                y.max = (this.ymax / 1.2) / (1 - 0.01 * 0.5),
                my.col = this.col,
                fit.surface = TRUE
            )
            mtext(side = 1 ,
                  text = 'Genetic Liability',
                  line = 1,
                  cex = 1.2)
            mtext(
                side = 2 ,
                text = 'Expected fitness',
                line = 1,
                las = 0
            )
            axis(
                side = 2,
                at = seq(0, 2 * this.ymax / 1.2, length.out = 3),
                labels = seq(0, 2, by = 1),
                las = 1
            )
            # mtext(
            #     side = 3 ,
            #     text = 'B' ,
            #     line = -0.4,
            #     at = -3,
            #     cex = 1.5
            # )
        }
        par(op4)
    }
    dev.off()
}

#### Figure 2 ####
{
    rm(list=ls())
    png(
        'figures/paperFiguresForRealForRealThisTime/SchematicFigure2.png',
        height = 18 ,
        width = 22,
        units = 'cm',
        res = 500
    )
    {
        my.cols <- c(
            wes_palette('Rushmore1')[4] ,
            wes_palette('FantasticFox1')[c(4,5)]
        )
        
        op1 <- par(mar = c(1.3, 1.7, 1, 0) + 0.1)
        source('scripts/figureFuncs.R')
        ## y.max <- 0.6
        t.pos <- qnorm(0.99)
        alphaS <- 0.09
        my.xlim <- c(1.5, 2.6)
        v.offset <- 0.007
        my.out <- makePhenWEffect(
            xlim = my.xlim,
            alphaS = alphaS,
            t.pos = t.pos,
            v.offset = v.offset ,
            xaxt = 'n',
            return.stuff = TRUE ,
            cols = my.cols
        )
        y.max <- my.out$y.max
        effect.x <- my.out[[3]]
        effect.y <- my.out[[4]]
    }
    par(op1)
    dev.off()
}

#### Figure 3 ####
{
    rm(list=ls())
    library(wesanderson)
    png(
        'figures/paperFiguresForRealForRealThisTime/Figure3.png',
        height = 7 ,
        width = 17.75,
        units = 'cm',
        res = 600
    )
    layout(matrix(c(1, 2, 3),nrow=1))
    
    {
#### 3A ####
        {
            op3 <- par(mar = c(3.1, 2.6, 0, 0.1),
                       mgp = c(3, 0.6, 0))
            my.gamma <-
                10 ^ seq(log(0.01, 10), log(100, 10), length.out = 1000)
            bias <-
                ifelse(my.gamma < 700, (exp(my.gamma) - 1) / (exp(my.gamma) + 1), 1)
            plot(
                x = my.gamma ,
                y = bias ,
                log = 'x' ,
                type = 'n' ,
                bty = 'n' ,
                lwd = 2 ,
                xlab = '' ,
                ylab = '',
                xaxt = 'n',
                yaxt = 'n'
            )
            mtext(
                text = 'Fixation bias',
                side = 2,
                line = 1.4,
                cex = 0.7
            )
            mtext(
                text = expression(paste('Scaled selection coefficient, ', gamma, sep = '')),
                side = 1,
                line = 1.6,
                cex = 0.7,
                at = 0.833666
            )
            axis(
                side = 1 ,
                at = c(0.01, 0.1, 1, 10, 100),
                labels = c('0.01', '0.1', '1', '10', '100'),
                cex.axis = 0.7
            )
            axis(
                side = 2 ,
                at = 0:5 / 5,
                labels = c('0', '0.2', '0.4', '0.6', '0.8', '1'),
                cex.axis = 0.7,
                las = 2
            )
            
            
                                        # regime markings and labels
            {
                polygon(
                    x = c(0.01, 0.2, 0.2, 0.01) ,
                    y = c(1, 1, 0, 0) ,
                    col = adjustcolor('grey89') ,
                    border = NA
                )
                polygon(
                    x = c(5, 100, 100, 5) ,
                    y = c(1, 1, 0, 0) ,
                    col = adjustcolor('grey89') ,
                    border = NA
                )
                text(
                    labels = 'effectively',
                    x = 0.04,
                    y = 0.2 ,
                    cex = 1
                )
                text(
                    labels = 'neutral',
                    x = 0.04,
                    y = 0.15,
                    cex = 1
                )
                text(
                    labels = 'weakly',
                    x = 0.7,
                    y = 0.9 ,
                    cex = 1
                )
                text(
                    labels = 'selected',
                    x = 0.7,
                    y = 0.85,
                    cex = 1
                )
                text(
                    labels = 'strongly',
                    x = 20,
                    y = 0.9 ,
                    cex = 1
                )
                text(
                    labels = 'selected',
                    x = 20,
                    y = 0.85,
                    cex = 1
                )
                
            }
            
            lines(x = my.gamma,
                  y = bias,
                  lwd = 2)
            
            #mtext(text = 'A',side = 2, line = 0.4 , las = 2 , at = 0.9)
            
            par(op3)
        }
        

#### 3B ####
        {
            freqSpec <- function(x, y) {
                2 * exp(-2 * y * x) / ((1 + exp(-2 * y)) * x * (1 - x))
            }
            my.x <- seq(1/40000,1-1/40000,length.out=1000000)
            my.y <- c(0.1, 1, 10,100)
            tmp.freq.specs <- lapply(my.y, function(Y)
                freqSpec(my.x, Y))
            prob.dens <- TRUE
            prob.freq.specs <- lapply(tmp.freq.specs , function (X ) X / sum(X))
            my.cdfs <- lapply(prob.freq.specs , cumsum )
            plot.these <- lapply ( my.cdfs , function(X) X < 0.999999)
            if(prob.dens){
                ymax <- freqSpec(0.01,1)/(7*sum(tmp.freq.specs[[2]]))
                
                my.freq.specs <- prob.freq.specs
            } else {
                ymax <- freqSpec(0.01,1)
                my.freq.specs <- tmp.freq.specs
            }
            
            
            my.cols <- wes_palette("AsteroidCity3")
            op4 <- par(mar = c(3.1, 2.1, 0, 0.4),
                       mgp = c(3, 0.5, 0))
            plot(
                NA,
                type = 'n',
                xlim = c(0,1),
                ylim = c(0,ymax) ,
                xaxt = 'n' ,
                yaxt = 'n' ,
                bty = 'n' #,log = 'y'
            )
            for (i in 1:length(my.freq.specs)) {
                lines(my.x[plot.these[[i]]] ,
                      my.freq.specs[[i]][plot.these[[i]]],
                      col = my.cols[i],
                      lwd = 1.7)
            }
            axis(1, cex.axis = 0.75)
            axis(
                2,
                at = c(0,ymax) ,
                labels = FALSE ,
                cex.axis = 0.75,
                las = 2
            )
            mtext(
                text = 'Frequency of liability increasing allele',
                side = 1,
                line = 1.5,
                cex = 0.7
            )
            mtext(
                text = 'Density',
                side = 2,
                line = 0.5,
                cex = 0.7
            )
            legend(
                x = 0.2135 ,
                y = 1.1e-6 ,
                lwd = 2.7 ,
                lty = 1,
                col = my.cols ,
                legend = c('0.1', '1', '10', '100'),
                title = expression(paste('Scaled coefficient (', gamma, ')' , sep = '')) ,
                bty = 'n'
            )
            # mtext(text = 'B',side = 2, line = 0.4 , las = 2, at = 0.9*ymax)
            par(op4)
        }
        
#### 3C ####
        {
            op5 <- par(mar = c(3.1, 2.5, 0, 0.4),
                       mgp = c(3, 0.8, 0))
            my.gamma <- exp(seq(log(0.01), log(100), length.out = 10000))
            bias <-
                ifelse(my.gamma < 700, (exp(my.gamma) - 1) / (exp(my.gamma) +
                                                              1), 1)
            het <- 2 * bias / my.gamma
            col1 <- wes_palette('FantasticFox1')[1]
            col2 <- wes_palette('FantasticFox1')[3]
            
            
            plot(
                my.gamma ,
                het ,
                log = 'xy' ,
                type = 'n' ,
                ylim = c(1 / 20, 1.4) ,
                xlim = c(0.01, 100) ,
                xaxt = 'n' ,
                yaxt = 'n' ,
                bty = 'n' ,
                xlab = '',
                ylab = ''
            )
            polygon(
                x = c(0.01, 0.2, 0.2, 0.01) ,
                y = c(1.4, 1.4, 1/20, 1/20) ,
                col = adjustcolor('grey89') ,
                border = NA
            )
            polygon(
                x = c(5, 100, 100, 5) ,
                y = c(1.4, 1.4, 1/20, 1/20) ,
                col = adjustcolor('grey89') ,
                border = NA
            )
            lines(
                x = c(0.01,100) ,
                y = c(1,1) ,
                lty = 3 ,
                lwd = 2.3 ,
                col = my.cols[1]
            )
            plot.these <- (2 / my.gamma) < 1.4
            lines(
                x = my.gamma[plot.these],
                y = 2 / my.gamma[plot.these],
                lty = 3 ,
                lwd = 2.3,
                col = my.cols[4]
            )
            lines(my.gamma ,
                  het ,
                  lwd = 3)
            text(
                x = 0.045,
                y = 1.14 ,
                cex = 1.2 ,
                labels = expression("E[h(x)]"%~~% theta ) ,
                col = my.cols[1]
            )
            text(
                x = 16.5,
                y = 0.195 ,
                cex = 1.2 ,
                srt = 360 - 74 ,
                labels = expression("E[h(x)]"%~~% theta/gamma == 2*u/s) ,
                col = my.cols[4]
            )
            axis(
                side = 1 ,
                at = c(0.01, 0.1, 1, 10, 100) ,
                labels = c('0.01', '0.1', '1', '10', '100') ,
                cex.axis = 0.7
            )
            axis(
                side = 2 ,
                at = c(1 / 20, 1.4) ,
                labels = FALSE ,
                las = 2
            )
            mtext(
                side = 2,
                text = '(log) Heterozygosity',
                line = 0.5,
                cex = 0.7
            )
            mtext(
                text = expression(paste('Scaled selection coefficient, ', gamma, sep = '')),
                side = 1,
                line = 1.6,
                cex = 0.7 ,
                at = 0.833666
            )
            # mtext(text = 'C',side = 2, line = 0.4 , las = 2, at = 0.9*1.2)
            par(op5)
        }
    }
    dev.off()
}



#### Figure 4 ####
{
    rm(list=ls())
    library(wesanderson)
    sim.results <- get(load('figures/smallEffectInsensitivityResultsTable.Rdata'))
    png(
        'figures/paperFiguresForRealForRealThisTime/Figure4.png',
        height = 7 ,
        width = 17.75,
        units = 'cm',
        res = 600
    )
    layout(matrix(c(1, 2, 3),nrow=1),widths = c(0.8,1,1))
    
#### 4A ####
    {
        sim.results$costFactor <- as.factor(sim.results$cost)
        my.cols <- wes_palette("Darjeeling1")[c(2, 4)]
        my.bt <- seq(1e-6, 1 - 1e-6 , length.out = 1e5)
        my.gamma <- 0.5 * log((1 + my.bt) / (1 - my.bt))
        these.costs <- unique(sim.results$cost)
        op3 <- par(mar = c(2.9, 2.6, 0, 0.1),
                   mgp = c(3, 0.6, 0))
        plot(
            NA,
            xlim = c(0, 1) ,
            ylim = c(0, 6.5) ,
            bty = 'n' ,
            ylab = '' ,
            xlab = ''
        )
        lines(x = my.bt ,
              y = my.gamma)
        mtext(
            side = 1 ,
            text = expression(paste('Relative threshold position, ' , b[T] , sep = '')) ,
            line = 1.5,
            cex = 0.7,
            )
        mtext(
            side = 2 ,
            text = expression(paste('Scaled selection coefficient, ', gamma, (a), sep = '')) ,
            line = 1.5,
            cex = 0.7,
            )
        points(
            x = sim.results$b ,
            y = 2 * sim.results$Ne * sim.results$cost * sim.results$sim.deltaR ,
            xlim = c(0, 1) ,
            pch = c(4, 20)[sim.results$costFactor] ,
            col = my.cols[sim.results$costFactor] ,
            cex = 2
        )
        legend(
            x = 0.5,
            y = 6,
            pch = c(4, 20) ,
            col = my.cols ,
            bty = 'n' ,
            title = 'Fitness cost' ,
            legend = these.costs ,
            pt.cex = 2
        )
        ##text(x = 0.02, y = 6.4,labels='A) The single effect model',cex = 2)
        par(op3)
    }
    
#### 4B ####
    {
### plotting params ###
        {
            op4 <- par(mar = c(2.9, 0, 0, 0),
                       mgp = c(3, 0.6, 0))
            x.min <- -4
            x.max <- 8
            y.min <- 0
            y.max <- 1.4 * dnorm(0)
            lia.max <- 4
            my.x <- seq(x.min, lia.max , length.out = 1e5)
            y1 <- dnorm(my.x)
            thr <- qnorm(0.99)
            y1ft <- dnorm(thr)
            mean.shift <-
                uniroot(function(X)
                    y1ft * these.costs[1] / these.costs[2] - dnorm(thr + X) ,
                    interval = c(0, 1))$root
            y2 <- dnorm(my.x , mean = -mean.shift)
            aspect.ratio <- (y.max - y.min) / (x.max - x.min)
            box.hw <- 0.4
        }
        
### main plot ###
        {
            plot(
                NA,
                xlim = c(x.min, x.max) ,
                ylim = c(y.min, y.max) ,
                bty = 'n' ,
                xaxt = 'n' ,
                yaxt = 'n'
            )
            # text(x = x.min + (x.max - x.min)*0.05, y = y.min + (y.max - y.min)*0.99,labels='B',cex = 2)
            axis(side = 1,
                 at = seq(x.min, lia.max),
                 labels = FALSE)
            mtext(
                side = 1 ,
                text = 'Liability' ,
                line = 1.2,
                cex = 0.7,
                at = 0
            )
            lines(x = my.x [y1 > y2] ,
                  y = y1[y1 > y2] ,
                  col = my.cols[1])
            lines(x = my.x [y1 <= y2] ,
                  y = y1[y1 <= y2] ,
                  col = adjustcolor(my.cols[1] , alpha.f = 0.3))
            polygon(
                x = c(my.x [y1 > y2], x.max , rev(my.x [y1 > y2])) ,
                y = c(y1[y1 > y2], 0 , rev(y2[y1 > y2])) ,
                col = adjustcolor(my.cols[1], alpha.f = 0.08),
                border = NA
            )
            lines(
                x = c(0, 0),
                y = c(0, dnorm(0,-mean.shift)),
                col = adjustcolor(my.cols[1], alpha.f = 0.3),
                lty = 2
            )
            lines(
                x = c(0, 0),
                y = c(dnorm(0,-mean.shift), dnorm(0) * 1.18),
                col = my.cols[1],
                lty = 2
            )
            lines(x = my.x ,
                  y = y2 ,
                  col = my.cols[2])
            polygon(
                x = c(x.min, my.x, x.max) ,
                y = c(0, y2, 0) ,
                col = adjustcolor(my.cols[2], alpha.f = 0.08),
                border = NA
            )
            lines(
                x = c(-mean.shift,-mean.shift),
                y = c(0, dnorm(0) * 1.2),
                col = my.cols[2],
                lty = 2
            )
            lines(x = rep(thr, 2),
                  y = c(0, dnorm(0) * 0.8))
        }
        
### small box ###
        {
            lines(x = thr + c(-box.hw, box.hw) ,
                  y = rep(2 * box.hw, 2) * aspect.ratio)
            lines(x = thr + rep(-box.hw, 2) ,
                  y = c(0, 2 * box.hw) * aspect.ratio)
            lines(x = thr + rep(box.hw, 2) ,
                  y = c(0, 2 * box.hw) * aspect.ratio)
            lines(x = thr + c(-box.hw, box.hw) ,
                  y = rep(0, 2))
        }
        
### parameters for bigger box ###
        {
            vo <- 0.11
            ho <- 3.4
            mult <- 6.6
            inset.ymin <- vo
            inset.ymax <- 2 * box.hw * aspect.ratio * mult + vo
            inset.xmin <- ho + thr - box.hw * mult
            inset.xmax <- ho + thr + box.hw * mult
        }
        
### zoom in lines ###
        {
            lines(
                x = c(thr - box.hw, inset.xmin) ,
                y = c(2 * box.hw * aspect.ratio , inset.ymax) ,
                lty = 2 ,
                lwd = 0.8 ,
                col = adjustcolor('black' , alpha.f = 0.4)
            )
            lines(
                x = c(thr + box.hw, inset.xmax) ,
                y = c(0 , inset.ymin) ,
                lty = 2 ,
                lwd = 0.8 ,
                col = adjustcolor('black' , alpha.f = 0.4)
            )
        }
        
### plot inside bigger box ###
        {
            n.pt <- 100000
            plot.x <- seq(inset.xmin, inset.xmax, length.out = n.pt)
            silent.x1 <- seq(thr - box.hw, thr + box.hw, length.out = n.pt)
            silent.x2 <-
                seq(thr - box.hw + mean.shift,
                    thr + box.hw + mean.shift,
                    length.out = n.pt)
            new.y1 <- inset.ymin + dnorm(silent.x1) * mult
            new.y2 <- inset.ymin + dnorm(silent.x2) * mult
            x.in.max <- max(plot.x[new.y1 < inset.ymax])
            
            
            ## shading inside box
            polygon(
                x = c(plot.x[new.y1 < inset.ymax] ,
                      rep(inset.xmax, 2) ,
                      rev(plot.x[new.y2 < inset.ymax]),
                      inset.xmin) ,
                y = c(new.y1[new.y1 < inset.ymax] ,
                      c(tail(new.y1, 1), tail(new.y2, 1)),
                      rev(new.y2[new.y2 < inset.ymax]),
                      inset.ymax) ,
                col = adjustcolor(my.cols[1], alpha.f = 0.08) ,
                border = NA
            )
            
            polygon(
                x = c(plot.x ,
                      rep(inset.xmax, 2) ,
                      rev(plot.x)) ,
                y = c(new.y2 ,
                      c(tail(new.y2, 1), inset.ymin),
                      rep(inset.ymin, length(plot.x))) ,
                col = adjustcolor(my.cols[2], alpha.f = 0.08) ,
                border = NA
            )
            
            ## density lines inside box
            lines(x = plot.x[new.y1 < inset.ymax] ,
                  y = new.y1[new.y1 < inset.ymax] ,
                  col = my.cols[1])
            lines(x = plot.x[new.y2 < inset.ymax] ,
                  y = new.y2[new.y2 < inset.ymax] ,
                  col = my.cols[2])
            
            ## threshold inside box
            lines(
                x = rep(ho + thr, 2),
                y = c(inset.ymin, inset.ymax),
                lwd = 1
            )
        }
        
### draw bigger box ###
        {
            ## top
            lines(x = thr + c(-box.hw, box.hw) * mult + ho ,
                  y = rep(inset.ymax, 2))
            ## left
            lines(x = thr + rep(-box.hw, 2) * mult + ho ,
                  y = c(vo, inset.ymax))
            ## right
            lines(x = thr + rep(box.hw, 2) * mult + ho ,
                  y = c(vo, inset.ymax))
            ## bottom
            lines(x = thr + c(-box.hw, box.hw) * mult + ho ,
                  y = rep(vo, 2))
        }
        
### annotations inside bigger box ###
        {
            ## text inside bigger box
            {
                text(
                    x = ho + thr - 1.4 ,
                    y = inset.ymin + (inset.ymax - inset.ymin) * 0.135 ,
                    labels = 'f(T|C=0.75)' ,
                    col = my.cols[2] ,
                    cex = 0.76
                )
                text(
                    x = ho + thr - 1.25 ,
                    y = inset.ymin + (inset.ymax - inset.ymin) * 0.67 ,
                    labels = 'f(T|C=0.25)' ,
                    col = my.cols[1] ,
                    cex = 0.76
                )
            }
            
            ## lines marking threshold density heights
            {
                                        # high cost
                arrows(
                    x0 = ho + thr - 0.2 ,
                    x1 = ho + thr - 0.2 ,
                    y0 = inset.ymin + 0.002 ,
                    y1 = inset.ymin + dnorm(thr,-mean.shift) * mult ,
                    col = my.cols[2] ,
                    angle = 90 ,
                    length = 0.022
                )
                arrows(
                    x0 = ho + thr - 0.2 ,
                    x1 = ho + thr - 0.2 ,
                    y1 = inset.ymin + 0.002 ,
                    y0 = inset.ymin + dnorm(thr,-mean.shift) * mult ,
                    col = my.cols[2] ,
                    angle = 90 ,
                    length = 0.022
                )
                
                                        # low cost
                arrows(
                    x0 = ho + thr + 0.2 ,
                    x1 = ho + thr + 0.2 ,
                    y0 = inset.ymin + 0.002 ,
                    y1 = inset.ymin + dnorm(thr) * mult ,
                    col = my.cols[1] ,
                    angle = 90 ,
                    length = 0.022
                )
                arrows(
                    x0 = ho + thr + 0.2 ,
                    x1 = ho + thr + 0.2 ,
                    y1 = inset.ymin + 0.002 ,
                    y0 = inset.ymin + dnorm(thr,-mean.shift) * mult ,
                    col = my.cols[1] ,
                    angle = 90 ,
                    length = 0.022
                )
            }
        }
        
### text labels in main plot ###
        {
            text.label.ho <- 0.6
            text(
                x = ho + thr - 1.7 + text.label.ho ,
                y = inset.ymax + 0.1 ,
                labels = 'f(T|C=0.75)' ,
                col = my.cols[2] ,
                cex = 0.9
            )
            text(
                x = ho + thr - 1.7 + text.label.ho,
                y = inset.ymax + 0.07 ,
                labels = 'f(T|C=0.25)' ,
                col = my.cols[1] ,
                cex = 0.9
            )
            lines(
                x = ho + thr + c(-3.15,-0.3) + text.label.ho ,
                y = rep(inset.ymax + 0.086, 2)
            )
            
            text(
                x = ho + thr + text.label.ho,
                y = inset.ymax + 0.086 ,
                labels = '=' ,
                cex = 1
            )
            text(
                x = ho + thr + 0.9 + text.label.ho,
                y = inset.ymax + 0.1 ,
                labels = '0.25' ,
                col = my.cols[1] ,
                cex = 0.9
            )
            text(
                x = ho + thr + 0.9 + text.label.ho,
                y = inset.ymax + 0.07 ,
                labels = '0.75' ,
                col = my.cols[2] ,
                cex = 0.9
            )
            lines(
                x = ho + thr + 0.3 + c(0, 1.15) + text.label.ho,
                y = rep(inset.ymax + 0.086, 2)
            )
        }
        par(op4)
    }
    
    
#### 4C ####
    solve.model.for.4C <- T
    op4 <- par(mar = c(2.9, 3.1, 0.1, 0.5),
               mgp = c(3, 0.6, 0))
    if (solve.model.for.4C) {
        source('scripts/solveTwoEffect.R')
        getSolnList <-
            function(Y, Z)
                sapply(Y, function(W)
                    sapply(W, function(X)
                        X[Z]))
        ba <- function(Y)
            ifelse(Y>50,1,(exp(2 * Y) - 1) / (exp(2 * Y) + 1))
        ya <- function(B) 0.5 * log((1 + B) / (1 - B))
        #yt <- c(1/2,1,3)
        bt <- c(0.5,0.9) #ba(yt)
        Ne = 20000
        L <- 1.5e7
        u <- 1e-8
        thetaL <- 4 * Ne * L * u
        C = c(0.25,0.75)
        yt <- ya(bt)
        
        ft <- sapply(yt, function(Y) Y / (2 * Ne * C) )
        min.ya <- 4#ya(1 - 100/(L))
        as <- 1
        soln <- list()
        prev <- list()
        lambda <- list()
        deltal <- list()
        Vt <- list()
        al <- list()
        bs <- list()
        raw.ft <- list()
        h2l <- list()
        n.pts <- 500
        #pl <- 0.000001
        pl <- 0.0005
        pl*L
        ## max.nl <- 3
        
        rf <- F
        ## fixed h2
        for(k in seq_along(C)){
          al[[k]] <- list()
          soln[[k]] <- list()
            for (j in seq_along(bt)) {
              se.vg <- thetaL * as ^ 2 * bt[j] / yt[j]
              min.al <- min.ya / (2 * Ne * ft[k,j] * C[k])
              al[[k]][[j]] <-
                  exp(seq(log(min.al), log(1100), length.out =
                                                      500))
              soln[[k]][[j]] <- list()
              last.bs <- NULL
              for (i in 1:length(al[[k]][[j]])) {
                  if (i > 1){
                      last.bs <- soln[[k]][[j]][[i - 1]]['bs']
                  }
                  soln[[k]][[j]][[i]] <- solveTwoEffect3D(
                      bt = bt[j],
                      init.bs = ifelse(i == 1, bt[j], last.bs),
                      Ne = 20000,
                      as = as,
                      al = al[[k]][[j]][i],
                      L = L,
                      gs = 1 - pl,
                      h2 = 1/2 ,
                      Ve = NULL ,
                      u = u,
                      C = C[k],
                      Bval = 1,
                      init.deltal = NULL,
                      init.tstar = NULL,
                      norm.deltal = NULL,
                      norm.init.tstar = NULL
                  )
              }
          }
        }
        al.mat <- lapply(al, function(AL) do.call(cbind, AL))
        prev.mat <- lapply(soln,function(SOLN) getSolnList(SOLN, 'prev'))
        deltal.mat <- lapply(soln,function(SOLN) getSolnList(SOLN,  'deltal'))
        ### change "lambda" to "mean.nl"
        lambda.mat <- lapply(soln,function(SOLN) getSolnList(SOLN,  'mean.nl'))
        Vt.mat <- lapply(soln,function(SOLN) getSolnList(SOLN,  'raw.Vt'))
        raw.ft.mat <- lapply(soln,function(SOLN) getSolnList(SOLN,  'raw.ft'))
        h2l.mat <- lapply(soln,function(SOLN) getSolnList(SOLN,  'h2l'))
        bs.mat <- lapply(soln,function(SOLN) getSolnList(SOLN,  'bs'))
        yl.mat <- lapply(soln,function(SOLN) getSolnList(SOLN,  'yl'))
        raw.ft.mat <- lapply(soln,function(SOLN) getSolnList(SOLN,  'raw.ft'))
        code.mat <- lapply(soln,function(SOLN) getSolnList(SOLN,  'code'))
        se.approx.deltal.mat <- mapply(AL=al.mat ,FT= raw.ft.mat,function(AL,FT) AL * FT )
        scaled.al.mat <- list()
        for( k in seq_along(C)){
          scaled.al.mat[[k]] <- 2*Ne*al.mat[[k]]*raw.ft.mat[[k]]*C[[k]]
        }
        
        ## two small effects solve
        ##small.al <- seq(1, 1000, length.out = 3000)
        ##small.deltal <- small.al*ft
        
    }
    {

        #### gamma plot ####
        my.cols <- wes_palette("Darjeeling1")[c(2, 4)] ##my.cols <- wes_palette('Darjeeling2')[c(1,2,5)]
        small.xlims <- c(0,40)
        big.xlims <- c(0,800)
        small.ylims <- c(0, 0.01)
        big.ylims <- c(0, 33500)
        {
          plot(
            NA,
            xlim = big.xlims,
            ylim = big.ylims,
            bty = 'n',
            lty = 1 ,
            col = my.cols ,
            bty = 'n',
            xlab = '',
            ylab = '',
            xaxt = 'n',
            yaxt = 'n'
          )
          abline(h = 2*Ne*C[1], lty = 3 , col = my.cols[1])
          abline(h = 2*Ne*C[2], lty = 3 , col = my.cols[2])
          text(x = 770 , y = 2*Ne*C[1]+950 , col = my.cols[1] , labels = '2NC' , cex = 1.2)
          text(x = 770 , y = 2*Ne*C[2]+950 , col = my.cols[2] , labels = '2NC' , cex = 1.2)
          #text(x = big.xlims[1] + (big.xlims[2] - big.xlims[1])*0.03, y = big.ylims[1] + (big.ylims[2] - big.ylims[1])*1,labels='C',cex = 2)
          axis(side = 1, at = c(0,400,800), labels = c('0','400','800') )
          axis(side = 2, at = c(0,10000,20000,30000), labels = c('0','10000','20000','30000') ,cex.axis = 0.64)
          mtext(side = 1,
                text = expression(paste('Liability effect of large effect sites, ', a[L], sep = '')),
                line = 1.6,
                cex = 0.7
          )
          mtext(side = 2,
                text = expression(paste('Scaled selection coefficient, ', gamma, (a[L]), sep = '')),
                line = 1.4,
                cex = 0.7)
          for(k in 1:length(C)){
            for (i in 1:ncol(al.mat[[k]])) {
              plot.these <- code.mat[[k]][, i] == 1 & al.mat[[k]][,i] < 800 & lambda.mat[[k]][,i] < 3
              lines(al.mat[[k]][plot.these, i],
                    2*Ne*C[k]*deltal.mat[[k]][plot.these, i],
                    col = my.cols[k],
                    lwd = 1,
                    lty = i)
            }
          }
        }
        

        #### 4C inset ####
        {
          #### draw little box ####
          aspect.ratio <- big.ylims[2]/big.xlims[2]
          box.width <- 15
          box.height <- aspect.ratio*box.width
          # lines(x=c(0,box.width),y=rep(0,2),lty = 1 )
          # lines(x=rep(0,2),y=c(0,box.height),lty = 1 )
          # lines(x=c(0,box.width),y=rep(box.height,2),lty = 1 )
          # lines(x=rep(box.width,2),y=c(0,box.height),lty = 1 )
          
          
          #### draw big box ####
          ho <- 41
          vo <- 17500
          mult <- 20
          big.box.width <- box.width*mult
          big.box.height <- box.width*aspect.ratio*0.9*mult
          big.box.xmax <- ho + big.box.width
          big.box.ymax <- vo + big.box.height
          lines(x=ho+c(0,big.box.width),y=rep(vo,2),lty = 1 )
          lines(x=rep(ho,2),y=vo+c(0,big.box.height),lty = 1 )
          lines(x=ho+c(0,big.box.width),y=vo+rep(big.box.height,2),lty = 1 )
          lines(x=ho+rep(big.box.width,2),y=vo+c(0,big.box.height),lty = 1 )
          
          ####plot lines inside inset ####
          plot.max.small.gamma <- 30
          x.rescale <- big.box.width / plot.max.small.gamma
          y.rescale <- big.box.height / plot.max.small.gamma
          for (k in seq_along(C)){
            for (i in 1:ncol(al.mat[[k]])) {
                plot.these1 <- code.mat[[k]][, i] == 1 & scaled.al.mat[[k]][,i] < plot.max.small.gamma
                tmp.x <- ho + x.rescale*scaled.al.mat[[k]][plot.these1, i]
                tmp.y <- vo + y.rescale*yl.mat[[k]][plot.these1, i]
                plot.these2 <- tmp.x < big.box.xmax & tmp.y < big.box.ymax
                plot.x <- tmp.x[plot.these2]
                plot.y <- tmp.y[plot.these2]
                lines(x = plot.x,
                      y = plot.y,
                      col = my.cols[k],
                      lty = i)
            }
          }
          
          
          #### draw lines connecting boxes ####
          lines(
              x = c(0, ho/3) ,
              y = c(0, vo*0.93),
              lty = 3,
              col = adjustcolor('black' , alpha.f = 0.3)
          )
          lines(
              x = c(box.width, ho+big.box.width) ,
              y = c(0, vo*0.92),
              lty = 3,
              col = adjustcolor('black' , alpha.f = 0.3)
          )
          
          ## x = y line
          lines(
            x = ho + c(0,big.box.width) ,
            y = vo + c(0,big.box.height) ,
            lty = 3,
            col = 'black'
          )
          
          #### ticks ####
          y.tick.label.offset <- 0.3
          my.ticks <- c(0,10,20,30)
          x.ticks <- ho + x.rescale*my.ticks
          for(i in 2:length(my.ticks)){
            lines(
              x = rep(x.ticks[i],2) ,
              y = vo*c(0.98,1)
            )
          }
          text(
            x = 2*ho*y.tick.label.offset ,
            y = vo*0.965 ,
            cex = 0.5 ,
            labels = '0'
          )
          text(
            x = ho + x.rescale*10,
            y = vo*0.95 ,
            cex = 0.5 ,
            labels = '10'
          )
          text(
            x = ho + x.rescale*20,
            y = vo*0.95 ,
            cex = 0.5 ,
            labels = '20'
          )
          text(
            x = ho + x.rescale*30,
            y = vo*0.95 ,
            cex = 0.5 ,
            labels = '30'
          )
          
          y.ticks <- vo + y.rescale*my.ticks
          for(i in 2:length(my.ticks)){
            lines(
              x = ho*c(0.75,1) ,
              y = rep(y.ticks[i],2)
            )
          }
          lines(
            x = c(ho*(1-0.4*0.4),ho) ,
            y = vo*c((1-0.02*0.75),1)
          )
          text(
            x = ho*y.tick.label.offset ,
            y = y.ticks[2] ,
            cex = 0.5 ,
            labels = '10'
          )
          text(
            x = ho*y.tick.label.offset ,
            y = y.ticks[3] ,
            cex = 0.5 ,
            labels = '20'
          )
          text(
            x = ho*y.tick.label.offset ,
            y = y.ticks[4] ,
            cex = 0.5 ,
            labels = '30'
          )
          
          #### inset axis labels ####
          {
            text(
              x= ho + (ho+big.box.width)/2 ,
              y = vo*0.90 ,
              labels=expression(paste("2", N, f(T), C, "\u00B7", a[L])),#expression(2 * N * f(T) * C %*% a[L]),
              cex = 0.5
            )
            text(
              x= -ho*y.tick.label.offset*3/4 ,
              y = vo + big.box.height/2 ,
              labels=expression(gamma(a[L])),
              cex = 0.5 ,
              srt = 90
            )
          }
          
        }
        
        #### legend ####
        {
          legend('bottomright',
                 legend = these.costs ,
                 title = 'C' ,
                 col = my.cols ,
                 lty = 1,
                 bty = 'n',
                 cex = 0.8
                 )
          legend(x = 375,
                 y = 5135,
                 legend = bt ,
                 title = expression(b[T]) ,
                 lty = c(1,2),
                 bty = 'n',
                 cex = 0.8
          )
        }
        
        #### draw main lines in big plot ####
        for(k in 1:length(C)){
          for (i in 1:ncol(al.mat[[k]])) {
            plot.these <- code.mat[[k]][, i] == 1 & al.mat[[k]][,i] < 800 & lambda.mat[[k]][,i] < 3
            lines(al.mat[[k]][plot.these, i],
                  2*Ne*C[k]*deltal.mat[[k]][plot.these, i],
                  col = my.cols[k],
                  lwd = 1,
                  lty = i)
          }
        }
    }
    dev.off()
  
## exploratory plot ##
    exploratory.plots <- F
    if(exploratory.plots) {
        my.cols <- wes_palette("FantasticFox1")
        std.al.mat <- al.mat / sqrt(Vt.mat)
        op4 <- par(mar = c(3.1, 3.1, 0, 0),
                   mgp = c(3, 0.6, 0))
        layout(t(matrix(1:4,nrow=2)))
        plot(
            NA ,
            xlim = c(0, 4) ,
            ylim = c(0, 1) ,
            type = 'l',
            lty = 1 ,
            col = my.cols ,
            bty = 'n' ,
            xlab = '',
            ylab = ''
        )
        mtext(side=1,text=expression(a[L]/sqrt(V[T])), line = 2)
        mtext(side=2,text=expression(delta[L]), line = 2)
        for(i in 1:ncol(std.al.mat)) {
            plot.these <- code.mat[,i]==1 & lambda.mat[,i]<max.nl
            lines(
                std.al.mat[plot.these,i],
                deltal.mat[plot.these==1,i], 
                col = my.cols[i]
            )
        }
        for(i in 1:ncol(std.al.mat)) {
            plot.these <- code.mat[,i]==1 & lambda.mat[,i]<max.nl
            lines(
                std.al.mat[plot.these,i],
                se.approx.deltal.mat[plot.these,i], 
                col = my.cols[i] ,
                lty = 2
            )
        }
        
        plot(
            NA ,
            xlim = c(0, 1500) ,
            ylim = c(0, 1) ,
            type = 'l',
            lty = 1 ,
            col = my.cols ,
            bty = 'n',
            xlab = '',
            ylab = ''
        )
        mtext(side=1,text=expression(a[L]), line = 2)
        mtext(side=2,text=expression(delta[L]), line = 2)
        for(i in 1:ncol(al.mat)) {
            plot.these <- code.mat[, i] == 1 & lambda.mat[, i] < max.nl
            lines(
                al.mat[code.mat[,i]==1,i],
                deltal.mat[code.mat[,i]==1,i], 
                col = my.cols[i]
            )
        }
        
        for(i in 1:ncol(al.mat)) {
            plot.these <- code.mat[,i]==1 & lambda.mat[,i]<max.nl
            lines(
                al.mat[code.mat[,i]==1,i],
                se.approx.deltal.mat[code.mat[,i]==1,i], 
                col = my.cols[i] ,
                lty = 2
            )
        }
        
        
        plot(
            NA ,
            xlim = c(0, 1000) ,
            ylim = c(0, 0.2) ,
            type = 'l',
            lty = 1 ,
            col = my.cols ,
            bty = 'n',
            xlab = '',
            ylab = ''
        )
        mtext(side=1,text=expression(a[L]), line = 2)
        mtext(side=2,text='Prevalence', line = 2)
        for(i in 1:ncol(al.mat)) {
            plot.these <- code.mat[,i]==1 & lambda.mat[,i]<max.nl
            lines(
                al.mat[code.mat[,i]==1,i],
                prev.mat[code.mat[,i]==1,i], 
                col = my.cols[i]
            )
        }
        
        plot(
            NA ,
            xlim = c(0, 1500) ,
            ylim = c(0, 1) ,
            type = 'l',
            lty = 1 ,
            col = my.cols ,
            bty = 'n',
            xlab = '',
            ylab = ''
        )
        mtext(side=1,text=expression(a[L]), line = 2)
        mtext(side=2,text=expression(h[L]^2), line = 2)
        for(i in 1:ncol(al.mat)) {
            plot.these <- code.mat[,i]==1 & lambda.mat[,i]<max.nl
            lines(
                al.mat[code.mat[,i]==1,i],
                h2l.mat[code.mat[,i]==1,i], 
                col = my.cols[i]
            )
        }
        legend(
            'bottomright' ,
            col = my.cols ,
            lty = 1 ,
            legend = L ,
            bty = 'n'
        )
    }
}


#### Figure 5 ####
{
    rm(list = ls())
    {
        source('scripts/figureFuncs.R')
        getSolnList <-
            function(Y, Z)
                sapply(Y, function(W)
                    sapply(W, function(X)
                        X[Z]))
        ba <- function(Y)
            ifelse(Y > 50, 1, (exp(2 * Y) - 1) / (exp(2 * Y) + 1))
        ya <- function(B)
            0.5 * log((1 + B) / (1 - B))
        library(wesanderson)
        my.cols <- wes_palette('Darjeeling2')
        sim.Ne = 1000
        filename <- paste('resultsFiles/singleEffectPrevalenceResultsTableN', sim.Ne, '.Rdata', sep ='' )
        # sim.results <- get(load(filename))
        # L <- unique(sim.results$target.size)
        # nl <- length(L)
        # nb <- length(unique(sim.results$b))
        # sim.prevs <- matrix(NA, nrow = nb, ncol = nl )
        # for(l in 1:nl){
        #     sim.prevs[,l] <- sim.results[sim.results$target.size == L[l],'sim.prev']
        # }
        png(
            'figures/paperFiguresForRealForRealThisTime/Figure5.png',
            height = 7 ,
            width = 17.75,
            units = 'cm',
            res = 800
        )
        layout(matrix(c(1, 2, 3), nrow = 1))
#### 5A ####
        {
            op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
                       mgp = c(3, 0.6, 0))
            my.bt <- seq(0.01, 0.999, by = 0.001)
            Ne <- 2e4
            u <- 1e-8
            h2 <- 0.5
            C <- 0.1
            L <- c(1.5e6, 5e6, 1.5e7, 5e7, 1.5e8)
            cparam.vals <- sqrt(L * u / (2 * Ne * h2)) / C
            bt.term <- sqrt(my.bt * log((1 + my.bt) / (1 - my.bt)))
            std.dens <- sapply(X = cparam.vals , function(X)
                X * bt.term)
            prevs <- 1 - pnorm(dnorminv(std.dens))
            save(prevs, file = 'output/Fig5Aprevs.Robj')
            save(my.bt, file = 'output/Fig5Abt.Robj')
            sim.results <- get(load("~/Dropbox/msdbPaper/resultsFiles/singleEffectPrevalenceResultsTableN1000.Rdata"))
            split.sims <- split(sim.results,sim.results$target.size)
            plot(
                NA,
                xlim=c(0,1),
                ylim = c(1e-4, 1e-1) ,
                yaxt = 'n',
                log = 'y' ,
                bty = 'n' ,
                )
            abline(v = 1/2 , col = 'seashell2',lty = 2)
            matplot(
                x = my.bt ,
                y = prevs ,
                type = 'l' ,
                lty = 1 ,
                col = my.cols ,
                add = TRUE,
                )
            for(j in 1:length(split.sims)){
              points(
                x = split.sims[[j]]$b ,
                y = split.sims[[j]]$sim.prev,
                col = my.cols[j] ,
                pch = 4 
              )
            }
            axis(side = 2 ,
                 at = 10 ^ seq(-4,-1, 1))
            mtext(
                side = 1 ,
                text = expression(paste(
                    'Relative threshold position, ' , b[T] , sep = ''
                )) ,
                line = 1.5,
                cex = 0.7,
                )
            mtext(
                side = 2 ,
                text = 'Prevalence' ,
                line = 1.5,
                cex = 0.7,
                )
            legend(
                'bottomright' ,
                lty = 1,
                col = my.cols ,
                legend = c(
                    expression(1.5 %*% 10 ^ 6) ,
                    expression(5 %*% 10 ^ 6) ,
                    expression(1.5 %*% 10 ^ 7) ,
                    expression(5 %*% 10 ^ 7) ,
                    expression(1.5 %*% 10 ^ 8)
                ) ,
                bty = 'n' ,
                title = 'Target size, L'
            )
            legend(
              'topleft' ,
              lty = c(NA,1),
              pch = c(4,NA),
              legend = c(
                'Simulations' ,
                'Theory'
              ) ,
              bty = 'n'
            )
            par(op3)
        }
        
#### 5B ####
        
        {
            source('scripts/newSolveTwoEffect.R')
            bt <- 0.5
            Ne = 20000
            C = 0.1
            L <- c(1.5e6, 5e6, 1.5e7, 5e7, 1.5e8)
            u <- 1e-8
            as <- 1
            al <- 100
            yt <- 0.5 * log((1 + bt) / (1 - bt))
            soln <- list()
            prev <- list()
            lambda <- list()
            bs <- list()
            raw.ft <- list()
            raw.Vt <- list()
            std.ft <- list()
            n.pts <- 600
            my.gs <- matrix(NA, nrow = n.pts + 1, ncol = length(L))
            max.pl <- 0.008
            max.nl <- 3
            
            rf <- FALSE
            ## fixed h2
            fixed.h2.prev <- list()
            fixed.h2.bs <- list()
            fixed.h2.deltal <- list()
            fixed.h2.lambda <- list()
            fixed.h2 = 1 / 2
            fixed.h2.soln <- list()
            for (j in seq_along(L)) {
                my.gs[, j] <-
                    sort(unique(c(
                        1 - max.pl / 2, seq (1 - 1 / L[j] , 1 - max.pl , length.out = n.pts)
                    )),decreasing=TRUE)  
                thetaL <- 4 * Ne * L[j] * u
                se.vg <- thetaL * as ^ 2 * bt / yt
                fixed.h2.soln[[j]] <- list()
                ## get initial dl
                {
                    this.bt <- bt # my.grid[k,1]
                    this.al <- al # my.grid[k,2]
                    init.ft <- (4*Ne*C*as)^(-1) * log ((1+this.bt)/(1-this.bt))
                    init.Vas <- 2*L[j]*my.gs[1, j]*u*as*this.bt / (C*init.ft)
                    init.Vtot <- init.Vas/fixed.h2
                    init.tstar <- dnorminv(init.ft,s=sqrt(init.Vtot))
                    init.Ft <- pnorm(init.tstar,sd=sqrt(init.Vtot),lower.tail=FALSE)
                    init.Fta <- pnorm(init.tstar - this.al,sd=sqrt(init.Vtot),lower.tail=FALSE)
                    init.dl <- init.Fta - init.Ft
                }
                ## i = 1
                for (i in 1:nrow(my.gs)) {
                    fixed.h2.soln[[j]][[i]] <- solveTwoEffect1D(
                        bt = bt,
                        Ne = 20000,
                        as = as,
                        al = al,
                        L = L[j],
                        gs = my.gs[i, j],
                        h2 = fixed.h2 ,
                        Ve = NULL ,
                        u = u,
                        C = C,
                        Bval = 1,
                        init.dl = ifelse(i==1,init.dl,prev.dl),
                        init.tstar = ifelse(i==1,init.tstar,prev.tstar)
                    )
                    if(is.na(fixed.h2.soln[[j]][[i]]['dl'])) break
                    prev.dl <- fixed.h2.soln[[j]][[i]]['dl']
                    prev.tstar <- fixed.h2.soln[[j]][[i]]['tstar']
                }
                fixed.h2.prev[[j]] <- sapply(fixed.h2.soln[[j]], function(X)
                    X['prev'])
                fixed.h2.bs[[j]] <- sapply(fixed.h2.soln[[j]], function(X)
                    X['bs'])
                fixed.h2.deltal[[j]] <- sapply(fixed.h2.soln[[j]], function(X)
                    X['deltal'])
                fixed.h2.lambda[[j]] <-
                    sapply(fixed.h2.soln[[j]], function(X)
                        X['lambda'])
            }
        }
        ## fixed h2 plot
        {
            ## my.cols <- wes_palette('Darjeeling2')
            op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
                       mgp = c(3, 0.6, 0))
            y.max <-
                max(fixed.h2.prev[[length(L)]][fixed.h2.lambda[[length(L)]] < max.nl])
            plot(
                NA ,
                type = 'l',
                bty = 'n' ,
                xlim = c(0, max.pl) ,
                ylim = c(1e-3, 0.1) ,
                log = 'y' ,
                yaxt = 'n' ,
                xaxt = 'n'
            )
            axis(2,
                 at = c(0.001, 0.01, 0.1),
                 labels = c('0.001', '0.01', '0.1'))
            axis(1,
                 at = c(0, 0.004, 0.008),
                 labels = c('0', '0.004', '0.008'))
            mtext(
                side = 2 ,
                text = 'Prevalence' ,
                line = 1.5,
                cex = 0.7,
                )
            mtext(
                side = 1 ,
                text = expression(paste(
                    'Fraction of sites with large effects, ' , p[L] , sep = ''
                )) ,
                line = 1.5,
                cex = 0.7,
                )
            for (j in seq_along(fixed.h2.prev)) {
                lines(1 - my.gs[fixed.h2.lambda[[j]] < max.nl , j] ,
                      fixed.h2.prev[[j]][fixed.h2.lambda[[j]] < max.nl] ,
                      col = my.cols[j] ,
                      lty = 1)
            }
            
            ## add points for 5C
            j <- 1
            next.cols <- wes_palette('Darjeeling1')[c(4,2)]#wes_palette('Rushmore1')[c(1, 3, 4)]
            my.pl <- max(1-my.gs [fixed.h2.lambda[[j]] < max.nl,j])
            these <- sapply(1-my.pl,function(X) which.min(abs(my.gs[,j] - X) ) )
            points(
                x = c(0,my.pl) ,
                y = c(prevs[my.bt==0.5,j],fixed.h2.prev[[j]][these]) , 
                pch = 21 , 
                col = next.cols,#[c(1,3)] ,
                bg = adjustcolor(next.cols,alpha.f=0.3) ,
                cex = 1.15
            )
        }
        {
            ## get points for sims
            these.ones <- round(seq(1,max(which(fixed.h2.prev[[1]] > 0 & fixed.h2.lambda[[1]] < 3 )),length.out=18),0)
            tmp.gs <- my.gs[these.ones,1]
            sim.gs <- list()
            for( j in 1:length(L)) {
                smallest.gs <- min(my.gs[fixed.h2.lambda[[j]] < 3,j])
                sim.gs[[j]] <- tmp.gs[tmp.gs >= smallest.gs]
            }
            setwd('~/Documents/academics/LTM_simulation')
            save(sim.gs, file = 'paramFiles/twoEffect.gs.Robj')
            setwd('~/Dropbox/msdbPaper')
        }
        par(op3)
#### 5C ####
        {
            j <- 1
            my.pls <- c(1)*max(1-my.gs [fixed.h2.lambda[[j]] < max.nl,j])
            these <- sapply(1-my.pls,function(X) which.min(abs(my.gs[,j] - X) ) )
            ## normal
            norm.var <- 4 * L[j] * Ne * u * bt / yt * (1/fixed.h2)
            ft.norm <- yt / (2 * Ne * C)
            std.ft.norm <- ft.norm * sqrt (norm.var)
            t.star.norm <- dnorminv(std.ft.norm) * sqrt (norm.var)
            my.xlims <- c(-600, 100)
            my.x <- seq(my.xlims[1], my.xlims[2], length.out = 1500)
            norm.dens <- dnorm(my.x, -t.star.norm, sqrt(norm.var))
            these.solns <- fixed.h2.soln[[j]][these]
            pois.dens <- list()
            op3 <- par(mar = c(3.1, 2.6, 0.4, 0.3),
                       mgp = c(3, 0.6, 0))
            my.cols <- wes_palette('Darjeeling1')[c(4,2)]#wes_palette('Rushmore1')[c(1, 3, 4)]
            y.max <- max(c(unlist(pois.dens)),norm.dens)
            plot(
                NA,
                xlim = my.xlims ,
                ylim = c(0, y.max),
                bty = 'n',
                xaxt = 'n',
                yaxt = 'n'
            )
            axis(
              side = 1,
              at = seq(my.xlims[1],my.xlims[2] , length.out = 12),
              labels=FALSE
            )
            axis(
              side = 2,
              at = seq(0, y.max, length.out = 12),
              labels = FALSE
            )
            mtext(
                side = 2 ,
                text = 'Density' ,
                line = 1.5,
                cex = 0.7,
                )
            mtext(
                side = 1 ,
                text = 'Liability' ,
                line = 1.5,
                cex = 0.7,
                )
            legend(
                x = -600,
                y = y.max*0.8,
                pch = 22,
                title = expression(p[L]) ,
                legend = round(c(0,my.pls),4),
                col = my.cols ,
                pt.bg = adjustcolor(my.cols , alpha.f = 0.3) ,
                bty = 'n' ,
                cex = 1
            )
            t.star.pois <- these.solns[[1]]['tstar'] ##* sqrt(these.solns[[1]]['Vt'])
            this.lambda <- these.solns[[1]]['lambda']
            al <- these.solns[[1]]['al']
            raw.Vt <- these.solns[[1]]['Vt']
            norm.sd <- sqrt(these.solns[[1]]['Vas'] + these.solns[[1]]['Ve'])   
            pois.dens <-
                sapply(my.x, function(X)
                    dPoisConv(X + t.star.pois, this.lambda, norm.sd, al))
            
            ## bottom normal layer
            lines(
              x = my.x [norm.dens > pois.dens ] ,
              y = norm.dens [norm.dens > pois.dens ],
              col = my.cols[1]
            )
            first.switch <- which(diff(norm.dens <= pois.dens)==-1)
            second.switch <- which(diff(norm.dens <= pois.dens)==1)
            lines(
                x = my.x [ 1:first.switch ] ,
                y = norm.dens [1:first.switch],
                col = adjustcolor(my.cols[1],alpha.f=0.3))
            # lines(
            #     x = my.x [ 1:first.switch ] ,
            #     y = norm.dens [1:first.switch],
            #     col = adjustcolor(my.cols[1],alpha.f=0.3))
            lines(
                x = my.x [ second.switch:length(my.x) ] ,
                y = norm.dens [second.switch:length(my.x)],
                col = adjustcolor(my.cols[1],alpha.f=0.3))
            alt.norm.dens <- dnorm(my.x , mean = -t.star.pois,sd=sqrt(raw.Vt))
            lines(
              x = my.x ,
              y = alt.norm.dens,
              lty = 2 ,
              col = my.cols[2]
            )
            lines(my.x,
                  pois.dens,
                  col = my.cols[2])
            ## polygon for top Poisson layer
            pois.polygon1.x <-
                c(my.x, rev(my.x))
            pois.polygon1.y <-
                c(pois.dens, rep(0,length(pois.dens)))
            polygon(
                pois.polygon1.x ,
                pois.polygon1.y ,
                col = adjustcolor(my.cols[2] , alpha.f = 0.08),
                border = NA
            )
            ## polygon for middle layer
            norm.polygon.x <-
                c(my.x [norm.dens > pois.dens], rev(my.x[norm.dens > pois.dens]))
            norm.polygon.y <-
                c(norm.dens[norm.dens > pois.dens], rev(pois.dens[norm.dens > pois.dens] ) )
            polygon(
                norm.polygon.x ,
                norm.polygon.y ,
                col = adjustcolor(my.cols[1] , alpha.f = 0.08),
                border = NA
            )
            ## polygon for top normal layer
            norm.polygon.x <-
                c(my.xlims[1], my.x, my.xlims[2], my.xlims[1])
            norm.polygon.y <- c(y.max, norm.dens, 0, 0)
            polygon(
                norm.polygon.x ,
                norm.polygon.y ,
                col = adjustcolor(my.cols[1] , alpha.f = 0.08),
                border = NA
            )
            lines(x = rep(0,2),y = c(0,y.max*0.96),lty = 1)
            #abline(v = 0 , lwd = 1)
            ##abline(v = -t.star.pois,lty = 2)
        }
        dev.off()
        
    }
    
#### Figure 5-1 ####
  if (FALSE) {
   png(
        'figures/paperFiguresForRealForRealThisTime/SuppFigure5-1.png',
        height = 7 ,
        width = 17.75,
        units = 'cm',
        res = 600
    )
    my.cols <- wes_palette('Darjeeling2')
    op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
               mgp = c(3, 0.6, 0))
    layout(t(matrix(1:2,nrow=2)))
    my.pl <- 1 - my.gs
    my.deltal <- getSolnList(soln,'deltal')
    my.h2l <- getSolnList(soln,'h2l')
    plot(
        NA ,
        type = 'l',
        lty = 1 ,
        col = my.cols,
        xlim = c(0,max.pl) ,
        ylim = c(0,max(my.deltal)) ,
        xlab = '' ,
        ylab = '' ,
        bty = 'n'
    )
    mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=0.7)
    mtext(side=2,text=expression(paste('Risk effect of large effect sites, ', delta^{R}, (a[L]), sep = '')),line=1.5,cex=0.7)
    for( j in 1:ncol(my.deltal)){
        plot.these <- fixed.h2.lambda[[j]] < max.nl
        lines(
            x = my.pl[plot.these,j] ,
            y = my.deltal[plot.these,j] ,
            col = my.cols[j]
        )
    }
    plot(
        NA ,
        type = 'l',
        lty = 1 ,
        col = my.cols,
        xlim = c(0,max.pl) ,
        ylim = c(0,1) ,
        xlab = '' ,
        ylab = '' ,
        bty = 'n'
    )
    mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=0.7)
    mtext(side=2,text='Heritability due to large effect sites',line=1.5,cex=0.7)
    for( j in 1:ncol(my.deltal)){
        plot.these <- fixed.h2.lambda[[j]] < max.nl
        lines(
            x = my.pl[plot.these,j] ,
            y = my.h2l[plot.these,j] ,
            col = my.cols[j]
        )
    }
    legend(
        'topleft',
        legend = L, 
        col = my.cols,
        lty = 1,
        title = 'Target size, L' ,
        bty = 'n' ,
        cex = 0.7
    )
    dev.off()
  }
    
#### Figure 5-2 ####
    if(FALSE){
        my.cols <- wes_palette('Darjeeling2')
        ## 1 - bs
        my.bs <- getSolnList(soln,'bs')
        ## 2 - raw.ft
        my.ft <- getSolnList(soln,'raw.ft')
        ## 3 - raw.Vt
        my.Vt <- getSolnList(soln,'raw.Vt')
        ## 4 - std.ft
        my.std.ft <- getSolnList(soln,'std.ft')
        
        png(
            'figures/paperFiguresForRealForRealThisTime/SuppFigure5-2.png',
            height = 17.75 ,
            width = 26.625,
            units = 'cm',
            res = 600
        )
        
        layout(matrix(1:6,nrow=2))
        op3 <- par(mar = c(3.1, 2.9, 1.4, 0.5),
                   mgp = c(3, 0.6, 0))
        ## 5-2 A
        plot(
            NA ,
            type = 'l',
            lty = 1 ,
            col = my.cols,
            xlim = c(0,max.pl) ,
            ylim = c(0,0.6) ,
            xlab = '' ,
            ylab = '' ,
            bty = 'n'
        )
        text(
            x = 0,
            y = 0.6,
            labels = 'A',
            cex = 1.6
        )
        mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=1)
        mtext(side=2,text='Fraction of small effect sites fixed for + allele',line=1.5,cex=1)
        lines(x = my.pl ,
              y = my.bs[, 1])
        
        ## 5-2 B
        plot(
            NA ,
            type = 'l',
            lty = 1 ,
            col = my.cols,
            xlim = c(0,max.pl) ,
            ylim = c(0,max(unlist(my.ft))*1.2) ,
            xlab = '' ,
            ylab = '' ,
            bty = 'n'
        )
        text(
            x = 0,
            y = max(unlist(my.ft))*1.2,
            labels = 'B',
            cex = 1.6
        )
        mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=1)
        mtext(side=2,text='Raw threshold density',line=1.5,cex=1)
        for( j in 1:ncol(my.ft)){
            ##  plot.these <- lambda[[j]] < max.nl
            lines(
                x = my.pl[,j] ,
                y = my.ft[,j]
            )
        }
        
        ## 5-2 C
        plot(
            NA ,
            type = 'l',
            lty = 1 ,
            col = my.cols,
            xlim = c(0,max.pl) ,
            ylim = c(0,max(sqrt(unlist(my.Vt)))) ,
            xlab = '' ,
            ylab = '' ,
            bty = 'n'
        )
        text(
            x = 0,
            y = max(sqrt(unlist(my.Vt))),
            labels = 'C',
            cex = 1.6
        )
        mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=1)
        mtext(side=2,text='Std. dev. of liability',line=1.5,cex=1)
        for( j in 1:ncol(my.Vt)){
            ##  plot.these <- lambda[[j]] < max.nl
            lines(
                x = my.pl[,j] ,
                y = sqrt(my.Vt[,j]) ,
                col = my.cols[j]
            )
        }
        
        ## 5-2 D
        plot(
            NA ,
            type = 'l',
            lty = 1 ,
            col = my.cols,
            xlim = c(0,max.pl) ,
            ylim = c(0,max(unlist(my.std.ft))*1.2) ,
            xlab = '' ,
            ylab = '' ,
            bty = 'n'
        )
        text(
            x = 0,
            y = max(unlist(my.std.ft))*1.2,
            labels = 'D',
            cex = 1.6
        )
        mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=1)
        mtext(side=2,text='Standardized threshold density',line=1.5,cex=1)
        for( j in 1:ncol(my.Vt)){
            ##  plot.these <- lambda[[j]] < max.nl
            lines(
                x = my.pl[,j] ,
                y = my.std.ft[,j] ,
                col = my.cols[j]
            )
        }
        



        ## 5-2 E/F
        rho3 <- list()
        rho4 <- list()
        rho5 <- list()
        for(i in 1:length(soln)){
            rho3[[i]] <- numeric()
            rho4[[i]] <- numeric()
            rho5[[i]] <- numeric()
            for( j in 1:length(soln[[i]])){
                t.star.pois <- soln[[i]][[j]]['tstar']
                Vt <- soln[[i]][[j]]['raw.Vt']
                this.lambda <- soln[[i]][[j]]['lambda']
                if(this.lambda>3) next
                std.al <- soln[[1]][[j]]['al']/sqrt(Vt)
                norm.sd <- sqrt((soln[[i]][[j]]['raw.Vas'] + soln[[i]][[j]]['raw.Ve'])/Vt)
                rho3[[i]][j] <- integrate( f=function(X) X^3*dPoisConv(X, this.lambda, norm.sd, std.al),lower=-Inf,upper=Inf)$value
                rho4[[i]][j] <- integrate( f=function(X) X^4*dPoisConv(X, this.lambda, norm.sd, std.al),lower=-Inf,upper=Inf)$value
                rho5[[i]][j] <- integrate( f=function(X) X^5*dPoisConv(X, this.lambda, norm.sd, std.al),lower=-Inf,upper=Inf)$value
            }
        }

        max.rho3 <- max(unlist(rho3))
        plot(
            NA ,
            type = 'l',
            lty = 1 ,
            col = my.cols,
            xlim = c(0,max.pl) ,
            ylim = c(0,max.rho3) ,
            xlab = '' ,
            ylab = '' ,
            bty = 'n'
        )
        text(
            x = 0,
            y = max.rho3,
            labels = 'E',
            cex = 1.6
        )
        abline(h = 0 , lty = 2 )
        mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=1)
        mtext(side=2,text='Standardized third central moment',line=1.5,cex=1)
        for( i in 1:length(rho3)){
            ##  plot.these <- lambda[[j]] < max.nl
            lines(
                x = my.pl[1:length(rho3[[i]]),i] ,
                y = rho3[[i]] ,
                col = my.cols[i]
            )
        }

        max.rho5 <- max(unlist(rho5))
        plot(
            NA ,
            type = 'l',
            lty = 1 ,
            col = my.cols,
            xlim = c(0,max.pl) ,
            ylim = c(0,max.rho5) ,
            xlab = '' ,
            ylab = '' ,
            bty = 'n'
        )
        text(
            x = 0,
            y = max.rho5,
            labels = 'F',
            cex = 1.6
        )
        legend(
            x=max.pl/3,
            y=max.rho5,
            lty = 1,
            col = my.cols ,
            legend = c(
                expression(1.5 %*% 10 ^ 6) ,
                expression(5 %*% 10 ^ 6) ,
                expression(1.5 %*% 10 ^ 7) ,
                expression(5 %*% 10 ^ 7) ,
                expression(1.5 %*% 10 ^ 8)
            ) ,
            bty = 'n' ,
            title = 'Target size, L'
        )
        abline(h = 0 , lty = 2 )
        mtext(side=1,text='Fraction of sites with large effects',line=1.5,cex=1)
        mtext(side=2,text='Standardized fifth central moment',line=1.5,cex=1)
        for( i in 1:length(rho3)){
            ##  plot.these <- lambda[[j]] < max.nl
            lines(
                x = my.pl[1:length(rho5[[i]]),i] ,
                y = rho5[[i]] ,
                col = my.cols[i]
            )
        }
        dev.off()
    }
}


#### Figure 6 ####
my.dir <- 'figures/smilePlotData'
my.filenames <- list.files(my.dir)
my.data <-
    lapply(my.filenames, function(x)
        read.delim(paste(my.dir, x, sep = '/'), sep = '\t', header = T))
library(wesanderson)
my.cols <- wes_palette("AsteroidCity3")[c(2,4)]

{
    png(
        'figures/paperFiguresForRealForRealThisTime/Figure6.png',
        height = 17.75 ,
        width = 17.75,
        units = 'cm',
        res = 600
    )
    layout(matrix(c(1,1,1, 2, 3, 6,4,5,6), nrow = 3),widths = c(1,12,12),heights = c(
                                                                             12,12,1))
    
    par(mar=c(0,0,0,0))
    plot.new()
    text(x=0.5,y=0.5,labels='Effect Size',srt=90,cex = 2.5)
    op6 <- par(mar = c(2.9, 3, 2, 0),
               mgp = c(3, 0.6, 0))
    
    #### CVD ####
    plot(
        x = my.data[[1]][, 1] ,
        y = my.data[[1]][, 2],
        ylim = c(0.005,0.05),
        xlim = c(0,1),
        pch = 21 ,
        bty = 'n',
        xlab = '',
        ylab = '' ,
        log = 'y' ,
        #main = 'CVD' ,
        bg = my.cols[1] ,
        col = my.cols[2] ,
        yaxt = 'n'
    )
    axis(side = 2, at = c(0.005,0.015,0.05), labels = c('0.005','0.015','0.05'))
    
    
    #### IBD ####
    plot(
        x = my.data[[3]][, 1] ,
        y = my.data[[3]][, 2],
        ylim = c(0.05,0.5),
        xlim = c(0,1),
        pch = 21 ,
        bty = 'n',
        xlab = '',
        ylab = '' ,
        log = 'y',
        #main = 'IBD',
        bg = my.cols[1] ,
        col = my.cols[2] ,
        yaxt = 'n'
    )
    axis(side = 2, at = c(0.05,0.15,0.5),labels = c('0.05','0.15','0.5'))
    
    
    #### Schizophrenia ####
    plot(
        x = my.data[[4]][, 1] ,
        y = my.data[[4]][, 2],
        ylim = c(0.045,0.25),
        xlim = c(0,1),
        pch = 21 ,
        bty = 'n',
        xlab = '',
        ylab = '' ,
        log = 'y',
        #main = 'SCZ',
        bg = my.cols[1] ,
        col = my.cols[2] ,
        yaxt = 'n'
    )
    axis(side = 2, at = c(0.05,0.125,0.25), labels = c('0.05','0.125','0.25'))
    
    #### T2D ####
    plot(
        x = my.data[[5]][, 1] ,
        y = my.data[[5]][, 2],
        ylim = c(0.03,0.3), #c(0.035,0.3),
        xlim = c(0,1),
        pch = 21 ,
        bty = 'n',
        xlab = '',
        ylab = '' ,
        log = 'y',
        yaxt = 'n',
        #main = 'T2D',
        bg = my.cols[1] ,
        col = my.cols[2]
    )
    axis(side = 2, at = c(0.03,0.1,0.3))#c(0.04,0.1,0.25))
    par(op6)
    par(mar=c(0,0,0,0))
    plot.new()
    text(x=0.5,y=0.5,labels='Frequency of risk increasing allele',cex = 2.5)
    dev.off()
}


#### Supplementary Figures ####


#### Figure S1 ####
{
  my.as <- seq(0,6,length.out = 500)
  prevs <- c(0.001,0.003,0.01,0.03,0.1)
  thrs <- qnorm(1-prevs)
  deltas <- lapply(thrs , function(THR) 
    pnorm(THR-my.as,lower.tail = FALSE) - pnorm(THR,lower.tail = FALSE)
  )
  library(wesanderson)
  my.cols <- wes_palette("Darjeeling1")
  png(
    'figures/paperFiguresForRealForRealThisTime/suppFigures/FigureS1.png',
    height = 10.75 ,
    width = 10.75,
    units = 'cm',
    res = 600
  )
  op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
             mgp = c(3, 0.6, 0))
  plot(
    NA, 
    xlim = c(0,6),
    ylim = c(0,1),
    bty = 'n',
    xlab = '',
    ylab = ''
  )
  for ( i in 1:length(deltas)){
    lines(
      x = my.as ,
      y = deltas[[i]] ,
      lty = 1,
      col = my.cols[i]
    )
  }
  mtext(
    text = expression(paste('Standardized liability scale effect size, ' , a/sqrt(V[T]) , sep = '')),
    side = 1,
    line = 2
  )
  mtext(
    text = expression(paste('Risk scale effect size, ' , delta[R](a) , sep = '')),
    side = 2,
    line = 1.4
  )
  legend(
    'topleft',
    legend=prevs,
    col = my.cols,
    lty = 1 ,
    bty = 'n',
    title = 'Prevalence'
  )
  dev.off()
}

#### Figure S_PrevDensityRelationship ####
{
  prevs <- seq(0.1,1e-4,length.out=10000)
  thr <- qnorm(prevs,lower.tail=F)
  dens <- dnorm(thr)

  png(
    'figures/paperFiguresForRealForRealThisTime/suppFigures/prevalenceDensityRelationshipSFigure.png',
    height = 7 ,
    width = 17.75*2/3,
    units = 'cm',
    res = 800
  )
  op3 <- par(mar = c(3.1, 3.3, 0.8, 0.1),
             mgp = c(3, 0.6, 0))
  plot(
    x = dens, 
    y = prevs, 
    xlim = c(0,0.2),
    ylim = c(0,0.1),
    bty = 'n',
    type = 'l' ,
    xlab = '',
    ylab = ''
  )
  mtext('Threshold density of a standard Normal r.v.', side = 1 , line = 1.6, cex = 1.2)
  mtext('Upper tail probability of', side = 2 , line = 2.3, cex = 1)
  mtext('a standard Normal r.v.', side = 2 , line = 1.4, cex = 1)
  abline(a = 0 , b = 0.5 , lty =2 )
  legend('topleft',
         legend = c('Truth',expression(1-Phi(x) %~~% 0.4 * phi(x)) ),
         lty = c(1,2),
         bty = 'n')
  dev.off()
}


#### Figure S_Flattening ####
{
  rm(list=ls())
  getB <- function(y) ifelse(y>50,1,(exp(2*y)-1) / (exp(2*y)+1))
  getDeltaR <- function(a,thr) pnorm(thr-a,lower.tail=FALSE) - pnorm(thr,lower.tail=FALSE)
  prev <- c(0.001,0.003,0.01,0.03,0.1)
  thr <- qnorm(prev,lower.tail=FALSE)
  my.as <- exp(seq(log(1e-5),log(6),length.out = 600))
  Ne <- 20000
  u <- 1e-8
  cost <- 0.1
  theta <- 4*Ne*u
  my.deltas <- lapply(thr,function(THR) getDeltaR(my.as,THR))
  my.gammas <- lapply(my.deltas,function(DELTA) 2*Ne*DELTA*cost )
  my.bs <- lapply(my.gammas,function(GAMMA) getB(GAMMA) )
  
  my.vars <- mapply(function(A,B,GAMMA) my.as^2 * theta * B / GAMMA , B = my.bs , GAMMA = my.gammas )
  y.max <- max(my.vars)
  my.cols <- wes_palette("Darjeeling1")
  
  TAU <- 1
  theta*integrate(function(x) 4*exp(-2*TAU*x*(1-x)),lower=0,upper=1/2)$value
  
  
  png(
    'figures/paperFiguresForRealForRealThisTime/suppFigures/flattening.png',
    height = 7 ,
    width = 17.75*2/3,
    units = 'cm',
    res = 800
  )
  layout(matrix(c(1, 2), nrow = 1))
  op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
             mgp = c(3, 0.6, 0))
  plot(
    NA,
    xlim = c(0,6),
    ylim = c(0,y.max),
    bty = 'n' ,
    xaxt = 'n',
    yaxt = 'n'
  )
  for(i in 1:ncol(my.vars)) {
    lines(
      my.as,
      my.vars[,i] ,
      col = my.cols[i]
    )
  }
  mtext(
    text = expression(paste('Standardized liability scale effect size, ' , a/sqrt(V[T]) , sep = '')),
    side = 1,
    line = 1.6,
    cex = 0.6
  )
  mtext(
    text = expression(paste('Per-site contribution to heritability' , sep = '')),
    side = 2,
    line = 1.6 ,
    cex = 0.6
  )
  axis(side=1,cex.axis=0.6)
  axis(side=2,cex.axis=0.6)
  legend(
    'topright',
    legend=prev,
    col = my.cols,
    lty = 1 ,
    bty = 'n',
    title = 'Prevalence',
    cex = 0.6
  )
     
  x.min2 <- min(unlist(my.gammas))
  x.max2 <- max(unlist(my.gammas))
  plot(
    NA,
    xlim = c(1e-2,x.max2),
    ylim = c(0,y.max),
    bty = 'n' ,
    log = 'x',
    xaxt = 'n',
    yaxt = 'n'
  )
  abline(v = 1 , lty = 2)
  for(i in 1:ncol(my.vars)) {
    lines(
      my.gammas[[i]],
      my.vars[,i] ,
      col = my.cols[i]
    )
  }
  mtext(
    text = expression(paste('Population scaled selection coefficient, ' , gamma(a) , sep = '')),
    side = 1,
    line = 1.6,
    cex = 0.6 ,
    at = 2.8
  )
  mtext(
    text = expression(paste('Per-site contribution to heritability' , sep = '')),
    side = 2,
    line = 1.6 ,
    cex = 0.6
  )
  axis(side=1,cex.axis=0.6)
  axis(side=2,cex.axis=0.6)
  dev.off()

  
  my.delta.vars <- mapply(function(DELTA,B,GAMMA) DELTA^2 * theta * B / GAMMA, DELTA = my.deltas , B = my.bs , GAMMA = my.gammas )
  delta.y.max <- max(my.delta.vars)
  png(
    'figures/paperFiguresForRealForRealThisTime/suppFigures/riskFlattening.png',
    height = 7 ,
    width = 17.75*2/3,
    units = 'cm',
    res = 800
  )
  layout(matrix(c(1, 2), nrow = 1))
  op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
             mgp = c(3, 0.6, 0))
  plot(
    NA,
    xlim = c(0,6),
    ylim = c(0,delta.y.max),
    bty = 'n' ,
    xaxt = 'n',
    yaxt = 'n'
  )
  for(i in 1:ncol(my.delta.vars)) {
    lines(
      my.as,
      my.delta.vars[,i] ,
      col = my.cols[i]
    )
  }
  mtext(
    text = expression(paste('Standardized liability scale effect size, ' , a/sqrt(V[T]) , sep = '')),
    side = 1,
    line = 1.6,
    cex = 0.6
  )
  mtext(
    text = expression(paste('Per-site additive contribution to risk variance' , sep = '')),
    side = 2,
    line = 1.6 ,
    cex = 0.6
  )
  axis(side=1,cex.axis=0.6)
  axis(side=2,cex.axis=0.6)
  legend(
    'bottomright',
    legend=prev,
    col = my.cols,
    lty = 1 ,
    bty = 'n',
    title = 'Prevalence',
    cex = 0.6
  )
  
  x.min2 <- min(unlist(my.gammas))
  x.max2 <- max(unlist(my.gammas))
  plot(
    NA,
    xlim = c(1e-2,x.max2),
    ylim = c(0,delta.y.max),
    bty = 'n' ,
    log = 'x',
    xaxt = 'n',
    yaxt = 'n'
  )
  abline(v = 1 , lty = 2)
  for(i in 1:ncol(my.delta.vars)) {
    lines(
      my.gammas[[i]],
      my.delta.vars[,i] ,
      col = my.cols[i]
    )
  }
  mtext(
    text = expression(paste('Population scaled selection coefficient, ' , gamma(a) , sep = '')),
    side = 1,
    line = 1.6,
    cex = 0.6 ,
    at = 2.8
  )
  mtext(
    text = expression(paste('Per-site additive contribution to risk variance' , sep = '')),
    side = 2,
    line = 1.6 ,
    cex = 0.6
  )
  axis(side=1,cex.axis=0.6)
  axis(side=2,cex.axis=0.6)
  dev.off()
}

#### Figure S_MeanFrequency ####
{
  my.cols <- wes_palette("Darjeeling1")
  my.gamma <- exp(seq(log(0.1),log(20),length.out=500))
  my.N <- c(2e3,2e4,2e5,2e6,2e7)
  freqSpec <- function(x,y) 2 / (1+exp(-2*y)) * exp(-2*y*x) / (x*(1-x) )
  
  realMean <- function(y,pstar) {
    num <- integrate(function(x) x*freqSpec(x,y),lower=pstar,upper=1-pstar )$value
    denom <- integrate(function(x) freqSpec(x,y),lower=pstar,upper=1-pstar )$value
    num / denom
  }
  approxMean <- function(y,Ne) 1 / ( 1 + exp(2*y*(1-2/log(4*Ne))))
  fixMean <- function(y) 1 / (1+exp(2*y))
  
  real.means <- do.call(rbind,lapply(my.gamma,function(Y) sapply(my.N,function(NE) realMean(Y,1/(4*NE)))))
  fix.means <- fixMean(my.gamma)
  approx.means <- t(sapply(my.gamma,function(Y) approxMean(Y,my.N) ))
  
  png(
    'figures/paperFiguresForRealForRealThisTime/suppFigures/figureSMean.png',
    height = 7 ,
    width = 17.75*2/3,
    units = 'cm',
    res = 800
  )
  op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
             mgp = c(3, 0.6, 0))
  plot(
    NA,
    xlim = c(0.1,10),
    ylim = c(3e-3,1/2) , 
    bty  = 'n' ,
    log = 'xy' ,
    yaxt = 'n' ,
    xaxt = 'n' ,
    xlab = '',
    ylab = ''
  )
  mtext(text='Scaled selection coefficient',side =1 , line =1.7)
  mtext(text='Mean risk allele frequency',side =2 , line =1.7)
  axis(side = 1, at = c(0.1,0.3,1,3,10) , labels = c('0.1','0.3','1','3','10'))
  axis(side = 2, at = 5*10^c(-3:-1) , labels = c('0.005','0.05','0.5'))
  for ( i in 1:ncol(real.means)){
    lines( 
      x = my.gamma ,
      y = real.means[,i] ,
      col = my.cols[i])
  }
  for ( i in 1:ncol(real.means)){
    lines( 
      x = my.gamma ,
      y = approx.means[,i] ,
      col = my.cols[i] ,
      lty = 3)
  }
  lines(
    x = my.gamma ,
    y = fix.means ,
    lty = 2
  )
  legend(
    'topright' ,
    legend= c('Exact','Approximate','Fixed') ,
    lty = c(1,3,2) ,
    bty = 'n'
  )
  legend(
    'bottomleft' ,
    legend= c('2,000','20,000','200,000','2,000,000','20,000,000') ,
    col = my.cols ,
    title = 'Population Size' ,
    bty = 'n' ,
    lty = 1
  )
  dev.off()
}



#### Figure S_Derived ####
{
  rm(list=ls())
  library(wesanderson)
  pder <- function(x,y) (exp(2*y) - exp(2*y*x)) / (exp(2*y) - 1)  
  
  my.y <- 1:5
  my.x <- seq(0,1,by = 0.001)
  my.pder <- sapply(my.y,function(Y)pder(my.x,Y))
  
  png(
    'figures/paperFiguresForRealForRealThisTime/suppFigures/pDerived.png',
    height = 7 ,
    width = 17.75*2/3,
    units = 'cm',
    res = 800
  )
  op3 <- par(mar = c(3.1, 2.6, 0.4, 0.1),
             mgp = c(3, 0.6, 0))
  my.cols <- wes_palette('Darjeeling1')
  plot(
    NA,
    xlim=c(0,1),
    ylim=c(0,1),
    bty = 'n',
    xlab = '',
    ylab = ''
  )
  mtext(side = 1 , line = 1.7 , text = 'Risk allele frequency, x')
  mtext(side = 2 , line = 1.7 , text = 'Probability risk allele is derived', at = 0.45)
  abline(a = 1 , b = -1 , lty = 2)
  for(i in 1:ncol(my.pder)){
    lines(
      x = my.x ,
      y = my.pder[,i] ,
      col = my.cols[i]
    )
  }
  legend('bottomleft',
         lty = c(2,rep(1,5)) ,
         col = c('black',my.cols) ,
         legend = c(0,my.y) ,
         title = 'Scaled coefficient' ,
         cex = 0.7, 
         bty = 'n')
  dev.off()
}

#### Figure S_MeanCoefficient ####
{
  rm(list=ls())
  setwd('~/Dropbox/msdbPaper')
  library(wesanderson)
  balpha = function(a){
    ifelse(a < 100, (exp(2 * a) - 1) / (exp(2 * a) + 1), 1)
  }
  bt.diff2 = function(my.mean, my.cv, my.bt,my.q=1e-5) {
    my.shape = 1 / my.cv ^ 2
    my.rate = 1 / (my.cv ^ 2 * my.mean)
    my.lower = qgamma(my.q, my.shape, my.rate)
    my.upper = qgamma(1 - my.q, my.shape, my.rate)
    ## balpha(my.mean)
    myFunc <- function(a)
      dgamma(a, my.shape, my.rate) * a * balpha(a)
    my.bt - (1 / my.mean) * integrate(myFunc, lower = my.lower, upper = my.upper)$value
  }
  my.bts = seq(1e-3, 0.999, length.out = 1000)
  my.cvs = sqrt(c(0.1, 1, 3, 10))
  se.gammas = 0.5 * log((1 + my.bts) / (1 - my.bts))
  Ne = 1e4
  U = 0.2
  h2 = 1 / 2
  cost = 1 / 2
  tmp = list()
  my.gammas = list()
  for (j in 1:length(my.cvs)) {
    tmp[[j]] = list()
    my.gammas[[j]] = rep(NA, length(my.bts))
    for (i in 1:length(my.bts)) {
      my.lower = 0.01 * se.gammas[i]
      my.upper =  2000 * se.gammas[i]
      tmp[[j]][[i]] = uniroot(
        function(X)
          bt.diff2(X, my.cv = my.cvs[j], my.bts[i]),
        lower = my.lower,
        upper = my.upper
      )
      if (tmp[[j]][[i]]$root < 0.002)
        next
      my.gammas[[j]][i] = tmp[[j]][[i]]$root
    }
  }
  
  {
    cex.lab = 1.3
    cex.axis = 1.2
    png(
      'figures/paperFiguresForRealForRealThisTime/suppFigures/MeanCoefficientSFigure.png',
      height = 10 ,
      width = 22,
      units = 'cm',
      res = 500
    )
    layout(matrix(c(1, 2), nrow = 1))
    op = par(mar = c(3.6, 4.6, 1, 0.5) + 0.1)
    ## 1
    plot(
      x = my.bts,
      y = se.gammas,
      type = 'l' ,
      xlab = '',
      ylab = '',
      cex.axis = cex.axis,
      cex.lab = cex.lab,
      ylim = c(0, 10),
      lwd = 2,
      lty = 1,
      bty = 'n'
    )
    mtext(
      side = 2,
      text = 'Mean population scaled',
      cex = cex.lab,
      line = 3.5
    )
    mtext(
      side = 2,
      text = expression(paste(
        'selection coefficient, ', bar(gamma(alpha)), sep = ''
      )),
      cex = cex.lab,
      line = 1.9
    )
    mtext(
      side = 1,
      text = expression(paste(
        'Relative threshold position, ', b[T], sep = ''
      )),
      cex = cex.lab,
      line = 2.4
    )
    my.cols <- wes_palette('GrandBudapest1', length(my.gammas))
    for (j in 1:length(my.gammas)) {
      plot.these = !is.na(my.gammas[[j]]) & my.gammas[[j]] < 10
      lines(
        x = my.bts[plot.these],
        y = my.gammas[[j]][plot.these],
        lty = 2,
        lwd = 2,
        col = my.cols[j]
      )
    }
    legend(
      'topleft',
      legend = c(0, my.cvs ^ 2),
      col = c('black', my.cols),
      lty = c(1, rep(2, length(my.cvs))),
      lwd = c(2, rep(2, length(my.cvs))),
      bty = 'n',
      title = expression(CV[alpha] ^ 2),
      ##'Effect dist\'n coef. of var.',
      cex = 1
    )
    ##abline(v = bt.inflect, lty = 2 , lwd = 2)
    par(op)
    
    
    op = par(mar = c(3.6, 4.6, 1, 0.5) + 0.1)
    plot(
      NA,
      xlab = '',
      ylab = '',
      cex.axis = cex.axis,
      cex.lab = cex.lab,
      xlim = c(0, 1),
      ylim = c(0, 1.2),
      lwd = 2,
      lty = 1,
      bty = 'n'
    )
    mtext(
      side = 2,
      text = expression(paste(
        'Fold change in ', bar(gamma(alpha)) , ' relative', sep = ''
      )),
      cex = cex.lab,
      line = 3.3
    )
    mtext(
      side = 2,
      text = 'to single effect model',
      cex = cex.lab,
      line = 2.1
    )
    mtext(
      side = 1,
      text = expression(paste(
        'Relative threshold position, ', b[T], sep = ''
      )),
      cex = cex.lab,
      line = 2.4
    )
    for (i in seq_along(my.cvs)) {
      abline(
        a = 1 / (1 + my.cvs[i] ^ 2) ,
        b = 0,
        col = my.cols[i],
        lwd = 0.3,
        lty = 3
      )
    }
    abline(a = 1 , b = 0 , lty = 2)
    for (j in 1:length(my.gammas)) {
      plot.these = !is.na(my.gammas[[j]]) & my.gammas[[j]] < 10
      lines(
        x = my.bts[plot.these],
        y = my.gammas[[j]][plot.these] / se.gammas[plot.these],
        lty = 2,
        lwd = 2,
        col = my.cols[j]
      )
    }
    par(op)
    dev.off()
  }
  {
    cex.lab = 1.3
    cex.axis = 1.2
    png(
      'figures/paperFiguresForRealForRealThisTime/suppFigures/EffectDistPrevSFigure.png',
      height = 12 ,
      width = 18,
      units = 'cm',
      res = 500
    )
    op = par(mar = c(3.6, 4.6, 1, 0.5) + 0.1)
    plot(
      NA,
      xlab = '',
      ylab = '',
      cex.axis = cex.axis,
      cex.lab = cex.lab,
      xlim = c(0, 1),
      ylim = c(0, 1.2),
      lwd = 2,
      lty = 1,
      bty = 'n'
    )
    mtext(
      side = 2,
      text = expression(paste(
        'Approximate fold change in prevalence', sep = ''
      )),
      cex = cex.lab,
      line = 3.3
    )
    mtext(
      side = 2,
      text = 'relative to single effect model',
      cex = cex.lab,
      line = 2.1
    )
    mtext(
      side = 1,
      text = expression(paste(
        'Relative threshold position, ', b[T], sep = ''
      )),
      cex = cex.lab,
      line = 2.4
    )
    for (i in seq_along(my.cvs)) {
      abline(
        a = sqrt(1 / (1 + my.cvs[i] ^ 2)) ,
        b = 0,
        col = my.cols[i],
        lwd = 0.3,
        lty = 3
      )
    }
    abline(a = 1 , b = 0 , lty = 2)
    for (j in 1:length(my.gammas)) {
      plot.these = !is.na(my.gammas[[j]]) & my.gammas[[j]] < 10
      lines(
        x = my.bts[plot.these],
        y = sqrt(my.gammas[[j]][plot.these] / se.gammas[plot.these]),
        lty = 2,
        lwd = 2,
        col = my.cols[j]
      )
    }
    legend(
      'bottomright',
      legend = c(0, my.cvs ^ 2),
      col = c('black', my.cols),
      lty = c(1, rep(2, length(my.cvs))),
      lwd = c(2, rep(2, length(my.cvs))),
      bty = 'n',
      title = expression(CV[alpha] ^ 2),
      ##'Effect dist\'n coef. of var.',
      cex = 1
    )
    par(op)
    dev.off()
  }
}


#### Figure S large effect grid ####
{
  rm(list=ls())
  setwd('~/Dropbox/msdbPaper')
  library(wesanderson)
  source('scripts/newSolveTwoEffect.R')
  bt <- c(0.1,0.5,0.9)
  al <- c(50,100,150)
  my.grid <- expand.grid(bt,al)
  
  Ne = 20000
  C = 0.1
  #ft <- log((1 + bt) / (1 - bt)) / (2 * Ne * C)
  L <- c(1.5e6, 5e6, 1.5e7, 5e7, 1.5e8)
  u <- 1e-8
  as <- 1
  ##deltal <- 0.2
  yt <- 0.5 * log((1 + bt) / (1 - bt))
  n.pts <- 301
  my.gs <- matrix(NA, nrow = n.pts + 1, ncol = length(L))
  max.pl <- 0.008
  max.nl <- 3
  
  #### prevalence ####
  png(
    'figures/paperFiguresForRealForRealThisTime/suppFigures/twoEffectPrevSFigure.png',
    height = 17.75 ,
    width = 17.75,
    units = 'cm',
    res = 800
  )
  layout(matrix(1:9,nrow=3))
  for(k in 1:nrow(my.grid)) {
    rf <- FALSE
    ## fixed h2
    prev <- list()
    h2 = 1 / 2
    tmp.soln <- list()
    soln <- list()
    
    for (j in seq_along(L)) {
      my.gs[, j] <-
        sort(unique(c(
          1 - max.pl / 2, seq (1 - 1 / L[j] , 1 - max.pl , length.out = n.pts)
        )),decreasing=TRUE)  
      tmp.soln[[j]] <- list()
      soln[[j]] <- list()
      
      ## get initial dl
      {
        this.bt <- my.grid[k,1]
        this.al <- my.grid[k,2]
        init.ft <- (4*Ne*C*as)^(-1) * log ((1+this.bt)/(1-this.bt))
        init.Vas <- 2*L[j]*my.gs[1, j]*u*as*this.bt / (C*init.ft)
        init.Vtot <- init.Vas/h2
        init.tstar <- dnorminv(init.ft,s=sqrt(init.Vtot))
        init.Ft <- pnorm(init.tstar,sd=sqrt(init.Vtot),lower.tail=FALSE)
        init.Fta <- pnorm(init.tstar - this.al,sd=sqrt(init.Vtot),lower.tail=FALSE)
        init.dl <- init.Fta - init.Ft
      }
        
      ## source('scripts/newSolveTwoEffect.R')
      for (i in 1:nrow(my.gs)) {
        print(paste(k,j,i,sep=';'))
        tmp.soln[[j]][[i]] <- solveTwoEffect1D(
          bt = my.grid[k,1],
          Ne = 20000,
          as = as,
          al = my.grid[k,2],
          L = L[j],
          gs = my.gs[i, j],
          h2 = h2 ,
          Ve = NULL ,
          u = u,
          C = C,
          Bval = 1,
          init.dl = ifelse(i==1,init.dl,prev.dl),
          init.tstar = ifelse(i==1,init.tstar,prev.tstar)
        )
        if(is.na(tmp.soln[[j]][[i]]['dl'])) break
        prev.dl <- tmp.soln[[j]][[i]]['dl']
        prev.tstar <- tmp.soln[[j]][[i]]['tstar']
      }
      soln[[j]] <- do.call(rbind,tmp.soln[[j]])
      soln[[j]] <- soln[[j]][!is.na(soln[[j]][,1]),]
    }
    dl <- lapply(soln,function(X)X[,'dl'])
    tstar <- lapply(soln,function(X)X[,'tstar'])
    lambda <- lapply(soln,function(X)X[,'lambda'])
    prev <- lapply(soln,function(X)X[,'prev'])
    
    my.cols <- wes_palette('Darjeeling2')
    if(k == 1){
      marg <- 1.6
      op3 <- par(mar = c(3.1, 2.6 + marg, 0.4 + marg, 0.6),
                 mgp = c(3, 0.6, 0))
    } else if(k %in% 2:3){
      op3 <- par(mar = c(3.1, 2.6 + marg, 0.4, 0.6),
                 mgp = c(3, 0.6, 0))
    } else if(k %in% c(4,7)){
      op3 <- par(mar = c(3.1, 2.6, 0.4 + marg, 0.6),
                 mgp = c(3, 0.6, 0))
    }
    # y.max <- max(unlist(prev))
    
    plot(
      NA ,
      type = 'l',
      bty = 'n' ,
      xlim = c(0, max.pl) ,
      ylim = c(1e-4, 0.1) ,
      log = 'y' ,
      yaxt = 'n' ,
      xaxt = 'n',
      xlab = '',
      ylab = ''
    )
    axis(2,
         at = c(1e-4,0.001, 0.01, 0.1),
         labels = c('0.0001','0.001', '0.01', '0.1'))
    axis(1,
         at = c(0, 0.004, 0.008),
         labels = c('0', '0.004', '0.008'))
    mtext(
      side = 2 ,
      text = 'Prevalence' ,
      line = 1.5,
      cex = 0.7,
    )
    mtext(
      side = 1 ,
      text = expression(paste(
        'Fraction of sites with large effects, ' , p[L] , sep = ''
      )) ,
      line = 1.5,
      cex = 0.7,
    )
    if(k == 1){
      mtext(
        side = 2 ,
        text = expression(
          paste(b[T], ' = 0.1', sep ='')) ,
          line = 2.5 ,
          cex = 0.7
      )
      mtext(
        side = 3 ,
        text = expression(
          paste(a[L]/a[S], ' = 50', sep ='')) ,
        line = 0.2 ,
        cex = 0.7
      )
    }
    if(k == 2){
      mtext(
        side = 2 ,
        text = expression(
          paste(b[T], ' = 0.5', sep ='')) ,
        line = 2.5 ,
        cex = 0.7
      )
    }
    if(k == 3){
      mtext(
        side = 2 ,
        text = expression(
          paste(b[T], ' = 0.9', sep ='')) ,
        line = 2.5 ,
        cex = 0.7
      )
    }
    if(k == 4){
      mtext(
        side = 3 ,
        text = expression(
          paste(a[L]/a[S], ' = 100', sep ='')) ,
        line = 0.2 ,
        cex = 0.7
      )
      legend(
        'topright' ,
        lty = 1,
        col = my.cols ,
        legend = c(
          expression(1.5 %*% 10 ^ 6) ,
          expression(5 %*% 10 ^ 6) ,
          expression(1.5 %*% 10 ^ 7) ,
          expression(5 %*% 10 ^ 7) ,
          expression(1.5 %*% 10 ^ 8)
        ) ,
        bty = 'n' ,
        title = 'Target size, L'
      )
    }
    if(k == 7){
      mtext(
        side = 3 ,
        text = expression(
          paste(a[L]/a[S], ' = 200', sep ='')) ,
        line = 0.2 ,
        cex = 0.7
      )
    }
    for (j in 1:length(prev)) {
      plot.these <- which(prev[[j]] > 0 & lambda[[j]] < 3 )
      lines(1 - my.gs[plot.these,j] ,
            prev[[j]][plot.these] ,
            col = my.cols[j] ,
            lty = 1)
    }
  }
  dev.off()
  
  
  #### variance ####
  {
    png(
      'figures/paperFiguresForRealForRealThisTime/suppFigures/twoEffectVarSFigure.png',
      height = 17.75 ,
      width = 17.75,
      units = 'cm',
      res = 800
    )
    layout(matrix(1:9,nrow=3))
    for(k in 1:nrow(my.grid)) {
      rf <- FALSE
      ## fixed h2
      prev <- list()
      h2 = 1 / 2
      tmp.soln <- list()
      soln <- list()
      
      for (j in seq_along(L)) {
        my.gs[, j] <-
          sort(unique(c(
            1 - max.pl / 2, seq (1 - 1 / L[j] , 1 - max.pl , length.out = n.pts)
          )),decreasing=TRUE)  
        tmp.soln[[j]] <- list()
        soln[[j]] <- list()
        
        ## get initial dl
        {
          this.bt <- my.grid[k,1]
          this.al <- my.grid[k,2]
          init.ft <- (4*Ne*C*as)^(-1) * log ((1+this.bt)/(1-this.bt))
          init.Vas <- 2*L[j]*my.gs[1, j]*u*as*this.bt / (C*init.ft)
          init.Vtot <- init.Vas/h2
          init.tstar <- dnorminv(init.ft,s=sqrt(init.Vtot))
          init.Ft <- pnorm(init.tstar,sd=sqrt(init.Vtot),lower.tail=FALSE)
          init.Fta <- pnorm(init.tstar - this.al,sd=sqrt(init.Vtot),lower.tail=FALSE)
          init.dl <- init.Fta - init.Ft
        }
        
        ## source('scripts/newSolveTwoEffect.R')
        for (i in 1:nrow(my.gs)) {
          print(paste(k,j,i,sep=';'))
          tmp.soln[[j]][[i]] <- solveTwoEffect1D(
            bt = my.grid[k,1],
            Ne = 20000,
            as = as,
            al = my.grid[k,2],
            L = L[j],
            gs = my.gs[i, j],
            h2 = h2 ,
            Ve = NULL ,
            u = u,
            C = C,
            Bval = 1,
            init.dl = ifelse(i==1,init.dl,prev.dl),
            init.tstar = ifelse(i==1,init.tstar,prev.tstar)
          )
          if(is.na(tmp.soln[[j]][[i]]['dl'])) break
          prev.dl <- tmp.soln[[j]][[i]]['dl']
          prev.tstar <- tmp.soln[[j]][[i]]['tstar']
        }
        soln[[j]] <- do.call(rbind,tmp.soln[[j]])
        soln[[j]] <- soln[[j]][!is.na(soln[[j]][,1]),]
      }
      dl <- lapply(soln,function(X)X[,'dl'])
      tstar <- lapply(soln,function(X)X[,'tstar'])
      lambda <- lapply(soln,function(X)X[,'lambda'])
      prev <- lapply(soln,function(X)X[,'prev'])
      Vt <- lapply(soln,function(X)X[,'Vt'])
      
      my.cols <- wes_palette('Darjeeling2')
      if(k == 1){
        marg <- 1.6
        op3 <- par(mar = c(3.1, 2.6 + marg, 0.4 + marg, 0.6),
                   mgp = c(3, 0.6, 0))
      } else if(k %in% 2:3){
        op3 <- par(mar = c(3.1, 2.6 + marg, 0.4, 0.6),
                   mgp = c(3, 0.6, 0))
      } else if(k %in% c(4,7)){
        op3 <- par(mar = c(3.1, 2.6, 0.4 + marg, 0.6),
                   mgp = c(3, 0.6, 0))
      }
      # y.max <- max(unlist(prev))
      
      keep.Vts <- list()
      for (j in seq_along(L)) {
        plot.these <- which(prev[[j]] > 0 & lambda[[j]] < 3 )
        keep.Vts[[j]] <- Vt[[j]][plot.these]
      }
      log10.min.vt <- floor(log10(min(unlist(keep.Vts))))
      log10.max.vt <- ceiling(log10(max(unlist(keep.Vts))))
      min.vt <- 10^log10.min.vt
      max.vt <- 10^log10.max.vt
      my.tick.pos <- 10^(log10.min.vt:log10.max.vt)
      plot(
        NA ,
        type = 'l',
        bty = 'n' ,
        xlim = c(0, max.pl) ,
        ylim = c(min.vt, max.vt) ,
        log = 'y' ,
        yaxt = 'n' ,
        xaxt = 'n',
        xlab = '',
        ylab = ''
      )
      axis(2,
           at = my.tick.pos)
      axis(1,
           at = c(0, 0.004, 0.008),
           labels = c('0', '0.004', '0.008'))
      mtext(
        side = 2 ,
        text = 'Total Liability Variance' ,
        line = 1.5,
        cex = 0.7,
      )
      mtext(
        side = 1 ,
        text = expression(paste(
          'Fraction of sites with large effects, ' , p[L] , sep = ''
        )) ,
        line = 1.5,
        cex = 0.7,
      )
      if(k == 1){
        mtext(
          side = 2 ,
          text = expression(
            paste(b[T], ' = 0.1', sep ='')) ,
          line = 2.5 ,
          cex = 0.7
        )
        mtext(
          side = 3 ,
          text = expression(
            paste(a[L]/a[S], ' = 50', sep ='')) ,
          line = 0.2 ,
          cex = 0.7
        )
      }
      if(k == 2){
        mtext(
          side = 2 ,
          text = expression(
            paste(b[T], ' = 0.5', sep ='')) ,
          line = 2.5 ,
          cex = 0.7
        )
      }
      if(k == 3){
        mtext(
          side = 2 ,
          text = expression(
            paste(b[T], ' = 0.9', sep ='')) ,
          line = 2.5 ,
          cex = 0.7
        )
      }
      if(k == 4){
        mtext(
          side = 3 ,
          text = expression(
            paste(a[L]/a[S], ' = 100', sep ='')) ,
          line = 0.2 ,
          cex = 0.7
        )
        legend(
          'topright' ,
          lty = 1,
          col = my.cols ,
          legend = c(
            expression(1.5 %*% 10 ^ 6) ,
            expression(5 %*% 10 ^ 6) ,
            expression(1.5 %*% 10 ^ 7) ,
            expression(5 %*% 10 ^ 7) ,
            expression(1.5 %*% 10 ^ 8)
          ) ,
          bty = 'n' ,
          title = 'Target size, L'
        )
      }
      if(k == 7){
        mtext(
          side = 3 ,
          text = expression(
            paste(a[L]/a[S], ' = 200', sep ='')) ,
          line = 0.2 ,
          cex = 0.7
        )
      }
      for (j in 1:length(prev)) {
        plot.these <- which(prev[[j]] > 0 & lambda[[j]] < 3 )
        lines(1 - my.gs[plot.these,j] ,
              Vt[[j]][plot.these] ,
              col = my.cols[j] ,
              lty = 1)
      }
    }
    dev.off()
  }
  
  #### small effect scaled coefficient ####
  {
    png(
      'figures/paperFiguresForRealForRealThisTime/suppFigures/twoEffectGammaSSFigure.png',
      height = 17.75 ,
      width = 17.75,
      units = 'cm',
      res = 800
    )
    layout(matrix(1:9,nrow=3))
    
    for(k in 1:nrow(my.grid)) {
      rf <- FALSE
      ## fixed h2
      h2 = 1 / 2
      
      my.gs <- sort(seq (1 - 1e-7 , 1 - max.pl , length.out = n.pts),decreasing=TRUE)  
      this.bt <- my.grid[k,1]
      this.al <- my.grid[k,2]
      amean <- as*my.gs + this.al*(1-my.gs)
      bs <- (this.bt*amean - (1-my.gs)*this.al) / (my.gs*as)
      keep.these <- bs > 0
      bs <- bs[bs > 0]
      gammaS <- 0.5*log ((1+bs)/(1-bs))
      
      
      my.cols <- wes_palette('Darjeeling2')
      if(k == 1){
        marg <- 1.6
        op3 <- par(mar = c(3.1, 2.6 + marg, 0.4 + marg, 0.6),
                   mgp = c(3, 0.6, 0))
      } else if(k %in% 2:3){
        op3 <- par(mar = c(3.1, 2.6 + marg, 0.4, 0.6),
                   mgp = c(3, 0.6, 0))
      } else if(k %in% c(4,7)){
        op3 <- par(mar = c(3.1, 2.6, 0.4 + marg, 0.6),
                   mgp = c(3, 0.6, 0))
      }
      # y.max <- max(unlist(prev))
      
      min.log10.gamma <- floor(log10(min(gammaS)))
      max.log10.gamma <- ceiling(log10(max(gammaS)))
      min.gamma <- 10^min.log10.gamma
      max.gamma <- 10^max.log10.gamma
      my.tick.pos <- 10^(min.log10.gamma:max.log10.gamma)
      plot(
        NA ,
        type = 'l',
        bty = 'n' ,
        xlim = c(0, max.pl) ,
        ylim = c(min.gamma, max.gamma) ,
        log = 'y' ,
        yaxt = 'n' ,
        xaxt = 'n',
        xlab = '',
        ylab = ''
      )
      axis(2,
           at = my.tick.pos)
      axis(1,
           at = c(0, 0.004, 0.008),
           labels = c('0', '0.004', '0.008'))
      mtext(
        side = 2 ,
        text = expression(paste('Small Effect ' , gamma[S] , sep ='')),
        line = 1.5,
        cex = 0.7,
      )
      mtext(
        side = 1 ,
        text = expression(paste(
          'Fraction of sites with large effects, ' , p[L] , sep = ''
        )) ,
        line = 1.5,
        cex = 0.7,
      )
      if(k == 1){
        mtext(
          side = 2 ,
          text = expression(
            paste(b[T], ' = 0.1', sep ='')) ,
          line = 2.5 ,
          cex = 0.7
        )
        mtext(
          side = 3 ,
          text = expression(
            paste(a[L]/a[S], ' = 50', sep ='')) ,
          line = 0.2 ,
          cex = 0.7
        )
      }
      if(k == 2){
        mtext(
          side = 2 ,
          text = expression(
            paste(b[T], ' = 0.5', sep ='')) ,
          line = 2.5 ,
          cex = 0.7
        )
      }
      if(k == 3){
        mtext(
          side = 2 ,
          text = expression(
            paste(b[T], ' = 0.9', sep ='')) ,
          line = 2.5 ,
          cex = 0.7
        )
      }
      if(k == 4){
        mtext(
          side = 3 ,
          text = expression(
            paste(a[L]/a[S], ' = 100', sep ='')) ,
          line = 0.2 ,
          cex = 0.7
        )
      }
      if(k == 7){
        mtext(
          side = 3 ,
          text = expression(
            paste(a[L]/a[S], ' = 200', sep ='')) ,
          line = 0.2 ,
          cex = 0.7
        )
      }
      
      ##plot.these <- which(prev[[j]] > 0 & lambda[[j]] < 3 )
      lines(1 - my.gs[keep.these] ,
            gammaS ,
            lty = 1)
    }
    dev.off()
  }
  
  
  #### skew inflation factor ####
  {
    png(
      'figures/paperFiguresForRealForRealThisTime/suppFigures/twoEffectVarSFigure.png',
      height = 17.75 ,
      width = 17.75,
      units = 'cm',
      res = 800
    )
    layout(matrix(1:9,nrow=3))
    for(k in 1:nrow(my.grid)) {
      rf <- FALSE
      ## fixed h2
      prev <- list()
      h2 = 1 / 2
      tmp.soln <- list()
      soln <- list()
      
      for (j in seq_along(L)) {
        my.gs[, j] <-
          sort(unique(c(
            1 - max.pl / 2, seq (1 - 1 / L[j] , 1 - max.pl , length.out = n.pts)
          )),decreasing=TRUE)  
        tmp.soln[[j]] <- list()
        soln[[j]] <- list()
        
        ## get initial dl
        {
          this.bt <- my.grid[k,1]
          this.al <- my.grid[k,2]
          init.ft <- (4*Ne*C*as)^(-1) * log ((1+this.bt)/(1-this.bt))
          init.Vas <- 2*L[j]*my.gs[1, j]*u*as*this.bt / (C*init.ft)
          init.Vtot <- init.Vas/h2
          init.tstar <- dnorminv(init.ft,s=sqrt(init.Vtot))
          init.Ft <- pnorm(init.tstar,sd=sqrt(init.Vtot),lower.tail=FALSE)
          init.Fta <- pnorm(init.tstar - this.al,sd=sqrt(init.Vtot),lower.tail=FALSE)
          init.dl <- init.Fta - init.Ft
        }
        
        ## source('scripts/newSolveTwoEffect.R')
        for (i in 1:nrow(my.gs)) {
          print(paste(k,j,i,sep=';'))
          tmp.soln[[j]][[i]] <- solveTwoEffect1D(
            bt = my.grid[k,1],
            Ne = 20000,
            as = as,
            al = my.grid[k,2],
            L = L[j],
            gs = my.gs[i, j],
            h2 = h2 ,
            Ve = NULL ,
            u = u,
            C = C,
            Bval = 1,
            init.dl = ifelse(i==1,init.dl,prev.dl),
            init.tstar = ifelse(i==1,init.tstar,prev.tstar)
          )
          if(is.na(tmp.soln[[j]][[i]]['dl'])) break
          prev.dl <- tmp.soln[[j]][[i]]['dl']
          prev.tstar <- tmp.soln[[j]][[i]]['tstar']
        }
        soln[[j]] <- do.call(rbind,tmp.soln[[j]])
        soln[[j]] <- soln[[j]][!is.na(soln[[j]][,1]),]
      }
      dl <- lapply(soln,function(X)X[,'dl'])
      tstar <- lapply(soln,function(X)X[,'tstar'])
      lambda <- lapply(soln,function(X)X[,'lambda'])
      prev <- lapply(soln,function(X)X[,'prev'])
      Vt <- lapply(soln,function(X)X[,'Vt'])
      norm.prev <- lapply(soln,function(X)X[,'norm.prev'])
      
      my.cols <- wes_palette('Darjeeling2')
      if(k == 1){
        marg <- 1.6
        op3 <- par(mar = c(3.1, 2.6 + marg, 0.4 + marg, 0.6),
                   mgp = c(3, 0.6, 0))
      } else if(k %in% 2:3){
        op3 <- par(mar = c(3.1, 2.6 + marg, 0.4, 0.6),
                   mgp = c(3, 0.6, 0))
      } else if(k %in% c(4,7)){
        op3 <- par(mar = c(3.1, 2.6, 0.4 + marg, 0.6),
                   mgp = c(3, 0.6, 0))
      }
      # y.max <- max(unlist(prev))
      
      keep.Vts <- list()
      for (j in seq_along(L)) {
        plot.these <- which(prev[[j]] > 0 & lambda[[j]] < 3 )
        keep.Vts[[j]] <- Vt[[j]][plot.these]
      }
      log10.min.vt <- floor(log10(min(unlist(keep.Vts))))
      log10.max.vt <- ceiling(log10(max(unlist(keep.Vts))))
      min.vt <- 10^log10.min.vt
      max.vt <- 10^log10.max.vt
      my.tick.pos <- 10^(log10.min.vt:log10.max.vt)
      plot(
        NA ,
        type = 'l',
        bty = 'n' ,
        xlim = c(0, max.pl) ,
        ylim = c(min.vt, max.vt) ,
        log = 'y' ,
        yaxt = 'n' ,
        xaxt = 'n',
        xlab = '',
        ylab = ''
      )
      axis(2,
           at = my.tick.pos)
      axis(1,
           at = c(0, 0.004, 0.008),
           labels = c('0', '0.004', '0.008'))
      mtext(
        side = 2 ,
        text = 'Total Liability Variance' ,
        line = 1.5,
        cex = 0.7,
      )
      mtext(
        side = 1 ,
        text = expression(paste(
          'Fraction of sites with large effects, ' , p[L] , sep = ''
        )) ,
        line = 1.5,
        cex = 0.7,
      )
      if(k == 1){
        mtext(
          side = 2 ,
          text = expression(
            paste(b[T], ' = 0.1', sep ='')) ,
          line = 2.5 ,
          cex = 0.7
        )
        mtext(
          side = 3 ,
          text = expression(
            paste(a[L]/a[S], ' = 50', sep ='')) ,
          line = 0.2 ,
          cex = 0.7
        )
      }
      if(k == 2){
        mtext(
          side = 2 ,
          text = expression(
            paste(b[T], ' = 0.5', sep ='')) ,
          line = 2.5 ,
          cex = 0.7
        )
      }
      if(k == 3){
        mtext(
          side = 2 ,
          text = expression(
            paste(b[T], ' = 0.9', sep ='')) ,
          line = 2.5 ,
          cex = 0.7
        )
      }
      if(k == 4){
        mtext(
          side = 3 ,
          text = expression(
            paste(a[L]/a[S], ' = 100', sep ='')) ,
          line = 0.2 ,
          cex = 0.7
        )
        legend(
          'topright' ,
          lty = 1,
          col = my.cols ,
          legend = c(
            expression(1.5 %*% 10 ^ 6) ,
            expression(5 %*% 10 ^ 6) ,
            expression(1.5 %*% 10 ^ 7) ,
            expression(5 %*% 10 ^ 7) ,
            expression(1.5 %*% 10 ^ 8)
          ) ,
          bty = 'n' ,
          title = 'Target size, L'
        )
      }
      if(k == 7){
        mtext(
          side = 3 ,
          text = expression(
            paste(a[L]/a[S], ' = 200', sep ='')) ,
          line = 0.2 ,
          cex = 0.7
        )
      }
      for (j in 1:length(prev)) {
        plot.these <- which(prev[[j]] > 0 & lambda[[j]] < 3 )
        lines(1 - my.gs[plot.these,j] ,
              Vt[[j]][plot.these] ,
              col = my.cols[j] ,
              lty = 1)
      }
    }
    dev.off()
  }
  
  
  
  
}

#### Figure S large effect sims ####
#### probably not a real S figure as I am currently working it into 5B (4/15/25)
{
  my.list <- get(load('output/twoEffectSuppFigureSolns.Robj'))
  gs.theory <- my.list[[1]]
  prev.theory <- my.list[[2]]
  my.sims <-
    get(load('resultsFiles/twoEffectPrevResultsTableN1000_old.Rdata'))
  my.sims$L <- round(my.sims[, 'Ls'] + my.sims[, 'Ll'], -1)
  my.sims$gs <- my.sims[, 'Ls'] / my.sims$L
  split.sims <- split(my.sims, my.sims$L)
  
  
  max.pl <- 0.008
  png(
    'figures/paperFiguresForRealForRealThisTime/suppFigures/twoEffectPrevSimsSFigure.png',
    height = 10.75 ,
    width = 10.75,
    units = 'cm',
    res = 600
  )
  op3 <- par(mar = c(3.1, 2.6, 0.4, 0.4),
             mgp = c(3, 0.6, 0))
  plot(
    NA ,
    type = 'l',
    bty = 'n' ,
    xlim = c(0, max.pl) ,
    ylim = c(1e-3, 0.1) ,
    log = 'y' ,
    yaxt = 'n' ,
    xaxt = 'n',
    xlab = '',
    ylab = ''
  )
  axis(2,
       at = c(0.001, 0.01, 0.1),
       labels = c('0.001', '0.01', '0.1'))
  axis(1,
       at = c(0, 0.004, 0.008),
       labels = c('0', '0.004', '0.008'))
  mtext(
    side = 2 ,
    text = 'Prevalence' ,
    line = 1.7,
    cex = 1,
  )
  mtext(
    side = 1 ,
    text = expression(paste(
      'Fraction of sites with large effects, ' , p[L] , sep = ''
    )) ,
    line = 2,
    cex = ,
  )
  my.cols <- wes_palette('Darjeeling2')
  for (j in 1:length(prev.theory)) {
    lines(1 - gs.theory[[j]] ,
          prev.theory[[j]] ,
          col = my.cols[j] ,
          lty = 1)
    points(1 - split.sims[[j]]$gs ,
           split.sims[[j]]$sim.prev ,
           col = my.cols[j],
           pch = 4)
  }
  legend(
    'topright' ,
    lty = 1,
    col = my.cols ,
    legend = c(
      expression(1.5 %*% 10 ^ 6) ,
      expression(5 %*% 10 ^ 6) ,
      expression(1.5 %*% 10 ^ 7) ,
      expression(5 %*% 10 ^ 7) ,
      expression(1.5 %*% 10 ^ 8)
    ) ,
    bty = 'n' ,
    title = 'Target size, L'
  )
  legend(
    'topleft' ,
    lty = c(NA,1),
    pch = c(4,NA),
    legend = c(
      'Simulations' ,
      'Numerical Solution'
    ) ,
    bty = 'n'
  )
  dev.off()
}
