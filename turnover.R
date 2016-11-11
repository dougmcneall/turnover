# Analysis of regional carbon-sink-to-source timing
#
# D.McNeall 27th April 2011

# Take timeseries of regional carbon sinks
# smooth
# Find point where they turn from gaining to losing carbon

# --------------------------------------------------------
# 0. source scripts, packages, functions
#
# --------------------------------------------------------

print.pdf <- TRUE

source('/home/h01/hadda/code/R/useful/dougpalettes.R')
source('/home/h01/hadda/code/R/useful/dougplotpars.R')

source('/home/h01/hadda/avoid/kitchenSink/RCP26analysis_work/ensembleFunctions.R')

RCP.col = 'tomato2'
A1B.col = 'skyblue3'

# legend vectors
colvec = c(A1B.col, RCP.col)
legvec = c('A1B','RCP2.6')


library(fields)
library(RColorBrewer)
byr <- brewer.pal(11, "RdYlBu")[11:1]
ryb <- brewer.pal(11, "RdYlBu")
yor <-  brewer.pal(11, "YlOrRd")

ryg <- brewer.pal(11, "RdYlGn")
g <-brewer.pal(9, "Greens")

rb <- brewer.pal(11, "RdBu")
br <- brewer.pal(11, "RdBu")[11:1]

text.colour <- 'black'

years <- 1860:2099

# Carbon budget obsevations/modelling studies from Tyndall, Le Quere etc.
#carbonBudget <- read.table('data/global_carbon_budget.txt', header = TRUE) #old data
carbonBudget <- read.table('data/global_carbon_budget2013.txt',skip = 19, header = TRUE)

# carbon emissions in RCP26 and A1B
emissions <- read.table(file = 'CO2emissGTc.txt', header = TRUE)



carbonChange <- function(soilfile,vegfile,rlfile,kscale = 1,f = 0.4, ...){
  # function which loads and processes soil and veg carbon data,
  # finds turnover, stores everything as an R object

  
  # soil carbon
  soilc <-  processTS(soilfile,
                      rlfile = rlfile,
                      year = 1860:2099,
                      exclude = 43
                      ) / kscale
  # veg carbon
  vegc <-  processTS(vegfile,
                     rlfile = rlfile,
                     year = 1860:2099,
                     exclude = 43
                     ) / kscale
  # total carbon
  totc <- soilc+vegc

  # anomalies of carbon
  soilcAnom <- anomalizeTS(soilc, 1:30)
  vegcAnom  <- anomalizeTS(vegc, 1:30)
  totcAnom  <- anomalizeTS(totc, 1:30)

  # lowess smoothed carbon
  soilcAnomSmooth <- loessSmoothdf(soilcAnom, f = f)
  vegcAnomSmooth  <- loessSmoothdf(vegcAnom, f = f)
  totcAnomSmooth  <- loessSmoothdf(totcAnom, f = f)


  # define the years that we want to look in
  years.i <- 120:240
  
  # find turnover year in soil, veg, and total carbon 
  soilcTurnoverYear <-  findTurnover(soilcAnomSmooth[ ,years.i]) 
  vegcTurnoverYear  <-  findTurnover(vegcAnomSmooth[ ,years.i])
  totcTurnoverYear  <-  findTurnover(totcAnomSmooth[ ,years.i])
  
  # find the year indices
  soilc.ix <- match(soilcTurnoverYear, colnames(soilcAnomSmooth))
  vegc.ix  <- match(vegcTurnoverYear, colnames(vegcAnomSmooth))
  totc.ix  <- match(totcTurnoverYear, colnames(totcAnomSmooth))

 
  # how many don't turn over?

  soilcNoTurn <- sum(is.na(soilc.ix))
  vegcNoTurn  <- sum(is.na(vegc.ix))
  totcNoTurn  <- sum(is.na(totc.ix))
    

  return(list(soilc = soilc, vegc = vegc, totc = totc,
              soilcAnom = soilcAnom, vegcAnom = vegcAnom, totcAnom = totcAnom,
              soilcAnomSmooth = soilcAnomSmooth, vegcAnomSmooth = vegcAnomSmooth,
              totcAnomSmooth = totcAnomSmooth,
              soilcTurnoverYear = soilcTurnoverYear,
              vegcTurnoverYear = vegcTurnoverYear,
              totcTurnoverYear = totcTurnoverYear,
              soilc.ix = soilc.ix, vegc.ix = vegc.ix,
              totc.ix  = totc.ix,
              soilcNoTurn = soilcNoTurn,vegcNoTurn = vegcNoTurn,
              totcNoTurn = totcNoTurn))
  
}


# Global (per m^2)
#GlobalcA1B   <- carbonChange ('A1BsoilcGlobal.txt', 'A1BvegcGlobal.txt', 'RCP26runlist.txt')
#GlobalcRCP26 <- carbonChange ('RCP26soilcGlobal.txt', 'RCP26vegcGlobal.txt', 'RCP26runlist.txt')

 
# Global Total Carbon
GlobalcA1B   <- carbonChange('A1BsoilcTotalGlobal.txt',
                             'A1BvegcTotalGlobal.txt',
                             'RCP26runlist.txt',
                             kscale = 1e12
                             )

GlobalcRCP26 <- carbonChange('RCP26soilcTotalGlobal.txt',
                             'RCP26vegcTotalGlobal.txt',
                             'RCP26runlist.txt',
                             kscale = 1e12
                             )

# Tropics total carbon
TropicscA1B <- carbonChange(soilfile = 'soilcTropicsA1B.txt',
                          vegfile = 'vegcTropicsA1B.txt',
                          rlfile = 'RCP26runlist.txt',
                          kscale = 1e12
                          )
TropicscRCP26 <- carbonChange(soilfile = 'soilcTropicsRCP26.txt',
                          vegfile = 'vegcTropicsRCP26.txt',
                          rlfile = 'RCP26runlist.txt',
                          kscale = 1e12
                          )

# Boreal(Taiga) total carbon
TaigacA1B <- carbonChange(soilfile = 'soilcTaigaA1B.txt',
                          vegfile = 'vegcTaigaA1B.txt',
                          rlfile = 'RCP26runlist.txt',
                          kscale = 1e12
                          )
TaigacRCP26 <- carbonChange(soilfile = 'soilcTaigaRCP26.txt',
                          vegfile = 'vegcTaigaRCP26.txt',
                          rlfile = 'RCP26runlist.txt',
                          kscale = 1e12
                          )


carbPlot1 <- function(carb){
  # plot up a single 

  par(dougpar_web)
  
  matplot(colnames(carb$totcAnom), t(carb$totcAnom),
        type = 'l',
        col = 'grey',
        ylim = c(-9,7),
        xlab = '',
        ylab = expression(paste('carbon store change from 1860-90 (kg m'^-2,')')) 
        )
  
  matlines(colnames(carb$totcAnom),t(carb$totcAnomSmooth),
           lty = 1,
           col = A1B.col)
  
  points( as.numeric(colnames(carb$totcAnom)[carb$totc.ix]),carb$totcAnomSmooth[cbind(1:57,carb$totc.ix)],
         pch = 20,
         col = 'dodgerblue3')
  
  legend('topleft', legvec,
         col = c(colvec),
         text.col = colvec,
         lwd = 1.5,
         pch = 20,
         lty = 1
         )
  
  mtext(side = 3, line = 1, adj = 0,'Amazon Sink to Source Timing', cex = 2, col = 'black')
  
  
}




carbPlotAll <- function(carbA1B, carbRCP26, ctype = 'tot', ylim = c(-700,700), maintext ='Sink to Source Timing', leg = TRUE, xax = FALSE, yax = FALSE){

  # Plot both scenarios on a single diagram

  par(dougpar_web)
  par(mgp = c(3,1,0))
  par(mar = c(1,1,1,1))
  
  matplot(colnames(carbA1B[[paste(ctype,'c', sep = '')]]),t(carbA1B[[paste(ctype,'cAnomSmooth', sep = '')]]),
          axes = FALSE,
          ylim = ylim,
          type = 'l',
          lty = 1,
          col = A1B.col,
          xlab = '',
          ylab = '',
          pty = 'n',
          xaxs = 'i', yaxs = 'i'
          )
  
  polygon(c(1860, 2100,20100, 1860), c(-700, -700, 700, 700), col = 'grey90', border = NA)
  abline(h = seq(from = -600, to = 600, by = 200), col = 'white')
  abline(v = seq(from = 1900, to = 2050, by = 50), col = 'white')
  
  
  matlines(colnames(carbA1B[[paste(ctype,'c', sep = '')]]),t(carbA1B[[paste(ctype,'cAnomSmooth', sep = '')]]),
          axes = FALSE,
          ylim = ylim,
          type = 'l',
          lty = 1,
          col = A1B.col,
          xlab = '',
          #ylab = expression(paste('carbon store change from 1860-90 (kg m'^-2,')'))
          ylab = 'GtC' 
        )
  
  matlines( colnames(carbRCP26[[paste(ctype,'c', sep = '')]]),t(carbRCP26[[paste(ctype,'cAnomSmooth', sep = '')]]),
           lty = 1,
           col = RCP.col
           )
  
  if(xax){axis(1)}
  if(yax){axis(2)}
  
  
  points(as.numeric(colnames(carbA1B[[paste(ctype,'c', sep = '')]])) [carbA1B[[paste(ctype,'c.ix',sep = '')]]],
         carbA1B[[paste(ctype,'cAnomSmooth', sep = '')]][cbind(1:57,carbA1B[[paste(ctype,'c.ix',sep = '')]])],
         pch = 21,
         col = 'dodgerblue3', bg = 'skyblue')

    points(as.numeric(colnames(carbRCP26[[paste(ctype,'c', sep = '')]])) [carbRCP26[[paste(ctype,'c.ix',sep = '')]]],
         carbRCP26[[paste(ctype,'cAnomSmooth', sep = '')]][cbind(1:57,carbRCP26[[paste(ctype,'c.ix',sep = '')]])],
         pch = 21,
         col = 'firebrick3', bg = 'tomato')
  

  rug(as.numeric(colnames(carbA1B[[paste(ctype,'c', sep = '')]])) [carbA1B[[paste(ctype,'c.ix',sep = '')]]]  , col = 'dodgerblue3', line = -2, lwd =1)
  rug(as.numeric(colnames(carbRCP26[[paste(ctype,'c', sep = '')]])) [carbRCP26[[paste(ctype,'c.ix',sep = '')]]], col = 'firebrick3', line = -1, lwd = 1)

  if(leg){
    legend('bottomleft', legvec,
           col = c(colvec),
           text.col = colvec,
           lwd = 1.5,
           pch = 21,
           pt.bg = c('skyblue', 'tomato'),
           lty = 1
           )
  }
  mtext(text = maintext, side = 3, line = -1.2, adj = 0, cex = 1, col = 'black')
  
}


# ----------------------------------------------------------------------
# Straight-up global land carbon stores through the runs
# ----------------------------------------------------------------------

dev.new(width = 7, height = 7)
par(mfrow = c(1,2), las = 1, mar = c(5,5,3,2), mgp = c(3.5,1,0), fg = 'darkgrey')
matplot(years, t(GlobalcA1B$totc),
        type = 'l',
        col = A1B.col,
        lty = 'solid',
        ylim = c(0,3000),
        xlab = 'Year',
        ylab = 'Land carbon store (GtC)',
        bty = 'l'
        )

matlines(years, t(GlobalcRCP26$totc),col = RCP.col, lty = 'solid')
mtext('a)',adj = 0, cex = 2, line = 0.7, col = text.colour)

matplot(years, t(GlobalcA1B$totcAnom),
        type = 'l',
        col = A1B.col,
        lty = 'solid',
        ylim = c(-650,650),
        xlab = 'Year',
        ylab = 'Land carbon store anomaly (GtC)',
        bty = 'l')

matlines(years, t(GlobalcRCP26$totcAnom),col = RCP.col, lty = 'solid')
legend('topleft',c('A1B', 'RCP2.6'), col = colvec, text.col = colvec,lty = 'solid', bty = 'n')

mtext('b)',adj = 0, cex = 2, line = 0.7, col = text.colour)

if(print.pdf){dev.print(dev = pdf, file = 'graphics/turnoverGlobalc.pdf', width = 7, height = 7)}

# ------------------------------------------------------------------------
# Global Vegetation and soil carbon store changes
# ------------------------------------------------------------------------
dev.new(width = 7, height = 7)
par(mfrow = c(1,2), las = 1, mar = c(5,5,3,2), mgp = c(3.5,1,0), fg = 'grey')

matplot(years, t(GlobalcA1B$vegcAnom),
        type = 'l',
        col = A1B.col,
        lty = 'solid',
        ylim = c(-600, 400),
        xlab = 'Year',
        ylab = 'Vegetation carbon store (GtC)',
        bty = 'l'
        )

matlines(years, t(GlobalcRCP26$vegcAnom),col = RCP.col, lty = 'solid')
mtext('a)',adj = 0, cex = 2, line = 0.7, col = text.colour)

matplot(years, t(GlobalcA1B$soilcAnom),
        type = 'l',
        col = A1B.col,
        lty = 'solid',
        ylim = c(-600,400),
        xlab = 'Year',
        ylab = 'Soil carbon store anomaly (GtC)',
        bty = 'l'
        )

matlines(years, t(GlobalcRCP26$soilcAnom),col = RCP.col, lty = 'solid')
legend('topleft',c('A1B', 'RCP2.6'), col = colvec, text.col = colvec,lty = 'solid', bty = 'n')

mtext('b)',adj = 0, cex = 2, line = 0.7, col = text.colour)

if(print.pdf){dev.print(dev = pdf, file = 'graphics/turnoverGlobalcvegsoil.pdf', width = 7, height = 7)}




# ----------------------------------------------------------------------
# Visualise all of the global and regional carbon store changes
# ----------------------------------------------------------------------
dev.new(width = 7, height = 10)
par(mfrow = c(3,3))
par(oma = c(1,4,1,0))
carbPlotAll(GlobalcA1B, GlobalcRCP26,
            ctype = 'tot', maintext = 'A. Global Total',ylim = c(-700,700), yax = TRUE)

carbPlotAll(GlobalcA1B, GlobalcRCP26,
            ctype = 'veg', maintext = 'B. Global Vegetation',ylim = c(-700,700), leg = FALSE)

carbPlotAll(GlobalcA1B, GlobalcRCP26,
            ctype = 'soil', maintext = 'C. Global Soil',ylim = c(-700,700),leg = FALSE)


carbPlotAll(TaigacA1B, TaigacRCP26,
            ctype = 'tot', maintext = 'D. Boreal Total',c(-700,700),leg = FALSE, yax = TRUE)

carbPlotAll(TaigacA1B, TaigacRCP26,
            ctype = 'veg', maintext = 'E. Boreal Vegetation',c(-700,700), leg = FALSE)

carbPlotAll(TaigacA1B, TaigacRCP26,
            ctype = 'soil', maintext = 'F. Boreal Soil',c(-700,700), leg = FALSE)

carbPlotAll(TropicscA1B, TropicscRCP26,
            ctype = 'tot', maintext = 'G. Tropics Total',c(-700,700), leg = FALSE, yax = TRUE, xax = TRUE)

carbPlotAll(TropicscA1B, TropicscRCP26,
            ctype = 'veg', maintext = 'H. Tropics Vegetation',c(-700,700), leg = FALSE, xax = TRUE)

carbPlotAll(TropicscA1B, TropicscRCP26,
            ctype = 'soil', maintext = 'I. Tropics Soil',c(-700,700), leg = FALSE,xax = TRUE)

if(print.pdf){dev.print(dev = pdf, file = 'graphics/turnoverAll.pdf', width = 7, height = 10)}

# --------------------------------------------------------------------
# paired differences in global carbon uptake (timeseries)
# --------------------------------------------------------------------

GlobalcDiffs <-GlobalcRCP26$totc -  GlobalcA1B$totc

dev.new(width = 7, height = 7)
matplot(years, t(GlobalcDiffs), type = 'l', lty = 'solid', lwd = 2,
        xlim = c(1950,2100),
        ylim = c(-600, 600),
        col = dougpal2,
        ylab = 'Terrestrial Carbon Difference (GtC) RCP2.6 - A1B ',
        )
abline(h = 0, lty = 'dashed')


# --------------------------------------------------------------------
# Timings of turnover
# --------------------------------------------------------------------
dev.new(width = 7, height = 7)
par(las = 1, mar = c(6,6,4,1), mgp = c(4,1,0), fg = 'grey')

plot(GlobalcRCP26$totcTurnoverYear,GlobalcA1B$totcTurnoverYear,
     pch = 21,
     xlab = 'RCP26 turnover year',
     ylab = 'A1B turnover year',
     col = 'black',
     xlim = c(1980,2100),
     ylim = c(1980,2100),
     cex.lab = 1.5,
     bg = 'lightgrey'
     )
abline(0,1, col = 'darkgrey')
 
points(TropicscRCP26$totcTurnoverYear,TropicscA1B$totcTurnoverYear,
       col = 'darkorange',bg = 'orange', pch = 25)

points(TaigacRCP26$totcTurnoverYear,TaigacA1B$totcTurnoverYear,
       col = 'black', pch = 24, bg = 'white')

legend('bottomright',c('Global','Tropics','Boreal'),
       col = c('black', 'darkorange', 'black'), pch = c(21,25,24),pt.bg = c('lightgrey', 'orange', 'white'), text.col = text.colour)

text(1980, 2095,'A1B later',col = A1B.col, pos = 4, cex = 3 )
text(2080, 1980,'RCP2.6 later',col = RCP.col, pos = 2, cex = 3)

if(print.pdf){dev.print(dev = pdf, file = 'graphics/turnoverPairs2.pdf', width = 7, height = 7)}

# -----------------------------------------------------------
# Pairwise comparison of turnover dates 
# -----------------------------------------------------------

# globally, the turnover point is 17 years later under
# A1B
globalDiff.totc <- GlobalcA1B$totcTurnoverYear - GlobalcRCP26$totcTurnoverYear
globalDiff.soilc <- GlobalcA1B$soilcTurnoverYear - GlobalcRCP26$soilcTurnoverYear
globalDiff.vegc <- GlobalcA1B$vegcTurnoverYear - GlobalcRCP26$vegcTurnoverYear


mean(globalDiff.totc, na.rm = TRUE)
sd(globalDiff.totc, na.rm = TRUE)

# How hard is the change in global land surface driven by
# a reduction in emissions?
dev.new()
par(mfrow = c(2,1))
hist(GlobalcA1B$totcTurnoverYear - 2051)
hist(GlobalcRCP26$totcTurnoverYear - 2020)


# histograms of turnover point

h.xlim <- c(1970, 2100)
h.ylim = c(0,28)

breaks <- seq(from = 1960,to = 2100, by = 10)

dev.new() 
par(mfrow = c(3,3))

hist(GlobalcRCP26$soilcTurnoverYear,
     xlim = h.xlim,ylim = h.ylim, col = RCP.col, main = '', breaks = breaks)
hist(GlobalcRCP26$vegcTurnoverYear,
     xlim = h.xlim,ylim = h.ylim, col = RCP.col, main = '', breaks = breaks)
hist(GlobalcRCP26$totcTurnoverYear,
     xlim = h.xlim,ylim = h.ylim, col = RCP.col, main = '', breaks = breaks)

hist(GlobalcA1B$soilcTurnoverYear,
     xlim = h.xlim,ylim = h.ylim, col = A1B.col, main = '', breaks = breaks)
hist(GlobalcA1B$vegcTurnoverYear,
     xlim = h.xlim,ylim = h.ylim, col = A1B.col, main = '', breaks = breaks)
hist(GlobalcA1B$totcTurnoverYear,
     xlim = h.xlim,ylim = h.ylim, col = A1B.col, main = '', breaks = breaks)

hist(globalDiff.soilc, col = 'grey', main = 'Pairwise difference', xlim = c(-20,60))
hist(globalDiff.vegc, col = 'grey', main = 'Pairwise difference', xlim = c(-20,60))
hist(globalDiff.totc, col = 'grey', main = 'Pairwise difference', xlim = c(-20,60))


dev.new()
par(mfrow = c(3,1))
hist(globalDiff.totc, col = 'grey', xlim = c(-20, 60), main = 'Total C Turnover year A1B - RCP2.6')
hist(globalDiff.soilc, col = 'grey',xlim = c(-20, 60), main = 'Soil C Turnover year A1B - RCP2.6')
hist(globalDiff.vegc, col = 'grey',xlim = c(-20, 60), main = 'Vegetation C Turnover year A1B - RCP2.6')


# Dot plots of turnover dates

totc.ix <- order(GlobalcRCP26$totcTurnoverYear, na.last = TRUE)
vegc.ix <- order(GlobalcRCP26$vegcTurnoverYear, na.last = TRUE)
soilc.ix <- order(GlobalcRCP26$soilcTurnoverYear, na.last = TRUE)


# This set ordered by the RCP2.6 ordering

dev.new(width = 9, height = 5)
par(las = 1, mfrow = c(1,3), mar = c(3,2,3,1), cex.axis = 1.5, cex.main = 1.5)
                                        # vegetation carbon
stripchart(as.list(GlobalcRCP26$vegcTurnoverYear[vegc.ix]),
           col = RCP.col,
           pch = 19,
           xlim = c(1970, 2100),
           main = 'Vegetation',
           axes = FALSE
           )
legend('topleft',legend = legvec, col = colvec, text.col = colvec, bty = 'n', pch = 19, cex = 1.3)
abline(h = 57, col = 'grey')


stripchart(as.list(GlobalcA1B$vegcTurnoverYear[vegc.ix]),
           col = A1B.col,
           pch = 19,
           add = TRUE
           )
abline(h = 57, col = 'grey')
axis(1)


# soil carbon
stripchart(as.list(GlobalcRCP26$soilcTurnoverYear[soilc.ix]),
           col = RCP.col,
           pch = 19,
           xlim = c(1970, 2100),
           main = 'Soil',
           axes = FALSE
           )
abline(h = 57, col = 'grey')
axis(1)
stripchart(as.list(GlobalcA1B$soilcTurnoverYear[soilc.ix]),
           col = A1B.col, pch = 19, add = TRUE)

# total carbon
stripchart(as.list(GlobalcRCP26$totcTurnoverYear[totc.ix]),
           col = RCP.col, pch = 19, xlim = c(1970, 2100),
           main = 'Land surface',
           axes = FALSE
           )
abline(h = 57, col = 'grey')
axis(1)
stripchart(as.list(GlobalcA1B$totcTurnoverYear[totc.ix]),
           add = TRUE,
           col = A1B.col,
           pch = 19
           )
if(print.pdf){dev.print(dev = pdf, file = 'graphics/turnoverDotPlotPaired.pdf', width = 9, height = 5)}



# This set ordered by their own, internal ordering
dev.new(width = 9, height = 5)
par(las = 1, mfrow = c(1,3), mar = c(3,2,3,1), cex.axis = 1.5, cex.main = 1.5)

totc.ix.RCP26 <- order(GlobalcRCP26$totcTurnoverYear, na.last = TRUE)
totc.ix.A1B <- order(GlobalcA1B$totcTurnoverYear, na.last = TRUE)

vegc.ix.RCP26 <- order(GlobalcRCP26$vegcTurnoverYear, na.last = TRUE)
vegc.ix.A1B <- order(GlobalcA1B$vegcTurnoverYear, na.last = TRUE)


soilc.ix.RCP26 <- order(GlobalcRCP26$soilcTurnoverYear, na.last = TRUE)
soilc.ix.A1B <- order(GlobalcA1B$soilcTurnoverYear, na.last = TRUE)

# vegetation carbon
stripchart(as.list(GlobalcRCP26$vegcTurnoverYear[vegc.ix.RCP26]),
           col = RCP.col,
           pch = 19,
           xlim = c(1970, 2100),
           main = 'Vegetation',
           axes = FALSE
           )
abline(h = 57, col = 'grey', lty = 'dashed')
legend('left',legend = legvec, col = colvec, text.col = colvec, bty = 'n', pch = 19, cex = 1.3)  

stripchart(as.list(GlobalcA1B$vegcTurnoverYear[vegc.ix.A1B]),
           col = A1B.col,
           pch = 19,
           add = TRUE
           )
axis(1)


# soil carbon
stripchart(as.list(GlobalcRCP26$soilcTurnoverYear[soilc.ix.RCP26]),
           col = RCP.col,
           pch = 19,
           xlim = c(1970, 2100),
           main = 'Soil',
           axes = FALSE
           )
abline(h = 57, col = 'grey', lty = 'dashed')
axis(1)
stripchart(as.list(GlobalcA1B$soilcTurnoverYear[soilc.ix.A1B]),
           col = A1B.col, pch = 19, add = TRUE)

# total carbon
stripchart(as.list(GlobalcRCP26$totcTurnoverYear[totc.ix.RCP26]),
           col = RCP.col, pch = 19, xlim = c(1970, 2100),
           main = 'Land surface',
           axes = FALSE
           )
abline(h = 57, col = 'grey', lty = 'dashed')
axis(1)
stripchart(as.list(GlobalcA1B$totcTurnoverYear[totc.ix.A1B]),
           add = TRUE,
           col = A1B.col,
           pch = 19
           )

if(print.pdf){dev.print(dev = pdf, file = 'graphics/turnoverDotPlotOrdered.pdf', width = 9, height = 5)}


# ------------------------------------------------------------------
# Output some data on turnover dates
# ------------------------------------------------------------------

# How many ensemble members transition from sink to source?

count.turnover <- function(x) {
  
  sum(is.finite(x))

}

counts.vec.A1B <- sapply(list(GlobalcA1B$totcTurnoverYear,
                          GlobalcA1B$soilcTurnoverYear,
                          GlobalcA1B$vegcTurnoverYear,
                          
                          TaigacA1B$totcTurnoverYear,
                          TaigacA1B$soilcTurnoverYear,
                          TaigacA1B$vegcTurnoverYear,
                          
                          TropicscA1B$totcTurnoverYear,
                          TropicscA1B$soilcTurnoverYear,
                          TropicscA1B$vegcTurnoverYear
                          ),
                     count.turnover
                     )

counts.vec.RCP26 <- sapply(list(GlobalcRCP26$totcTurnoverYear,
                          GlobalcRCP26$soilcTurnoverYear,
                          GlobalcRCP26$vegcTurnoverYear,
                          
                          TaigacRCP26$totcTurnoverYear,
                          TaigacRCP26$soilcTurnoverYear,
                          TaigacRCP26$vegcTurnoverYear,
                          
                          TropicscRCP26$totcTurnoverYear,
                          TropicscRCP26$soilcTurnoverYear,
                          TropicscRCP26$vegcTurnoverYear                                      
                          ),
                     count.turnover
                     )

counts.df.A1B <- data.frame(matrix(counts.vec.A1B, ncol = 3, byrow = TRUE), row.names = c('global', 'boreal', 'tropics') )
colnames(counts.df.A1B) <- c('total', 'soil', 'veg')


counts.df.RCP26 <- data.frame(matrix(counts.vec.RCP26, ncol = 3, byrow = TRUE), row.names = c('global', 'boreal', 'tropics') )
colnames(counts.df.RCP26) <- c('total', 'soil', 'veg')

# write out data files
library(MASS)

write.matrix(file = 'data/turnoverCountsA1B.txt', counts.df.A1B)
write.matrix(file = 'data/turnoverCountsRCP26.txt', counts.df.RCP26)

later <- sum(GlobalcA1B$totcTurnoverYear > GlobalcRCP26$totcTurnoverYear, na.rm = TRUE)

same <- sum(GlobalcA1B$totcTurnoverYear == GlobalcRCP26$totcTurnoverYear, na.rm = TRUE)

turnovers <- cbind(GlobalcA1B$totcTurnoverYear, GlobalcRCP26$totcTurnoverYear)

notinscope <- sum(is.finite(turnovers[,2]) & is.na(turnovers[,1])  )

notinscope + later

earlier <-  sum(GlobalcA1B$totcTurnoverYear < GlobalcRCP26$totcTurnoverYear, na.rm = TRUE)


# ------------------------------------------------------------
# Is there a bias/response relationship in turnover dates?
#
# ------------------------------------------------------------

# doesn't look like it from the equilibrium (although there might be
# a relationship if you throw out some which don't respond much)
x11()
plot(GlobalcA1B$totc[ ,1], GlobalcA1B$totcTurnoverYear)
x11()
plot(GlobalcRCP26$totc[ ,1], GlobalcRCP26$totcTurnoverYear)


# There is a strong relationship between trends over the 20th century,
# and the turnover point

years <- 1860:2099
ix.1900 <- match(1900, years)
ix.1999 <- match(1999, years)

diff20C.A1B <- GlobalcA1B$totcAnomSmooth[ ,ix.1999] -  GlobalcA1B$totcAnomSmooth[ ,ix.1900]
diff20C.RCP26 <- GlobalcRCP26$totcAnomSmooth[ ,ix.1999] -  GlobalcRCP26$totcAnomSmooth[ ,ix.1900]

dev.new()
plot(diff20C.A1B,  GlobalcA1B$totcTurnoverYear, col = A1B.col, pch = 19,
     xlab = "Land surface 20C change (GtC)",
     ylab = "Turnover year"
     )
points(diff20C.RCP26,  GlobalcRCP26$totcTurnoverYear, col = RCP.col, pch = 19)




# --------------------------------------------------------------
# A relationship between 20th century change and
# overall uptake at end of 21st century?
# --------------------------------------------------------------


# Is there a relationship between 19th century carbon store, and
# 21st century carbon store?

Cstore1860to90 <- apply(yearSelect(GlobalcA1B$totc,1860:1890),1,mean, na.rm = TRUE)

CstoreAnom2069to99A1B <- apply(yearSelect(GlobalcA1B$totcAnom,2069:2099),1,mean, na.rm = TRUE)
CstoreAnom2069to99RCP26 <- apply(yearSelect(GlobalcRCP26$totcAnom,2069:2099),1,mean, na.rm = TRUE)


# Final change in percentage terms
CstoreAnomPerc2069to99A1B   <- (CstoreAnom2069to99A1B / Cstore1860to90) * 100
CstoreAnomPerc2069to99RCP26 <- (CstoreAnom2069to99RCP26 / Cstore1860to90) * 100

plot(Cstore1860to90,CstoreAnom2069to99A1B, col = A1B.col, pch = 19)
points(Cstore1860to90,CstoreAnom2069to99RCP26, col = RCP.col, pch = 19)

dev.new(width = 7, height = 7)
par(las = 1, fg = 'grey')
plot(Cstore1860to90,CstoreAnomPerc2069to99A1B,
     col = A1B.col,
     pch = 19,
     xlab = 'Land carbon store 1860-90',
     ylab = 'Land carbon store change by 2100 (%)',
     bty = 'l'
     )
points(Cstore1860to90,CstoreAnomPerc2069to99RCP26, col = RCP.col, pch = 19)
abline(h = 0, lty = 'dashed')
legend('topleft', legvec, col = colvec, text.col = colvec, pch = 19)

if(print.pdf){dev.print(dev = pdf, file = 'graphics/LandCarbonStoreVsResponse.pdf', width = 7, height = 7)}


# is this repeated?
# there is a strong, linear relationship between 20th century land
# carbon uptake, and eventual land carbon uptake
dev.new()
par(las = 1)
plot(diff20C.A1B, GlobalcA1B$totcAnomSmooth[,240],
     col = A1B.col,
     pch = 19,
     #bg = 'blue',
     xlab = '20th Century terrestrial carbon change (GtC)',
     ylab = 'Total terrestrial carbon change by 2100 (GtC)'
     )
points(diff20C.RCP26, GlobalcRCP26$totcAnomSmooth[,240],
       col = RCP.col,
       #bg = 'tomato3',
       pch = 19
       )
abline(h = 0, v = 0, lty = 'dashed')
#dev.off()


# tyndall estimates the land carbon sink from 1960 - 2010


# Sum of all fluxes to the 
obs.sink <- sum(carbonBudget$landSink)

ix.1959 <- match(1959, years)
ix.2012 <- match(2012, years)

LandSink59_12.A1B <- GlobalcA1B$totcAnomSmooth[ ,ix.2012] -  GlobalcA1B$totcAnomSmooth[ ,ix.1959]
LandSink59_12.RCP26 <- GlobalcRCP26$totcAnomSmooth[ ,ix.2012] -  GlobalcRCP26$totcAnomSmooth[ ,ix.1959]


# uncertainty on the global sum of the land sink is approx:

# 1 sd = ~ 5% of fossil fuel
ffc.sd <- 0.05 * carbonBudget$fossilFuelCement
luc.sd <- rep(0.7, length(ffc.sd))
ag.sd <- rep(0.25, length(ffc.sd))
os.sd <- rep(0.5, length(ffc.sd))

# standard deviation of the uncertainty through time
unc.tot.mat <- cbind(ffc.sd,luc.sd,ag.sd,os.sd)

# total uncertainty on landSink through time
#unc.tot.ts <- sqrt(apply(unc.tot.mat^2,1,sum))

# newer data estimate sigma as 0.8 on average
unc.tot.ts <- rep(0.8, nrow(carbonBudget))

# total uncertainty of land sink change
unc.tot <- sqrt(sum(unc.tot.ts^2))

# total land carbon change at end of 21st century
futureSinkA1B <- GlobalcA1B$totcAnomSmooth[,240]
futureSinkRCP26 <- GlobalcRCP26$totcAnomSmooth[,240]


# there is a strong, linear relationship between 20th century land
# carbon uptake, and eventual land carbon uptake
dev.new(width = 7, height = 7)
par(las = 1, fg = 'grey')
plot(LandSink59_12.A1B,futureSinkA1B ,
     col = A1B.col,
     pch = 19,
     xlab = 'Terrestrial carbon change 1960 to 2012 (GtC)',
     ylab = 'Total terrestrial carbon change 1860 to 2100 (GtC)'
     )
points(LandSink59_12.RCP26,futureSinkRCP26 ,
       col = RCP.col,
       pch = 19
       )
abline(h = 0, v = 0, lty = 'dashed')
abline(v = obs.sink, col = 'darkgrey')
abline(v = obs.sink+ (2*unc.tot), col = 'darkgrey', lty = 'dashed')
abline(v = obs.sink- (2*unc.tot), col = 'darkgrey', lty = 'dashed')

legend('topleft',legvec, col = colvec, pch = 19, bty = 'n',text.col = colvec, inset = 0.07)

text(109, -390, expression(paste('Observations '%+-%2,sigma)), col = 'darkgrey', srt = 90)

if(print.pdf){dev.print(dev = pdf, file = 'graphics/LandUptake20Cvs21C.pdf', width = 7, height = 7)}


# -------------------------------------------------------------------------
# Are there observational constraints on turnover time?
#
# -------------------------------------------------------------------------
dev.new(width = 7, height = 7)
par(las = 1, fg = 'grey')
plot(LandSink59_12.A1B, GlobalcA1B$totcTurnoverYear, pch = 19,col = A1B.col,
     xlab = 'Terrestrial carbon change 1960 to 2012 (GtC)',
     ylab = 'Sink to source transition year')
points(LandSink59_12.RCP26, GlobalcRCP26$totcTurnoverYear, pch = 19,col = RCP.col)


A1B.turnover.naix <- is.na(GlobalcA1B$totcTurnoverYear)
RCP26.turnover.naix <- is.na(GlobalcRCP26$totcTurnoverYear)

abline(h = 0, v = 0, lty = 'dashed')
abline(v = obs.sink, col = 'darkgrey')
abline(v = obs.sink+ (2*unc.tot), col = 'darkgrey', lty = 'dashed')
abline(v = obs.sink- (2*unc.tot), col = 'darkgrey', lty = 'dashed')

legend('topleft',legvec, col = colvec, pch = 19, bty = 'n',text.col = colvec, inset = 0.07)

text(109, 2010, expression(paste('Observations '%+-%2,sigma)), col = 'darkgrey', srt = 90)
rug(LandSink59_12.A1B[A1B.turnover.naix], col = A1B.col, side = 3, lwd = 1.5)
rug(LandSink59_12.RCP26[RCP26.turnover.naix], col = RCP.col, side = 3, lwd = 1.5)

if(print.pdf){dev.print(dev = pdf, file = 'graphics/LandUptake20Cvsturnover.pdf', width = 7, height = 7)}

# find the places where there is no turnover point, and plot them in a rug plot
# (or similar)



#points(LandSink59_12.A1B, GlobalcA1B$vegcTurnoverYear, col = A1B.col, pch = 19)
#points(LandSink59_12.A1B, GlobalcA1B$soilcTurnoverYear, pch = 20,col = A1B.col)

# -------------------------------------------------------------------------




# -------------------------------------------------------------------------
# build a linear model for the land response in the future
# -------------------------------------------------------------------------


# get everything into data frames
sink.A1B.df <- data.frame(histSink = LandSink60_10.A1B,futureSink = futureSinkA1B)
sink.RCP26.df <- data.frame(histSink = LandSink60_10.RCP26,futureSink = futureSinkRCP26)

lm.futureSink.A1B <- lm(futureSink ~ histSink, data = sink.A1B.df) 
lm.futureSink.RCP26 <- lm(futureSink ~ histSink, data = sink.RCP26.df)

# build a prediction from the linear models
pred.dat <- data.frame(histSink = seq(from = -10, to = 150, by = 1))

pred.A1B <- predict(lm.futureSink.A1B,newdata = pred.dat, se.fit = TRUE)
pred.RCP26 <- predict(lm.futureSink.RCP26,newdata = pred.dat, se.fit = TRUE)

lines(pred.dat$histSink, pred.A1B$fit, col = A1B.col )
lines(pred.dat$histSink, pred.A1B$fit+pred.A1B$se.fit, col = A1B.col, lty = 'dashed' )
lines(pred.dat$histSink, pred.A1B$fit-pred.A1B$se.fit, col = A1B.col, lty = 'dashed' )

lines(pred.dat$histSink, pred.RCP26$fit, col = RCP.col )
lines(pred.dat$histSink, pred.RCP26$fit+pred.RCP26$se.fit, col = RCP.col, lty = 'dashed' )
lines(pred.dat$histSink, pred.RCP26$fit-pred.RCP26$se.fit, col = RCP.col, lty = 'dashed' )


# ----------------------------------------------------------------
# Uncertainty in carbon cycle uptake rate 
# ----------------------------------------------------------------
library(zoo)

totcGlobalRollA1B <- t(rollapply(ts(t(GlobalcA1B$totc)),
                                 30,
                                 FUN = mean,
                                 by.column = TRUE)
                       )

totcGlobalRollRCP26 <- t(rollapply(ts(t(GlobalcRCP26$totc)),
                                 30,
                                 FUN = mean,
                                 by.column = TRUE)
                       )

totcGlobalDiffA1B <- t(apply(totcGlobalRollA1B,1,diff))
totcGlobalDiffRCP26 <- t(apply(totcGlobalRollRCP26,1,diff))


dev.new(width = 7, height = 7)
par(las = 1, fg = 'grey')
matplot(1875:2084,t(totcGlobalDiffA1B),
        type = 'l',
        lty = 1,
        col = A1B.col,
        xlab = 'year',
        ylab = expression(paste('Land Carbon Uptake (GtC/year)')) 
        )

matlines(1875:2084,t(totcGlobalDiffRCP26),
        type = 'l',
        lty = 1,
        col = RCP.col,
        )
abline(h = 0, col = 'black', lty = 'dashed')

legend('bottomleft', legvec,
       col = c(colvec),
       text.col = colvec,
       lwd = 1,
       lty = 1
       )
if(print.pdf){dev.print(dev = pdf, file = 'graphics/GlobalLandCarbonUptake2.pdf', width = 7, height = 7)}


colnames(totcGlobalDiffA1B) <- 1875:2084
colnames(totcGlobalDiffRCP26) <- 1875:2084

# A different visualization 
x11()
set.panel(1,2)
image.plot(1980:2084,1:57, t(yearSelect(totcGlobalDiffA1B, 1980:2084)),
           col = ryb, zlim= c(-15,15))

image.plot(1980:2084,1:57, t(yearSelect(totcGlobalDiffRCP26, 1980:2084)),
           col = ryb, zlim= c(-5,5))


# Split out to vegetation and soil
soilcGlobalRollA1B <- t(rollapply(ts(t(GlobalcA1B$soilc)),
                                 30,
                                 FUN = mean,
                                 by.column = TRUE)
                       )

soilcGlobalRollRCP26 <- t(rollapply(ts(t(GlobalcRCP26$soilc)),
                                 30,
                                 FUN = mean,
                                 by.column = TRUE)
                       )

soilcGlobalDiffA1B <- t(apply(soilcGlobalRollA1B,1,diff))
soilcGlobalDiffRCP26 <- t(apply(soilcGlobalRollRCP26,1,diff))



vegcGlobalRollA1B <- t(rollapply(ts(t(GlobalcA1B$vegc)),
                                 30,
                                 FUN = mean,
                                 by.column = TRUE)
                       )

vegcGlobalRollRCP26 <- t(rollapply(ts(t(GlobalcRCP26$vegc)),
                                 30,
                                 FUN = mean,
                                 by.column = TRUE)
                       )

vegcGlobalDiffA1B <- t(apply(vegcGlobalRollA1B,1,diff))
vegcGlobalDiffRCP26 <- t(apply(vegcGlobalRollRCP26,1,diff))



test <- vegcGlobalDiffA1B[1,] + soilcGlobalDiffA1B[1,]


plot(vegcGlobalDiffA1B[1,], type = 'l', ylim = c(-7,2))
lines(soilcGlobalDiffA1B[1,], col = 'red')
lines(test, col = 'blue')

# here
# next: drivers of uptake in the land surface
# (CO2, temperature etc.)

# load temperature data
# ---------------------------------------------------------------
# Global Mean Temperature
# ---------------------------------------------------------------

GlobTempRCP26 <- processTS(file = 'RCP26temps.txt',
                        rlfile = 'RCP26runlist.txt',
                        year = 1860:2099,
                        exclude = 43
                        )

GlobTempAnomRCP26 <- anomalizeTS(GlobTempRCP26, 1:30)


GlobTempA1B <- processTS('A1Btemps.txt',
                      rlfile = 'RCP26runlist.txt',
                      year = 1860:2099,
                      exclude = 43
                      )

# decade labels
decs <- seq(from = 1865, to = 2095, by = 10)

GlobTempSmoothRCP26 <- loessSmoothdf(GlobTempRCP26)
GlobTempSmoothA1B <- loessSmoothdf(GlobTempA1B)

GlobTempAnomA1B <- anomalizeTS(GlobTempA1B, 1:30)

# Non-overlapping decadal smooth
GlobTempDecadalRCP26 <- matAggDec(GlobTempRCP26, ndeltat = 10,FUN = mean, na.rm = TRUE, cn = decs)  
GlobTempDecadalA1B   <- matAggDec(GlobTempA1B, ndeltat = 10,FUN = mean, na.rm = TRUE, cn = decs)

#matplot(decs, t(GlobTempDecadalRCP26), type = 'l')


# decadal Global land surface uptake rates
# Land carbon uptake rates
# *check this!*
decs.short <- seq(from = 1875,to = 2080, by = 10) # shorter timeseries decadal

LandCarbonUptakeDecadalRCP26 <- matAggDec(totcGlobalDiffRCP26,
                                          ndeltat = 10,
                                          FUN = mean,
                                          na.rm = TRUE,
                                          cn = decs.short
                                          )

LandCarbonUptakeDecadalA1B   <- matAggDec(totcGlobalDiffA1B,
                                          ndeltat = 10,
                                          FUN = mean,
                                          na.rm = TRUE,
                                          cn = decs.short
                                          )


# Land carbon uptake in the years 1990 -> 2000
mean(LandCarbonUptakeDecadalRCP26[,'1995'])
mean(LandCarbonUptakeDecadalA1B[,'1995'])


LandUptake90sRCP26 <- apply(totcGlobalDiffRCP26[,match(1990:2000,years)],1,mean)
LandUptake90sA1B <- apply(totcGlobalDiffA1B[,match(1990:2000,years)],1,mean)



x11()
par(mfrow = c(2,1), las = 1)
hist(LandUptake90sRCP26, col = RCP.col, main = 'RCP26')
hist(LandUptake90sA1B, col = A1B.col, main = 'A1B')


totcDiffA1B.nonsmooth <- apply(GlobalcA1B$totc,1,diff)
totcDiffRCP26.nonsmooth <- apply(GlobalcRCP26$totc,1,diff)

# difficult to see amongst the noise
x11()
matplot(1860:2098,totcDiffA1B.nonsmooth, type = 'l', lty = 'solid', col = A1B.col)
matlines(1860:2098,totcDiffRCP26.nonsmooth, type = 'l', lty = 'solid', col = RCP.col)
lines(carbonBudget$year, carbonBudget$landSink, col = 'black') 


# Place on previous land uptake plot
x11()
par(las = 1)
matplot(1875:2084,t(totcGlobalDiffA1B),
        type = 'l',
        lty = 1,
        col = A1B.col,
        xlab = 'year',
        ylab = expression(paste('Land Carbon Uptake (GtC/year)')) 
        )
matlines(1875:2084,t(totcGlobalDiffRCP26),
        type = 'l',
        lty = 1,
        col = RCP.col,
        xlab = 'year',
        ylab = expression(paste('Land Carbon Uptake (GtC/year)')) 
        )
lines(carbonBudget$year, carbonBudget$landSink, col = 'black') 




# decadal temperature against land surface uptake
GlobTempDecadalA1B.short <- yearSelect(GlobTempDecadalA1B, decs.short) - 273.15
GlobTempDecadalRCP26.short <- yearSelect(GlobTempDecadalRCP26, decs.short) - 273.15



dev.new(width = 7, height = 7)
par(las = 1)
plot(c(GlobTempDecadalA1B.short[1, ]), c(LandCarbonUptakeDecadalA1B[1, ]),
     col = A1B.col,
     pch = 21,
     bg = 'blue',
     type = 'o',
     xlim = range(GlobTempDecadalA1B.short, na.rm = TRUE),
     ylim = range(   LandCarbonUptakeDecadalA1B, na.rm = TRUE) ,
     xlab = expression(paste('Global Mean Temperature (',degree,'C)')),
     ylab = expression(paste('Land Carbon Uptake (GtC/year)')) 
     )

for(i in 2:nrow(GlobTempDecadalA1B.short)){
  
  points(c(GlobTempDecadalA1B.short[i, ]), c(LandCarbonUptakeDecadalA1B[i, ]),
         col = A1B.col,
         bg = 'blue',
         pch = 21,
         #bg = 'blue',
         type = 'o'
         )  
 
 
}

for(i in 1:nrow(GlobTempDecadalRCP26.short)){

 points(c(GlobTempDecadalRCP26.short[i, ]), c(LandCarbonUptakeDecadalRCP26[i, ]),
         col = 'tomato2',
         bg  = 'tomato3',
         pch = 21, type = 'o'
         )
}

legend('topright', legvec, col = colvec, lty = 'solid', text.col = colvec, pch = 21, pt.bg = c('blue','tomato3'), bty = 'n')

abline(h = 0, lty = 'dashed')

if(print.pdf){dev.print(dev = pdf, file = 'graphics/turnoverGlobTempVsLandCUptake.pdf', width = 7, height = 7)}



# Is there a relationship between temperature change and carbon uptake by
# the land surface in the 20th century?

# how much has the world warmed in the 20th century?

hadcrut3 <- read.table('/home/h01/hadda/code/R/useful/AnnualHadCRUT3.dat',
                    col.names = c('year','anom','usg','lsg','uc','lc','ub','lb','usgc','lsgc','usgcb','lsgcb')
                      )

# rollmean to get an estimate of the 20C temp change
library(zoo)
hadcrut3.rm.anom <- rollapply(hadcrut3$anom, 30, FUN = mean)

tempyear <- hadcrut3$year[16:(16 + length(hadcrut3.rm.anom) - 1)]

temp.1900 <- hadcrut3.rm.anom[tempyear==1900]
temp.1996 <- hadcrut3.rm.anom[tempyear==1996]

temp.change <- temp.1996 - temp.1900


# temperature change in the ESE
tempchange20c <-yearSelect(GlobTempSmoothRCP26,1996) -  yearSelect(GlobTempSmoothRCP26,1900)

dev.new(width = 5.5, height = 10)
par(mfrow = c(2,1))
quilt.plot(tempchange20c, diff20C.A1B,GlobalcA1B$totcAnomSmooth[,240],
           xlim = c(-0.1,1.75), ylim = c(0, 190),zlim = c(-650, 650),
           nx = 30, ny = 30 ,
           xlab = '20th Century temperature change (K)',
           ylab = '1959 - 2012 land carbon uptake (GtC)',
           legend.args = list('Carbon uptake\nby 2100\n(GtC)', las = 1,side = 3, line = 1),
           main = 'A1B',
           col = ryg
           )

points(temp.change, obs.sink, pch = 4, cex = 2, col = 'black',lwd = 2)
text(temp.change, obs.sink,'observations',pos = 4)



quilt.plot(tempchange20c, diff20C.RCP26,GlobalcA1B$totcAnomSmooth[,240],
           xlim = c(-0.1,1.75), ylim = c(0, 190),zlim = c(-650, 650),
           nx = 30, ny = 30 ,
           xlab = '20th Century temperature change (K)',
           ylab = '1959 - 2012 land carbon uptake (GtC)',
           legend.args = list('Carbon uptake\nby 2100\n(GtC)', las = 1,side = 3, line = 1),
           main = 'RCP2.6',
           col = ryg)

points(temp.change, obs.sink, pch = 4, cex = 2, col = 'black',lwd = 2)
text(temp.change, obs.sink,'observations',pos = 4)


if(print.pdf){dev.print(dev = pdf, file = 'graphics/tempCarbonObs.pdf', width = 5.5, height = 10)}



# -------------------------------------------------------------------------
# CO2 in the atmosphere
#
# -------------------------------------------------------------------------

# copied across from carbonBudget.R
# Co2 Mass mixing ratio
co2mmrRCP26 <- processTS(file = 'RCP26CO2mmr.txt',
                        rlfile = 'RCP26runlist.txt',
                        year = 1860:2099,
                        exclude = 43
                        )



co2mmrA1B <- processTS(file = 'A1BCO2mmr.txt',
                        rlfile = 'RCP26runlist.txt',
                        year = 1860:2099,
                        exclude = 43
                        )


ensCompTShist(co2mmrA1B,co2mmrRCP26,
          colvec = colvec,
          legvec = legvec,
          mainvec = 'CO2 Mass mixing ratio',
          xlim = c(1860, 2100),
          xlab = '',
          ylab = 'Carbon',
          )

# calculate the carbon in the atmosphere

#ppmv
co2ppmA1B <- co2mmrA1B * (290 / 4.4e-4)
co2ppmRCP26 <- co2mmrRCP26 * (290 / 4.4e-4)

ensCompTShist(co2ppmA1B,co2ppmRCP26,
          colvec = colvec,
          legvec = legvec,
          mainvec = '',
          xlim = c(1860, 2100),
          xlab = '',
          ylab = 'CO2 concentration (ppmv)',
          )


# Gt Carbon
co2GtCA1B <- co2ppmA1B * 2.136
co2GtCRCP26 <- co2ppmRCP26 * 2.136

co2ppmDecadalRCP26 <- matAggDec(co2ppmRCP26, ndeltat = 10, FUN = mean, cn = decs, na.rm = TRUE)
co2ppmDecadalA1B   <- matAggDec(co2ppmA1B, ndeltat = 10, FUN = mean, cn = decs, na.rm = TRUE)

co2ppmDecadalRCP26.short <- yearSelect(co2ppmDecadalRCP26, decs.short)
co2ppmDecadalA1B.short <- yearSelect(co2ppmDecadalA1B, decs.short)


dev.new(width = 7, height = 7)
par(las = 1)
plot(c(co2ppmDecadalA1B.short[1, ]), c(LandCarbonUptakeDecadalA1B[1, ]),
     col = A1B.col,
     pch = 21,
     bg = 'blue',
     type = 'o',
     xlim = range(co2ppmDecadalA1B.short, na.rm = TRUE),
     ylim = range(   LandCarbonUptakeDecadalA1B, na.rm = TRUE) ,
     xlab = 'Atmospheric CO2 concentration (PPMV)',
     ylab = expression(paste('Land Carbon Uptake (GtC/year)')) 
     )

for(i in 2:nrow(co2ppmDecadalA1B.short)){
  
  points(c(co2ppmDecadalA1B.short[i, ]), c(LandCarbonUptakeDecadalA1B[i, ]),
         col = A1B.col,
         bg = 'blue',
         pch = 21,
         type = 'o'
         )  
 
 
}

for(i in 1:nrow(co2ppmDecadalRCP26.short)){

 points(c(co2ppmDecadalRCP26.short[i, ]), c(LandCarbonUptakeDecadalRCP26[i, ]),
         col = 'tomato2',
         bg  = 'tomato3',
         pch = 21, type = 'o'
         )
}

legend('topright', legvec, col = colvec, lty = 'solid', text.col = colvec, pch = 21, pt.bg = c('blue','tomato3'), bty = 'n')

abline(h = 0, lty = 'dashed')

if(print.pdf){dev.print(dev = pdf, file = 'graphics/turnoverCO2ppmVsLandCUptake.pdf', width = 7, height = 7)}


# -------------------------------------------------------------------------------------
# Temperature and CO2 concentration at the point where the land carbon
# cycle turns form sink to source 
#
# -------------------------------------------------------------------------------------

# toy = 'turn over year'
toy.A1B   <- GlobalcA1B$totcTurnoverYear
toy.RCP26 <- GlobalcRCP26$totcTurnoverYear


GlobTempA1B.col.ix <- match(toy.A1B, colnames(GlobTempA1B))
GlobTempRCP26.col.ix <- match(toy.RCP26, colnames(GlobTempRCP26))


co2ppmA1B.col.ix <- match(toy.A1B, colnames(co2ppmA1B))
co2ppmRCP26.col.ix <- match(toy.RCP26, colnames(co2ppmRCP26))


turnovertemps.A1B <- as.matrix(GlobTempA1B)[cbind(1:57,GlobTempA1B.col.ix )]
turnovertemps.RCP26 <- as.matrix(GlobTempRCP26)[cbind(1:57,GlobTempRCP26.col.ix )]


turnoverppm.A1B <- as.matrix(co2ppmA1B)[cbind(1:57,co2ppmA1B.col.ix )]
turnoverppm.RCP26 <- as.matrix(co2ppmRCP26)[cbind(1:57,co2ppmRCP26.col.ix )]


# plot up
x11()
plot(c(co2ppmA1B, recursive = TRUE), c(GlobTempA1B, recursive = TRUE), col = 'grey')
points(c(co2ppmRCP26, recursive = TRUE), c(GlobTempRCP26, recursive = TRUE), col = 'darkgrey')

points(turnoverppm.A1B, turnovertemps.A1B, col = A1B.col, pch = 19)
points(turnoverppm.RCP26, turnovertemps.RCP26, col = RCP.col, pch = 19)


# --------------------------------------------------------------------
# Relationship between emissions and carbon stores
# --------------------------------------------------------------------

plot(emissions$A1B, c(GlobalcA1B$totcAnom[1,], recursive = TRUE),
     xlim = c(0, 17), ylim = c(-600, 600),
     col = A1B.col,
     type = 'l',
     #pch = 20

     )

 for (i in 2:57){
   
   points(emissions$A1B, c(GlobalcA1B$totcAnom[i,], recursive = TRUE),
          col = A1B.col, type = 'l')

 }
 for (i in 1:57){
   
   points(emissions$RCP26, c(GlobalcRCP26$totcAnom[i,], recursive = TRUE),
          col = RCP.col, type = 'l')

 }



plot(emissions$A1B[emissions$Year %in% 1875:2084], abs(totcGlobalDiffA1B[1,]), type = 'o',
     ylim = c(0,12), xlim = c(0,17))

 for (i in 2:57){
   
   points(emissions$A1B[emissions$Year %in% 1875:2084],abs(totcGlobalDiffA1B[i,]),
          col = A1B.col, type = 'o')

 }



 for (i in 1:57){
   
   points(emissions$RCP26, c(GlobalcRCP26$totcAnom[i,], recursive = TRUE),
          col = RCP.col, type = 'l')

 }





# -------------------------------------------------------------------
# Airborne fraction
#
#
# -------------------------------------------------------------------

# plot of global CO2 emissions
par(bty = 'l', mar = c(5,4,2,1), las = 1)
plot(emissions$Year, emissions$A1B,
     type = 'l',
     col = A1B.col,
     lwd = 2,
     xlab = 'Year',
     ylab = 'Emissions (GtC)'
     )
lines(emissions$Year, emissions$RCP26, type = 'l', col = RCP.col, lwd = 2)
text(2090, 12, 'A1B', col = A1B.col, cex = 1.5)
text(2060, 6, 'RCP2.6', col = RCP.col, cex = 1.5)



# what is the change in atmospheric carbon since 1860?
 
delco2GtCA1B <- sweep(co2GtCA1B, 1,co2GtCA1B[,1])
delco2GtCRCP26 <- sweep(co2GtCRCP26, 1,co2GtCRCP26[,1])

delco2GtCsmoothA1B <- loessSmoothdf(delco2GtCA1B )
delco2GtCsmoothRCP26 <-  loessSmoothdf(delco2GtCRCP26)


# divide change in CO2 content from 1860, by cumulative emissions from 1860
airborneFractionA1B <- sweep(x = delco2GtCA1B,
                             MARGIN =  2,
                             STATS = cumsum(emissions$A1B),
                             FUN = '/'
                             )

airborneFractionRCP26 <- sweep(x = delco2GtCRCP26,
                               MARGIN =  2,
                               STATS = cumsum(emissions$RCP26),
                               FUN = '/'
                               )

# Airborne fraction anomaly
airborneFractionAnomA1B <- anomalizeTS(airborneFractionA1B, 100:130)
airborneFractionAnomRCP26 <- anomalizeTS(airborneFractionRCP26, 100:130)


# Airborne fraction plot
dev.new(width = 7, height = 7)
par(mfrow = c(1,2), las = 1, fg = 'grey')

# AF as a percentage
matplot(years[120:240], t(airborneFractionA1B[,120:240] *100),
        ylim = c(0,100),
        type = 'l',
        lty = 'solid',
        col = A1B.col,
        xlab = 'Year',
        ylab = 'Airborne Fraction (%)'
        )
matlines(years[120:240],t(airborneFractionRCP26[,120:240] *100),
         type = 'l',
         lty = 'solid',
         col = RCP.col
         )
mtext('a)',adj = 0, cex = 2, line = 0.7, col = text.colour)
legend('topleft', legvec, col = colvec, lty = 'solid', text.col = colvec,bty = 'n')


# AF anomaly from 1960 - 1990
matplot(years[120:240],t(airborneFractionAnomA1B[,120:240] *100),
        ylim = c(-50,50),
        type = 'l',lty = 'solid', col = A1B.col,
        xlab = 'Year',
        ylab = 'Airborne Fraction change from 1960-90 (%)'
        )
matlines(years[120:240], t(airborneFractionAnomRCP26[,120:240] *100),
         type = 'l',
         lty = 'solid',
         col = RCP.col
         )

mtext('b)',adj = 0, cex = 2, line = 0.7, col = text.colour)

if(print.pdf){dev.print(dev = pdf, file = 'graphics/turnoverAirborneFraction.pdf', width = 7, height = 7)}


# Fraction of emissions that go into the land surface
# divide change in CO2 content from 1860, by cumulative emissions from 1860
LandborneFractionA1B <- sweep(x = GlobalcA1B$totcAnom,
                             MARGIN =  2,
                             STATS = cumsum(emissions$A1B),
                             FUN = '/'
                             )

LandborneFractionRCP26 <- sweep(x = GlobalcRCP26$totcAnom,
                               MARGIN =  2,
                               STATS = cumsum(emissions$RCP26),
                               FUN = '/'
                               )


# Landborne Fraction / Airborne Fraction Plot


dev.new(width = 7, height = 7)
par(las = 1)
plot(yearSelect(LandborneFractionA1B,2000), yearSelect(airborneFractionA1B,2000),
     xlim = c(-0.4,0.7), ylim = c(0,1),
     col = 'skyblue', pch = 21,
     xlab = 'Land Fraction', ylab = 'Airborne Fraction'
     )


af <- function(l, o = 0.1){

  out <- 1 - (l + o)
  out
}

l <- seq(from = -1, to = 1, by = 0.01)
os <- seq(from = -2, to = 2, by = 0.2)

test <- lapply(os, FUN = af, l = l)

## for(i in 1:length(test)){
  
## lines(l,test[[i]], col = 'grey')
## }
abline(h = 0, v = 0, col = 'grey', lty = 'dashed')

for (i in 1:length(os)){
  abline(os[i],-1, col = 'grey')
}

for (i in 1:length(os)){
  
  text(0,af(0,os[i]),as.character(os[i]))
}

points(yearSelect(LandborneFractionRCP26,2000), yearSelect(airborneFractionA1B,2000), col = 'tomato', pch = 21)


segments(x0 = yearSelect(LandborneFractionA1B,2000), y0 = yearSelect(airborneFractionA1B,2000),
         x1 = yearSelect(LandborneFractionA1B,2099), y1 = yearSelect(airborneFractionA1B,2099), col = 'grey')

segments(x0 = yearSelect(LandborneFractionRCP26,2000), y0 = yearSelect(airborneFractionA1B,2000),
         x1 = yearSelect(LandborneFractionRCP26,2099), y1 = yearSelect(airborneFractionRCP26,2099),
         col = 'grey')
points(yearSelect(LandborneFractionA1B,2099), yearSelect(airborneFractionA1B,2099), col = A1B.col, pch = 19)
points(yearSelect(LandborneFractionRCP26,2099), yearSelect(airborneFractionRCP26,2099), col = RCP.col, pch = 19)

text(0,1,'Ocean Fraction',pos = 4, srt = -45)

legend('topright', legend = c(legvec,'2000','2100'), col = c(colvec, colvec[1],colvec[1]), pch = c(19,19,21,19), text.col = c(colvec, colvec[1],colvec[1]), bg = 'white' )

if(print.pdf){dev.print(dev = pdf, file = 'graphics/ALOFractionplot.pdf', width = 7, height = 7)}



# Landborne Fraction / Airborne Fraction Anomaly Plot
LandborneFractionAnomA1B <- anomalizeTS(LandborneFractionA1B, 100:130)
LandborneFractionAnomRCP26 <- anomalizeTS(LandborneFractionRCP26, 100:130)

x11()
plot(yearSelect(LandborneFractionAnomA1B,2000), yearSelect(airborneFractionAnomA1B,2000),
     xlim = c(-0.8,0.1), ylim = c(-0.2,0.7),
     col = 'skyblue', pch = 21,
     xlab = 'Land Fraction Anomaly', ylab = 'Airborne Fraction Anomaly'
     )
abline(h = 0, v = 0, col = 'grey', lty = 'dashed')
points(yearSelect(LandborneFractionAnomRCP26,2000), yearSelect(airborneFractionAnomA1B,2000), col = 'tomato', pch = 21)


segments(x0 = yearSelect(LandborneFractionAnomA1B,2000), y0 = yearSelect(airborneFractionAnomA1B,2000),
         x1 = yearSelect(LandborneFractionAnomA1B,2099), y1 = yearSelect(airborneFractionAnomA1B,2099), col = 'grey')

segments(x0 = yearSelect(LandborneFractionAnomRCP26,2000), y0 = yearSelect(airborneFractionAnomA1B,2000),
         x1 = yearSelect(LandborneFractionAnomRCP26,2099), y1 = yearSelect(airborneFractionAnomRCP26,2099),
         col = 'grey')


points(yearSelect(LandborneFractionAnomA1B,2099), yearSelect(airborneFractionAnomA1B,2099), col = A1B.col, pch = 19)
points(yearSelect(LandborneFractionAnomRCP26,2099), yearSelect(airborneFractionAnomRCP26,2099), col = RCP.col, pch = 19)



stop()

# ------------------------------------------------------------------
# Woah there! Surely Q10 has a strong influence on the carbon cycle?
#
# ------------------------------------------------------------------

load('/home/h01/hadda/avoid/kitchenSink/makeEns/RCP26design.rda')


rlfile = 'RCP26runlist.txt'
runlist <- scan(file = rlfile, what = c('character'))[-43]

#rownames(RCP26design) %in% runlist

keep.ix <- match(runlist, rownames(RCP26design))
design <- RCP26design[keep.ix, ]

 
#x11(width = 8, height = 6)
pdf(width = 8, height = 6, file = 'graphics/UptakeVscarbonParamsMarginal.pdf')
par(mfrow = c(2,3), las = 1, pch = 19, cex.lab = 1.2)

plot(design$Q10, GlobalcA1B$totcAnomSmooth[,240], xlab = 'Q10', ylab = 'Carbon change (GtC)', col = A1B.col)
points(design$Q10, GlobalcRCP26$totcAnomSmooth[,240], xlab = 'Q10', ylab = 'Carbon change (GtC)', col = RCP.col)
legend('topright', legvec,col = colvec,text.col = colvec, pch = 19)

plot(design$NL0, GlobalcA1B$totcAnomSmooth[,240], xlab = 'NL0', ylab = 'Carbon change (GtC)', col = A1B.col)
points(design$NL0, GlobalcRCP26$totcAnomSmooth[,240], xlab = 'Q10', ylab = 'Carbon change (GtC)', col = RCP.col)

plot(design$TUPP, GlobalcA1B$totcAnomSmooth[,240], xlab = 'TUPP', ylab = 'Carbon change (GtC)', col = A1B.col)
points(design$TUPP, GlobalcRCP26$totcAnomSmooth[,240], xlab = 'Q10', ylab = 'Carbon change (GtC)', col = RCP.col)

plot(design$F0, GlobalcA1B$totcAnomSmooth[,240],xlab = 'F0', ylab = 'Carbon change (GtC)', col = A1B.col)
points(design$F0, GlobalcRCP26$totcAnomSmooth[,240], xlab = 'Q10', ylab = 'Carbon change (GtC)', col = RCP.col)

plot(design$LAI_MIN, GlobalcA1B$totcAnomSmooth[,240],xlab = 'LAI_MIN', ylab = 'Carbon change (GtC)', col = A1B.col)
points(design$LAI_MIN, GlobalcRCP26$totcAnomSmooth[,240], xlab = 'Q10', ylab = 'Carbon change (GtC)', col = RCP.col)

plot(design$V_CRIT_ALPHA, GlobalcA1B$totcAnomSmooth[,240],xlab = 'V_CRIT_ALPHA', ylab = 'Carbon change (GtC)', col = A1B.col)
points(design$V_CRIT_ALPHA, GlobalcRCP26$totcAnomSmooth[,240], xlab = 'Q10', ylab = 'Carbon change (GtC)', col = RCP.col)
dev.off()

x11(width = 8, height = 6)
par(mfrow = c(2,3), las = 1)

plot(design$Q10, GlobalcA1B$soilcAnomSmooth[,240], xlab = 'Q10', ylab = 'Soil C change')

plot(design$NL0, GlobalcA1B$soilcAnomSmooth[,240], xlab = 'NL0', ylab = 'Soil C change')

plot(design$TUPP, GlobalcA1B$soilcAnomSmooth[,240], xlab = 'TUPP', ylab = 'Soil C change')

plot(design$F0, GlobalcA1B$soilcAnomSmooth[,240],xlab = 'F0', ylab = 'Soil C change')

plot(design$LAI_MIN, GlobalcA1B$soilcAnomSmooth[,240],xlab = 'LAI_MIN', ylab = 'Soil C change')

plot(design$V_CRIT_ALPHA, GlobalcA1B$soilcAnomSmooth[,240],xlab = 'V_CRIT_ALPHA', ylab = 'Soil C change')



x11(width = 8, height = 6)
par(mfrow = c(2,3), las = 1)

plot(design$Q10, GlobalcA1B$vegcAnomSmooth[,240], xlab = 'Q10', ylab = 'Veg C change')

plot(design$NL0, GlobalcA1B$vegcAnomSmooth[,240], xlab = 'NL0', ylab = 'Veg C change')

plot(design$TUPP, GlobalcA1B$vegcAnomSmooth[,240], xlab = 'TUPP', ylab = 'Veg C change')

plot(design$F0, GlobalcA1B$vegcAnomSmooth[,240],xlab = 'F0', ylab = 'Veg C change')

plot(design$LAI_MIN, GlobalcA1B$vegcAnomSmooth[,240],xlab = 'LAI_MIN', ylab = 'Veg C change')

plot(design$V_CRIT_ALPHA, GlobalcA1B$vegcAnomSmooth[,240],xlab = 'V_CRIT_ALPHA', ylab = 'Veg C change')




# metaparameter design
meta <- read.table('/home/h01/hadda/avoid/kitchenSink/makeEns/metaparameters.txt', skip = 3, header=TRUE)[-43, ]


cparam <- subset(design, select = c(Q10, NL0, TUPP, F0, LAI_MIN, V_CRIT_ALPHA))


# could build a gp emulator, and do a formal sensitivity analysis here?

# how was the 'metaparameter' for the carbon cycle chosen?
plot(meta$c, GlobalcA1B$soilcAnomSmooth[,240])



# --------------------------------------------------------------
# Carbon flux behaviour by metaparameter
#
# --------------------------------------------------------------

decs.20Con <-  seq(from = 1895, to = 2075, by = 5)

tempsA1B <- yearSelect(GlobTempDecadalA1B - 273.15, decs.20Con )
cfluxA1B <- yearSelect(LandCarbonUptakeDecadalA1B, decs.20Con )
cconcA1B <- yearSelect(co2ppmDecadalA1B,decs.20Con) 

tempsRCP26 <- yearSelect(GlobTempDecadalRCP26 - 273.15,decs.20Con)
cfluxRCP26 <- yearSelect(LandCarbonUptakeDecadalRCP26, decs.20Con)
cconcRCP26 <- yearSelect(co2ppmDecadalRCP26, decs.20Con)


 
#x11(width = 7, height = 9)
pdf(file = 'graphics/metaSplitTempvCcycle.pdf', width = 7, height = 9)
par(mfrow = c(5,3), mar = c(1,1,1,1), oma = c(5,5,1,1), las = 1, cex.axis = 1.2)
for (i in (0:16)[-c(14,16)]) {
  
  print(i)                                  # ensemble rows with apropriate C-cycle metaparameter
  j <- which(meta[,2] == i)
 
  plot(c(tempsA1B[j[1] , ], recursive = TRUE),   c(cfluxA1B[j[1] , ] , recursive = TRUE)  , xlim = c(13,21), ylim = c(-12,7),
       type = 'o',
       pch = 21,
       bg = 'blue',
       xlab = '',#expression(paste('Global Mean Temperature ', degree, '(C)')),
       ylab = '',#'Carbon flux to land surface (GtC)',
       col = A1B.col,
       axes = FALSE
       )
  
  if (i %in% c(12,14,16)) axis(1)
  if (i %in% c(0,3,6,9,12)) axis(2)
  if (i %in% 0) legend('topright', legvec, col = colvec, pch = 21, pt.bg = c('blue', 'red'), lty = 'solid', bty = 'n')
  
  points(c(tempsRCP26[j[1] , ], recursive = TRUE),   c(cfluxRCP26[j[1] , ] , recursive = TRUE)  ,type = 'o',pch = 21, col = RCP.col, bg = 'red')
  abline(h = 0, lty = 'dashed')
  
  
  for ( r in 2:length(j)){
    
    points(c(tempsA1B[j[r] , ], recursive = TRUE),   c(cfluxA1B[j[r] , ] , recursive = TRUE), type = 'o', pch = 21,col = A1B.col , bg = 'blue'  )
    points(c(tempsRCP26[j[r] , ], recursive = TRUE),   c(cfluxRCP26[j[r] , ] , recursive = TRUE) ,type = 'o',pch = 21,  col = RCP.col, bg = 'red'  )
    
  }


  mtext(expression(paste('Global Mean Temperature ', degree, '(C)')), side = 1, outer = TRUE, cex = 1.5, line = 3)
  mtext('Land Carbon Uptake (GtC/year)', side = 2, outer = TRUE, cex = 1.5, line = 3, las =0)
}

dev.off()



x11(width = 8, height = 12)
#pdf(file = 'graphics/metaSplitCcycle.pdf', width = 8, height = 12)
par(mfrow = c(5,3), mar = c(4,4,2,1))

for (i in (0:16)[-c(14,16)]) {
                                        # ensemble rows with apropriate C-cycle metaparameter
  j <- which(meta[,2] == i)
  
  plot(c(cconcA1B[j[1] , ], recursive = TRUE),   c(cfluxA1B[j[1] , ] , recursive = TRUE)  , xlim = c(270,900), ylim = c(-12,7),
       xlab = expression(paste('CO2 concentration ', degree, '(C)')),
       ylab = 'Carbon flux to land surface (GtC)',
       col = A1B.col
       )
  
  points(c(cconcRCP26[j[1] , ], recursive = TRUE),   c(cfluxRCP26[j[1] , ] , recursive = TRUE)  , col = RCP.col)
  abline(h = 0, lty = 'dashed')
  
  
  for ( r in 2:length(j)){
    
    points(c(cconcA1B[j[r] , ], recursive = TRUE),   c(cfluxA1B[j[r] , ] , recursive = TRUE), col = A1B.col   )
    points(c(cconcRCP26[j[r] , ], recursive = TRUE),   c(cfluxRCP26[j[r] , ] , recursive = TRUE) ,col = RCP.col  )
    
  }
  
}








# ----------------------------------------------------------
# Maps of terrestrial carbon changes
#
#
#
# ----------------------------------------------------------






# HadCM3 longitudes centered on 0 longitude
HadCM3.lons.wrt <- seq(from = -180, to = 176.25, by = 3.75)

mapSwap <- function(mapmat){
  # shift map data matrix so that it is
  # centered on 0 degrees longitude
  
  nc <- ncol(mapmat)
  
  out <- cbind(mapmat[ , ((nc/2) + 1):nc ], mapmat[ , 1:(nc/2) ])
                    
}

mapEns <- function(dat, nrow = 73, ncol = 96){
  # create a summary map ensemble object, from the raw
  # ensemble data
  # outputs maps of mean, standard deviation min, max and
  # signal-to-noise ratio

  ens.sd <- apply(dat, 2, sd, na.rm = TRUE)
  sdmap <- matrix(ens.sd, nrow = nrow, ncol = ncol, byrow = TRUE)

  ens.mean <- apply(dat, 2, mean, na.rm = TRUE)
  meanmap <- matrix(ens.mean, nrow = nrow, ncol = ncol, byrow = TRUE)

  ens.max <- apply(dat, 2, max, na.rm = TRUE)
  maxmap  <- matrix(ens.max, nrow = nrow, ncol = ncol, byrow = TRUE)
  
  ens.min <- apply(dat, 2, min, na.rm = TRUE)
  minmap  <- matrix(ens.min, nrow = nrow, ncol = ncol, byrow = TRUE)
  
  snrmap <- meanmap / sdmap
  
  return(list(sdmap = sdmap ,meanmap = meanmap,maxmap = maxmap, minmap = minmap, snrmap = snrmap)) 

}


mapBGcol <- function(lons, lats, mapdat,BGcol = 'lightgrey', legend.args,axis.args, mapcol, zlim,titext,world = FALSE,axes = FALSE, ...){
  # a function to map HadCM3 data on a grey background

  nr <- nrow(mapdat)

  x0 <- min(lons)
  x1 <- max(lons)

  y0 <- min(lats)
  y1 <- max(lats)
  # first, a null plot

  plot(c(x0,x1), c(y0,y1),type = 'n',xaxs = 'i', yaxs = 'i' ,xlab = '', ylab = '', axes = FALSE)
  #plot(c(x0,x1), c(y0,y1),xaxs = 'i', yaxs = 'i' ,xlab = '', ylab = '')
  
 # plot(c(0,1), c(0,1),type = 'n',xaxs = 'i', yaxs = 'i' ,xlab = '', ylab = '', axes = FALSE)

  # and a polygon
  polygon(x = c(x0,x1,x1,x0), y = c(y0,y0,y1,y1), col = BGcol, border = NA)
  #polygon(x = c(0,1,1,0), y = c(0,0,1,1), col = BGcol, border = NA)

  # the map itself
  image(lons, lats, t(mapdat[nr:1, ]), col = mapcol , zlim = zlim, axes = axes,xlab="", ylab="", add = TRUE,...)

  if(world){
    world(add = TRUE)
  }

  mtext(titext,side = 3, adj = 0,line = 0.5,cex = 1.2)
  
  # add a legend
  image.plot(t(mapdat[nr:1, ]), zlim = zlim,
             xlab="", ylab="",
             col = mapcol,
             legend.only = TRUE,
             horizontal = TRUE,
             legend.args = legend.args,
             axis.args = axis.args
             )
}



# Carbon maps
VegCarbonPreIndTotal <-read.table(file = 'data/VegCarbonPreIndTotal.txt',
                                  na.strings = c('-999.0000000','-1073741824.0000000'))[-43, ]

VegCarbonEOCTotalA1B <- read.table(file = 'data/VegCarbonEOCTotalA1B.txt',
                                  na.strings = c('-999.0000000','-1073741824.0000000'))[-43, ]


VegCarbonEOCTotalRCP26 <- read.table(file = 'data/VegCarbonEOCTotalRCP26.txt',
                                  na.strings = c('-999.0000000','-1073741824.0000000'))[-43, ]

SoilCarbonPreIndTotal <-read.table(file = 'data/SoilCarbonPreIndTotal.txt',
                                  na.strings = c('-999.0000000','-1073741824.0000000'))[-43, ]

SoilCarbonEOCTotalA1B <- read.table(file = 'data/SoilCarbonEOCTotalA1B.txt',
                                  na.strings = c('-999.0000000','-1073741824.0000000'))[-43, ]


SoilCarbonEOCTotalRCP26 <- read.table(file = 'data/SoilCarbonEOCTotalRCP26.txt',
                                  na.strings = c('-999.0000000','-1073741824.0000000'))[-43, ]


# ensemble summaries
VegmapsPreInd <- mapEns(VegCarbonPreIndTotal)
VegmapsEOCA1B <- mapEns(VegCarbonEOCTotalA1B)
VegmapsEOCRCP26 <- mapEns(VegCarbonEOCTotalRCP26)

SoilmapsPreInd <- mapEns(SoilCarbonPreIndTotal)
SoilmapsEOCA1B <- mapEns(SoilCarbonEOCTotalA1B)
SoilmapsEOCRCP26 <- mapEns(SoilCarbonEOCTotalRCP26)


x11(width = 8, height = 10)
set.panel(3,2)

z <- c(0,30)

par(las = 1)
image.plot(HadCM3.lons.wrt, HadCM3.lats, t(mapSwap(VegmapsPreInd$meanmap)[73:1, ]),
           col = byr,
           xlab = '',
           ylab = '',
           zlim = z,
           main = 'PreInd Veg C'
           )
world(add = TRUE)

par(las = 1)
image.plot(HadCM3.lons.wrt, HadCM3.lats, t(mapSwap(SoilmapsPreInd$meanmap)[73:1, ]),
           col = byr,
           xlab = '',
           ylab = '',
           zlim = z,
           main = 'PreInd Soil C'
           )
world(add = TRUE)


par(las = 1)
image.plot(HadCM3.lons.wrt, HadCM3.lats, t(mapSwap(VegmapsEOCA1B$meanmap)[73:1, ]),
           col = byr,
           xlab = '',
           ylab = '',
           zlim = z,
           main = 'EOC Veg C A1B'
           )
world(add = TRUE)

par(las = 1)
image.plot(HadCM3.lons.wrt, HadCM3.lats, t(mapSwap(SoilmapsEOCA1B$meanmap)[73:1, ]),
           col = byr,
           xlab = '',
           ylab = '',
           zlim = z,
           main = 'EOC Soil C A1B'
           )
world(add = TRUE)

par(las = 1)
image.plot(HadCM3.lons.wrt, HadCM3.lats, t(mapSwap(VegmapsEOCRCP26$meanmap)[73:1, ]),
           col = byr,
           xlab = '',
           ylab = '',
           zlim = z,
           main = 'EOC Veg C RCP26'
           )
world(add = TRUE)


par(las = 1)
image.plot(HadCM3.lons.wrt, HadCM3.lats, t(mapSwap(SoilmapsEOCRCP26$meanmap)[73:1, ]),
           col = byr,
           xlab = '',
           ylab = '',
           zlim = z,
           main = 'EOC Soil C RCP26'
           )
world(add = TRUE)


# Total terrestrial carbon
# pre-industrial
LandCarbonPreInd <- VegCarbonPreIndTotal + SoilCarbonPreIndTotal

# End of century
LandCarbonEOCA1B <- VegCarbonEOCTotalA1B + SoilCarbonEOCTotalA1B
LandCarbonEOCRCP26 <- VegCarbonEOCTotalRCP26 + SoilCarbonEOCTotalRCP26

LandCarbonEOCA1BMaps <- mapEns(LandCarbonEOCA1B)
LandCarbonEOCRCP26Maps <- mapEns(LandCarbonEOCRCP26)


# Plots of land surface carbon at end of century
x11(width = 7, height = 10)
par(las = 1)
set.panel(2,1)
image.plot(HadCM3.lons.wrt, HadCM3.lats, t(mapSwap(LandCarbonEOCA1BMaps$meanmap)[73:1, ]),
           col = g,
           xlab = '',
           ylab = '',
           main = 'EOC A1B Land carbon',
           zlim = c(0,40)
           )
world(add = TRUE)

image.plot(HadCM3.lons.wrt, HadCM3.lats, t(mapSwap(LandCarbonEOCRCP26Maps$meanmap)[73:1, ]),
           col = g,
           xlab = '',
           ylab = '',
           main = 'EOC RCP26 Land carbon',
           zlim = c(0,40)
           )
world(add = TRUE)



# Paired difference between land carbon in RCP2.6 and A1B
LandCarbonEOCDiff <- LandCarbonEOCRCP26 - LandCarbonEOCA1B

LandCarbonEOCDiffMaps <- mapEns(LandCarbonEOCDiff)


# Mean paired difference between RCP26 and A1B
x11(width = 7, height = 5)
par(las = 1)
image.plot(HadCM3.lons.wrt, HadCM3.lats, t(mapSwap(LandCarbonEOCDiffMaps$meanmap)[73:1, ]),
           col = br,
           xlab = '',
           ylab = '',
           main = 'EOC RCP26 - A1B',
           zlim = c(-3,3)
           )
world(add = TRUE)


# ------------------------------------------------------------
# Maps of ensemble Land carbon changes by end of century
# (postage stamp style)
# ------------------------------------------------------------

LandCarbonChangeA1B <- LandCarbonEOCA1B - LandCarbonPreInd
LandCarbonChangeRCP26 <- LandCarbonEOCRCP26 - LandCarbonPreInd

LandCarbonChangePercA1B <- (LandCarbonChangeA1B / LandCarbonPreInd) *100
LandCarbonChangePercRCP26 <- (LandCarbonChangeRCP26 / LandCarbonPreInd) *100



sumEOC <- apply(LandCarbonChangeA1B, 1,sum, na.rm = TRUE)
ix <- sort(sumEOC, index.return = TRUE)$ix




pdf(file = 'graphics/LandCEOCchangeA1B.pdf', width = 8, height = 10)
#x11(width = 8, height = 9)
par(mfrow = c(10,6), mar = c(0.3,0.3,0.3,0.3), oma = c(0.5, 0.5, 0.5, 0.5))

zl <- c(-10,10)

for(i in ix){

  test <-  t(mapSwap(matrix(c(LandCarbonChangeA1B[i,], recursive = TRUE), nrow = 73, ncol = 96, byrow = TRUE))[73:1, ])
  
  test[ test < zl[1]  ] <- zl[1]
  test[ test > zl[2] ] <- zl[2]
  # chop off antarctica in the plotting
  image(HadCM3.lons.wrt, HadCM3.lats[13:73], test[,13:73],
             col = rb,
             xlab = '',
             ylab = '',
             zlim = zl,
        axes = FALSE
             )
  world(add = TRUE)
  
}

plot(1:10, type = 'n', axes = FALSE, xlab = '', ylab = '')
image.plot(t(test[73:1,]), zlim = zl,
        xlab="", ylab="", col = rb,
           axes = FALSE,
           legend.width = 5,
           legend.mar = 14.1,
           legend.only = TRUE,
           horizontal = TRUE,
           legend.args = list(text = expression('kg m'^-2), line = 0.5)
           )

dev.off() 


sumEOC <- apply(LandCarbonChangeRCP26, 1,sum, na.rm = TRUE)
ix <- sort(sumEOC, index.return = TRUE)$ix

pdf(file = 'graphics/LandCEOCchangeRCP26.pdf', width = 8, height = 10)
#x11(width = 8, height = 9)
par(mfrow = c(10,6), mar = c(0.3,0.3,0.3,0.3), oma = c(0.5, 0.5, 0.5, 0.5))

zl <- c(-10,10)

for(i in ix){

  test <-  t(mapSwap(matrix(c(LandCarbonChangeRCP26[i,], recursive = TRUE), nrow = 73, ncol = 96, byrow = TRUE))[73:1, ])
  
  test[ test < zl[1]  ] <- zl[1]
  test[ test > zl[2] ] <- zl[2]
  # chop off antarctica in the plotting
  image(HadCM3.lons.wrt, HadCM3.lats[13:73], test[,13:73],
             col = rb,
             xlab = '',
             ylab = '',
             zlim = zl,
        axes = FALSE,
        bg = 'grey'
             )
  world(add = TRUE)
  
}

plot(1:10, type = 'n', axes = FALSE, xlab = '', ylab = '')
image.plot(t(test[73:1,]), zlim = zl,
        xlab="", ylab="", col = rb,
           axes = FALSE,
           legend.width = 5,
           legend.mar = 14.1,
           legend.only = TRUE,
           horizontal = TRUE,
           legend.args = list(text = expression('kg m'^-2), line = 0.5)
           )

dev.off() 


# ----------------------------------------------------------------
# Postage stamp diagrams of the paired differences between the 
# ensembles
#
# ----------------------------------------------------------------

# first, order the ensemble by the sum of
# the difference between the scenarios
diff.sum <- apply(LandCarbonEOCDiff, 1,sum, na.rm = TRUE)
ix <- sort(diff.sum, index.return = TRUE)$ix

pdf(file = 'graphics/LandCEOCpairDiffs.pdf', width = 8, height = 10)
par(mfrow = c(10,6), mar = c(0.3,0.3,0.3,0.3), oma = c(0.5, 0.5, 0.5, 0.5))

zl <- c(-4,4)

for(i in ix){

  test <-  t(mapSwap(matrix(c(LandCarbonEOCDiff[i,], recursive = TRUE), nrow = 73, ncol = 96, byrow = TRUE))[73:1, ])

  # chop the data for the image colorbar
 test[ test < zl[1]  ] <- zl[1]
  test[ test > zl[2] ] <- zl[2]

  # chop off antarctica in the plotting
  image(HadCM3.lons.wrt, HadCM3.lats[13:73], test[,13:73],
             col = br,
             xlab = '',
             ylab = '',
             zlim = zl,
        axes = FALSE
             )
  world(add = TRUE)
  
}

plot(1:10, type = 'n', axes = FALSE, xlab = '', ylab = '')
image.plot(t(test[73:1,]), zlim = zl,
        xlab="", ylab="", col = rb,
           axes = FALSE,
           legend.width = 5,
           legend.mar = 14.1,
           legend.only = TRUE,
           horizontal = TRUE,
           legend.args = list(text = expression('kg m'^-2), line = 0.5)
           )

dev.off() 




















PercMapsA1B <- mapEns(LandCarbonChangePercA1B)
PercMapsRCP26 <- mapEns(LandCarbonChangePercRCP26)

x11()
z = c(-50,50)
set.panel(2,1)
par(las = 1)
image.plot(HadCM3.lons.wrt, HadCM3.lats, t(mapSwap(PercMapsA1B$meanmap)[73:1, ]),
           col = byr,
           xlab = '',
           ylab = '',
           main = 'Land carbon mean change A1B',
           zlim = z
           )
world(add = TRUE)

image.plot(HadCM3.lons.wrt, HadCM3.lats, t(mapSwap(PercMapsRCP26$meanmap)[73:1, ]),
           col = byr,
           xlab = '',
           ylab = '',
           main = 'Land carbon mean change A1B',
           zlim = z
           )
world(add = TRUE)






LandCarbonPreIndTotal <- VegCarbonPreIndTotal + SoilCarbonPreIndTotal

LandmapsPreInd <- mapEns(LandCarbonPreIndTotal)

par(las = 1)
image(HadCM3.lons.wrt, HadCM3.lats, t(mapSwap(LandmapsPreInd$meanmap)[73:1, ]),
      col = byr,
      xlab = '',
      ylab = '')
world(add = TRUE)


zr <- range(LandmapsPreInd$meanmap, na.rm = TRUE)
x11(width = 7, height = 5)
mapBGcol(HadCM3.lons.wrt,HadCM3.lats, mapSwap(LandmapsPreInd$meanmap),
         BGcol = 'lightgrey',
         mapcol = byr,
         axis.args = list(cex.axis = 1.2),
         legend.args = list(text = 'GtC',cex.axis = 2, line = 0.5),
         zlim = zr,
         titext = 'test',
         world = TRUE,
         axes = TRUE)

# need to sort legend width

x11( )
image.plot(HadCM3.lons.wrt,HadCM3.lats, t(mapSwap(SoilmapsPreInd$meanmap)[73:1, ]), col = byr)
world(add = TRUE)

x11()
image.plot(HadCM3.lons.wrt,HadCM3.lats, t(mapSwap(VegmapsEOCA1B$meanmap)[73:1, ]), col = byr)
world(add = TRUE)

x11()
image.plot(HadCM3.lons.wrt,HadCM3.lats, t(mapSwap(VegmapsEOCRCP26$meanmap)[73:1, ]), col = byr)
world(add = TRUE)



## -----------------------------------------------------------
# Possible constraints on the land surface carbon uptake
#
# ------------------------------------------------------------

# building a linear model?


co2ppmDecadalA1B.test <- yearSelect(co2ppmDecadalA1B, seq(from = 1965, to = 2075, by = 5))

GlobTempDecadalA1B.test <- yearSelect(GlobTempDecadalA1B, seq(from = 1965, to = 2075, by = 5))

LandCarbonUptakeDecadalA1B.test <- yearSelect(LandCarbonUptakeDecadalA1B, seq(from = 1965, to = 2075, by = 5))

dat <- data.frame(co2 = c(co2ppmDecadalA1B.test, recursive = TRUE), temp = c(GlobTempDecadalA1B.test, recursive = TRUE), lu = c(LandCarbonUptakeDecadalA1B.test, recursive = TRUE) )


lmtest <- lm(lu ~ temp + co2, data = dat)




coefficients(lmtest)

# subtract out the co2 driver

coef.co2 <- lmtest$coefficients['co2']
coef.temp <- lmtest$coefficients['temp']
intercept <- lmtest$coefficients['(Intercept)']

plot( (dat$co2 * coef.co2) intercept )


plot(dat$temp, dat$lu - (dat$co2 * coef.co2)  )
plot(dat$co2, dat$lu - (dat$temp * coef.temp)  )
