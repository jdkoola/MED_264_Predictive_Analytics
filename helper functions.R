ApplyFigureTheme<-function(p,markingscolor = "black") {
  
  p +
    theme(
      panel.background = element_rect(fill = "transparent") # bg of the panel
      , plot.background = element_rect(fill = "transparent") # bg of the plot
      #, panel.grid.major = element_blank() # get rid of major grid
      , panel.grid.minor = element_blank() # get rid of minor grid
      , legend.background = element_rect(fill = "transparent", colour = "transparent") # get rid of legend bg
      , legend.box.background = element_rect(fill = "transparent", colour = "transparent") # get rid of legend panel bg
      , panel.border = element_blank()
    )  +
    theme(
      rect = element_rect(fill = "transparent") # bg of the panel
    ) +
    theme(axis.text=element_text(size=18, colour=markingscolor),
          axis.title=element_text(size=18,face="bold",colour=markingscolor)
          ,axis.text.x= element_text(colour=markingscolor)
          ,axis.ticks = element_line(colour=markingscolor)
          ,legend.text = element_text(colour=markingscolor)
          ,legend.title = element_text(colour=markingscolor)
    ) +
    theme(axis.line = element_line(color=markingscolor))
  
  
}
ApplyFigureThemeWithGrids<-function(p,markingscolor = "black") {
  
  p +
    theme(
     panel.background = element_rect(fill = "transparent") # bg of the panel
      , plot.background = element_rect(fill = "transparent") # bg of the plot
      #, panel.grid.major = element_blank() # get rid of major grid
     # , panel.grid.minor = element_blank() # get rid of minor grid
      , legend.background = element_rect(fill = "transparent", colour = "transparent") # get rid of legend bg
      , legend.box.background = element_rect(fill = "transparent", colour = "transparent") # get rid of legend panel bg
      , panel.border = element_blank()
    )  +
  #  theme(
  #    rect = element_rect(fill = "transparent") # bg of the panel
  #  ) +
    theme(axis.text=element_text(size=18, colour=markingscolor),
          axis.title=element_text(size=18,face="bold",colour=markingscolor)
          ,axis.text.x= element_text(colour=markingscolor)
          ,axis.ticks = element_line(colour=markingscolor)
          ,legend.text = element_text(colour=markingscolor)
          ,legend.title = element_text(colour=markingscolor)
    ) +
    theme(axis.line = element_line(color=markingscolor))
  
  
}
ApplyFigureThemeLargeFontOnly<-function(p,markingscolor = "black") {
  p +
    theme(axis.text=element_text(size=18, colour=markingscolor),
          axis.title=element_text(size=18,face="bold",colour=markingscolor)
          ,axis.text.x= element_text(colour=markingscolor)
    ) 
}
ApplyFigureThemeCalCurve <- function(p) {
  p + 
    geom_abline(intercept=0,slope=1,size=1, alpha=0.7) + theme_bw() + 
    xlab('Predicted')+ylab('Observed')+
    xlim(0,1) + ylim(0,1) + 
    theme(legend.title=element_blank())+
    theme(legend.justification=c(1,0),legend.position=c(1,0))
}
ApplyFigureThemeCalCurvePoints <- function(p) {
  ApplyFigureThemeCalCurve(p  + geom_point(size=3, color='blue') )
   
}
ApplyFigureThemeCalCurveSmooth <- function(p) {
  ApplyFigureThemeCalCurve(p  + geom_line(color='blue', size=1.3) + geom_ribbon(color=NA, fill='blue',alpha=0.3))
}

SaveStdSquareFigure <- function(p, filename) {
  ggsave(filename = filename, 
         plot=p, width=6, height=6, units='in', dpi=300)
}

SaveMedRectFigure <- function(p, filename) {
  ggsave(filename = filename, 
         plot=p, width=8, height=6, units='in', dpi=300)
}

# The following cal curve plotting functions are taken from Dr. Sharon Davis (Sharon.Davis@vanderbilt.edu)
CreateSmoothedCalCurveSmoothData <- function(predicted, observed, ci=0.95) {
  dl_adv_plot = data.frame(o = observed, p = predicted)
  if (missing(ci)) ci=0.95
  
  attach(list(), name="design.options")
  assign('dd', datadist(dl_adv_plot), pos='design.options')
  options(datadist="dd")
  
  
  glm.trend <- glm('o ~ rcs(p, 5)', dl_adv_plot, family='binomial')
  
  X_plot <- seq(0.01, .99, by=.01)
  fitVal<-predict(glm.trend, newdata=data.frame(p=X_plot),
                  se.fit=TRUE)
  
  fv<-data.frame(X=X_plot, Y=plogis(fitVal$fit),
                 LCL=plogis(fitVal$fit + qnorm((1-ci)/2)* fitVal$se.fit),
                 UCL=plogis(fitVal$fit + qnorm(ci+(1-ci)/2)* fitVal$se.fit))
  
  byObs<-predict(glm.trend, se.fit=TRUE)
  
  
  
  byObsFV<-data.frame(X=dl_adv_plot$p, Y=plogis(byObs$fit),
                      LCL=plogis(byObs$fit + qnorm((1-ci)/2)* byObs$se.fit),
                      UCL=plogis(byObs$fit + qnorm(ci+(1-ci)/2)* byObs$se.fit))
  byObsFV %>% unique()
}

CreateSmoothedCalCurvePlot <- function(predicted, observed, ci=0.95) {
  byObsFV<-CreateSmoothedCalCurveSmoothData(predicted, observed, ci)
  #ggplot(fv, aes(x = X, y = Y, ymin = LCL, ymax=UCL)) + geom_line() + geom_ribbon(alpha=0.3)
  p<-ApplyFigureThemeCalCurveSmooth(ggplot(byObsFV, aes(x = X, y = Y, ymin = LCL, ymax=UCL)))
  p
}

generateSyntheticDataFactory <- function(baseDF, method = 'syn', random_seed) {
  if (missing(random_seed)) random_seed = 12347691
  set.seed(random_seed)
  if (method=='syn') {
    brinati.syn <- syn.strata(baseDF, strata=c('SWAB'), minstratumsize = 100, m=1, k=100000)
    function(prev, sampSize) {
      brinati.syn$syn %>% 
        group_by(SWAB) %>% 
        group_modify(function(x,y) {
          posCases = round(0.1*sampSize)
          negCases = sampSize - posCases
          if (y$SWAB=='1') {
            x %>% sample_n(posCases)
          } else {
            x %>% sample_n(negCases)
          }
        })}
  } else if (method=='boot') {
    function(prev, sampSize) {
      if (missing(sampSize)) sampSize = nrow(baseDF)
      baseDF %>% 
        mutate(weight = ifelse(SWAB=='1', prev, 1-prev)) %>% 
        sample_n(size=sampSize, weight=weight, replace = TRUE)
    }
  }
}

confusionMatrix<-function(predicted, observed, cutoff=0.5) {
  cm=table(predicted > cutoff, observed) 
  obsCounts = colSums(cm)
  tstCounts = rowSums(cm)
  list(
    cm= cm, sens = cm[2,2]/obsCounts['1'],
    spec = cm[1,1]/obsCounts['0'],
    ppv = cm[2,2]/tstCounts['TRUE'],
    npv = cm[1,1]/tstCounts['FALSE']
  )
}

# The following functions taken from Dr. Sharon Davis (shar.davis@vanderbilt.edu)
assessPerf<-function(predicted, observed) {
  tryCatch( expr = {
    # adjust predictions that are exactly 0 or 1
    predicted<-as.numeric(predicted)
    predicted[predicted==0]<-.00000000001
    predicted[predicted==1]<-1-.00000000001
    if(is.factor(observed)) observed<-as.numeric(as.character(observed))
    val.prob(p = predicted, y=observed, pl=FALSE, smooth=FALSE, logistic.cal=FALSE)
  }, error = function(e) {
    print(e)
    return(rep(NA, 18))
  })
}

assessSingleTrainTest <- function(trainIdx, testIdx, df, glm.model, outcome.var.name) {
  tryCatch(
    expr = {
      glm.fit <- glm(glm.model, data = df[trainIdx, ], family='binomial')
      y.train.pred <- plogis(predict(glm.fit))
      y.test.pred <- plogis(predict(glm.fit, newdata=df[testIdx,]))
      train.perf <- assessPerf(predicted = y.train.pred, observed=df[trainIdx, outcome.var.name]) 
      test.perf <- assessPerf(predicted = y.test.pred, observed=df[testIdx, outcome.var.name]) 
      ret=list("y.train" = df[trainIdx, outcome.var.name], "y.train.pred" = y.train.pred, train.perf = train.perf,
               "y.test" = df[testIdx, outcome.var.name], "y.test.pred" = y.test.pred, test.perf = test.perf)
      ret$completd = TRUE
      ret
    },
    error = function(e) {
      print(e); 
      list("y.train" = NA, "y.train.pred" = NA, train.perf = NA,
           "y.test" = NA, "y.test.pred" = NA, test.perf = NA, completed=FALSE)
    }
  )
}
