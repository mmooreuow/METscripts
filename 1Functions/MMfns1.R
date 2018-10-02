##### MM's functions ---------------------------------------------------------------------
# Version 1, 24th September 2018



##### Asreml Call Writing Function -------------------------------------------------------

# Input 
# - mt: mt is a list with a component for each term to be fitted in the model. 
#       Each component is a vector of Experiment names for which the term will be fitted.


# covtest.call
# 
# not.cov <- c("lrow", "lcol", "crep", "rrep", "rrow", "rcol", "resid")
# mt.cov <- paste0('mt$',names(mt)[!names(mt)%in%not.cov])
# mt.cov <- mt.cov[-1] # removes mt$YrCon 

lintest.call <- function(mt.tmp){
  mt.tmp.fixed <- paste0('mt.tmp$',names(mt.tmp)[names(mt.tmp)%in%c('lrow','lcol')])
  fixed.call <- c('Experiment')
  if(sum(mt.tmp.fixed=="mt.tmp$lrow")==1) fixed.call <- c(fixed.call,'at(Experiment, mt.tmp$lrow):lin(Row)')
  if(sum(mt.tmp.fixed=="mt.tmp$lcol")==1) fixed.call <- c(fixed.call,'at(Experiment, mt.tmp$lcol):lin(Range)')
  
  mt.tmp.residual <- paste0('mt.tmp$resid$',names(mt.tmp$resid)[names(mt.tmp$resid)%in%c('aa','ai','ia','ii')])
  resid.call <- NULL
  if(sum(mt.tmp.residual=="mt.tmp$resid$aa")==1) resid.call <- c(resid.call,'dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt.tmp$resid$aa)')
  if(sum(mt.tmp.residual=="mt.tmp$resid$ia")==1) resid.call <- c(resid.call,'dsum(~id(Range):ar1(Row)| Experiment, levels = mt.tmp$resid$ia)')
  if(sum(mt.tmp.residual=="mt.tmp$resid$ai")==1) resid.call <- c(resid.call,'dsum(~ar1(Range):id(Row)| Experiment, levels = mt.tmp$resid$ai)')
  if(sum(mt.tmp.residual=="mt.tmp$resid$ii")==1) resid.call <- c(resid.call,'dsum(~id(Range):id(Row)| Experiment, levels = mt.tmp$resid$ii)')
  
  paste0("asreml(yield ~ ", paste0(fixed.call,collapse=' + '),",",
  "sparse=~ at(Experiment):Variety,",
  "random = ~ at(Experiment, mt.tmp$crep):ColRep + at(Experiment, mt.tmp$rrep):RowRep + at(Experiment, mt.tmp$rrow):Row + at(Experiment, mt.tmp$rcol):Range,",
  "residual = ~",  paste0(resid.call,collapse=' + '),',',
  "na.action = na.method(y='include', x='include'), data=tmpdata)")
}

# require(asreml) 
# temp.diag <- eval(parse(text = lintest.call(mt.tmp)))
  