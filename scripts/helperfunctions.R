library(bayesplot)
library(gridExtra)
fit_model_to_df <- function(df, condition_col, id_cols, model_path=NA, model_code=NA, 
                            value_col, variable_col,
                            ignore_cols = c(), nchains=5,DROP=F, impute=F, mask=NA,
                            debug=FALSE){
  # does the dirty work of getting a tall dataframe ready for analysis
  # assumes that the values in the `value_col` are numeric (or can properly be coerced)
  # mask is used to subset the dataframe
  # if not "drop", NA values will be imputed with Rstan
  if (length(mask) == nrow(df)){
    masked_df <- df[mask,]
  }
  masked_df <- masked_df[, !colnames(masked_df) %in% ignore_cols]
  # pad numbers by X_0 to make character vector easier to access 
  # masked_df[,variable_col] <- sprintf("X_%04d", masked_df[, variable_col])
  
  # we now cast to wide.  If there is replication, here we take a simple mean aggregation, as opposed to dcast's defaulting to length
  # you dont want to see the posteriors when you model after defaulting to length; some things you cant unsee
  wide_df <- dcast(
    fun.aggregate = mean,
    masked_df, 
    formula = paste(paste(id_cols, "+", collapse=" "), condition_col,  "~", variable_col), 
    value.var = value_col)
  if (DROP){  # for some reason, I cant get the cast call to drop NAs...
    print("dropping missing data rows")
    wide_df <- na.omit(wide_df)
  }
  missing_wide_df <- missing_data.frame(wide_df)
  if(any(is.na(wide_df))){
    show(missing_wide_df)
    image(missing_wide_df)
    if (impute){
      print("'mi' imputation depreciated")
      #wide_df <- auto_impute()
    }   
  }
  datacols <- colnames(wide_df)[!colnames(wide_df) %in% c(condition_col, id_cols)]
  x <- as.numeric(datacols)
  yobs <- as.matrix(wide_df[, datacols])  
  ymiss <- ifelse(is.na(yobs), 1, 0)
  ymiss_vector <- as.numeric(is.na(t(yobs))) # we traverse rowwise later so get the tranpose of the matrix
  # now we know where the NAs r, so we can replace them with 9999 and ignore them so it doesnt break stans data transformation
  yobs[is.na(yobs)] <- 9999
  conditions = factor(wide_df[, condition_col])
  C = length(unique(conditions))
  N = nrow(yobs)
  nT = ncol(yobs)
  data <- list(N=N, nT=nT, C=C, conditions=as.numeric(conditions), x=x, yobs=yobs, ymiss=ymiss, ymiss_vector=ymiss_vector)
   if (debug) {
    fit = list()
  }else{
    if (is.na(model_path)){
      fit = stan(model_code=model_code, data=data)
    } else{
      fit = stan(model_path, data=data)
    }
  }
  return(list("data"=data, 
              "fit"=fit, 
              "missing"=missing_wide_df,
              "conditions"=conditions))
  
}
basic_fit_plots <- function(fit, condition_names=c(), ignored_params = c("y_hat")){
  
  # ignore y_hat, sigma and log posterior for now
  fit_df <- as.data.frame(summary(fit)$summary)
  fit_df$param <- row.names(fit_df)
  ignored_cols <- rep(TRUE, length(fit_df$param))
  for (param in ignored_params){
    ignored_cols <- ifelse(grepl(param, fit_df$param), FALSE, ignored_cols)
  }
  fit_df <- fit_df[ignored_cols, ]
  
  global_pars = fit@model_pars[fit@model_pars %in% names(fit)]
  sample_pars = fit@model_pars[!fit@model_pars %in% names(fit)]
  sample_pars <- sample_pars[!sample_pars %in% ignored_params]
  rhats <- fit_df$Rhat
  names(rhats) <- fit_df$param
  plot_params = fit@model_pars[!fit@model_pars %in% ignored_params]
  nparams <- length(plot_params)
  for (param in plot_params){
    print(paste0("making plot ", match(param, fit@model_pars), " of ", nparams))
    these_rhats <- rhats[grepl(paste0("^",param),  names(rhats))]

    traces <- plot(fit, plotfun = "trace", pars=param)
    rhatplot <- mcmc_rhat(these_rhats) +
      yaxis_text() 
    p <- stan_plot(fit, show_density = TRUE, ci_level = 0.8, pars = c(param)) +
      labs(title=paste("Parameter:", param))
    
    if (length(condition_names) == 0 | param %in% global_pars){
      print("skipping renaming parameter")
    } else{
      p <- p + scale_y_discrete(limits=sort(unique(as.character(condition_names)),decreasing = T))
    }
    grid.arrange(traces, rhatplot, p, 
                 layout_matrix = rbind(c(1,3),
                                       c(1,3),
                                       c(2,3)))
  }
}
basic_fit_plots_pp <- function(fit, pooled_params=c(), condition_names=c()){
  # ignore sigma and log posterior for now
  global_pars = fit@model_pars[fit@model_pars %in% names(fit)]
  # get the params that we have one per condition by ignoring the pooled params (ie, biorep or techrep level)
  sample_pars = fit@model_pars[!fit@model_pars %in% c(pooled_params, names(fit))]
  rhats <- fit %>%
    rhat(pars = fit@model_pars)
  names(rhats) <- names(fit)
  
  for (param in fit@model_pars){
    # traces
    traces <- plot(fit, plotfun = "trace", pars=param, inc_warmup = TRUE)
    rhats <- fit %>%
      rhat(pars = param)
    rhatplot <- mcmc_rhat(rhats) +
      yaxis_text()
    p <- stan_plot(fit, show_density = TRUE, ci_level = 0.5, pars = c(param)) +
      labs(title=paste("Parameter:", param))
    
    if (length(condition_names) == 0 | !param %in% sample_pars){
      print(paste("skipping renaming for parameter", param))
    } else{
#      combined <- rep(unique(condition_names), length(sample_pars))
      # names(fit) <- c(
      #   rep(unique(condition_names), length(sample_pars),
      #   rep(unique(as.character(condition_names)), length(pooled_params))
      #   global_pars
      # )
      p <- p + scale_y_discrete(limits=sort(unique(condition_names),decreasing = T))
    }
    grid.arrange(traces, rhatplot, p, 
                 layout_matrix = rbind(c(1,3),
                                       c(1,3),
                                       c(2,3)))
  }
}


fit_model_to_df_pp <- function(df, condition_col, id_cols, model_path=NA, model_code=NA, 
                            value_col, variable_col, nchains=5,DROP=F, impute=F, mask=NA,
                            debug=FALSE){
  # Does the same as above, but doesnt bother casting to wide. "pp" means partial pooling,
  # as here we get estimates for the parameters per biologial replicate, but also estimate  
  #  error for both the biological replication and the technical replication 
  # 
  # assumes that the values in the `value_col` are numeric (or can properly be coerced)
  # mask is used to subset the dataframe
  # if not "drop", NA values will be imputed with Rstan
  if (length(mask) == nrow(df)){
    masked_df <- df[mask,]
  } else {
    masked_df <- df
  }
    
  # we often filter out entire treatments;  we we like factors, but dont want factor 
  # levels that are no longer presesnt, lets get rid of the now-unused levels
  masked_df[, condition_col] <- droplevels(masked_df[, condition_col])
  # we drop missing values here; no sense keeping them in
  masked_df <-masked_df[!is.na(masked_df[,value_col]),]
  data <- list(N=nrow(masked_df),
               x=masked_df[, variable_col],
               y=masked_df[, value_col],
               biorep=as.numeric(masked_df[, "biorep"]), 
               techrep=as.numeric(masked_df[, "techrep"]),
               nB=length(unique(masked_df[, "biorep"])),
               nT=length(unique(masked_df[, "techrep"])),
               nC=length(unique(masked_df[, condition_col])), 
               condition=as.numeric(masked_df[, condition_col]))
  if (debug) {
    fit = list()
  }else{
    if (is.na(model_path)){
      fit = stan(model_code=model_code, data=data)
    } else{
      fit = stan(model_path, data=data)
    }
  }
  return(list("data"=data, 
              "fit"=fit, 
              "conditions"=masked_df[, condition_col]))
  
}
plot_all_fits <- function(fits, param, condition_names){
  
  plots <- list()
  for (f in 1:length(fits)){
    plotcol = ifelse("m" %in% strsplit(param, split = "")[[1]] ,"#B2001D" , "#0083B2")
    p <- plot(fits[[f]], show_density = TRUE, ci_level = 0.8, pars = c(param), fill_color=plotcol)
    if (f == 1){
       p <- p + labs(title=paste(names(fits)[f]), subtitle=paste("parameter:", param)) + 
        scale_y_discrete(limits=sort(unique(as.character(condition_names)),decreasing = T)) + theme(panel.grid=element_line(color="gray10"))
        
    } else{
      p <- p + labs(title=paste(names(fits)[f]), subtitle=" ") + theme(axis.text.y = element_blank(), panel.grid=element_line(color="gray10"))
    }
    plots <- c(plots, list(p))
  }
  do.call(grid.arrange, c(plots, ncol = length(fits)))
  # layout_matrix = rbind(c(1,3),
  #                       c(1,3),
  #                       c(2,3)))
  
}

