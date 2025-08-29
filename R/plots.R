#' Plot all model parameters and proliferation effects
#'
#' Internal generic function to plot proliferation and dispersion model parameters.
#'
#' @param village A Village object.
#' @return An updated village object. Also saves a diagnostic plot.
#' @keywords internal
plot_params <- function(village) {
  UseMethod("plot_params")
}
#' @exportS3Method plot_params Village
plot_params.Village <- function(village) {

  df_all <- extract_eta(village)
  if (length(village$treatcol) != 0) {
    df_treat <- df_all[[2]]
    df_cont <- df_all[[3]]
    df_all <- df_all[[1]]
  } else {
    df_cont <- df_all[[2]]
    df_all <- df_all[[1]]
  }

  df_donor <- village$data_noT0 |>
    ungroup() |>
    select(donor, donorid, !!!syms(village$donorcols)) |>
    distinct(donor, donorid, !!!syms(village$donorcols))

  if (village$num_donorcov > 0) {

    model <- str_remove_all(village$model, "\\s*\\+?\\s*\\w*:treatment_scaled|treatment_scaled:\\w*|treatment_scaled")

    # Step 2: Clean up redundancies
    model <- str_replace_all(model, "\\+\\s*\\+", "+") # Replace double '+' with a single '+'
    model <- str_remove_all(model, "\\+\\s*$")         # Remove trailing '+'
    model <- str_remove_all(model, "^~\\s*\\+?")       # Remove leading '+' after '~'
    model <- str_trim(model)                           # Trim whitespace

    # Step 3: Ensure the formula starts with '~'
    model <- paste0("~ ", model)

    df_predictors <- model.matrix(as.formula(model), data=df_donor)
    df_predictors <- df_predictors[, -1, drop = FALSE]

    df_donor <- cbind(df_donor, df_predictors)
  }

  # df_donor <- village$data_noT0[,c('donor', 'donorid')] |> distinct(donor, donorid)
  df_all <- merge(df_all, df_donor, by = 'donorid')

  lowci <- (1-village$credinterval)/2
  highci <- 1-lowci

  # Plot betas
  plot_eta(village, lowci, highci, df_all, df_donor)

  # Plot donor covariates and return samples df for proliferation effect plot
  if(length(village$donorcov) > 0){
    df_dcov <- plot_donorcov(village, lowci, highci)[[2]]
    df_dcov$chain <- df_dcov$name |>
      str_remove("^.*:") |>
      str_remove("\\..*$")
    df_dcov <- df_dcov |>
      group_by(cov.id, cov, chain) |>
      mutate(iter = row_number())  |>
      ungroup()

    df_dcov <- df_dcov[, c('iter', 'chain', 'value', 'cov')] |>
      pivot_wider(names_from = cov, values_from = 'value')
  }

  # Plot technical variation parameters
  plot_overdispersion(village, lowci, highci)

  # Plot proliferation effects
  if (length(village$treatcol) != 0) {
    data <- plot_growthmetric(village, lowci, highci, df_treat, df_cont, df_donor, df_dcov, returnboth=TRUE)
    if (village$num_donorcov > 0) {
      village$df_lfsr_treat <- data[[3]]
    } else {
      village$df_lfsr_treat <- data[[2]]
    }
  } else {
    data <- plot_growthmetric(village, lowci, highci, df_treat=NULL, df_cont, df_donor, df_dcov, returnboth=TRUE)
  }

  if(village$num_donorcov > 0){
    village$df_lfsr <- data[[2]]
  }

  return(village)
}

#' Plot posterior predictive checks
#'
#' Internal generic function to plot posterior predictive checks.
#'
#' @param village A Village object.
#' @return An updated village object. Also saves a diagnostic plot.
#' @keywords internal
## Plot posterior predictive check
plot_ppck <- function(village) {
  UseMethod("plot_ppck")
}
#' @exportS3Method plot_ppck Village
plot_ppck.Village <- function(village) {
  if(village$ppck == 0) {
    stop('No posterior predictive checks (ppck) performed. To include ppcks set ppck = TRUE when running model')
  }

  df_ppck <- rstan::extract(village$fit, pars= village$summary_ppck$paramname[startsWith(village$summary_ppck$paramname, 'post')],
                            permuted=FALSE) |> as.data.frame()
  df_ppck <- pivot_longer(df_ppck, cols=1:ncol(df_ppck))
  df_ppck$sample <- stringr::str_remove(df_ppck$name, ".*,") |> str_remove('\\]') |> as.factor()
  df_ppck$donorid <- stringr::str_remove(df_ppck$name, '.*\\[') |> str_remove(",.*") |> as.numeric()
  df_ppck$name <- NULL

  df_donor <- village$data_noT0[,c('donor', 'donorid')] |> distinct(donor, donorid)
  df_ppck <- merge(df_ppck, df_donor, by.x='donorid')

  if (length(village$treatcol) != 0) {
    df_time <- village$data[village$data$time !=0, c('sample', 'time', 'replicate', village$treatcol)] |> distinct()
  } else {
    df_time <- village$data[village$data$time !=0, c('sample', 'time','replicate')] |> distinct()
  }

  if(min(df_time$sample) != 1) {
    df_time$sample <- df_time$sample - min(df_time$sample) +1
  }
  df_ppck <- merge(df_ppck, df_time, by = 'sample')
  df_ppck$time <- as.factor(df_ppck$time)

  lowci <- (1- village$credinterval)/2
  highci <- 1-lowci

  # Calculate samples in ppck cred interval
  if (length(village$treatcol) != 0) {
    df_ppck_stats <- df_ppck |> group_by(!!sym(village$treatcol), time, replicate, donorid, donor) |>
      summarise(mean=mean(value), median=median(value), low_ci=quantile(value, probs = lowci),
                high_ci=quantile(value, probs = highci), .groups = "keep")
  } else{
    df_ppck_stats <- df_ppck |> group_by(time, replicate, donorid, donor) |>
      summarise(mean=mean(value), median=median(value), low_ci=quantile(value, probs = lowci),
                high_ci=quantile(value, probs = highci), .groups = "keep")
  }

  df_ppckpt <- village$data[village$data$time !=0,]
  df_ppckpt$donorid <- as.factor(df_ppckpt$donorid)
  df_ppckpt$time <- as.factor(df_ppckpt$time)
  df_ppckpt$donor <- as.factor(df_ppckpt$donor)


  mergecols <- colnames(df_ppckpt)[colnames(df_ppckpt) %in% colnames(df_ppck_stats)]

  if (length(village$treatcol) != 0) {
    df_ppckpt <- merge(df_ppckpt[,c("donorid", "donor","time", 'replicate', village$treatcol, 'representation')],
                       df_ppck_stats[,c('time', 'replicate', 'donorid', "donor", 'low_ci', 'high_ci', 'mean', village$treatcol)],
                       by = mergecols)
  } else {
    df_ppckpt <- merge(df_ppckpt[,c("donorid", "donor","time", 'replicate', 'representation')],
                       df_ppck_stats[,c('time', 'replicate', 'donorid', "donor", 'low_ci', 'high_ci', 'mean')],
                       by = mergecols)
  }

  df_ppckpt$ppck <- ifelse(df_ppckpt$representation < df_ppckpt$high_ci &
                             df_ppckpt$representation > df_ppckpt$low_ci,
                           1, 0)

  village$ppck_ratio <- sum(df_ppckpt$ppck)/nrow(df_ppckpt)

  # Make plot
  if(village$sim == FALSE) {
    # Stats for plot
    if (length(village$treatcol) != 0) {
      df_ppck_stats <- df_ppck |> group_by(!!sym(village$treatcol), time, donorid, donor) |>
        summarise(mean=mean(value), median=median(value), low_ci=quantile(value, probs = lowci),
                  high_ci=quantile(value, probs = highci), .groups = "keep")
    } else{
      df_ppck_stats <- df_ppck |> group_by(time, donorid, donor) |>
        summarise(mean=mean(value), median=median(value), low_ci=quantile(value, probs = lowci),
                  high_ci=quantile(value, probs = highci), .groups = "keep")
    }
    df_ppckpt <- village$data[village$data$time !=0,]
    df_ppckpt$donorid <- as.factor(df_ppckpt$donorid)
    df_ppckpt$time <- as.factor(df_ppckpt$time)
    df_ppckpt$donor <- as.factor(df_ppckpt$donor)


    mergecols <- colnames(df_ppckpt)[colnames(df_ppckpt) %in% colnames(df_ppck_stats)]

    if (length(village$treatcol) != 0) {
      df_ppckpt <- merge(df_ppckpt[,c("donorid", "donor","time", 'replicate', village$treatcol, 'representation')],
                         df_ppck_stats[,c('time', 'donorid', "donor", 'low_ci', 'high_ci', 'mean', village$treatcol)],
                         by = mergecols)
    } else {
      df_ppckpt <- merge(df_ppckpt[,c("donorid", "donor","time", 'replicate', 'representation')],
                         df_ppck_stats[,c('time', 'donorid', "donor", 'low_ci', 'high_ci', 'mean')],
                         by = mergecols)
    }
    df_ppckpt$point <- ifelse(df_ppckpt$representation > df_ppckpt$low_ci & df_ppckpt$representation < df_ppckpt$high_ci, 19, 23)
    df_ppckpt$color[df_ppckpt$point == 19] <- 'black'
    df_ppckpt$color[df_ppckpt$point == 23] <- 'darkred'

    df_ppck_stats$donorid <- as.factor(df_ppck_stats$donorid)
    df_ppck_stats$donor <- as.factor(df_ppck_stats$donor)

    if (length(village$treatcol) != 0) {
      df_ppck_stats <- df_ppck_stats |> rename(Dose = !!sym(village$treatcol))
      df_ppckpt <- df_ppckpt |> rename(Dose = !!sym(village$treatcol))
      pltwidth <- 3 * village$num_doses
      p <- ggplot(df_ppck_stats) +
        geom_col(aes(x=mean, y=time, fill = donor, alpha = time)) +
        geom_errorbar(aes(y= time, xmin=low_ci, xmax=high_ci), width =0, lwd=0.5) +
        geom_point(data = df_ppckpt, aes(x=representation, y=time, color = color), color=df_ppckpt$color,shape = '|') +
        geom_point(data = df_ppckpt[df_ppckpt$point == 23,], aes(x=high_ci + 0.05, y=time), color = 'darkred', size=2, shape =4) +
        scale_x_continuous(labels = function(x) ifelse(x == 0, "0", x)) +
        ylab(paste0('Time (', village$timeunit, ')')) +
        theme_classic() +
        theme(axis.text = element_text(size=13),
              axis.title = element_text(size = 13),
              axis.title.x = element_blank(),
              legend.position = 'none') +
        facet_grid(rows = vars(donor), cols = vars(Dose),
                   labeller = labeller(.cols = label_both, .rows = label_parsed)) +
        theme(strip.text.y = element_text(angle = 0, margin = margin(t =5, b = 5, l= 5, r=5)),
              strip.text.x = element_text(margin = margin(r = 5, t=5, b=5)),
              strip.text = element_text(size = 13),
              panel.spacing = unit(0.005, "cm"),
              strip.background = element_rect(color = "darkgrey", linewidth = 2))
    } else {
      pltwidth <- 8

      p <- ggplot(df_ppck_stats) +
        geom_col(aes(x=mean, y=time, fill = donor, alpha = time)) +
        geom_errorbar(aes(y= time, xmin=low_ci, xmax=high_ci), width =0, lwd=0.5) +
        geom_point(data = df_ppckpt, aes(x=representation, y=time, color = color), color=df_ppckpt$color,shape = '|') +
        geom_point(data = df_ppckpt[df_ppckpt$point == 23,], aes(x=high_ci + 0.05, y=time), color = 'darkred', size=2, shape =4) +
        scale_x_continuous(labels = function(x) ifelse(x == 0, "0", x)) +
        ylab(paste0('Time (', village$timeunit, ')')) +
        theme_classic() +
        theme(axis.text = element_text(size=13),
              axis.title = element_text(size = 13),
              axis.title.x = element_blank(),
              legend.position = 'none') +
        facet_grid(rows = vars(donor),
                   labeller = labeller(.cols = label_both, .rows = label_parsed)) +
        theme(strip.text.y = element_text(angle = 0, margin = margin(t =5, b = 5, l= 5, r=5)),
              strip.text.x = element_text(margin = margin(r = 5, t=5, b=5)),
              strip.text = element_text(size = 13),
              panel.spacing = unit(0.005, "cm"),
              strip.background = element_rect(color = "darkgrey", linewidth = 2))
    }

    pltheight <- 0.2 * village$num_donors * village$num_timepts

    if(!is.null(village$color)) {
      p <- p + scale_fill_manual(values=village$color)
    }

    ggsave(paste0(village$outdir, village$name, '_ppck_barplot.png'), bg = 'white', width = pltwidth, height = pltheight, units = 'in')
  }

  return(village)
}

#' Plot effective sample size
#'
#' Internal generic function to plot effective sample size of all parameters.
#'
#' @param village A Village object.
#' @return None.
#' @keywords internal
plot_ESS <- function(village) {
  UseMethod("plot_ESS")
}
#' @exportS3Method plot_ESS Village
plot_ESS.Village <- function(village) {

  ESS <- village$summary_params[!startsWith(village$summary_params$paramname, 'beta_raw') &
                                  !startsWith(village$summary_params$paramname, 'v'),]
  pltwidth <- 0.1*length(ESS$paramname)

  p <- ggplot(ESS, aes(paramname, n_eff)) +
    geom_point() +
    theme_classic() +
    ylim(0, max(ESS$n_eff)) +
    theme(axis.text.x = element_text(angle=90),
          axis.title.x = element_blank())

  ggsave(paste0(village$outdir, village$name, '_ESS.png'), width = pltwidth, height =4, units = 'in')
}


#' Plot raw data
#'
#' Internal generic function to plot scatterplot and stacked bar plot of raw representation data.
#'
#' @param village A Village object.
#' @return None.
#' @keywords internal
plot_rawdata <- function(village) {
  UseMethod("plot_rawdata")
}
#' @exportS3Method plot_rawdata Village
plot_rawdata.Village <- function(village) {
  timeunit = village$timeunit
  data <- village$data
  if (length(village$treatcol) == 1){
    data <- data |> rename(Dose = !!sym(village$treatcol))
    plt_height <- 4 *village$num_doses
  } else {
    plt_height <- 4
    data$Dose <- 'Control'
  }

  data$time <- as.factor(data$time)
  plt_width  <- 0.2 * village$num_timepts * village$num_donors
  if (plt_width > 15) {
    plt_width <- 15
  }

  shapes <- c(15, 16, 17, 18, 19, 25, 9, 12)
  shapes <- shapes[1:(village$num_timepts +1)]

  plt <- list()
  plt[[1]] <- ggplot(data) +
    geom_point(aes(x=time, y=representation, shape=time, fill=donor), na.rm = TRUE) +
    scale_shape_manual(values= shapes) +
    scale_fill_manual(values= village$color) +
    theme_classic() +
    facet_grid(cols=vars(donor), rows=vars(Dose),
               labeller = labeller(.cols = label_parsed, .rows = label_both)) +
    xlab(paste0('Time (', village$timeunit, ')')) +
    theme(legend.position = 'none',
          axis.text = element_text(size=15),
          axis.title = element_text(size = 15),
          strip.text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          strip.text.x = element_text(angle = 90, margin = margin(t =5, b = 5, l= 100, r=100)),
          strip.background = element_rect(color = "darkgrey", linewidth = 1),
          panel.spacing = unit(0, "cm")) +
    scale_x_discrete(expand = expansion(mult = 0.3))


  ggsave(paste0(village$outdir, village$name, '_rawdata_scatter.png'), bg = 'white', width = plt_width, height = plt_height, units = 'in')

  data <- data |> group_by(Dose) |> mutate(diff= min(replicate) - 1) |> ungroup()
  data$replicate <- data$replicate - data$diff

  data <- data |> rename(!!sym(timeunit):=time)

  plt_width <- 1 * village$num_timepts * village$num_reps
  if (plt_width > 15) {
    plt_width <- 15
  }

  plt[[2]] <- ggplot(data, aes(x = replicate, y = representation, fill = donor)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values= village$color) +
    theme_classic() +
    theme(
      axis.text = element_text(size=10),
      axis.title = element_text(size = 10),
      strip.text = element_text(size= 10),
      legend.text = element_text(size=8),
      legend.title = element_text(size = 8),
      legend.key = element_blank()) +
    facet_grid(cols = vars(!!sym(timeunit)), rows= vars(Dose), labeller = label_both, scales ='free_x')

  ggsave(paste0(village$outdir, village$name, '_rawdata_stackbar.png'), bg = 'transparent', width = plt_width, height = plt_height, units = 'in')
}


#' Plot donor covariates
#'
#' Internal generic function to plot violin plot of donor covariates.
#'
#' @param village A Village object.
#' @param lowci Lower bound of credible interval
#' @param highci Higher bound of credible interval
#' @return List containing plot and data.
#' @keywords internal
plot_donorcov <- function(village, lowci, highci) {
  UseMethod("plot_donorcov")
}
#' @exportS3Method plot_donorcov Village
plot_donorcov.Village <- function(village, lowci, highci) {
  df <- rstan::extract(village$fit,
                       pars= village$summary_params$paramname[startsWith(village$summary_params$paramname, 'tau_g[')],
                       permuted=FALSE) |> as.data.frame()
  df <- pivot_longer(df, cols=1:ncol(df))
  df$cov.id <- df$name |> stringr::str_remove('.*\\[') |>  str_remove(",.*") |> as.numeric()
  df$cov <- unlist(village$donorcov[df$cov.id]) |> as.factor()

  df <- df |> group_by(cov) |> mutate(mean= mean(value),
                                      low_ci=quantile(value, probs = c(lowci)),
                                      high_ci=quantile(value, probs = c(highci)))

  pltwidth <- 0.5 + length(village$predictors)
  p <- ggplot(df) +
    geom_violin(aes(x=cov, y=value, fill=cov)) +
    geom_point(aes(x=cov, y=mean)) +
    geom_errorbar(aes(x=cov, ymin=low_ci, ymax=high_ci), width =0, lwd=0.4, , color='black') +
    ylab('Proliferation effect') +
    xlab(NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'none')

  ggsave(paste0(village$outdir, village$name, '_regressors.png'), bg = 'white', width = pltwidth, height = 4, units = 'in')

  returnlst <- list(plot = p, data = df)
  return(returnlst)
}

#' Plot overdispersion parameters
#'
#' Internal generic function to plot violin of theta, phi, and phi_r.
#'
#' @param village A Village object.
#' @param lowci Lower bound of credible interval
#' @param highci Higher bound of credible interval
#' @return None
#' @keywords internal
plot_overdispersion <- function(village, lowci, highci) {
  UseMethod("plot_overdispersion")
}
#' @exportS3Method plot_overdispersion Village
plot_overdispersion.Village <- function(village, lowci, highci) {
  df <- rstan::extract(village$fit,
                       pars= village$summary_params$paramname[startsWith(village$summary_params$paramname, 'theta')],
                       permuted=FALSE) |> as.data.frame()
  df <- pivot_longer(df, cols=1:ncol(df))
  df$name <- 'theta'
  df$mean <- mean(df$value)
  df$low_ci <- quantile(df$value, probs = c(lowci))
  df$high_ci <- quantile(df$value, probs = c(highci))


  ggplot(df) +
    geom_violin(aes(x=name, y=value, fill=name)) +
    geom_point(aes(x=name, y=mean)) +
    geom_errorbar(aes(x=name, ymin=low_ci, ymax=high_ci), width =0, lwd=0.4, , color='black') +
    scale_fill_manual(values = c('royalblue')) +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    ylab(NULL) +
    xlab(NULL) +
    theme_classic() +
    theme(axis.text = element_text(size = 15),
          axis.text.x = element_text(size=15),
          strip.text = element_blank(),
          legend.position = 'none')
  ggsave(paste0(village$outdir, village$name, '_theta.png'), bg = 'white', width = 1.5, height = 4, units = 'in')

  df <- rstan::extract(village$fit,
                       pars= village$summary_params$paramname[startsWith(village$summary_params$paramname, 'phi_r')],
                       permuted=FALSE) |> as.data.frame()
  df <- pivot_longer(df, cols=1:ncol(df))
  df$replicate <- df$name |> stringr::str_remove('.*\\[') |>  str_remove("]") |> as.numeric()
  df$name <- paste0('phi[r', df$replicate, ']')
  df <- df |> group_by(name) |> mutate(mean= mean(value),
                                       low_ci=quantile(value, probs = c(lowci)),
                                       high_ci=quantile(df$value, probs = c(highci)),
                                       id='phi_r')

  df_phi <- rstan::extract(village$fit,
                           pars= 'phi',
                           permuted=FALSE) |> as.data.frame()
  df_phi <- pivot_longer(df_phi, cols=1:ncol(df_phi))
  df_phi$name <- 'phi'

  df_phi <- df_phi |> mutate(mean=mean(df_phi$value),
                             low_ci=quantile(df_phi$value, probs = c(lowci)),
                             high_ci=quantile(df_phi$value, probs = c(highci)),
                             replicate=NA,
                             id= 'phi')

  df <- rbind(df, df_phi)

  p <- ggplot(df) +
    geom_violin(aes(x=name, y=value, fill=name, alpha = replicate), fill= 'seagreen3') +
    geom_point(aes(x=name, y=mean), size=1) +
    geom_errorbar(aes(x=name, ymin=low_ci, ymax=high_ci), width =0, lwd=0.4, , color='black') +
    labs(x = NULL) +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    ylab(NULL) +
    xlab(NULL) +
    ggh4x::facet_manual(
      vars(id),
      design = matrix(c(1, 2, 2, 2, 2), nrow = 1),
      scales = "free_x",
      labeller = 'label_parsed'
    ) +
    theme_classic() +
    theme(axis.text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_blank(),
          legend.position = 'none')

  pltwidth <- village$num_reps*0.4
  if(pltwidth < 3) {
    pltwidth <- 3
  }
  ggsave(paste0(village$outdir, village$name, '_phi_r.png'), bg = 'white', width = village$num_reps*0.5, height = 4, units = 'in')
}

## Plot violin of donor proliferation effect per treatment dose
#' Plot total proliferation effect
#'
#' Internal generic function to plot violin of total donor proliferation effect per treatment dose.
#'
#' @param village A Village object.
#' @param lowci Lower bound of credible interval
#' @param highci Higher bound of credible interval
#' @param df_treat Optional data frame containing treatment parameter draws.
#' @param df_cont Data frame containing control intercept draws.
#' @param df_donor Data frame with distinct donor ids/names.
#' @param df_dcov Donor covariate parameter draws.
#' @param returnboth If returnboth = TRUE A list of data and plots. Data includes lfsr calculations for treatment/donor covariate parameters.
#' @return List of data (lfsr) and plot.
#' @keywords internal
plot_growthmetric <- function(village, lowci, highci, df_treat=NULL, df_cont, df_donor, df_dcov, returnboth=FALSE) {
  UseMethod("plot_growthmetric")
}
#' @exportS3Method plot_growthmetric Village
plot_growthmetric.Village <- function(village, lowci, highci, df_treat=NULL, df_cont, df_donor, df_dcov, returnboth=FALSE) {
  pltwidth <- 0.2 * village$num_donors
  if (!is.null(df_treat)) {

    df_treat$control <- df_cont$value
    df_treat$chain <- rep(1:village$mcmc_chains, times=(village$mcmc_samples-village$mcmc_warmup)*village$num_donors)
    df_treat <- df_treat |>
      group_by(donorid, chain) |>
      mutate(iter = row_number())  |>
      ungroup()

    doses_scaled <- unique(village$data_noT0$treatment_scaled)
    doses <- village$doses[order(village$doses)]


    if(village$num_donorcov > 0){

      df_treat <- merge(df_treat, df_dcov, by=c('iter', 'chain'))

      df_lfsr <- df_treat |>
        distinct(iter, across(all_of(village$donorcov))) |>
        summarise(across(all_of(village$donorcov),
                         ~min(mean(. >= 0), mean(. <= 0)),
                         .names = "lfsr_{.col}"))

      colnames(df_treat) [7:ncol(df_treat)] <- paste0(colnames(df_treat)[7:ncol(df_treat)], '_pred')

      dcov <- colnames(df_donor[, -c(1:2, which(colnames(df_donor) %in% village$donorcols)), drop = FALSE])

      # check for treatment interaction term and extract donor cov column
      interaction_term <- village$predictors[str_detect(village$predictors, ':.*treatment_|treatment_.*:')]

      if (length(interaction_term) > 0){
        remove_term <- str_extract(interaction_term, "(?<!treatment_)\\w+(?=:)|(?<=:)(?!treatment_)\\w+")
        dinteract <- interaction_term |> str_remove(remove_term) |> str_remove(':')
      }
    }
    df_treat <- merge(df_treat, df_donor)

    df_lfsr_treat <- df_treat |> group_by(donor, donorid) |> summarise(lfsr=min(mean(value >= 0), mean(value <= 0)))

    # calculate proliferation effect by treatment dose
    for(d in 1:village$num_doses) {
      df_treat[, as.character(doses[d])] <- df_treat$control + df_treat$value * doses_scaled[d]
      if(village$num_donorcov > 0){
        for(dc in 1:length(dcov)) {
          df_treat[, as.character(doses[d])] <- df_treat[, as.character(doses[d])] + df_treat[,dcov[dc]] * df_treat[,paste0(village$donorcols[dc], '_pred')]
        }
      if (length(interaction_term) > 0){
        df_treat[, as.character(doses[d])] <- df_treat[, as.character(doses[d])] + df_treat[,dcov[str_detect(dcov, dinteract)]] * df_treat[,paste0(interaction_term, '_pred')] * doses_scaled[d]
      }
      }
    }
    df_treat$Treatment <- village$treatment
    df_treat <- df_treat |> pivot_longer(cols = as.character(village$doses), names_to = 'Dose', values_to = 'growthrate')
    df_treat$Dose <- as.numeric(df_treat$Dose)
    df_treat <- df_treat[order(df_treat$Dose),]

    df_beta_stats <- df_treat |> group_by(donor, donorid, Treatment, Dose) |> summarise(mean=mean(growthrate), median=median(growthrate),
                                                                                          low_ci=quantile(growthrate, probs = c(lowci)),
                                                                                          high_ci=quantile(growthrate, probs = c(highci)),
                                                                                          .groups = "keep")
    pltheight <- 2 * village$num_doses

    p <- ggplot(df_treat) +
      geom_violin(aes(x=donor, y=growthrate, fill=donor), ) +
      geom_point(data=df_beta_stats, aes(x=donor, y=mean), size=1) +
      geom_errorbar(data=df_beta_stats, aes(x=donor, ymin=low_ci, ymax=high_ci), width =0, lwd=0.4, , color='black') +
      geom_hline(yintercept = 0, linetype = 'FF', lwd = 0.1, color='black') +
      ylab('Total proliferation effect') +
      theme_classic() +
      facet_grid(rows = vars(Dose), cols=vars(Treatment), labeller= label_both) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')

    df <- df_treat

    if(is.null(village$color) == FALSE) {
      p <- p + scale_fill_manual(values= village$color)
    }
  } else {
    if(length(village$predictors) > 0){

      df_cont$chain <- rep(1:village$mcmc_chains, times=(village$mcmc_samples-village$mcmc_warmup)*village$num_donors)
      df_cont <- df_cont |>
        group_by(donorid, chain) |>
        mutate(iter = row_number())  |>
        ungroup()
      df_cont <- merge(df_cont, df_dcov, by=c('iter', 'chain'))

      df_lfsr <- df_cont |>
        distinct(iter, across(all_of(village$predictors))) |>
        summarise(across(all_of(village$predictors),
                         ~min(mean(. >= 0), mean(. < 0)),
                         .names = "lfsr_{.col}"))

      colnames(df_cont) [6:ncol(df_cont)] <- paste0(colnames(df_cont)[6:ncol(df_cont)], '_pred')

    }

    df_cont <- merge(df_cont, df_donor, by = 'donorid')
    df_cont$growthrate <- df_cont$value
    if(length(village$predictors) > 0){
      predcols <- colnames(df_cont)[endsWith(colnames(df_cont), '_pred')]
      designcols <- colnames(df_donor)[-1:-(ncol(df_donor) - village$num_donorcov)]
      for(p in 1:length(village$predictors)) {
        df_cont$growthrate <- df_cont$growthrate + df_cont[,designcols[p]] * df_cont[,predcols[p]]
      }
    }

    df_beta_stats <- df_cont |> group_by(donor, donorid) |> summarise(mean=mean(growthrate),
                                                                        median=median(growthrate),
                                                                        low_ci=quantile(growthrate, probs = c(lowci)),
                                                                        high_ci=quantile(growthrate, probs = c(highci)),
                                                                        .groups = "keep")
    donor_order <- unique(df_beta_stats$donor[order(df_beta_stats$mean)])
    df_cont$donor <- factor(df_cont$donor, levels=donor_order)

    df_beta_stats$donor <- factor(df_beta_stats$donor, levels=donor_order)

    pltheight <- 4
    p <- ggplot(df_cont) +
      geom_violin(aes(x=donor, y=growthrate, fill=donor), ) +
      geom_point(data=df_beta_stats, aes(x=donor, y=mean), size=1) +
      geom_errorbar(data=df_beta_stats, aes(x=donor, ymin=low_ci, ymax=high_ci), width =0, lwd=0.4, , color='black') +
      geom_hline(yintercept = 0, linetype = 'FF', lwd = 0.1, color='black') +
      ylab('Total proliferation effect') +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')

    df <- df_cont

    if(is.null(village$color) == FALSE) {
      p <- p + scale_fill_manual(values= village$color)
    }
  }
  ggsave(paste0(village$outdir, village$name, '_growthrates.png'), plot = p, bg = 'white', width = pltwidth, height = pltheight, units = 'in')

  if(returnboth == TRUE) {
    if(exists('df_lfsr_treat')) {
      if(village$num_donorcov > 0) {
        return(list(p, df_lfsr, df_lfsr_treat))
      } else {
        return(list(p, df_lfsr_treat))
      }
    } else {
      if(village$num_donorcov > 0) {
        return(list(p, df_lfsr))
      } else {
        return(list(p))
      }
    }
  } else {
    return(p)
  }

}

#' Plot donor beta parameters per treatment
#'
#' Internal generic function to plot violin plot of beta parameters per treatment.
#'
#' @param village A Village object.
#' @param lowci Lower bound of credible interval
#' @param highci Higher bound of credible interval
#' @param df_all Data frame of draws.
#' @param df_donor Donor identification info
#' @return None
#' @keywords internal
plot_eta <- function(village, lowci, highci, df_all, df_donor) {
  UseMethod("plot_eta")
}
#' @exportS3Method plot_eta Village
plot_eta.Village <- function(village, lowci, highci, df_all, df_donor) {
  df_beta_stats <- df_all |> group_by(donor, donorid, name) |> summarise(mean=mean(value), median=median(value),
                                                                           low_ci=quantile(value, probs = c(lowci)),
                                                                           high_ci=quantile(value, probs = c(highci)),
                                                                           .groups = "keep")
  donor_order <- df_beta_stats$donor[df_beta_stats$name == 'beta' & order(df_beta_stats$mean)]
  donor_order <- df_beta_stats$donor[df_beta_stats$name == 'beta'][order(df_beta_stats$mean[df_beta_stats$name == 'beta'])]


  df_all$donor <- factor(df_all$donor, levels=donor_order)

  df_beta_stats$donor <- factor(df_beta_stats$donor, levels=donor_order)

  df_beta_stats$point <- 23
  df_beta_stats$point[df_beta_stats$low_ci < 0 & df_beta_stats$high_ci >0 | df_beta_stats$name == 'beta'] <- 16

  pltwidth <- 0.2*village$num_donors

  p <- ggplot(df_all) +
    geom_violin(aes(x=donor, y=value, fill=donor)) +
    geom_point(data=df_beta_stats, aes(x=donor, y=mean, shape=as.factor(point)), shape=df_beta_stats$point, size=1) +
    geom_errorbar(data=df_beta_stats, aes(x=donor, ymin=low_ci, ymax=high_ci), width =0, lwd=0.4, , color='black') +
    geom_hline(yintercept = 0, linetype = 'FF', lwd = 0.1, color='black') +
    ylab('Proliferation effect') +
    theme_classic() +
    facet_grid(rows = vars(name), labeller= label_parsed) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')

  if(is.null(village$color) == FALSE) {
    p <- p + scale_fill_manual(values= village$color)
  }
  ggsave(paste0(village$outdir, village$name, '_eta_params.png'), bg = 'white', width = pltwidth, height = 4, units = 'in')
}
