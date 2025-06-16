# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @importFrom stringr str_remove str_remove_all str_replace_all str_trim str_detect str_extract
#' @importFrom rlang .data sym syms
#' @importFrom stats median quantile as.formula
#' @importFrom utils read.csv
#' @importFrom grid unit
#' @importFrom methods representation
#' @importFrom stats cov model.matrix time
#' @importFrom data.table :=

## usethis namespace: end
NULL


#' Initiate village experiment data
#' @description
#' Formats users village composition data and checks that required columns () are present to run townlet analysis. Plots raw data for inspection and sets baseline donor.
#'
#'
#' @param datapath File path to census composition data (see github tutorial for df format)
#' @param model Regression model including all optional donor covariates and treatments (e.g. ~sex + treatment_lead + sex:treatment_lead). If no covariates/treatments included in analysis then model=NULL.
#' @param normalize Force normalization of representation so that it sums to 1.
#' @param outdir Out directory to save results
#' @param name Name for output directory/files
#' @param cls Village donor colors to use in plots
#' @param T0_cutoff Threshold for excluding under represented donors at T0.
#' @param alldonors Include all donors? If False, normalization will be forced.
#' @param sim If using townlet to run simulations model fit will be saved, but no plots will be generated.
#' @param timeunit Time unit for samples (e.g., Days, Cell passages, etc.)
#' @param baseline Choose baseline donor (warning: see github tutorial before overriding default)
#' @param ebayes Use empirical bayes method to set dispersion prior variance
#'
#' @returns Village object with correctly formatted data for running townlet donor proliferation analysis
#' @export
#'
#' @examples
#' \dontrun{
#' test_village <- init_village(datapath='./village_census.csv',
#'   model= '~sex + treatment_lead + sex:treatment_lead',
#'   name='test_townlet', outdir='./')
#'}
init_village <- function(datapath, model=NULL, normalize=FALSE, outdir= './', name='test_townlet', cls=NULL, T0_cutoff=NULL, alldonors=TRUE, sim=FALSE, timeunit='days', baseline=NULL, ebayes=TRUE) {
  print('Initiating village object')
  village <- Village(datapath, model, normalize, outdir, name, cls, T0_cutoff, alldonors, sim, timeunit, baseline, ebayes)

  print('Checking data structure')
  village <- validate_village(village)



  if(is.null(village$color)) {
    cls <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
    village$color <- rep(cls, length.out = village$num_donors)
  }

  if(village$sim == FALSE) {
    print('Plotting raw data')
    plot_rawdata(village)
  }
  return(village)
}

#' Helper function that initiates village object
#' @keywords internal
Village <- function(datapath, model, normalize=FALSE, outdir= './', name='testmodel', cls, T0_cutoff, alldonors=TRUE, sim=FALSE, timeunit, baseline=NULL, ebayes) {
  structure(list(
    datapath=datapath,
    model=model,
    normalize=normalize,
    outdir=outdir,
    name=name,
    color=cls,
    T0_cutoff=T0_cutoff,
    alldonors = alldonors,
    sim = sim,
    timeunit = timeunit,
    baseline = baseline,
    ebayes = ebayes
  ), class = 'Village')
}

#' Check donor time-0 representation
#'
#' Internal generic function to check donor representation at time 0.
#'
#' @param village A Village object.
#' @return An updated village object. Also saves a diagnostic plot.
#' @keywords internal
t0_checkrep <- function(village){
  UseMethod('t0_checkrep')
}
#' @exportS3Method t0_checkrep Village
t0_checkrep.Village <- function(village) {

  t0 <- village$data[village$data$time == 0,]

  village$donors <- unique(village$data$donor)
  village$num_donors <- length(village$donors)
  print(paste('Number donors: ', village$num_donors))

  equalrep <- 1/village$num_donors
  if (is.null(village$T0_cutoff)) {
    village$T0_cutoff <- equalrep/2
  }

  if(village$sim == FALSE) {
    ggplot(t0) +
      geom_histogram(aes(x=representation), bins=30) +
      geom_vline(aes(xintercept=village$T0_cutoff), color = 'blue',  linetype ='dashed') +
      annotate("text", x = village$T0_cutoff - village$T0_cutoff/10, y = Inf, label = "Low representation cutoff", angle = 90, vjust = 0, hjust=1) +
      geom_vline(aes(xintercept=equalrep), color = 'purple',  linetype ='dashed') +
      annotate("label", x = equalrep, y = 0.5, label = paste("1/", village$num_donors, 'donors =', round(1/village$num_donors,4)), fill='snow2', angle = 0, vjust = 0, hjust=0.5) +
      theme_classic() +
      ylab(NULL) +
      xlab('T0 representation')

    ggsave(paste0(village$outdir, village$name, '_T0donor_representation.png'), bg = 'white', width = 5, height = 4, units = 'in')
  }


  village$num_T0lowdonors <- length(village$T0_lowdonors)

  if (village$num_T0lowdonors > 0) {
    village$T0_lowdonors <- unique(t0$donor[t0$representation < village$T0_cutoff])
    if(village$alldonors == FALSE) {
      donors <- donors[!(donors %in% village$T0_lowdonors)]
      village$num_donors <- length(donors)
      village$data <- village$data[village$data$donor %in% donors,]
      village$normalize <- TRUE

      warning(paste(village$num_T0lowdonors, 'donors have representation below the cutoff at T0:', paste(village$T0_lowdonors, collapse = ', '),
                    "\nRemoving these donors from village analysis (New donor number: ", village$num_donors, "). \nYou can set a new cutoff threshold (T0_cutoff = 1/2*(1/# donors)) during village initiation.
                  \nSee plots: ",
                    '\n',paste0(village$outdir, village$name, '_T0donor_representation.png'),
                    '\n',paste0(village$outdir, village$name, '_T0donor_zscore.png')))
    }
  } else {
    if (village$num_T0lowdonors > 0 & village$alldonors == TRUE) {
      warning(paste(village$num_T0lowdonors, 'donors have representation below the cutoff at T0:', paste(village$T0_lowdonors, collapse = ', '),
                    "\nTo remove these low representation donors set option 'alldonors = FALSE' when initiating village. \nYou can also set a new cutoff threshold (T0_cutoff = 1/2*(1/# donors)) during village initiation.
                \nSee plots: ",
                    '\n',paste0(village$outdir, village$name, '_T0donor_representation.png'),
                    '\n',paste0(village$outdir, village$name, '_T0donor_zscore.png')))
    } else{
      if (village$num_T0lowdonors == 0) {
        print('All donors met the T0 representation threshold.')
      }
    }
  }
  return(village)
}

#' Checks data and model structure
#'
#' Internal generic function that checks that all init_village() inputs are in correct format.
#'
#' @param village A Village object.
#' @return An updated village object. Also saves a diagnostic plot.
#' @keywords internal
validate_village <- function(village) {
  UseMethod("validate_village")
}
#' @exportS3Method validate_village Village
validate_village.Village <- function(village) {

  ## Check outdir and name are type character
  ## Create outdir
  if (!is.character(village$outdir) || length(village$outdir) != 1) {
    stop("Outdir attribute must be a single string", call. = FALSE)
  } else {
    village$outdir <- paste0(village$outdir, village$name, '/')
    if (dir.exists(village$outdir) == TRUE) {
      warning(paste0('Output files will be saved to existing directory (all existing files will be overwritten): ', village$outdir))
    } else {
      dir.create(village$outdir, recursive = TRUE)
      warning(paste('Creating out file directory:', paste0(village$outdir, village$name, '/')))
    }
  }
  if (!is.character(village$name) || length(village$name) != 1) {
    stop("Name attribute must be a single string", call. = FALSE)
  }

  ## Read in df
  village$data <- read.csv(village$datapath)

  ## check if user defined donor covariate(s) or treatment column to include in model
  if(!is.null(village$model)) {
    required_cols <- unlist(strsplit(village$model, "[~+]"))
    required_cols <- gsub("\\s", "", required_cols)
    required_cols <- unique(required_cols[required_cols != ""])
    village$predictors <- required_cols
    required_cols <- required_cols[!grepl(":", required_cols)]
  } else {
    required_cols <- c()
  }

  # check if all required columns are in df input
  required_cols <- c(c('donor', 'time', 'replicate', 'representation'), required_cols)

  missing_cols <- setdiff(required_cols, colnames(village$data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required column(s):",
               paste(missing_cols, collapse = ", "),
               "Double check dataframe columns and input model structure"),
         call. = FALSE)
  }

  ## check for treatments in village
  if(!is.null(village$predictors)) {
    village$treatcol <- village$predictors[startsWith(village$predictors, 'treatment_') & !grepl(':', village$predictors)]
  }

  if (length(village$treatcol) > 1) {
    stop(paste("Please provide only one treatment column per analysis with 'treatment_' column name prefix and treatment dose numeric row values (e.g. 0, 1, 5, 7.5).
               Current", length(village$treatcol), "columns provided:",
               paste(village$treatcol, collapse = ", ")),
         call. = FALSE)
  } else {
    if(length(village$treatcol) == 1) {
      if(!(is.numeric(village$data[[village$treatcol]]) || is.integer(village$data[[village$treatcol]]))) {
        stop(paste("The following column(s) need to be numeric or integer:",
                   paste(village$treatcol)),
             call. = FALSE)
      } else {
        #%%%% Add number doses, doses and treatment name attributes %%%%
        ### Add treatment_scaled columnn
        village$num_doses <- length(unique(village$data[[village$treatcol]]))
        village$doses <- unique(village$data[[village$treatcol]])
        village$treatment <- str_remove(village$treatcol, '.*\\_')
        treatment_scaled <- village$data[, village$treatcol] / max(village$data[, village$treatcol])
        village$data$treatment_scaled <- treatment_scaled
        warning(paste(village$treatment, 'treatment included in village experiment'))
      }
    } else {
      if(length(village$treatcol) == 0) {
        warning("No treatment included in village experiment. To include a treatment in analysis, \n provide only one treatment column containing treatment doses (e.g. 0, 1, 5, 7.5) with 'treatment_' \n column name prefix and add to defined model with optional interactions (e.g. '~sex + treatment_lead + sex:treatment_lead').")
      }
    }
  }

  ## check that all donors meet minimum representation threshold at T0
  village <- t0_checkrep(village)

  ## donor name column and add donorid column
  if (!is.character(village$data$donor)) {
    stop("The 'donor' column must be type character", call. = FALSE)
  }

  ## check for donor covariates included in data
  #%%%% Add donor covariate(s) names and total number to attributes %%%%
  if(!is.null(village$predictors)) {
    village$donorcov <- village$predictors[!(village$predictors %in% village$treatcol)]
    village$donorcols <- village$donorcov[!grepl(":", village$donorcov)]
  }

  village$num_donorcov <- length(village$donorcov)

  if(village$sim == FALSE) {
    print('Setting baseline donor')
    village <- baseline_donor(village)
  }

  if (village$num_donorcov > 0) {
    warning(paste('Donor covariate(s) included in analysis:', paste(village$donorcov, collapse = ', ')))
    cols <- c('donor', 'donorid', village$donorcols)
    village$df_dcov <- village$data |> dplyr::distinct(across(all_of(cols)))
  } else {
    if(village$num_donorcov == 0) {
      warning("No donor covariates included in village analysis. To include donor covariates,
      provide a desired model with optional interaction terms (e.g. 'model = ~sex + ancestry + sex:ancestry")
    }
  }

  ### time column (scale for running model)
  if (!(is.numeric(village$data$time) || is.integer(village$data$time))) {
    stop("The 'time' column must be type numeric", call. = FALSE)
  }
  village$data$time_scaled <- village$data$time/(max(village$data$time)/length(unique(village$data$time[village$data$time > 0])))
  warning(paste("New scaled time values for modeling growth rates:",
                paste(sort(unique(village$data$time_scaled)), collapse = ", ")))
  #%%%% Add number sample time points (excluding T0) attribute %%%%
  village$num_timepts <- length(unique(village$data$time_scaled[village$data$time_scaled > 0]))

  ### replicate column
  if (!all(village$data$replicate %% 1 == 0)) {
    stop("The 'replicate' column must be type integer", call. = FALSE)
  }
  #%%%% Add number replicates attribute %%%%
  village$num_reps <- max(village$data$replicate)

  ## Check if treatment time points have same number of replicates
  groupcols <- c('time', 'donor', village$treatcol)
  village$ckreps <- village$data |> group_by(!!!syms(groupcols)) |> summarise(reps = length(replicate), .groups = "keep")

  if(length(village$treatcol) == 0) {
    ckval <- village$num_reps
  } else{
    ckval <- village$num_reps/(village$num_doses)
    if(!(ckval %% 1 == 0)){
      warning('Not all treatment have the same number of replicates')
    }
  }
  if(any(village$ckreps$reps != ckval)){
    warning('Different number of replicates per treatment or sample time points. Inspect number of replicates for each sample by calling View(village$ckreps)')
  }

  ## arrange columns by time, replicate & treatment group
  cols <- c('time', 'replicate', 'donorid', village$treatcol)
  village$data <- village$data |> arrange(!!!syms(cols))
  rownames(village$data) <- NULL

  ### Add sample column
  groupcols <- c('time_scaled', 'replicate', village$treatcol)
  village$data <- village$data |> group_by(!!!syms(groupcols)) |> mutate(sample=cur_group_id())

  ### donor representation numeric and sums to 1
  if (!is.numeric(village$data$representation)) {
    stop("The 'representation' column must be type numeric", call. = FALSE)
  }

  village$data <- village$data |> group_by(sample) |> mutate(sum=sum(representation)) |> ungroup()
  if (village$normalize == TRUE | village$alldonors == FALSE) {
    village$data$representation <- village$data$representation/village$data$sum
    warning('Forcing normalization of donor representation data.')
  } else {
    if (any(village$data$sum != 1) & any(village$data$sum > 0.998)) {
      warning("Sample representation proportions need to sum to 1. Since samples sums are close to 1 (0.998 < sum < 1), performing sample normalization")
      village$data$representation <- village$data$representation/village$data$sum
    }
    if (any(village$data$sum != 1) & any(village$data$sum < 0.998)) {
      stop('Sample representation proportions do not all sum to 1. Double check data or set normalize option to TRUE to force normalization')
    }
  }
  village$data <- as.data.frame(village$data)
  village$data$sum <- NULL

  return(village)
}

#' Sets baseline donor
#'
#' Internal generic function that either sets a recommended default or user provided baseline donor. Function also recommends possible alternative donors that have median proliferation rates in the village experiment.
#'
#' @param village A Village object.
#' @return An updated village object.
#' @keywords internal
baseline_donor <- function(village) {
  UseMethod("baseline_donor")
}
#' @exportS3Method baseline_donor Village
baseline_donor.Village <- function(village) {

  # Define donors with median growth (alternative baseline donors)
  df_meddonors <- village$data

  ### Fix this later
  if (length(village$treatcol) == 1) {

    df_meddonors <- df_meddonors |>
      filter(.data[[village$treatcol]] == 0)

    result <- df_meddonors |>
      group_by(replicate) |>
      summarise(unique_times = list(unique(time))) |>
      ungroup() |>
      summarise(common_min = max(sapply(unique_times, min)))

    df_meddonors <- df_meddonors |> group_by(donor, replicate, !!!syms(village$donorcols)) |>
      summarise(gr=abs(representation[time == max(time)]-representation[time == result$common_min])) |>
      ungroup()

  } else {

    result <- df_meddonors |>
      group_by(replicate) |>
      summarise(unique_times = list(unique(time))) |>
      ungroup() |>
      summarise(common_min = max(sapply(unique_times, min)))

    df_meddonors <- df_meddonors |> group_by(donor, replicate, !!!syms(village$donorcols)) |>
      summarise(gr=abs(representation[time == max(time)]-representation[time == result$common_min])) |>
      ungroup()
  }


  df_meddonors <- df_meddonors |>
    group_by(donor, !!!syms(village$donorcols)) |>
    summarise(sumdiff= sum(gr)) |>
    arrange(sumdiff)

  med <- median(df_meddonors$sumdiff)

  df_meddonors$median_gr <- abs(df_meddonors$sumdiff-med)
  df_meddonors <- df_meddonors[order(df_meddonors$median_gr),]
  village$alt_baseline <- df_meddonors$donor[1:10]

  # Set default baseline donor and restructure data
  if(is.null(village$baseline)){
    village$baseline <- df_meddonors$donor[1]
  }
  data_reset <- village$data[village$data$donor != village$baseline, ]
  data_reset <- data_reset |> group_by(donor) |> mutate(donorid=cur_group_id()) |> ungroup()
  village$data <- village$data[village$data$donor == village$baseline,]
  village$data$donorid <- village$num_donors

  village$data <- rbind(data_reset, village$data)
  warning(paste('Set baseline donor to: ', village$baseline, '\n View alternative baseline donors run: village$alt_baseline'))

  return(village)
}

