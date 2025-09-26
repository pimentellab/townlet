#' Check Package Installation
#'
#' Runs test model to verify that Townlet work correctly after installation.
#' @return Fitted model results from test data
#' @export
test_townlet <- function() {
  # write temp data file
  data("testdata", package = "townlet")
  TAF::mkdir('./test_townlet')
  tmp_csv <- file.path("./test_townlet/df_temp.csv")
  write.csv(testdata, tmp_csv, row.names = FALSE)

  # test townlet
  village <- init_village(datapath = tmp_csv,
                          model = '~sex',
                          outdir = './',
  )
  village <- run_townlet(village = village,
                         cores = 4)

  # Check that all files created
  townlet_files <- files <- list.files("./test_townlet", full.names = TRUE)

  files <- c("./test_townlet/test_townlet_ESS.png",
             "./test_townlet/test_townlet_eta_params.png",
             "./test_townlet/test_townlet_growthrates.png",
             "./test_townlet/test_townlet_phi_r.png",
             "./test_townlet/test_townlet_ppck_barplot.png",
             "./test_townlet/test_townlet_rawdata_scatter.png",
             "./test_townlet/test_townlet_rawdata_stackbar.png",
             "./test_townlet/test_townlet_regressors.png",
             "./test_townlet/test_townlet_T0donor_representation.png",
             "./test_townlet/test_townlet_theta.png",
             "./test_townlet/test_townlet.RDS"
  )

  if(all(townlet_files %in% files)) {
    print('Townlet installation successful!')
  }

  on.exit(file.remove(tmp_csv), add = TRUE)

  return(village)
}


