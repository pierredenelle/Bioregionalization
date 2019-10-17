
run_oslom <- function(dat, n_runs = 10, t_param = 0.1, cp_param = 0.5,
                      saving_directory){
  # Controls: dat must be a data.frame with three columns containing id1, id2
  # and similarity metric
  if(!(is.data.frame(dat))){
    stop("dat must be a data.frame with three columns containing id1, id2
  and similarity metric")
  }

  if(ncol(dat) != 3){
    stop("dat must be a data.frame with three columns containing id1, id2
  and similarity metric")
  }

  if(!(is.numeric(dat[, 3]))){
    stop("dat must be a data.frame with three columns containing id1, id2
  and similarity metric")
  }

  if(!(abs(n_runs - round(n_runs)) < .Machine$double.eps^0.5)){
    stop("n_runs must be an integer setting the number of runs.")
  }

  if(!(is.numeric(t_param))){
    stop("t_param must be numeric.")
  }

  if(t_param > 1 | t_param < 0){
    stop("t_param must be comprised between 0 and 1.")
  }

  if(!(is.numeric(cp_param))){
    stop("cp_param must be numeric.")
  }

  if(cp_param > 1 | cp_param < 0){
    stop("cp_param must be comprised between 0 and 1.")
  }

  if(!(is.character(saving_directory))){
    stop("saving_directory must be a path where the OSLOM .tp file containing
         the bioregions identified will be saved.")
  }

  # Move OSLOM to saving_directory
  # path1 <- paste0(.libPaths(), "/Bioregionalization/OSLOM/")
  # current_files <- list.files(path1, full.names = TRUE)
  # path2 <- saving_directory
  # file.copy(from = current_files, to = path2, recursive = TRUE,
  #           copy.mode = TRUE)

  # Save input dataset as a .txt file into OSLOM folder
  write.table(dat,
              # paste0(saving_directory, "dataset.txt"),
              paste0(.libPaths(), "/Bioregionalization/OSLOM/dataset.txt"),
              row.names = FALSE)

  # Change working directory so the file is saved in the proper place
  current_path <- getwd()
  setwd(paste0(.libPaths(), "/Bioregionalization/"))

  # Set up the command with required parameters
  if(.Platform$OS.type == "windows"){
    # cmd <- paste0(
    #   "D:/PIERRE_DENELLE/CarHab/Bioregionalization/OSLOM/oslom_undir_win.exe -f OSLOM/dataset.txt -w",
    #   " -r ", n_runs, " -t ", t_param, " -cp ", cp_param)
    cmd <-
      paste0(.libPaths(),
             "/Bioregionalization/OSLOM/oslom_undir_win.exe -f OSLOM/dataset.txt -w",
             " -r ", n_runs, " -t ", t_param, " -cp ", cp_param)

    # cmd <-
    #   paste0(saving_directory,
    #          "oslom_undir_win.exe -f OSLOM/dataset.txt -w",
    #          " -r ", n_runs, " -t ", t_param, " -cp ", cp_param)

  } else if(.Platform$OS.type == "unix"){
    stop("To do")
  } else{
    stop("Windows or Unix distributions only.")
  }

  # Execute the command from R
  system(command = cmd)

  # Import tp file created
  tp_res <- readLines(
    paste0(.libPaths(),
           "/Bioregionalization/OSLOM/dataset.txt_oslo_files/tp"))
  # return("OSLOM/dataset.txt_oslo_files/tp")

  # Saving .tp file into chosen saving_directory
  # version = 2: files readable for R versions from 1.4.0 to 3.5.0
  save(tp_res, file = paste0(saving_directory, "/tp"), version = 2)

  # Remove .oslo_files created and the dataset
  # file.remove("OSLOM/dataset.txt")
  # file.remove("OSLOM/time_seed.dat")
  file.remove(paste0(.libPaths(), "/Bioregionalization/OSLOM/dataset.txt"))
  file.remove(paste0(.libPaths(), "/Bioregionalization/OSLOM/time_seed.dat"))

  # Remove all filed in .oslo_files folder
  # unlink("OSLOM/dataset.txt_oslo_files/*", recursive = FALSE)
  # unlink(paste0(.libPaths(),
  #               "/Bioregionalization/OSLOM/dataset.txt_oslo_files/*"),
  #        recursive = FALSE)
  file.remove(paste0(.libPaths(),
                     "/Bioregionalization/OSLOM/dataset.txt_oslo_files/",
                     dir(paste0(
                       .libPaths(),
                       "/Bioregionalization/OSLOM/dataset.txt_oslo_files/"),
                       pattern = "net")))
  file.remove(paste0(.libPaths(),
                     "/Bioregionalization/OSLOM/dataset.txt_oslo_files/",
                     dir(paste0(
                       .libPaths(),
                       "/Bioregionalization/OSLOM/dataset.txt_oslo_files/"),
                       pattern = "partitions_level")))
  file.remove(paste0(.libPaths(),
                     "/Bioregionalization/OSLOM/dataset.txt_oslo_files/",
                     dir(paste0(
                       .libPaths(),
                       "/Bioregionalization/OSLOM/dataset.txt_oslo_files/"),
                       pattern = "statistics_level")))
  file.remove(paste0(.libPaths(),
                     "/Bioregionalization/OSLOM/dataset.txt_oslo_files/",
                     dir(paste0(
                       .libPaths(),
                       "/Bioregionalization/OSLOM/dataset.txt_oslo_files/"),
                       pattern = "short_tp")))
  file.remove(paste0(.libPaths(),
                     "/Bioregionalization/OSLOM/dataset.txt_oslo_files/",
                     dir(paste0(
                       .libPaths(),
                       "/Bioregionalization/OSLOM/dataset.txt_oslo_files/"),
                       pattern = "tp")))

  # Reset previous working directory
  setwd(current_path)
}
