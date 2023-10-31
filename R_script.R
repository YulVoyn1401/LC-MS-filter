# This script filters a Thermo .raw LC-MS file for specified fragment ions and specified mass differences and outputs the result as a xlsx file.

cat("This script filters a Thermo .raw LC-MS file for specified fragment ions and specified mass differences and outputs the result as an .xlsx file")

cat("The mass spectra file should be transformed to an mgf file, refer to the 'msconvert' program:\n \nhttp://proteowizard.sourceforge.net/download.html\n")

#If you do not have these packages install them by entering:

	#install.packages(tidyverse) 	
	#install.packages(readr) 		
	#install.packages(dplyr) 		
	#install.packages(sos)			
	#install.packages(readxl)	
	#install.packages(data.table)
	#install.packages(openxlsx)	

cat("\nIf you do not have these packages, install them by entering: \n\ninstall.packages('tidyverse') \ninstall.packages('readr') \ninstall.packages('dplyr')  \ninstall.packages('sos') \ninstall.packages('readxl') \ninstall.packages('data.table') \ninstall.packages('openxlsx')\n")

library(tidyverse) 		#call the "tidyverse" library
library(readr) 			#call the "readr" library
library(dplyr) 			#call the "dplyr" library
library(sos)			#call the "sos" library
library(readxl)			#call the "readxl" library
library(data.table)		#call the "data.table" library
library(openxlsx)		#call the "openxlsx" library

invisible(readline(prompt="press [Enter] to continue, or [Esc] to exit"))
cat("In the working directory, you should have a file called 'settings_file', which contains infromation on the data filtering tasks. \n")
cat('\nEnter the path to the working directory:')

work_dir <- back2ForwardSlash()

invisible(readline(prompt="parsing the mgf files, press [Enter] to continue"))

data_filtering <- read_excel(file.path(work_dir, "settings_file.xlsx"), sheet = 1)
settings_info <- read_excel(file.path(work_dir, "settings_file.xlsx"), sheet = 2)

	min_int <- as.numeric(settings_info$value[1])
	mz_tol <- as.numeric(settings_info$value[2])
	polarity <- settings_info$value[3]

mgf_folder_path <- file.path(work_dir, "mgfs")
	fnames_mgf <- file.path(mgf_folder_path, list.files(mgf_folder_path))

# parsing the mgf files, it may take 1-10 min per mgf file
dat_split_list <- lapply(fnames_mgf, mgf_process_func)
names(dat_split_list) <- str_sub(list.files(mgf_folder_path), end = -5)

invisible(readline(prompt="mgf files parsed, press [Enter] to continue"))

cat("\n result files will be saved in your working directory, 'results' folder")

invisible(readline(prompt="\nPress [enter] to continue"))

# filtering the data

	LC_MS_tab_list <- LC_MS_tab_list_tot <- list()

	for (z in seq_along(fnames_mgf)) {

		dat_split <- dat_split_list[[z]]
		
		for (zz in seq_along(data_filtering$classes)) {
			
			combined_df_list <- list()
			
			# exclude_ions <- extract_MS2(data_filtering$exclude_MS2[zz])
			searched_frags <- as.numeric(str_extract_all(data_filtering$MS2[zz], boundary("word"), simplify = TRUE))			
			searched_diffs <- as.numeric(str_extract_all(data_filtering$differences[zz], boundary("word"), simplify = TRUE))	

			dat_split_filt <- dat_split[frag_diff_filt_func(dat_split, frags = searched_frags, diffs = searched_diffs, mz_tol = mz_tol, func = all)]
			
				dsc_cut <- do.call(rbind, lapply(dat_split_filt, function(x) {
					return(x[1, c(3, 7, 8, 9)])
				}))
				
				already_grouped <- vector()
				for (i in seq_along(dsc_cut$molec_ion)) {
					if (!(i %in% already_grouped)) {
						pos <- determine_central_func(dsc_cut$molec_ion, dsc_cut$ret_t, dsc_cut$FS_int, i, mz_tol)
						
						total_dfs <- do.call(rbind, dat_split_filt[pos])
						combined_df <- combine_scans_func(total_dfs)
						combined_df_list <- append(combined_df_list, list(combined_df))
						
						already_grouped <- append(already_grouped, pos)
					}
				}	
				
				names(combined_df_list) <- fix_dupl_str(sapply(combined_df_list, function(x) {
							return(format(x$molec_ion[1], nsmall = 4))
						}))
			LC_MS_tab <- tab_func(combined_df_list)
			LC_MS_tab_list[[zz]] <- LC_MS_tab
		}
		LC_MS_tab_list_tot[[z]] <- LC_MS_tab_list
	}
	
	LC_MS_tab_list_tot <- lapply(LC_MS_tab_list_tot, function(x) {
		names(x) <- data_filtering$classes
		return(x)
	})
	names(LC_MS_tab_list_tot) <- names(dat_split_list)
	# ----------------------

	# saving results
	if (!(dir.exists(file.path(work_dir, "results")))) {
		dir.create(file.path(work_dir, "results"), recursive = TRUE)
	}
	
	for (i in seq_along(LC_MS_tab_list_tot)) {
		write.xlsx(LC_MS_tab_list_tot[[i]], file.path(work_dir, "results", str_c(names(LC_MS_tab_list_tot)[i], ".xlsx")))
	}
	
	Readme <- file(file.path(work_dir, "results", "README.txt"), "w")
	writeLines("Every xlsx result file represents the results for a particular mgf file.\n", Readme, sep = "")
	writeLines("\nEvery xlsx file is divided in several sheets representing the class of searched compounds, defined in the 'settings.xlsx' file.\n", Readme, sep = "")
	writeLines("\n")
	writeLines("\n'FS_int' is the Full-scan intensity of the peak, ", Readme, sep = "")
	writeLines("'scan_numb' are the scans from which the substance is derived.", Readme, sep = "")
	close(Readme)

	cat("\n result files were saved in your working directory, 'results' folder")