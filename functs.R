	# options(future.globals.maxSize= 2097152000)
	# ------------------------------------------------
	
	# indicates the mz deviation in ppm
	calc_ppm_mz <- function (mz_tol, mz) {
		mz_ppm <- mz*(mz_tol/(10^6))
		return(max(mz_ppm, 0.00065))
	}

	#indicates the ret_t deviation
	calc_ret_t_tol <- function (ret_t, ret_tol = NULL) {
		
		if (!is.null(ret_tol)) {
			oper = str_extract(ret_tol, "[\\+\\-\\*\\/\\^]")
			num = as.numeric(str_extract(ret_tol, "\\d\\.?\\d*"))
			if (oper == "+") {
				return((0.10*ret_t^0.3)+num)
			} else if (oper == "-") {
				return((0.10*ret_t^0.3)-num)
			} else if (oper == "*") {
				return((0.10*ret_t^0.3)*num)
			} else if (oper == "/") {
				return((0.10*ret_t^0.3)/num)
			} else if (oper == "^") {
				return((0.10*ret_t^0.3)^num)
			}
		} else {
			return(0.10*ret_t^0.3)
		}
	}
	
	# ------------------------------------------------
	# parsing the mgf files
	mgf_split_func <- function(split_mgf, min_frags = 3) {
		split_mgf <- lapply(split_mgf, function(x) {
				
				x$precursor <- round(as.numeric(str_extract(x[2, 4], "\\d+\\.\\d{1,6}")), 5)

				if (any(pos <- str_detect(x$IONS, "CHARGE"))) {
					x$charge <- str_extract(x$IONS[pos], "\\d.+")
					if (as.numeric(str_extract(x$charge[1], "\\d")) > 1) {
						ioniz_mat <- str_extract_all(x$charge[1], c("\\d+", "[\\+\\-]"), simplify = TRUE)
						mult <- as.numeric(ioniz_mat[1, 1])
						ioniz <- ioniz_mat[2, 1]
						if (ioniz == "+") {
							x$molec_ion <- round((mult*x$precursor - (mult - 1)*1.007275), 5)
						} else {
							x$molec_ion <- round((mult*x$precursor + (mult - 1)*1.007275), 5)
						}
					} else {
						x$molec_ion <- x$precursor
					}
				} else {
					x$molec_ion <- x$precursor
				}
				x$ret_t <- round(as.numeric(str_extract(x[1, 4], "\\d+.\\d{1,}"))/60, 2)
				x$FS_int <- round(as.numeric(x$INTENSITY[which.max((x$INTENSITY))]))
				
				x <- x[str_detect(x$IONS, "^\\d"), ]
				x$IONS <- round(as.numeric(x$IONS), 5)
				x$INTENSITY <- round(as.numeric(x$INTENSITY), 1)
				x$INT_PERCENTAGE <- round((as.numeric(x$INTENSITY)/max(as.numeric(x$INTENSITY)))*100, 2)

				x <- x[order(x$IONS, decreasing = TRUE), ]

			return(x)	
			
		})
		split_mgf <- split_mgf[sapply(split_mgf, function(x) {(nrow(x) > min_frags)})]
		return(split_mgf)
	}
	# -------------------------------------
	# make a list of data.frames, every data.frame is data for a scan
	# min_int => minimum intensity, ions below this intensity are deleted, mz_tol the allowed ppm error 
	mgf_process_func <- function(fname_mgf, mz_tol = 15, min_int = 5000, pol = polarity) {
		
		start_time = Sys.time()
		
		dat_mgf <- read.delim(fname_mgf, sep = " ", strip.white = TRUE, header = FALSE)
		dat_mgf <- dat_mgf[, c(1, 2, ncol(dat_mgf))]

		pos_title_num <- str_which(dat_mgf[, 3], "^(NativeID)")
		
		scan_numb <- as.numeric(str_extract(dat_mgf[pos_title_num, 3], "(?<=scan=)\\d+"))
		
		dat_mgf <- cbind(charge = NA, precursor = NA, molec_ion = NA, dat_mgf, ret_t = NA, FS_int = NA, scan_numb = NA)

		colnames(dat_mgf) <- c("charge", "precursor", "molec_ion", "IONS", "INTENSITY", "INT_PERCENTAGE", "ret_t", "FS_int", "scan_numb")
		
		for (i in 1:(length(scan_numb)-1)) {
			dat_mgf$scan_numb[(pos_title_num[i] + 1):(pos_title_num[i + 1] - 3)] <- scan_numb[i]
		}
		
		dat_mgf$scan_numb[(pos_title_num[length(scan_numb)] + 1):(nrow(dat_mgf) - 1)] <- scan_numb[length(scan_numb)]

		split_mgf <- split(dat_mgf, dat_mgf$scan_numb)

			na.charge_nums <- sapply(split_mgf, function(x) {
				is.na(x$charge[1])
			})
			
			split_mgf[na.charge_nums] <- lapply(split_mgf[na.charge_nums], function(x) {
				x$charge = str_c(ifelse(pol == "pos", "+", "-"), "1")
				return(x)
			}) 
		# ------------------------------------------------------------
		split_mgf <- split_mgf[sapply(split_mgf, function(x) {nrow(x) > 6})]
	
		split_mgf <- mgf_split_func(split_mgf, min_frags = 3)
		
			split_mgf <- lapply(split_mgf, function(x) {
				x <- x[x$INTENSITY > min_int, ]
				x <- na.omit(x)
				return(x)
			})
			split_mgf <- split_mgf[sapply(split_mgf, function(x) {nrow(x) > 3})]			
		
		end_time = Sys.time()
		print(end_time - start_time)
		
		return(split_mgf)
	}
	
	# --------------------
	# filter the scans based on the presence of specific fragments and mass differences; 
	# func = all is to have all fragments present; func = any is to have at least 1 fragment ion present
	# it is mandatory to have at least 1 fragment ion for filtering
	frag_diff_filt_func <- function(dat_split, frags, diffs = NA, mz_tol = 15, coeff = 1.2, func = all) {

		if (length(dat_split) == 0) {
			return(FALSE)
		}
		
		logic_list_frags <- lapply(dat_split, function(x) {
				sapply(frags, function(y) {
						near(y, x$IONS, calc_ppm_mz(mz_tol, y)*coeff)
				})
			})
		
		logic_list_frags <- logic_list_frags[pos <- sapply(logic_list_frags, function(xx) {
								func(apply(xx, 2, sum))
							})]
							
						pos <- sapply(logic_list_frags, sum) >= length(frags)	
						pos <- which(names(dat_split) %in% names(pos))
						
		if (!is.na(diffs)) {
			logic_list_diffs <- lapply(dat_split[pos], function(x)
				{
					sapply(x$IONS, function(xx) {
						near(xx - x$IONS, diffs, calc_ppm_mz(mz_tol, xx))
					})
				})	
				
			pos <- sapply(logic_list_diffs, sum) >= length(diffs)
			pos <- which(names(dat_split) %in% names(pos))
		}
		
		return(pos)
	}
		
	# ---------------------------
	
	# combine the scans with similar features
	combine_scans_func <- function(total_dfs) {
	
		already_grouped <- to_del <- vector()
		i = 1
		while (i < length(total_dfs$IONS)) {

			if (!(i %in% already_grouped)) {

				pos1 <- determine_central_func(total_dfs$IONS, total_dfs$ret_t, total_dfs$INTENSITY, i, mz_tol)
				pos1 <- pos1[which(!(pos1 %in% already_grouped))]
				
				if (length(pos1) > 1) {
					total_dfs$IONS[pos1] <- round(weighted.mean(total_dfs$IONS[pos1], total_dfs$INTENSITY[pos1]), 4)
					total_dfs$INTENSITY[pos1] <- round(weighted.mean(total_dfs$INTENSITY[pos1], total_dfs$INTENSITY[pos1]), 2)
					to_del <- append(to_del, pos1[-1])
				}
				already_grouped <- append(already_grouped, pos1)
				
			}

			i = i + 1
		}

		total_dfs$ret_t <- round(total_dfs$ret_t[total_dfs$FS_int == max(total_dfs$FS_int)][1], 2)

		total_dfs$FS_int <- round(max(unique(total_dfs$FS_int)), 1) 

		scans_tot <- str_c(unique(total_dfs$scan_numb)[order(unique(total_dfs$scan_numb))], collapse = ", ")
		
		if (length(to_del) > 0) {
			total_dfs <- total_dfs[-to_del, ]
		}
		
		for (i in seq_along(total_dfs$INT_PERCENTAGE)) {
			max_int <- max(total_dfs$INTENSITY)
			total_dfs$INT_PERCENTAGE[i] <- round((total_dfs$INTENSITY[i]/max_int)*100, 2)
		}
		total_dfs$scan_numb <- str_c(unique(total_dfs$scan_numb)[order(unique(total_dfs$scan_numb))], collapse = ", ")
		total_dfs <- total_dfs[order(total_dfs$IONS, decreasing = TRUE), ]

		total_dfs$scan_numb <- scans_tot 
		
		return(total_dfs)
	}
	# ------------------------------------
	
	# determine which scan has the highest INTENSITY
		determine_central_func <- function(x, y, z, i, mz_tol) {
			
			pos1 <- which(near(x, x[i], calc_ppm_mz(mz_tol, x[i])) & near(y, y[i], calc_ret_t_tol(y[i])))

			if (length(pos1) > 1) {
				pos1.1 <- pos1[which.max(z[pos1])]

				return(which(near(x, mean(x[pos1.1]), calc_ppm_mz(mz_tol, x[pos1.1])[1]) & near(y, y[pos1.1], calc_ret_t_tol(y[pos1.1]))))
				
			} else {
				return(pos1)
			}
		}
		
		# ---------------------------------
		
	# create a publication ready table with MS data
	# top = n; n equals the number of ions to be included in the MS table
	# min_int_perc is the minimum intensity percentage for ion to be included
	tab_func <- function (combined_df_list, top = 30, min_int_perc = 0.5) {
	
		if (length(combined_df_list) == 0) {
			return(NA)
		}
	
		combined_df_list <- lapply(combined_df_list, function(x) 
			{
				x <- x[order(x$INT_PERCENTAGE, decreasing = TRUE), ]
				
			})
			
		combined_df_list <- combined_df_list[sapply(combined_df_list, function(x) {return(nrow(x)>0)}, simplify = TRUE)]
		
		
		if (length(combined_df_list) > 0) {
			
			pos_precursor <- which(str_detect(colnames(combined_df_list[[1]]), "precursor"))
			pos_IONS <- which(str_detect(colnames(combined_df_list[[1]]), "IONS"))
			pos_molec_ion <- which(str_detect(colnames(combined_df_list[[1]]), "molec_ion"))
			
			numeric_pos <- c(pos_precursor, pos_molec_ion, pos_IONS)
			
			combined_dfs_test <- lapply(combined_df_list, function(x) 
				{
					x[, numeric_pos] <- apply(x[, numeric_pos], 2, function(y) {format(round(y, 4), nsmall = 4)})
					return(x)
				})
				
			MS2 <- sapply(combined_dfs_test, function(x) {
					x <- str_c(apply(x, 1, function(y) {
						str_c(y[4], " (", y[6], ")")
					}), collapse = ", ")
					x <- str_replace_all(x, "\\(\\s{1,}", "\\(")
					return(x)
				}) 

			dfs_tab <- data.frame(code_num = names(combined_df_list),
								compound = NA,  
								chem_formula = NA, 
								charge = sapply(combined_df_list, function(x) {x$charge[1]}), 
								precursor = sapply(combined_df_list, function(x) {x$precursor[1]}), 
								molec_ion = sapply(combined_df_list, function(x) {x$molec_ion[1]}), 
								ppm = NA, 
								MS2 = MS2,
								rt_min = sapply(combined_df_list, function(x) {x$ret_t[1]}), 
								FS_int = sapply(combined_df_list, function(x) {x$FS_int[1]}), 
								scan_numb = sapply(combined_df_list, function(x) {x$scan_numb[1]}),
								row.names = NULL
							)	

			return(dfs_tab)
			
		} else {
			return(NULL)
		}
	}
		
	# ---------------------
			
	# fix duplicated STRINGS
		fix_dupl_str <- function (strng, str1 = "_") {
			for (i in seq_along(strng)) {
				if (!is.null(strng[i])) {

					if (sum(strng %in% strng[i]) > 1) {
						pos <- which(strng %in% strng[i])
						for (j in seq_along(2:length(pos)) + 1) {
							strng[pos][j] <- str_c(strng[pos][j], str1, j-1)
						}
					}
				}	
			}
			return(strng)
		}