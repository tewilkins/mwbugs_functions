load_taxon_all <- function(){
  if (!exists("taxon_all_db")){
    db_mwbugs <- RPostgres::dbConnect(RPostgres::Postgres(), "mwbugs", host = "localhost", port = 5432, user = "readonly", password = "reachcode_42")
    taxon_all_db <- DBI::dbGetQuery(db_mwbugs, "SELECT * FROM taxon_all;")
    taxon_all_db$taxon[taxon_all_db$shortcode == "M"] <- "Acarina" # Change "M" from arachnida (class of arachnids) to hydrachnidia (unranked water mites)
    assign("taxon_all_db", taxon_all_db, envir = .GlobalEnv)
  } else {
    warning("taxon_all_db is already loaded in the global environment.")
  }
}

taxonMatch <- function(x){
  # If the passed vector is empty:
  y <- x
  if (length(y) == 0){
    warning(paste0(x, " is an empty vector."))
    return(x)
  }
  # Loads taxon_all from Melbourne Water database if not loaded already
  if (!exists("taxon_all_db")){
    db_mwbugs <- RPostgres::dbConnect(RPostgres::Postgres(), "mwbugs", host = "localhost", port = 5432, user = "readonly", password = "reachcode_42")
    taxon_all_db <- DBI::dbGetQuery(db_mwbugs, "SELECT * FROM taxon_all;")
    taxon_all_db$taxon[taxon_all_db$shortcode == "M"] <- "Acarina" # Change "M" from arachnida (class of arachnids) to hydrachnidia (unranked water mites)
    assign("taxon_all_db", taxon_all_db, envir = .GlobalEnv)
  }
  
  # Taxonomy dataframe no longer needed with full taxon_all_db 
  #  if (!exists("taxonomy")){
  #    taxonomy <- readxl::read_excel("~/uomShare/wergStaff/ChrisW/git-data/barcoding_data_mgt/spring_2018_bug_data_2023.xlsx", sheet = "taxonomy") 
  #  }
  
  taxon_all_db$taxon <- gsub("\\s*\\(Unitent\\.\\)", "", taxon_all_db$taxon)
  
  for (i in 1:length(x)){
    # taxon_all_db dataframe does not use "00" or "99" as trailing taxon codes. These should be removed from supplied taxoncodes.
    if (substr(x[i], nchar(x[i])-1, nchar(x[i])) %in% c("00","99")) {
      x[i] <- substr(x[i], 1, nchar(x[i])-2)
      warning("Taxonomic name(s) truncated to ", toString(x[i]), ".")
    }
    
    # Matches the taxoncode to taxon in taxon_all:
    y[i] <- taxon_all_db$taxon[match(x[i], taxon_all_db$shortcode)]
    
    #     # If a match isn't found in taxon_all:
    #     if (is.na(y[i]) == TRUE){
    #       
    #       # Check for order-level taxonomic levels that aren't explicit in taxonomy dataframe:
    #       if (nchar(x[i]) == 2) {
    #         y[i] <- taxonomy$order[match(x[i], substr(taxonomy$taxoncode, 1, 2))]
    #       } else {
    # #      y[i] <- taxonomy$species[match(x[i], taxonomy$taxoncode)]
    #         # Match to highest level known taxonomy, given taxoncode:
    #         y[i] <- taxonomy[max(which(!is.na(taxonomy[taxonomy$taxoncode == x[i],]))), match(x[i], taxonomy$taxoncode)]
    #       }
    #    }
  }
  # Returns a warning in case of NAs
  if (any(is.na(y))) warning("No taxonomic name(s) found for ", toString(x[which(is.na(y))]), ".")
  return(y)
  return(taxon_all_db)
}

# This function should only be applied to the species-level DNA metabarcoding data. The Morph-ID datasets follow a slightly different set of rules.
speciesToFamily <- function(x, rm = TRUE){ # rm = TRUE removes families that have ambiguously identified families, but the species is distinct from others in the dataset.
  if (length(x) == 0){
    return(x)
  } else {
    for (i in 1:length(x)){
      x[i] <- substr(x[i], 1, 4)
      # Combine all Acarina (there are only MP and MC taxa, no spiders (MA))
      if (substr(x[i], 1, 1) == "M") x[i] <- "M"
      
      # Combine all Oligochaeta
      if (substr(x[i], 1, 2) == "LO") x[i] <- "LO"
      
      # Combine all Nemertea
      if (substr(x[i], 1, 2) == "IH") x[i] <- "IH"
      
      # Combine Limoniidae (QD02) with Tipulidae
      if (substr(x[i], 1, 4) == "QD02") x[i] <- "QD01"
      
      # Combine Chironomini, Tanytarsini and Pseudochironomoni into Chironiminae
      if (substr(x[i], 1, 4) %in% c("QDAG","QDAH","QDAI")) x[i] <- "QDAJ"
      
      # Remove chironomids, odonates and isopods of uncertain family, given the morph and dna 
      # data both include specimens of these groups identified to family or lower
      if (rm == TRUE){
        if (substr(x[i], 1, 4) %in% c("QDAZ","OR99","QO99")) x[i] <- 0
      } else {
        if (substr(x[i], 1, 4) == "QDAZ") x[i] <- "QD"
        if (substr(x[i], 1, 4) == "OR99") x[i] <- "OR"
        if (substr(x[i], 1, 4) == "QO99") x[i] <- "QO"
      }
    }
  }
  return (x[!x == 0])
}

# Function for formatting tables when knitting to word from RMarkdown
table_for_word <- function(input_table, font = "Helvetica", pgwidth = 6.69){
  ft <- flextable::regulartable(input_table)
  ft <- flextable::separate_header(ft, split = "_")
  ft <- flextable::align(ft, align = "center", part = "all")
  ft <- flextable::valign(ft, valign = "top", part = "all")
  ft <- flextable::font(ft,fontname = font, part = "all")
  ft <- flextable::fontsize(ft, size = 8, part = "all")
  ft <- flextable::padding(ft, padding.top = 2,  padding.bottom = 2, part = "all")
  ft <- flextable::autofit(ft)
  # fit to window for MS Word (adapted from
  # https://stackoverflow.com/questions/57175351/flextable-autofit-in-a-rmarkdown-to-word-doc-causes-table-to-go-outside-page-mar)
  ft <- flextable::width(ft, width = dim(ft)$widths*pgwidth /(flextable_dim(ft)$widths))
  ft 
}

counterfactPlot <- function(predset_list, ylim = c(min(c(predset_list[[1]]$low,predset_list[[2]]$low))*0.9, max(c(predset_list[[1]]$high, predset_list[[2]]$high))*1.1), var = "", cf = "", ylab = "", xlab = "", title = "", logx = FALSE, rm_legend = FALSE, col1 = rgb(0.2, 0.2, 0.8, 0.25), col2 = rgb(0.4, 0.4, 0.4, 0.25)){
  curve1 <- predset_list[[1]]
  xax_labs <- curve1[[var]]
  if(logx == TRUE){
    curve1[[var]] <- log(curve1[[var]] + 0.05)
  }
  plot(curve1[[var]], curve1$mean, ylim = ylim, xlim = c(min(curve1[[var]]), max(curve1[[var]])), type = "n", ylab = ylab, xlab = xlab, xaxt = "n")
  axis(1, at = curve1[[var]], labels = xax_labs)
  lines(curve1[[var]], curve1$mean)
  polygon(c(curve1[[var]], rev(curve1[[var]])), 
          c(curve1$low, rev(curve1$high)), 
          col=col1, border=NA)
  curve2 <- predset_list[[2]]
  if(logx == TRUE){
    curve2[[var]] <- log(curve2[[var]] + 0.05)
  }
  lines(curve2[[var]], curve2$mean)
  polygon(c(curve2[[var]], rev(curve2[[var]])), 
          c(curve2$low, rev(curve2$high)), 
          col=col2, border=NA)
  title(main = title, adj = 0, font.main = 1, cex.main = 0.75)
  if(rm_legend == FALSE){
    legend("topleft", legend = c(paste0(cf, " = ", unique(curve1[[cf]])), paste0(cf, " = ", unique(curve2[[cf]]))), bg="transparent", 
           fill=c(col1,col2))
  }
}

list2layout <- function(cf_list, ncol = 3){
  x <- vector()
  icounter = 1
  for (i in 1:(length(cf_list)*ncol)){
    if (i %% ncol == 1){
      x[i] <- (length(cf_list)*(ncol-1))+1
    } else {
      x[i] <- icounter
      icounter = icounter + 1
    }
  }
  x <- append(x, c(0, seq(from = icounter + 1, to = icounter + ncol - 1, by = 1)))
  y <- matrix(data = x, nrow = length(cf_list) + 1, ncol = ncol, byrow = TRUE)
  return(y)
}

ct <- function (rows, 
                cols, 
                values = NULL, 
                FUN = sum, 
                convertNAToZero = TRUE,...) 
{
  if(!is.vector(rows)) rows <- as.vector(rows)
  if(!is.vector(cols)) cols <- as.vector(cols)
  if(is.null(values)) values <- rep(1,length(rows))
  results <- tapply(values, list(rows, cols), FUN, ...)
  if(convertNAToZero)
    results[is.na(results)] <- 0
  results
}

# Helper function generated by ChatGPT:
save_current_chunk <- function(label, output_dir = "figures/", output_file = NULL) {
  if (!requireNamespace("rstudioapi", quietly = TRUE)) {
    stop("rstudioapi package is required for this function.")
  }
  
  ctx <- rstudioapi::getSourceEditorContext()
  doc_lines <- ctx$contents
  
  start <- grep(paste0("^#\\| label: ", label), doc_lines)
  if (length(start) == 0) stop("Label not found")
  
  # Remove control options from code
  knitr_end <- 0
  for (line in doc_lines){
    if (substr(line, 1, 2) == "#|"){
      knitr_end <- knitr_end + 1
    }
  }
  
  chunk <- doc_lines[(start + knitr_end):length(doc_lines)]
  next_label <- grep("^#\\| label:", chunk[-1])
  if (length(next_label)) {
    chunk <- chunk[1:(next_label[1])]
  }
  
  chunk <- gsub("save_current_chunk\\(.*", "", chunk)
  
  if (is.null(output_file)) {
    output_file <- paste0(label, ".R")
  }
  
  writeLines(chunk, paste0(output_dir, output_file))
}

# requiredPackages <- c("readxl","tidyverse", "flextable", "rio", "rstan", "spdep", "geosphere")
load_packages <- function(requiredPackages, lib.loc = .libPaths()){
  for (p in 1:length(requiredPackages)){
    if(require(requiredPackages[p], lib.loc = lib.loc, character.only = TRUE)){
      print(paste0(requiredPackages[p], " is loaded correctly"))
    } else {
      print(paste0("installing ", requiredPackages[p], "..."))
      install.packages(paste0(requiredPackages[p]), lib.loc = lib.loc)
      if(require(requiredPackages[p], lib.loc = lib.loc, character.only = TRUE)){
        print(paste0(requiredPackages[p], " is installed and loaded"))
      } else {
        stop(paste0("could not install ", requiredPackages[p]))
      }
    }
  }
}

# Calc EBFMI from stanfit object:
check_energy <- function(stanfit){
  sampler_params <- get_sampler_params(stanfit, inc_warmup=FALSE)
  EBFMI <- rep(0, times = length(sampler_params))
  for (n in 1:length(sampler_params)) {
    energies <- sampler_params[n][[1]][,'energy__']
    numer <- sum(diff(energies)**2) / length(energies)
    denom <- var(energies)
    EBFMI[n] <- numer / denom
  }
  return(EBFMI)
}

# Get number of divergent transitions from a stanfit obejct
num_divergent <- function(stanfit){
  sampler_params <- get_sampler_params(stanfit, inc_warmup=FALSE)
  num_divergent_per_chain <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
  num_divergent_all_chains <- sum(num_divergent_per_chain)
  return(num_divergent_all_chains)
}

# makes a layout of specified columns and rows. split_x allows the last row to be either column specific or column shared.
makeLayout <- function(cols, rows, split_x = FALSE, split_y = FALSE, nudge = 0, nudge_x = 0, nudge_y = 0){
  lo_vector <- integer()
  y_rows <- 0
  for (i in 1:rows-1){
    if (split_y == FALSE){
      lo_vector <- c(lo_vector, cols*rows+1, (1:cols)+(i*cols))
    } 
    if (split_y == TRUE){
      lo_vector <- c(lo_vector, cols*rows+1+i, (1:cols)+(i*cols))
      y_rows <- i
    }
  }
  if (split_x == FALSE){
    lo_vector <- c(lo_vector, 0, rep(cols*rows+2+y_rows, times = cols))
  }
  if (split_x == TRUE){
    lo_vector <- c(lo_vector, 0, seq(from = cols*rows+2+y_rows, by = 1, length.out = cols))
  }
  lo <- layout(matrix(lo_vector, ncol = cols+1, byrow = T), widths = c(1+(cols/1.5) - (nudge + nudge_y), rep(5, times = cols)), heights = c(rep(5, times = rows), 2 - (nudge + nudge_x)))
  return(lo)
}
