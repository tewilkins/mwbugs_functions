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
    # taxon_all_db dataframe does not use "00" as trailing taxon codes. These should be removed from supplied taxoncodes.
    if (substr(x[i], 7, 8) == "00") {
      x[i] <- substr(x[i], 1, 6)
    }
    # same with "99" as trailing taxon codes. These should be removed from supplied taxoncodes.
    if (substr(x[i], nchar(x[i])-1, nchar(x[i])) == "99") {
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
