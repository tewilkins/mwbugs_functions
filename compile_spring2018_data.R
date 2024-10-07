source("https://tools.thewerg.unimelb.edu.au/documents/misc/bugDatabaseFunctions.R")

# Import sites, samples, and biota tables as a list
data_list <- rio::import_list("~/uomShare/wergStaff/ChrisW/git-data/species_response_46_sites/data/biota_for_models.xlsx")
# Convert each table in data_list into a data.frame.
list2env(data_list,globalenv())

# sort samples list by sitecode (for ease of indexing/checking)
samples <- samples[order(samples$sitecode),]
# assign site_nos to sitecodes, and match to sample table
sites$site_no <- 1:nrow(sites)
samples$site_no <- sites$site_no[match(samples$sitecode,sites$sitecode)]

# scale the predictor variables
# keeping them as scaled objects for ease of back-transforming later
I <- scale(sites$lai)
F <- scale(sites$AF)
Q <- scale(sites$lmeanq)
CA <- scale(sites$lcarea)

#### samples table for morph samples 
samples_morph <- samples[samples$smpcode %in% unique(biota_morph_ab$smpcode),]

#### predictor matrix u for morph samples
I_morph <- as.vector(I)[match(samples_morph$sitecode, sites$sitecode)]
F_morph <- as.vector(F)[match(samples_morph$sitecode, sites$sitecode)]
Q_morph <- as.vector(Q)[match(samples_morph$sitecode, sites$sitecode)]
CA_morph <- as.vector(CA)[match(samples_morph$sitecode, sites$sitecode)]
hab <- as.numeric(samples_morph$collection_method == "RBA edge (sweep)") 
u_morph <- cbind(rep(1,nrow(samples_morph)),  # General intercept
                 I_morph,                 # Effective imperviousness
                 F_morph,                 # Attenuated Forest cover
                 I_morph * F_morph,       # Interaction betwen I and F
                 Q_morph,                 # mean annual runoff depth
                 Q_morph^2,               # square of Q for modal response
                 CA_morph,                # catchment area
                 hab)                     # habitat (1 = edge, 0 = riffle)

#### count matrix for morph_ab model
# check all counts are integers sum(biota_morph_ab$count %% 1 > 0) # 0, all good
c_morph_ab <- with(biota_morph_ab, ct(smpcode, taxoncode, count))
# ensure c and u rows are ordered the same way
c_morph_ab <- c_morph_ab[match(samples_morph$smpcode,row.names(c_morph_ab)),]

#### presence_absence matrix for morph samples 
pa_morph <- c_morph_ab
pa_morph[pa_morph > 0] <- 1

#### subsample fraction matrix (applies to both c_morph_ab and pa_morph)
s_morph <- pa_morph
for(i in 1:nrow(s_morph)){
  s_morph[i,] <- samples_morph$subsample_perc[samples_morph$smpcode == 
                                                row.names(s_morph)[i]]/100
}
# adjust coarsepick subsample fractions to 1
coarsepicks <- biota_morph_ab[biota_morph_ab$coarsepick == 1,]
for(i in 1:nrow(coarsepicks)){
  s_morph[row.names(s_morph) == coarsepicks$smpcode[i],
          coarsepicks$taxoncode[i]] <- 1
}

#### samples table for DNA samples 
samples_dna <- samples[samples$smpcode %in% unique(biota_dna_sp$smpcode),]

#### predictor matrix u for dna samples
I_dna <- as.vector(I)[match(samples_dna$sitecode, sites$sitecode)]
F_dna <- as.vector(F)[match(samples_dna$sitecode, sites$sitecode)]
Q_dna <- as.vector(Q)[match(samples_dna$sitecode, sites$sitecode)]
CA_dna <- as.vector(CA)[match(samples_dna$sitecode, sites$sitecode)]
hab <- as.numeric(samples_dna$collection_method == "RBA edge (sweep)") 
u_dna <- cbind(rep(1,nrow(samples_dna)),  # General intercept
               I_dna,                 # Effective imperviousness
               F_dna,                 # Attenuated Forest cover
               I_dna * F_dna,       # Interaction betwen I and F
               Q_dna,                 # mean annual runoff depth
               Q_dna^2,               # square of Q for modal response
               CA_dna)                # catchment area

#### count matrix for dna_sp  model
pa_dna_sp <- with(biota_dna_sp, ct(smpcode, taxoncode, count))
# ensure c and u rows are ordered the same way
pa_dna_sp <- pa_dna_sp[match(samples_dna$smpcode,row.names(pa_dna_sp)),]

#### subsample fraction matrix for pa_dna_sp
s_dna_sp <- pa_dna_sp
for(i in 1:nrow(s_dna_sp)){
  s_dna_sp[i,] <- samples_dna$subsample_perc[samples_dna$smpcode == 
                                               row.names(s_dna_sp)[i]]/100
}
# adjust coarsepick subsample fractions to 1
coarsepicks <- biota_dna_sp[biota_dna_sp$coarsepick == 1,]
if (nrow(coarsepicks) > 0){
  for(i in 1:nrow(coarsepicks)){
    s_dna_sp[row.names(s_dna_sp) == coarsepicks$smpcode[i],
             coarsepicks$taxoncode[i]] <- 1
  }
}

#### count matrix for dna_fam  model
pa_dna_fam <- with(biota_dna_fam, ct(smpcode, taxoncode, count))
# ensure c and u rows are ordered the same way
pa_dna_fam <- pa_dna_fam[match(samples_dna$smpcode,row.names(pa_dna_fam)),]

#### subsample fraction matrix for pa_dna_fam
s_dna_fam <- pa_dna_fam
for(i in 1:nrow(s_dna_fam)){
  s_dna_fam[i,] <- samples_dna$subsample_perc[samples_dna$smpcode == 
                                                row.names(s_dna_fam)[i]]/100
}
# adjust coarsepick subsample fractions to 1
coarsepicks <- biota_dna_fam[biota_dna_fam$coarsepick == 1,]
if (nrow(coarsepicks) > 0){
  for(i in 1:nrow(coarsepicks)){
    s_dna_fam[row.names(s_dna_fam) == coarsepicks$smpcode[i],
              coarsepicks$taxoncode[i]] <- 1
  }
}

#### Data lists for the 4 Stan models

c_morph_ab_data <- list(n_obs = nrow(c_morph_ab), # no. samples
                        n_pred = ncol(u_morph),       # no. predictors
                        n_taxa = ncol(c_morph_ab),  # no. taxa
                        n_site = nrow(sites),       # no. sites
                        u = as.matrix(u_morph),       # predictor model matrix
                        c = as.matrix(c_morph_ab), # taxa pa matrix
                        s = as.matrix(s_morph),  # subsample proportion for each obs
                        site = samples_morph$site_no) # site number for each obs

pa_morph_data <- list(n_obs = nrow(pa_morph), # no. samples
                      n_pred = ncol(u_morph),       # no. predictors
                      n_taxa = ncol(pa_morph),  # no. taxa
                      n_site = nrow(sites),       # no. sites
                      u = as.matrix(u_morph),       # predictor model matrix
                      pa = as.matrix(pa_morph), # taxa pa matrix
                      s = as.matrix(s_morph),  # subsample proportion for each obs
                      site = samples_morph$site_no) # site number for each obs

pa_dna_sp_data <- list(n_obs = nrow(pa_dna_sp), # no. samples
                       n_pred = ncol(u_dna),       # no. predictors
                       n_taxa = ncol(pa_dna_sp),  # no. taxa
                       n_site = nrow(sites),       # no. sites
                       u = as.matrix(u_dna),       # predictor model matrix
                       pa = as.matrix(pa_dna_sp), # taxa pa matrix
                       s = as.matrix(s_dna_sp),  # subsample proportion for each obs
                       site = samples_dna$site_no) # site number for each obs

pa_dna_fam_data <- list(n_obs = nrow(pa_dna_fam), # no. samples
                        n_pred = ncol(u_dna),       # no. predictors
                        n_taxa = ncol(pa_dna_fam),  # no. taxa
                        n_site = nrow(sites),       # no. sites
                        u = as.matrix(u_dna),       # predictor model matrix
                        pa = as.matrix(pa_dna_fam), # taxa pa matrix
                        s = as.matrix(s_dna_fam),  # subsample proportion for each obs
                        site = samples_dna$site_no) # site number for each obs
