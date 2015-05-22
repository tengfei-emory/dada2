#include "dada.h"
#include <Rcpp.h>

using namespace Rcpp;
//' @useDynLib dadac
//' @importFrom Rcpp evalCpp

B *run_dada(Raw **raws, int nraw, double score[4][4], Rcpp::NumericMatrix errMat, double gap_pen, bool use_kmers, double kdist_cutoff, int band_size, double omegaA, bool use_singletons, double omegaS, int max_clust, double min_fold, int min_hamming, bool use_quals);

//------------------------------------------------------------------
//' Run DADA on the provided unique sequences/abundance pairs. 
//' 
//' @return List.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List dada_uniques(Rcpp::CharacterVector seqs,  Rcpp::IntegerVector abundances,
                        Rcpp::NumericMatrix err,
                        Rcpp::NumericMatrix quals,
                        Rcpp::NumericMatrix score, Rcpp::NumericVector gap,
                        Rcpp::LogicalVector use_kmers, Rcpp::NumericVector kdist_cutoff,
                        Rcpp::NumericVector band_size,
                        Rcpp::NumericVector omegaA, 
                        Rcpp::LogicalVector use_singletons, Rcpp::NumericVector omegaS,
                        Rcpp::NumericVector maxClust,
                        Rcpp::NumericVector minFold, Rcpp::NumericVector minHamming,
                        Rcpp::LogicalVector useQuals,
                        int qMin, int qMax) {

  int i, j, pos, len1, len2, nrow, ncol;
  int f, r, s, nzeros, nones;
  size_t index;
  double tote;
  
  len1 = seqs.size();
  len2 = abundances.size();
  if(len1 != len2) {
    Rcpp::Rcout << "C: Different input lengths:" << len1 << ", " << len2 << "\n";
    return R_NilValue;
  }
  int nraw = len1;
    
  std::vector<std::string> cpp_seqs(nraw);
  std::vector<int> cpp_abunds(nraw);
  int max_seq_len = 0;
  for(i=0;i<nraw;i++) {
    cpp_seqs[i] = std::string(seqs[i]);
    if(cpp_seqs[i].length() > max_seq_len) { max_seq_len = cpp_seqs[i].length(); }
    cpp_abunds[i] = abundances[i];
  }
  if(max_seq_len >= SEQLEN) { Rcpp::stop("Input sequences exceed the maximum allowed string length of", SEQLEN-1); }
  // Need one extra byte for the string termination character
//  std::vector<int> cpp_abunds = Rcpp::as<std::vector<int> >(abunds);
  
  // Each sequence is a COLUMN, each row is a POSITION
  bool has_quals = false;
  if(quals.nrow() > 0) { has_quals = true; }
  int qlen = quals.nrow();

  // Make the raws (uniques)
  char seq[SEQLEN];
  double qual[SEQLEN];
  Raw **raws = (Raw **) malloc(nraw * sizeof(Raw *)); //E
  if (raws == NULL)  Rcpp::stop("Memory allocation failed!\n");

  for (index = 0; index < nraw; index++) {
    strcpy(seq, cpp_seqs[index].c_str());
    nt2int(seq, seq);
    for(pos=0;pos<strlen(seq);pos++) {
      qual[pos] = quals(pos, index);
    }
    if(has_quals) {
      raws[index] = raw_qual_new(seq, qual, cpp_abunds[index]);
    } else {
      raws[index] = raw_new(seq, cpp_abunds[index]);
    }
    raws[index]->index = index;
  }
  
  // Copy score into a C style array
  nrow = score.nrow();
  ncol = score.ncol();
  if(nrow != 4 || ncol != 4) {
    Rcpp::Rcout << "C: Score matrix malformed:" << nrow << ", " << ncol << "\n";
    return R_NilValue;
  }
  double c_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = score(i,j);
    }
  }

  nrow = err.nrow();
  ncol = err.ncol();
  if(nrow != 16) {
    Rcpp::Rcout << "C: Error matrix malformed:" << nrow << ", " << ncol << "\n";
    return R_NilValue;
  }

  // Check rest of params for Length-One-ness and make into C versions
  if(gap.size() != 1) { Rcpp::Rcout << "C: Gap penalty not length 1\n"; return R_NilValue; }
  double c_gap = as<double>(gap);
  
  // Copy use kmers into a C++ bool
  if(use_kmers.size() != 1) { Rcpp::Rcout << "C: Use_kmers not length 1\n"; return R_NilValue; }
  bool c_use_kmers = as<bool>(use_kmers);

  // Copy kdist_cutoff into a C double
  if(kdist_cutoff.size() != 1) { Rcpp::Rcout << "C: Kdist cutoff not length 1\n"; return R_NilValue; }
  double c_kdist_cutoff = as<double>(kdist_cutoff);

  if(band_size.size() != 1) { Rcpp::Rcout << "C: Band_size not length 1\n"; return R_NilValue; }
  int c_band_size = as<int>(band_size);

  if(omegaA.size() != 1) { Rcpp::Rcout << "C: OmegaA not length 1\n"; return R_NilValue; }
  double c_omegaA = as<double>(omegaA);

  if(use_singletons.size() != 1) { Rcpp::Rcout << "C: use_singletons not length 1\n"; return R_NilValue; }
  bool c_use_singletons = as<bool>(use_singletons);

  if(omegaS.size() != 1) { Rcpp::Rcout << "C: OmegaS not length 1\n"; return R_NilValue; }
  double c_omegaS = as<double>(omegaS);

  if(maxClust.size() != 1) { Rcpp::Rcout << "C: MaxClust not length 1\n"; return R_NilValue; }
  int c_max_clust = as<int>(maxClust);

  if(minFold.size() != 1) { Rcpp::Rcout << "C: MinFold not length 1\n"; return R_NilValue; }
  double c_min_fold = as<double>(minFold);

  if(minHamming.size() != 1) { Rcpp::Rcout << "C: MinHamming not length 1\n"; return R_NilValue; }
  int c_min_hamming = as<int>(minHamming);

  if(useQuals.size() != 1) { Rcpp::Rcout << "C: UseQuals not length 1\n"; return R_NilValue; }
  bool c_use_quals = as<bool>(useQuals);

  // Run DADA
  B *bb = run_dada(raws, nraw, c_score, err, c_gap, c_use_kmers, c_kdist_cutoff, c_band_size, c_omegaA, c_use_singletons, c_omegaS, c_max_clust, c_min_fold, c_min_hamming, c_use_quals);

  // Extract output from Bi objects
  char **oseqs = (char **) malloc(bb->nclust * sizeof(char *)); //E
  if (oseqs == NULL)  Rcpp::stop("Memory allocation failed!\n");
  
  for(i=0;i<bb->nclust;i++) {
    oseqs[i] = (char *) malloc((strlen(bb->bi[i]->seq)+1) * sizeof(char)); //E
    if (oseqs[i] == NULL)  Rcpp::stop("Memory allocation failed!\n");
    ntcpy(oseqs[i], bb->bi[i]->seq);
  }

  // Convert to R objects and return
  Rcpp::CharacterVector Rseqs;
  Rcpp::NumericVector Rabunds(bb->nclust);
  Rcpp::NumericVector Rzeros(bb->nclust);
  Rcpp::NumericVector Rones(bb->nclust); 
  Rcpp::NumericVector Rraws(bb->nclust); 
  Rcpp::NumericVector Rfams(bb->nclust); 
  Rcpp::NumericVector Rbirth_pvals(bb->nclust);
  Rcpp::NumericVector Rbirth_folds(bb->nclust);
  Rcpp::NumericVector Rbirth_hams(bb->nclust);
  Rcpp::NumericVector Rbirth_es(bb->nclust);
  Rcpp::CharacterVector Rbirth_types;
  Rcpp::NumericVector Rbirth_qaves(bb->nclust);
  Rcpp::NumericVector Rpvals(bb->nclust);


  for(i=0;i<bb->nclust;i++) {
    // The basics
    Rseqs.push_back(std::string(oseqs[i]));
    Rabunds[i] = bb->bi[i]->reads;
    
    // Extra info on the clusters
    nzeros = 0; nones = 0;
    for(f=0;f<bb->bi[i]->nfam;f++) {
      if(bb->bi[i]->fam[f]->sub) {
        if(bb->bi[i]->fam[f]->sub->nsubs == 0) { nzeros += bb->bi[i]->fam[f]->reads; }
        if(bb->bi[i]->fam[f]->sub->nsubs == 1) { nones += bb->bi[i]->fam[f]->reads; }
      }
    }
    Rzeros[i] = nzeros;
    Rones[i] = nones;
    Rraws[i] = bb->bi[i]->nraw;
    Rfams[i] = bb->bi[i]->nfam;
    
    // Record information from the cluster's birth
    Rbirth_types.push_back(std::string(bb->bi[i]->birth_type));
    Rbirth_pvals[i] = bb->bi[i]->birth_pval;
    Rbirth_folds[i] = bb->bi[i]->birth_fold;
    if(bb->bi[i]->birth_sub) { Rbirth_hams[i] = bb->bi[i]->birth_sub->nsubs; }
    else { Rbirth_hams[i] = NA_REAL; }
    Rbirth_es[i] = bb->bi[i]->birth_e;

    // Make qave column
    Sub *sub;
    double qave;
    if(has_quals) {
      qave = 0.0;
      sub = bb->bi[i]->birth_sub;
      if(sub && sub->q1) {
        for(int s=0;s<sub->nsubs;s++) {
          qave += (sub->q1[s]);
        }
        qave = qave/((double)sub->nsubs);
      }
      Rbirth_qaves[i] = qave;
    }
    
    // Calculate post-hoc pval
    tote = 0.0;
    for(j=0;j<bb->nclust;j++) {
      if(i != j) {
        tote += bb->bi[j]->e[bb->bi[i]->center->index];
      }
    }
    Rpvals[i] = calc_pA(1+bb->bi[i]->reads, tote); // better to have expected number of center sequence here?
  }
  Rcpp::DataFrame df_clustering = Rcpp::DataFrame::create(_["sequence"] = Rseqs, _["abundance"] = Rabunds, _["n0"] = Rzeros, _["n1"] = Rones, _["nunq"] = Rraws, _["nfam"] = Rfams, _["pval"] = Rpvals, _["birth_type"] = Rbirth_types, _["birth_pval"] = Rbirth_pvals, _["birth_fold"] = Rbirth_folds, _["birth_ham"] = Rbirth_hams, _["birth_qave"] = Rbirth_qaves);

  // Get error (or substitution) statistics
  int32_t otrans[4][4];
  b_get_trans_matrix(bb, otrans);

  Rcpp::IntegerMatrix Rtrans(4, 4);  // R INTS ARE SIGNED 32 BIT
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      Rtrans(i,j) = otrans[i][j];
    }
  }
  
  Rcpp::DataFrame df_subpos = b_get_positional_subs(bb);
  Rcpp::IntegerMatrix transMat = b_get_quality_subs2(bb, has_quals, qMin, qMax);
  
  // Get output average qualities
  Rcpp::NumericMatrix Rquals(qlen, bb->nclust);
  double qq;
  int nreads, pos0, pos1;
  Sub *sub;
  Raw *raw;
  if(has_quals) {
    for(i=0;i<bb->nclust;i++) {
      len1 = strlen(bb->bi[i]->seq);
      if(len1 > qlen) {
        if(tVERBOSE) { printf("Warning: C%i sequence too long for output q matrix (%i vs. %i).\n", i, len1, qlen); }
        len1 = qlen;
      }
      for(pos0=0;pos0<len1;pos0++) {
        qq=0.0;
        nreads=0;
        for(f=0;f<bb->bi[i]->nfam;f++) {
          for(r=0;r<bb->bi[i]->fam[f]->nraw;r++) {
            raw = bb->bi[i]->fam[f]->raw[r];
            sub = bb->bi[i]->sub[raw->index];
            if(sub) {
              pos1 = sub->map[pos0];
              if(pos1 == GAP_GLYPH) { // Gap
                continue;
              }
              nreads += raw->reads;
              qq += (raw->qual[pos1] * raw->reads);
            } else {
//              if(tVERBOSE) { printf("Warning: No sub for R%i in C%i\n", r, i); }
              if(VERBOSE) { printf("Warning: No sub for R%i in C%i\n", r, i); }
            }
          }
        }
        Rquals(pos0,i) = qq/nreads;
      } // for(pos0=0;pos0<len1;pos0++)
    }
  }

  // Make map from uniques to cluster
  Rcpp::IntegerVector Rmap(nraw);
  for(i=0;i<bb->nclust;i++) {
    for(f=0;f<bb->bi[i]->nfam;f++) {
      for(r=0;r<bb->bi[i]->fam[f]->nraw;r++) {
        Rmap(bb->bi[i]->fam[f]->raw[r]->index) = i+1; // +1 for R 1-indexing
      }
    }
  }

  // Make output data.frame of birth_subs
  char buf[2];
  buf[0] = buf[1] = '\0';
  Rcpp::IntegerMatrix bsubs(16,max_seq_len);
  Rcpp::IntegerVector bs_pos;
  Rcpp::CharacterVector bs_nt0;
  Rcpp::CharacterVector bs_nt1;
  Rcpp::NumericVector bs_qual;
  Rcpp::IntegerVector bs_clust;
  for(i=0;i<bb->nclust;i++) {
    if(bb->bi[i]->birth_sub) {
      for(s=0;s<bb->bi[i]->birth_sub->nsubs;s++) {
        bs_pos.push_back(bb->bi[i]->birth_sub->pos[s]);
        buf[0] = bb->bi[i]->birth_sub->nt0[s];
        int2nt(buf, buf);
        bs_nt0.push_back(std::string(buf));
        buf[0] = bb->bi[i]->birth_sub->nt1[s];
        int2nt(buf, buf);
        bs_nt1.push_back(std::string(buf));
        if(has_quals) {
          bs_qual.push_back(bb->bi[i]->birth_sub->q1[s]);
        } else {
          bs_qual.push_back(Rcpp::NumericVector::get_na());
        }
        bs_clust.push_back(i+1); // R 1 indexing
      }
    }
  }
  Rcpp::DataFrame df_birth_subs = Rcpp::DataFrame::create(_["pos"] = bs_pos, _["ref"] = bs_nt0, _["sub"] = bs_nt1, _["qual"] = bs_qual, _["clust"] = bs_clust);

  // Free memory
  for(i=0;i<bb->nclust;i++) {
    free(oseqs[i]);
  }
  free(oseqs);
  b_free(bb);
  for(index=0;index<nraw;index++) {
    raw_free(raws[index]);
  }
  free(raws);
  
  // Organize return List  
  return Rcpp::List::create(_["clustering"] = df_clustering, _["subpos"] = df_birth_subs, _["subqual"] = transMat, _["trans"] = Rtrans, _["clusterquals"] = Rquals, _["map"] = Rmap);
}

B *run_dada(Raw **raws, int nraw, double score[4][4], Rcpp::NumericMatrix errMat, double gap_pen, bool use_kmers, double kdist_cutoff, int band_size, double omegaA, bool use_singletons, double omegaS, int max_clust, double min_fold, int min_hamming, bool use_quals) {
  int newi=0, nshuffle = 0;
  bool shuffled = false;
  double inflation = 1.0;
  size_t index;
  
  B *bb;
  bb = b_new(raws, nraw, score, gap_pen, omegaA, use_singletons, omegaS, band_size, use_quals); // New cluster with all sequences in 1 bi and 1 fam
  b_lambda_update(bb, FALSE, 1.0, errMat); // Everyone gets aligned within the initial cluster, no KMER screen
  b_fam_update(bb);     // Organizes raws into fams, makes fam consensus/lambda
  b_p_update(bb);       // Calculates abundance p-value for each fam in its cluster (consensuses)
  
  if(max_clust < 1) { max_clust = bb->nraw; }
  
  while( (bb->nclust < max_clust) && (newi = b_bud(bb, min_fold, min_hamming)) ) {
    if(tVERBOSE) printf("----------- New Cluster C%i -----------\n", newi);
    ///! TESTING TESTING TESTING - Centers locked at time of bud
    ///! b_consensus_update(bb);
    b_lambda_update(bb, use_kmers, kdist_cutoff, errMat);
    // SHOULD GET_SELF ALSO USE QUALS??
    
    // Temporarily inflate the E's for the new cluster based on the expected number of reads from its center
    if((int) (bb->bi[newi]->center->reads/bb->bi[newi]->self) > bb->bi[newi]->reads) {
      inflation = (bb->bi[newi]->center->reads/bb->bi[newi]->self)/bb->bi[newi]->reads;
      for(index=0;index<bb->nraw;index++) {
        bb->bi[newi]->e[index] = bb->bi[newi]->e[index] * inflation;
        ///! TESTING TESTING TESTING - REMOVING UNNECESSARY LAMBDA UPDATES
        ///! bb->bi[newi]->update_lambda = TRUE;
      }
    }
    
    // Keep shuffling and updating until no more shuffles
    nshuffle = 0;
    do {
      shuffled = b_shuffle(bb);
      ///! TESTING TESTING TESTING - Centers locked at time of bud
      ///! b_consensus_update(bb);
      ///! TESTING TESTING TESTING - REMOVING UNNECESSARY LAMBDA UPDATES
      ///! b_lambda_update(bb, use_kmers, kdist_cutoff, errMat);
      b_e_update(bb); ///!
      if(tVERBOSE) { printf("S"); }
    } while(shuffled && ++nshuffle < MAX_SHUFFLE);
    
    if(tVERBOSE && nshuffle >= MAX_SHUFFLE) { printf("\nWarning: Reached maximum (%i) shuffles.\n", MAX_SHUFFLE); }
    
    b_fam_update(bb); // If centers can move, must have lambda_update before fam_update
    b_p_update(bb);
  }
  
  printf("\nALIGN: %i aligns, %i shrouded (%i raw).\n", bb->nalign, bb->nshroud, bb->nraw);
  
  return bb;
}

