package_state <- new.env()

.onLoad <- function(libname, pkgname) {
    package_state$activity_scores <- readr::read_rds(fs::path_package("extdata","cfd_activity_scores.rds",package='CFDscoRe'))
}

#' CFD Score
#'
#' Calculates the CFD Score for a given guide/target site alignment.
#'
#' @param rna Character representation of the guide portion of a guide/target site alignment. Should be represented as DNA, valid characters include ['A','C','G','T','-']. Must be 20 nucleotides long.
#' @param rna Character representation of the target portion of a guide/target site alignment. Should be represented as DNA, valid characters include ['A','C','G','T','-']. No length requirement.
#' @param pam Character representation of the inferred PAM. A 'GG' PAM will receive no penalty.
#' @export
cfd_score <- function(rna,dna,pam){
  position_table <- alignment_to_position_table(toupper(rna),toupper(dna))
  classified_positions <- classify_positions(position_table)
  score_positions(classified_positions,pam)
}


alignment_to_position_table <- function(rna,dna){
  tibble(
    rna = strsplit(rna,"")[[1]],
    dna = strsplit(dna,"")[[1]]
  ) %>%
    mutate(index=row_number()) %>%
    mutate(is_dna_bulge = rna == '-') %>%
    mutate(tot_delete = cumsum(is_dna_bulge)) %>%
    mutate(tot_offset = c(0,tot_delete) %>% head(n=n())) %>%
    mutate(index=index-tot_offset) %>%
    select(-is_dna_bulge,-tot_delete,-tot_offset)
}

classify_positions <- function(position_table) {
  position_table %>%
    mutate(matching = case_when(
      rna == dna ~ 'match',
      rna == '-' ~ 'delete',
      dna == '-' ~ 'insert',
      rna != dna ~ 'mismatch'
    ))
}

score_positions <- function(classified_positions,PAM){
  pam_score <- package_state$activity_scores$pam %>%
    filter(pam==PAM) %>%
    pull(activity)

  mismatch_positions <- classified_positions %>%
    filter(matching=='mismatch')
  insert_positions <- classified_positions %>%
    filter(matching=='insert')
  delete_positions <- classified_positions %>%
    filter(matching=='delete')

  mismatch_scores <- score_mismatches(mismatch_positions)$activity
  insert_scores <- score_inserts(insert_positions)$activity
  delete_scores <- score_deletes(delete_positions)$activity

  prod(c(1,pam_score,mismatch_scores,insert_scores,delete_scores))
}

dna_complement <- c('A'='T','C'='G','G'='C','T'='A')

score_mismatches <- function(mismatched_positions){
  mismatched_positions %>%
    mutate(dna=dna_complement[dna]) %>%
    left_join(package_state$activity_scores$mismatch,by=join_by(rna, dna, index))
}

score_inserts <- function(insert_positions){
  insert_positions %>%
    dplyr::rename(insertion=rna) %>%
    left_join(package_state$activity_scores$rna_bulge,by=join_by(insertion, index))
}

score_deletes <- function(delete_positions){
  delete_positions %>%
    mutate(dna=dna_complement[dna]) %>%
    dplyr::rename(deletion=dna) %>%
    left_join(package_state$activity_scores$dna_bulge,by=join_by(deletion, index))
}
