package_state <- new.env()

.onLoad <- function(libname, pkgname) {
    package_state$activity_scores <- readr::read_rds(fs::path_package("extdata","cfd_activity_scores.rds",package='CFDscoRe'))
}

#' CFD Score
#'
#' Calculates the CFD Score for a given guide/target site alignment.
#'
#' @param rna Character representation of the guide portion of a guide/target site alignment. Should be represented as DNA, valid characters include ['A','C','G','T','-']. Must be 20 nucleotides long.
#' @param dna Character representation of the target portion of a guide/target site alignment. Should be represented as DNA, valid characters include ['A','C','G','T','-']. No length requirement.
#' @param pam Character representation of the inferred PAM. A 'GG' PAM will receive no penalty.
#' @return the CFD score, a numeric quantity between 0 and 1 that represents the fraction of guide activity
#' @export
cfd_score <- function(rna,dna,pam){
  if( nchar(rna) != nchar(dna) ){
      stop('Invalid alignment, guide and target sequences are not equal length')
  }
  position_table <- alignment_to_position_table(toupper(rna),toupper(dna))
  validate_input( position_table, pam )
  classified_positions <- classify_positions(position_table)
  score_positions(classified_positions,pam)
}

validate_input <- function( position_table, pam ){
    guide_length <- position_table %>%
        filter(rna != '-') %>%
        nrow()

    allowed_characters <- c( 'A','T','G','C','-' )

    if( guide_length != 20){ stop('CFD can only be calculated for guides that are 20nt long' ) }

    if( !grepl('^[ATGC]{2,2}$', pam) ){ stop('Invalid PAM sequence, must be length 2 and contain only (A,C,G,T)') }

    if( length( setdiff(position_table$rna, allowed_characters) ) > 0 ){
        stop('Invalid guide sequence')
    }

    if( length( setdiff(position_table$dna, allowed_characters) ) > 0 ){
        stop('Invalid target sequence')
    }

    if( nrow( filter(position_table, rna == '-' & dna == '-') > 1 ) ){
        stop('Invalid alignment, all positions must be either match, mismatch, DNA bulge, or RNA bulge.')
    }

    if( nrow( filter(position_table, ( rna == '-' | dna == '-' ) & index == 1 ) > 1 ) ){
        stop('Invalid alignment, CFD is undefined for a bulge at position 1')
    }
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

#' CFD-Optimal Alignment
#'
#' Uses a modified Needleman-Wunsch algorithm to calculate the alignment with the maximum CFD score.
#'
#' @param rna Character representation of the guide. Should be represented as DNA, valid characters include ['A','C','G','T']. Must be 20 nucleotides long.
#' @param dna Character representation of the target target sequences to search for an alignment. Should be represented as DNA, valid characters include ['A','C','G','T']. No length requirement.
#' @param allow_bulge Length-1 logical vector indicating whether or not the alignment algorith is allowed to generate alignments that contain bulges.
#' @param search_both_strands Length-1 logical vector indicating whether or not the reverse-complement sequence should be searched for an alignment with the guide.
#' @return A data.frame will be returned with one row for each genome sequence provided, containing the optimal alignment and CFD score, and information about the location of the alignment.
#' @name optimal_alignment
#' @export
optimal_alignment <- function(query, genome, allow_bulge = TRUE, search_both_strands = TRUE) {
    private_optimal_alignment(package_state$activity_scores, query, genome, allow_bulge[1], search_both_strands[1])
}

#' Fast CFD Score
#'
#' Calls a Cpp implementation of the bulged CFD score calculation.
#'
#' @param alignment Data Frame with Character Vector columns ['guide', 'genome', 'pam'].
#' @param strict A length-1 logical vector indicating whether the guide must be exactly 20 nt long. If not strict, the positions are labeled relative to the nucleotide closest to the PAM, which is identified as position 20.
#' @return A data.frame will be returned with one row for each genome sequence provided, containing the optimal alignment and CFD score, and information about the location of the alignment.
#' @name fast_cfd_score
#' @export
fast_cfd_score <- function( alignment, strict = TRUE ) {
    private_cfd_score( package_state$activity_scores, alignment$guide, alignment$genome, alignment$pam, strict )
}
