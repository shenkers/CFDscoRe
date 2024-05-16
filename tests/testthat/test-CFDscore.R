test_that("verify loading of activity scores", {
    expect_equal( length( package_state$activity_scores ) , 4 )
    expect_true( length( setdiff( names(package_state$activity_scores), c('dna_bulge', 'rna_bulge', 'mismatch', 'pam') ) ) == 0 )
})
