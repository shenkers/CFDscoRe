test_that("activity scores load correctly", {
    expect_equal( length( package_state$activity_scores ) , 4 )
    expect_true( length( setdiff( names(package_state$activity_scores), c('dna_bulge', 'rna_bulge', 'mismatch', 'pam') ) ) == 0 )
})

test_that("obtain same CFD scores as from Supp Table 20", {
    expect_equal( cfd_score('AGGGGCCGGAGTATTGGGAG','AGGGGCCCGAGTATTGGGAG','GG'), 0.615384615 )
    expect_equal( cfd_score('TCCACAGATACCTGAAGAAC','TCCGCAGATACCTGAAGAAC','GG'), 0.625 )
    expect_equal( cfd_score('CCAGGGCGGCCGCCAGCAGC','CCAGGGCGGCCGCCAACAGC','GG'), 1 )
    expect_equal( cfd_score('ATGTCTCTCCGAGATTGTAA','CTGTCTCTCCCAGATTGTAA','GG'), 0.367346939 )
    expect_equal( cfd_score('CATGCCGTGTGTACCATGAG','CATGCCATGTGTACCATCAG','GG'), 0.476190476 )
    expect_equal( cfd_score('CTACCTGGAGGGCGAGTGCG','CTACCTGGAGGGCACGTGCG','GG'), 0.204545455 )
    expect_equal( cfd_score('CGACGCAAGTGGGAGCAGAG','AAACACAAGTGGGAGCAGGC','GG'), 0.117857143 )
    expect_equal( cfd_score('GGCCATCATTGGAGCTGTGG','AATAGTCACTGGAGCTGTGG','GG'), 0.293022734 )
    expect_equal( cfd_score('CATTACAAGGCCTACCTGGA','AGACTCAGGGCCTACCTGGA','GG'), 0.09859944 )
})

test_that("validation catches bad input", {
    # guide too long (>20)
    expect_error( cfd_score('ACGGGGCCGGAGTATTGGGAG','AGGGGCCCGAGTATTGGGAG','GG') )
    # guide too short (<20)
    expect_error( cfd_score('GGGGCCGGAGTATTGGGAG','AGGGGCCCGAGTATTGGGAG','GG') )
    # non allowed guide character
    expect_error( cfd_score('GXGGGCCGGAGTATTGGGAG','AGGGGCCCGAGTATTGGGAG','GG') )
    # non allowed target character
    expect_error( cfd_score('CATTACAAGGCCTACCTGGA','AGZCTCAGGGCCTACCTGGA','GG') )
    # invalid alignment
    expect_error( cfd_score('CA-TTACAAGGCCTACCTGGA','AG-ACTCAGGGCCTACCTGGA','GG') )
    # undefined for gap at position 1
    expect_error( cfd_score('C-ATTACAAGGCCTACCTGGA','-AGACTCAGGGCCTACCTGGA','GG') )
    # undefined for gap at position 1
    expect_error( cfd_score('-CATTACAAGGCCTACCTGGA','A-GACTCAGGGCCTACCTGGA','GG') )
    # unequal length
    expect_error( cfd_score('CATTACAAGGCCTACCTGGA','A-GACTCAGGGCCTACCTGGA','GG') )
    # gap in pam
    expect_error( cfd_score('CATTACAAGGCCTACCTGGA','AGACTCAGGGCCTACCTGGA','G-G') )
    # PAM too short
    expect_error( cfd_score('CATTACAAGGCCTACCTGGA','AGACTCAGGGCCTACCTGGA','G') )
    # PAM too long
    expect_error( cfd_score('CATTACAAGGCCTACCTGGA','AGACTCAGGGCCTACCTGGA','TGG') )
    # PAM invalid characters
    expect_error( cfd_score('CATTACAAGGCCTACCTGGA','AGACTCAGGGCCTACCTGGA','NG') )
})

test_that("PAM proximal mismatches result in lower scores in context of DNA bulge", {
    expect_true(
        cfd_score(
            'CA-TGCCGTGTGTACCATGAG',
            'CAGTGCCATGTGTACCATCAG',
            'GG') >
        cfd_score(
            'CA-TGCCGTGTGTACCATGAG',
            'CAGTGCCATGTGTACCATCAC',
            'GG')
    )
})
