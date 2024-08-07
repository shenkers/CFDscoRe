#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <Rcpp.h>

using namespace std;

enum class Traceback { Match, Insert, Delete };

class Cas9Alignment {
public:

    Cas9Alignment(){}

    Cas9Alignment( string& GUIDE, string& TARGET, string& PAM, double LOG_SCORE, int OFFSET, int TARGET_LENGTH )
        : guide(GUIDE), target(TARGET), pam(PAM), log_score(LOG_SCORE), offset(OFFSET), target_length(TARGET_LENGTH) {}

    string guide;
    string target;
    string pam;
    double log_score;
    int offset;
    int target_length;
    string strand;
};

vector<map<pair<string,string>,double>> mismatch_table;
vector<map<string,double>> insert_table;
vector<map<string,double>> delete_table;
map<string,double> pam_table;

string complement( string dna ) {
    if( dna == "A" ) return "T";
    else if( dna == "T" ) return "A";
    else if( dna == "C" ) return "G";
    else return "C";
}

string reverse_complement( string dna ) {
    string rc = "";
    for( int i=dna.length(); i>0; i--){
        rc += complement(string(1,dna[i-1]));
    }
    return rc;
}

inline double cfd_mismatch_score( int i, string rna, string dna) {
    return mismatch_table[i-1][make_pair(rna, complement(dna) )];
}

inline double cfd_insert_score( int i, string rna) {
    if( i > 20 )
        return -1000;
    return insert_table[i-1][rna];
}

inline double cfd_delete_score( int i, string dna) {
    if( i >= 20 )
        return -1000;
    return delete_table[i][complement(dna)];
}

class Cas9Aligner {
    public:

    Cas9Aligner(string& GUIDE, string& TARGET)
        : RNA(GUIDE), FULL_DNA(TARGET) {
        DNA = FULL_DNA.substr(0, FULL_DNA.length() - 3);
        int n = RNA.length();
        int m = DNA.length();
        prefix_score = vector<vector<double>>( n+1, vector<double>( m+1, 0) );
        traceback = vector<vector<Traceback>>( n+1, vector<Traceback>( m+1 ) );
    }

inline double score_match_or_mismatch( int i, int j, string rna, string dna) {
    if ( rna == dna ) {
        return prefix_score[i-1][j-1]; // no penalty if matches ( equivalent to previous + 0 )
    } else {
        return prefix_score[i-1][j-1] + cfd_mismatch_score( i, rna, dna );
    }
}

// rna bulge
inline double score_insert_pos( int i, int j, string rna) {
    return prefix_score[i-1][j] + cfd_insert_score(i, rna);
}

// dna bulge
inline double score_delete_pos( int i, int j, string dna) {
    return prefix_score[i][j-1] + cfd_delete_score(i, dna);
}

inline double needleman_wunsch(bool allow_bulge)
{
    int n = RNA.length();
    int m = DNA.length();
    for (int i=1;i<=n;i++) {
        // if the score is constant on the margin (rather than increase
        // with the offset) it won't penalize insertions/deletions
        prefix_score[i][0] = -DBL_MAX;
        traceback[i][0] = Traceback::Insert;
    }
    if( allow_bulge ) {
        for (int i=1;i<=n;i++)
        {
            for (int j=1;j<=m;j++)
            {
                string rna = string(1, RNA[i-1]);
                string dna = string(1, DNA[j-1]);

                double score_match = score_match_or_mismatch(i, j, rna, dna);
                double score_insert = score_insert_pos(i, j, rna);
                double score_delete = score_delete_pos(i, j, dna);

                if( score_match >= score_insert && score_match >= score_delete ) {
                    prefix_score[i][j] = score_match;
                    traceback[i][j] = Traceback::Match;
                } else if( score_insert > score_delete ) {
                    prefix_score[i][j] = score_insert;
                    traceback[i][j] = Traceback::Insert;
                } else {
                    prefix_score[i][j] = score_delete;
                    traceback[i][j] = Traceback::Delete;
                }
            }
        }
    } else {
        for (int i=1;i<=n;i++)
        {
            for (int j=1;j<=m;j++)
            {
                string rna = string(1, RNA[i-1]);
                string dna = string(1, DNA[j-1]);

                double score_match = score_match_or_mismatch(i, j, rna, dna);

                prefix_score[i][j] = score_match;
                traceback[i][j] = Traceback::Match;
            }
        }
    }
    return prefix_score[n][m];
}

struct traceback_init {
    int rna_i;
    int dna_j;
    double log_score;
    string pam;
};

traceback_init get_traceback_start(){
    traceback_init start = { -1, -1, -DBL_MAX, "" };
    double max_score = -DBL_MAX;
    int n = RNA.length();
    for( int j=0; j <= DNA.length(); j++ ){
        string pam = FULL_DNA.substr(j+1,2);
        double pam_score = pam_table[pam];
        double score_j =  prefix_score[n][j] + pam_score;
        if( score_j > start.log_score ) {
            start = { n, j, score_j, pam };
        }
    }
    return start;
}

inline Cas9Alignment get_optimal_alignment() {
    int n = RNA.length();
    int m = DNA.length();
    string alignmentRNA, alignmentDNA;
    stack<char> tracebackRNA, tracebackDNA;
    traceback_init start = get_traceback_start();
    if( start.rna_i < 0 )
        throw runtime_error("No non-zero CFD alignment exists");
    int rna_i = start.rna_i, dna_j = start.dna_j;
    while (rna_i != 0)
    {
        if (dna_j == 0)
        {
            tracebackRNA.push(RNA[rna_i-1]);
            tracebackDNA.push('-');
            rna_i--;
        }
        else
        {
            if (traceback[rna_i][dna_j] == Traceback::Match)
            {
                tracebackRNA.push(RNA[rna_i-1]);
                tracebackDNA.push(DNA[dna_j-1]);
                rna_i--; dna_j--;
            }
            else if (traceback[rna_i][dna_j] == Traceback::Insert) {
                tracebackRNA.push(RNA[rna_i-1]);
                tracebackDNA.push('-');
                rna_i--;
            }
            else
            {
                tracebackRNA.push('-');
                tracebackDNA.push(DNA[dna_j-1]);
                dna_j--;
            }
        }
    }
    while (!tracebackRNA.empty())
    {
        alignmentRNA += tracebackRNA.top();
        alignmentDNA += tracebackDNA.top();
        tracebackRNA.pop();
        tracebackDNA.pop();
    }
    return Cas9Alignment( alignmentRNA, alignmentDNA, start.pam, start.log_score, dna_j, start.dna_j - dna_j );
}

    string RNA;
    string DNA;
    string FULL_DNA;
    vector<vector<double>> prefix_score;
    vector<vector<Traceback>> traceback;
};

vector<map<string,double>> load_indel_table( Rcpp::DataFrame data_frame ) {

    vector<map<string,double>> table(20);

    Rcpp::IntegerVector indices = data_frame["index"];
    Rcpp::CharacterVector xnas = data_frame[0];
    Rcpp::NumericVector activity = data_frame["activity"];

    for( int i=0; i< data_frame.nrows(); i++ ){

        int index = indices[i];
        string xna = Rcpp::as<string>(xnas[i]);
        float score = activity[i];

        table[index-1][xna] = log(score);
    }

    vector<string> nt = { "A", "C", "G", "T" };
    for( string dna : nt )
        table[0][dna] = -1000;

    return table;
}

vector<map<pair<string,string>,double>> load_mismatch_table( Rcpp::DataFrame data_frame ) {
    vector<map<pair<string,string>,double>> table(20);

    Rcpp::IntegerVector indices = data_frame["index"];
    Rcpp::CharacterVector rnas = data_frame["rna"];
    Rcpp::CharacterVector dnas = data_frame["dna"];
    Rcpp::NumericVector activity = data_frame["activity"];

    for( int i=0; i< data_frame.nrows(); i++ ){

        int index = indices[i];
        string rna = Rcpp::as<string>(rnas[i]);
        string dna = Rcpp::as<string>(dnas[i]);
        double score = activity[i];

        table[index-1][make_pair(rna,dna)] = log(score);
    }

    return table;
}

map<string,double> load_pam_table( Rcpp::DataFrame data_frame ) {

    map<string,double> table;

    Rcpp::CharacterVector pams = data_frame["pam"];
    Rcpp::NumericVector activity = data_frame["activity"];

    for( int i=0; i< data_frame.nrows(); i++ ){

        string pam = Rcpp::as<string>(pams[i]);
        double score = activity[i];

        table[pam] = log(score);
    }

    return table;
}

double score_alignment(string guide, string genome, string pam, bool strict) {
    int n = guide.length();
    int guide_position = 20;
    double log_score = pam_table[pam];
    for (int alignment_column=n-1; alignment_column>-1; alignment_column--) {
        string rna = string(1, guide[alignment_column]);
        string dna = string(1, genome[alignment_column]);

        if( !( rna == "-" || rna == "A" || rna == "C" || rna == "G" || rna == "T" ) )
            throw runtime_error("CFD undefined for this alignment");
        if( !( dna == "-" || dna == "A" || dna == "C" || dna == "G" || dna == "T" ) )
            throw runtime_error("CFD undefined for this alignment");
        if( rna == "-" && dna == "-" )
            throw runtime_error("CFD undefined for this alignment");
        if( guide_position == 1 && rna == "-" )
            throw runtime_error("CFD undefined for this alignment");

        if( rna != dna ) {
            if( dna == "-" ){
               log_score += cfd_insert_score(guide_position, rna);
            } else if( rna == "-" ) {
               log_score += cfd_delete_score(guide_position - 1, dna);
            } else {
               log_score += cfd_mismatch_score(guide_position, rna, dna);
            }
        }

        if( alignment_column > 0 && string(1,guide[alignment_column-1]) != "-" ){
            guide_position--;
            if( guide_position < 1 )
                throw runtime_error("CFD undefined for this alignment");
        }
    }

    if( strict && guide_position != 1 )
        throw runtime_error("CFD undefined for this alignment");

    return log_score;
}

Cas9Alignment optimal_target(string guide, string genome, bool allow_bulge){

    Cas9Aligner aligner = Cas9Aligner( guide, genome );

    double log_score = aligner.needleman_wunsch(allow_bulge);
    Cas9Alignment cas9alignment = aligner.get_optimal_alignment();

    return cas9alignment;
}

Cas9Alignment optimal_fwd_rev_target(string guide, string genome, bool allow_bulge ){
    double max_cfd = -DBL_MAX;
    Cas9Alignment optimal;

    try {
        Cas9Alignment fwd = optimal_target( guide, genome, allow_bulge );
        fwd.strand = "+";
        optimal = fwd;
        max_cfd = fwd.log_score;
    } catch( exception& e ) { }

    try {
        Cas9Alignment rev = optimal_target( guide, reverse_complement(genome), allow_bulge );
        rev.strand = "-";
        if( rev.log_score > max_cfd ) {
            optimal = rev;
            max_cfd = rev.log_score;
        }
    } catch( exception& e ) { }

    if( ! ( max_cfd > -DBL_MAX ) )
        throw runtime_error("No non-zero CFD alignment exists on either strand");

    return optimal;
}

Cas9Alignment optimal_fwd_target(string guide, string genome, bool allow_bulge ){
    double max_cfd = -DBL_MAX;
    Cas9Alignment optimal;

    try {
        Cas9Alignment fwd = optimal_target( guide, genome, allow_bulge );
        fwd.strand = "+";
        optimal = fwd;
        max_cfd = fwd.log_score;
    } catch( exception& e ) { }

    if( ! ( max_cfd > -DBL_MAX ) )
        throw runtime_error("No non-zero CFD alignment exists on either strand");

    return optimal;
}

struct EditDistance {
    int edit_distance;
    int n_mismatch;
    int n_rna_bulge;
    int n_dna_bulge;
    int n_pam_mismatch;
};

EditDistance get_edit_distance( string guide, string target, string pam ) {
    EditDistance distance = { 0, 0, 0, 0, 0 };
    for( int i=0; i < guide.length(); i++ ) {
        string rna = string(1,guide[i]);
        string dna = string(1,target[i]);
        if( rna != dna ) {
            if( dna == "-" ) distance.n_rna_bulge++;
            else if( rna == "-" ) distance.n_dna_bulge++;
            else distance.n_mismatch++;
        }
    }
    if( string(1,pam[0]) != "G" ) distance.n_pam_mismatch++;
    if( string(1,pam[1]) != "G" ) distance.n_pam_mismatch++;
    distance.edit_distance = distance.n_mismatch + distance.n_rna_bulge + distance.n_dna_bulge + distance.n_pam_mismatch;
    return distance;
}

//
//' CFD Score
//'
//' Cpp implementation of the bulged CFD score calculation
//'
//' @param rna Character representation of the guide. Should be represented as DNA, valid characters include ['A','C','G','T']. Must be 20 nucleotides long.
//' @param dna Character representation of the target target sequences to search for an alignment. Should be represented as DNA, valid characters include ['A','C','G','T']. No length requirement.
//' @return A data.frame will be returned with one row for each genome sequence provided, containing the optimal alignment and CFD score, and information about the location of the alignment.
//' @name private_cfd_score
// [[Rcpp::export]]
Rcpp::List private_cfd_score( Rcpp::List activity_scores, Rcpp::CharacterVector guide, Rcpp::CharacterVector genome, Rcpp::CharacterVector pam, Rcpp::LogicalVector strict ) {

    mismatch_table = load_mismatch_table( Rcpp::as<Rcpp::DataFrame>(activity_scores["mismatch"]) );
    insert_table = load_indel_table( Rcpp::as<Rcpp::DataFrame>(activity_scores["rna_bulge"]) );
    delete_table = load_indel_table( Rcpp::as<Rcpp::DataFrame>(activity_scores["dna_bulge"]) );
    pam_table = load_pam_table( Rcpp::as<Rcpp::DataFrame>(activity_scores["pam"]) );

    int l = genome.length();

    Rcpp::NumericVector score( l );
    Rcpp::IntegerVector edit_distance( l );
    Rcpp::IntegerVector n_mismatch( l );
    Rcpp::IntegerVector n_rna_bulge( l );
    Rcpp::IntegerVector n_dna_bulge( l );
    Rcpp::IntegerVector n_pam_mismatch( l );

    score.fill( score.get_na() );
    edit_distance.fill( edit_distance.get_na() );
    n_mismatch.fill( n_mismatch.get_na() );
    n_rna_bulge.fill( n_rna_bulge.get_na() );
    n_dna_bulge.fill( n_dna_bulge.get_na() );
    n_pam_mismatch.fill( n_pam_mismatch.get_na() );

    for( int i=0; i< l; i++ ) {
        try {
            double cfd = score_alignment( Rcpp::as<string>(guide[i]), Rcpp::as<string>(genome[i]), Rcpp::as<string>(pam[i]), Rcpp::is_true(Rcpp::all(strict)) );
            score[i] = exp(cfd);
            EditDistance distance = get_edit_distance( Rcpp::as<string>(guide[i]), Rcpp::as<string>(genome[i]), Rcpp::as<string>(pam[i]) );
            edit_distance[i] = distance.edit_distance;
            n_mismatch[i] = distance.n_mismatch;
            n_rna_bulge[i] = distance.n_rna_bulge;
            n_dna_bulge[i] = distance.n_dna_bulge;
            n_pam_mismatch[i] = distance.n_pam_mismatch;
        } catch( exception& e ) {
            score[i] = score.get_na();
        }
    }


    Rcpp::DataFrame result = Rcpp::DataFrame::create(
        Rcpp::Named("guide") = guide,
        Rcpp::Named("target") = genome,
        Rcpp::Named("pam") = pam,
        Rcpp::Named("score") = score,
        Rcpp::Named("edit_distance") = edit_distance,
        Rcpp::Named("n_mismatch") = n_mismatch,
        Rcpp::Named("n_rna_bulge") = n_rna_bulge,
        Rcpp::Named("n_dna_bulge") = n_dna_bulge,
        Rcpp::Named("n_pam_mismatch") = n_pam_mismatch
    );

    result.attr("class") = Rcpp::CharacterVector::create("tbl_df","tbl","data.frame");

    return result;
}

//' CFD-Optimal Alignment
//'
//' Uses a modified Needleman-Wunsch algorithm to calculate the alignment with the maximum CFD score.
//'
//' @param rna Character representation of the guide. Should be represented as DNA, valid characters include ['A','C','G','T']. Must be 20 nucleotides long.
//' @param dna Character representation of the target target sequences to search for an alignment. Should be represented as DNA, valid characters include ['A','C','G','T']. No length requirement.
//' @return A data.frame will be returned with one row for each genome sequence provided, containing the optimal alignment and CFD score, and information about the location of the alignment.
//' @name private_optimal_alignment
// [[Rcpp::export]]
Rcpp::List private_optimal_alignment( Rcpp::List activity_scores, Rcpp::CharacterVector query, Rcpp::CharacterVector genome, Rcpp::LogicalVector allow_bulge, Rcpp::LogicalVector search_both_strands ) {

    mismatch_table = load_mismatch_table( Rcpp::as<Rcpp::DataFrame>(activity_scores["mismatch"]) );
    insert_table = load_indel_table( Rcpp::as<Rcpp::DataFrame>(activity_scores["rna_bulge"]) );
    delete_table = load_indel_table( Rcpp::as<Rcpp::DataFrame>(activity_scores["dna_bulge"]) );
    pam_table = load_pam_table( Rcpp::as<Rcpp::DataFrame>(activity_scores["pam"]) );

    int l = genome.length();

    Rcpp::CharacterVector guide( l );
    Rcpp::CharacterVector target( l );
    Rcpp::CharacterVector pam( l );
    Rcpp::NumericVector score( l );
    Rcpp::IntegerVector edit_distance( l );
    Rcpp::IntegerVector n_mismatch( l );
    Rcpp::IntegerVector n_rna_bulge( l );
    Rcpp::IntegerVector n_dna_bulge( l );
    Rcpp::IntegerVector n_pam_mismatch( l );
    Rcpp::IntegerVector offset( l );
    Rcpp::IntegerVector target_length( l );
    Rcpp::CharacterVector strand( l );

    guide.fill( guide.get_na() );
    target.fill( target.get_na() );
    pam.fill( pam.get_na() );
    score.fill( score.get_na() );
    edit_distance.fill( edit_distance.get_na() );
    n_mismatch.fill( n_mismatch.get_na() );
    n_rna_bulge.fill( n_rna_bulge.get_na() );
    n_dna_bulge.fill( n_dna_bulge.get_na() );
    n_pam_mismatch.fill( n_pam_mismatch.get_na() );
    offset.fill( offset.get_na() );
    target_length.fill( target_length.get_na() );
    strand.fill( strand.get_na() );

    for( int i=0; i< l; i++ ) {
        try {
            if( Rcpp::as<string>(genome[i]).length() < 5 )
                throw runtime_error("CFD is undefined for sequences with length less than 5.");

            Cas9Alignment optimal;
            if( Rcpp::is_true(Rcpp::all(search_both_strands)) )
                optimal = optimal_fwd_rev_target( Rcpp::as<string>(query[0]), Rcpp::as<string>(genome[i]), Rcpp::is_true(Rcpp::all(allow_bulge)) );
            else
                optimal = optimal_fwd_target( Rcpp::as<string>(query[0]), Rcpp::as<string>(genome[i]), Rcpp::is_true(Rcpp::all(allow_bulge)) );
            guide[i] = optimal.guide;
            target[i] = optimal.target;
            pam[i] = optimal.pam;
            score[i] =  exp(optimal.log_score);
            EditDistance distance = get_edit_distance( optimal.guide, optimal.target, optimal.pam );
            edit_distance[i] = distance.edit_distance;
            n_mismatch[i] = distance.n_mismatch;
            n_rna_bulge[i] = distance.n_rna_bulge;
            n_dna_bulge[i] = distance.n_dna_bulge;
            n_pam_mismatch[i] = distance.n_pam_mismatch;
            offset[i] = optimal.offset;
            target_length[i] = optimal.target_length;
            strand[i] = optimal.strand;
        } catch( exception& e ) {
            guide[i] = guide.get_na();
            target[i] = target.get_na();
            pam[i] = pam.get_na();
            score[i] = score.get_na();
            offset[i] = offset.get_na();
            target_length[i] = target_length.get_na();
            strand[i] = strand.get_na();
        }
    }

    Rcpp::DataFrame result = Rcpp::DataFrame::create(
        Rcpp::Named("guide") = guide,
        Rcpp::Named("target") = target,
        Rcpp::Named("pam") = pam,
        Rcpp::Named("score") = score,
        Rcpp::Named("edit_distance") = edit_distance,
        Rcpp::Named("n_mismatch") = n_mismatch,
        Rcpp::Named("n_rna_bulge") = n_rna_bulge,
        Rcpp::Named("n_dna_bulge") = n_dna_bulge,
        Rcpp::Named("n_pam_mismatch") = n_pam_mismatch,
        Rcpp::Named("offset") = offset,
        Rcpp::Named("target_length") = target_length,
        Rcpp::Named("strand") = strand
    );

    result.attr("class") = Rcpp::CharacterVector::create("tbl_df","tbl","data.frame");

    return result;
}

int main()
{
    printf("lookup: %.2f\n", insert_table[1]["A"]);
    printf("lookup3: %.2f\n", delete_table[1]["A"]);
    printf("lookup2: %.2f\n", mismatch_table[1][make_pair("A","C")]);
//    n = 5, m = 7;
//    A = "ACACT";
//    DNA = "ACGACTG";
//    mismatch only
//    RNA = "CATGCCGTGTGTACCATGAC";
//    DNA = "CAGGCCATGTGTACCATCAG";
//    indel and mismatch
//    RNA = "CATGCCGTGTGTACCATGAC";
//    DNA = "CAGTGCCATGTGTACCATCAG";
    string RNA = "CGTGCCATGTGTACCATGAG";
    string DNA = "CGGCCATGTGTACCATCGAG";

    Cas9Aligner al1 = Cas9Aligner( RNA, DNA );
    printf("score: %f\n",exp(al1.needleman_wunsch(true)));
    return 0;
}
