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

    Cas9Alignment( string& GUIDE, string& TARGET, string& PAM, double LOG_SCORE )
        : guide(GUIDE), target(TARGET), pam(PAM), log_score(LOG_SCORE) {}

    string guide;
    string target;
    string pam;
    double log_score;
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
        : RNA(GUIDE), DNA(TARGET) {
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

inline double needleman_wunsch()
{
    int n = RNA.length();
    int m = DNA.length();
    for (int i=0;i<=n;i++) {
        // if the score is constant on the margin (rather than increase
        // with the offset) it won't penalize insertions/deletions
        prefix_score[i][0] = -1000*i;
        traceback[i][0] = Traceback::Insert;
    }
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
    return prefix_score[n][m];
}

inline pair<string, string> get_optimal_alignment() {
    int n = RNA.length();
    int m = DNA.length();
    string alignmentRNA, alignmentDNA;
    stack<char> tracebackRNA, tracebackDNA;
    int rna_i = n, dna_j = m;
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
    return make_pair(alignmentRNA, alignmentDNA);
}

    string RNA;
    string DNA;
    vector<vector<double>> prefix_score;
    vector<vector<Traceback>> traceback;
};

vector<map<string,double>> load_indel_table( string csv_path ) {
    std::ifstream file(csv_path);

    vector<map<string,double>> table(20);

    std::string line;
    // skip header
    std::getline(file, line);
    while(std::getline(file, line)) {
        std::istringstream ss(line);
        std::string token;

        // Read each column separated by comma
        int index;
        string xna;
        float score;

        if (std::getline(ss, token, ','))
            index = std::stoi(token); // Convert string to int
        if (std::getline(ss, token, ','))
            xna = token; // Extract first character
        if (std::getline(ss, token, ','))
            score = std::stof(token); // Convert string to float

        table[index-1][xna] = log(score);
        // Output the values
    }

    return table;
}

vector<map<pair<string,string>,double>> load_mismatch_table( string csv_path ) {
    std::ifstream file(csv_path);

    vector<map<pair<string,string>,double>> table(20);

    std::string line;
    // skip header
    std::getline(file, line);
    while(std::getline(file, line)) {
        std::istringstream ss(line);
        std::string token;

        // Read each column separated by comma
        int index;
        string rna;
        string dna;
        double score;

        if (std::getline(ss, token, ','))
            index = std::stoi(token); // Convert string to int
        if (std::getline(ss, token, ','))
            rna = token; // Extract first character
        if (std::getline(ss, token, ','))
            dna = token; // Extract first character
        if (std::getline(ss, token, ','))
            score = std::stod(token); // Convert string to float

        table[index-1][make_pair(rna,dna)] = log(score);
    }

    return table;
}

map<string,double> load_pam_table( string csv_path ) {
    std::ifstream file(csv_path);

    map<string,double> table;

    std::string line;
    // skip header
    std::getline(file, line);
    while(std::getline(file, line)) {
        std::istringstream ss(line);
        std::string token;

        // Read each column separated by comma
        string pam;
        double score;

        if (std::getline(ss, token, ','))
            pam = token; // Extract first character
        if (std::getline(ss, token, ','))
            score = std::stod(token); // Convert string to float

        table[pam] = log(score);
    }

    return table;
}

Cas9Alignment optimal_target(string guide, string genome){

    double max_cfd = -DBL_MAX;

    Cas9Alignment optimal_alignment;

    for(int i=0; i<genome.length()-22; i++){
        string target = genome.substr(i,20);
        string pam = genome.substr(i+21,2);

        Cas9Aligner aligner = Cas9Aligner( guide, target );

        double log_score = aligner.needleman_wunsch();
        pair<string, string> alignment = aligner.get_optimal_alignment();

        log_score += pam_table[pam];
        Cas9Alignment cas9alignment( alignment.first, alignment.second, pam, log_score );

        if( log_score > max_cfd ) {
            optimal_alignment = cas9alignment;
            max_cfd = log_score;
        }
    }

    return optimal_alignment;
}

Cas9Alignment optimal_fwd_rev_target(string guide, string genome){
    double max_cfd = -DBL_MAX;

    Cas9Alignment fwd = optimal_target( guide, genome );
    Cas9Alignment rev = optimal_target( guide, reverse_complement(genome) );

    if( fwd.log_score > rev.log_score ){
        return fwd;
    } else {
        return rev;
    }
}

// [[Rcpp::export]]
Rcpp::List optimal_alignment( Rcpp::CharacterVector vpenalty_dir, Rcpp::CharacterVector query, Rcpp::CharacterVector genome ) {
    string penalty_dir = Rcpp::as<string>(vpenalty_dir[0]);
    mismatch_table = load_mismatch_table( penalty_dir + "/mismatch.csv" );
    insert_table = load_indel_table(penalty_dir + "/rna_bulge.csv");
    delete_table = load_indel_table(penalty_dir + "/dna_bulge.csv");
    pam_table = load_pam_table(penalty_dir + "/pams.csv");
    Cas9Alignment optimal = optimal_fwd_rev_target( Rcpp::as<string>(query[0]), Rcpp::as<string>(genome[0]) );
    Rcpp::CharacterVector guide = { optimal.guide };
    Rcpp::CharacterVector target = { optimal.target };
    Rcpp::CharacterVector pam = { optimal.pam };
    Rcpp::NumericVector score = { exp(optimal.log_score) };
    return Rcpp::List::create(
        Rcpp::Named("guide") = guide,
        Rcpp::Named("target") = target,
        Rcpp::Named("pam") = pam,
        Rcpp::Named("score") = score
    );
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
    printf("score: %f\n",exp(al1.needleman_wunsch()));
    pair<string, string> alignment = al1.get_optimal_alignment();
    printf("rna: %s\ndna: %s\n", alignment.first.c_str(), alignment.second.c_str());
    printf("\ncfd_score('%s','%s','GG')\n",alignment.first.c_str(), alignment.second.c_str());
    Cas9Alignment forward = optimal_fwd_rev_target(RNA, "AGCCGATCGTGCCGTGTGTACCATGAGCGGCCATGTGTACCATCGAGGGAGCACTCATGGTACACATGGCACGTGACGAGCC");
    printf("forward\n");
    printf(" score: %f\n", exp(forward.log_score) );
    printf(" guide: %s\n", forward.guide.c_str() );
    printf("target: %s\n", forward.target.c_str() );
    printf("   pam: %s (%.2f) \n", forward.pam.c_str(), pam_table[forward.pam] );
    printf("\ncfd_score('%s','%s','GG')\n",forward.guide.c_str(), forward.target.c_str());
    return 0;
}
