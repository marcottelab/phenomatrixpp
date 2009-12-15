/* 
 * File:   knearest.cpp
 * Author: jwoods
 * 
 * Created on July 23, 2009, 11:41 AM
 */

#include "knearest.h"

knearest::knearest(list<marshall>& mar, marshall& predmar, size_t kval, const set<gene_id_t>& test_set, const multimap<gene_id_t, phene_id_t>& cell_test_set, const string& identifier, const string& distance_function)
: oracle(mar, predmar, test_set, cell_test_set, identifier, distance_function, true),
        sorted(dest->phene_vector().size(), dist_phene_multiset_t()),
        k_(kval)
{
       // Get the two items over which we will iterate
    const vector<phene_id_t>& phenes = dest->phene_vector();
    assert(phenes[0] != phenes[1]);

    cout << "Sorting distance matrices." << endl;
    // Store sortings in advance in a vector.
    for (size_t j = 0; j < sorted.size(); ++j) {
        cout << "Current position: " << j << "/" << phenes.size() << endl;
        //sorted[j] = sorting(p);
        add_to_sorting(sorted[j], phenes[j], kval, 1);
    }

    // Set the method and column strings, which may vary depending on k
    if (kval == 0) {
        method_str = "all non-infinite-distance neighbors";
        column_str = "nb-all";
    } else {
        ostringstream meth;
        meth << "k = " << kval << " nearest non-infinite distance neighbors";
        method_str = meth.str();
        
        ostringstream col;
        col << "nb-knn" << kval;
        column_str = col.str();
    }
}


// Saves the sorted distance matrix.
bool knearest::save_sorted(const vector<phene_id_t>& j_to_phene, const string& filename, size_t k) const {
/*    ostringstream kfilename;
    kfilename << filename << "." << k;
    string kfilename_str = kfilename.str();

    if (exists(kfilename_str))
        return false;

    ofstream fout(kfilename_str.c_str());

    // Save with full precision and in scientific notation.
    fout.precision(dbprec1::digits10);
    fout.setf(std::ios::scientific,std::ios::floatfield);
    
    for (size_t j = 0; j < sorted.size(); ++j) {
        // Convert to phene
        phene_id_t p = j_to_phene[j];
        // Now write to file.

        // Write phenotype (destination) and the size of the multiset.
        fout << p << '\t' << sorted[j].size() << endl;

        for (multiset<dist_phene_t>::const_iterator dp = sorted[j].begin(); dp != sorted[j].end(); ++dp)
            // Now write phenotype (source)\tdistance
            fout << dp->second() << '\t' << dp->first() << endl;

    }

    fout.close(); */
    cerr << "Remember that you commented out save_sorted." << endl;

    return true;
}


bool knearest::load_sorted(const vector<phene_id_t>& j_to_phene, const string& filename, size_t k) {
    //ostringstream kfilename;
    //kfilename << filename << "." << k;
    string kfilename_str = filename; //kfilename.str();

    if (!exists(kfilename_str))
        return false;
    
    // Build a map to reverse the vector (which is itself an int to phene map)
    map<phene_id_t,size_t> phene_to_j;
    for (size_t j = 0; j < sorted.size(); ++j) {
        phene_to_j[j_to_phene[j]] = j;
    }

    // Read from the file
    ifstream fin(kfilename_str.c_str());

    phene_id_t to_phene = 0;
    size_t num_items = 0;
    fin >> to_phene;
    fin >> num_items;
    while (fin) {

        size_t j = phene_to_j[to_phene];

        multiset<dist_phene_t>::iterator sorted_j_it = sorted[j].begin();

        // Read in one line per item
        for (size_t n = 0; n < num_items; ++n) {
            phene_id_t from_phene = 0;
            dist_t distance = 0;

            fin >> from_phene;
            fin >> distance;

            // Insert and get an iterator to help next time.
            sorted_j_it = sorted[j].insert(sorted_j_it, dist_phene_t(distance, from_phene));
        }

        // Priming read.
        fin >> to_phene;
        fin >> num_items;
    }

    fin.close();

    cout << "Loaded sorted distance matrix for knearest." << endl;
    
    return true;
}


/* knearest::knearest(const knearest& orig) {
} */


knearest::~knearest() {
}


// Predict for all phenotypes in the dest matrix.
void knearest::predict() {

    // Get the two items over which we will iterate
    const vector<phene_id_t>& phenes = dest->phene_vector();
    const vector<gene_id_t>&  genes  = dest->gene_vector();

    cout << "Beginning Naive Bayes calculations." << endl;
    // FORMULA FOR NAIVEBAYES:
    // P(gene is involved in disease X | all phenologs involving disease X) =
    //   1 - PRODUCT(for phenologs involving disease X)
    //              ( 1 - P(gene involved in disease X | phenolog is correct) * P(phenolog is correct) )
    cout << "Total genes: " << genes.size() << endl << endl;
    progress_display show_progress(genes.size());

    for (size_t i = 0; i < genes.size(); ++i) {
        const gene_id_t& g = genes[i];
        //cout << "Looking at gene " << genes[i] << " (" << i << " of " << genes.size() << ")" << endl;

        for (size_t j = 0; j < phenes.size(); ++j) {
            //cout << "Looking at phene " << phenes[j] << " (" << j << " of " << phenes.size() << ")" << endl;
            const phene_id_t& to_phene = phenes[j];

            double cumul = 1.0;
            
            // the closest phenotypes to dest_p
            // get them:
            const multiset<dist_phene_t>& p_sorting = sorted[j];
            size_t kpos = 0;
            dist_t last_dist = 0.0;
            for (multiset<dist_phene_t>::const_iterator kt = p_sorting.begin(); kt != p_sorting.end(); ++kt) {
                const dist_t& dist              = kt->first();
                const phene_id_t& from_phene    = kt->second();
                phene_sp_map_t::const_iterator ps = phene_sp.find(from_phene);
                if (ps == phene_sp.end()) {
                    cerr << "Error: Could not find species for phenotype " << from_phene << endl; throw;
                }
                sp_genephene_map_t::const_iterator sg = genephenes.find(ps->second);
                if (sg == genephenes.end()) {
                    cerr << "Error: Could not find genephene for species " << ps->second << endl; throw;
                }

                if (sg->second->has_gene(g) && sg->second->has_phene(from_phene) && sg->second->has_association(g, from_phene)) {
                    // Probability that some gene is involved in disease X given that the phenolog is correct.
                    // In the equation, this is the k/n term.
                    double k_over_n = sg->second->common_items(to_phene,from_phene) / double(sg->second->num_genes_with_phene(from_phene));

                    cumul *= (1.0 - k_over_n * (1.0 - dist));
                }

                if (k_ > 0 && kpos >= k_ && dist > last_dist)
                    break; // Done looping, reached k
                else {
                    last_dist = dist;
                    kpos++;
                }
            }
            // Set association to be 1.0 - cumul
            dest->add_to_association(g, to_phene, 1.0 - cumul);

        }

        ++show_progress;
    }
}
