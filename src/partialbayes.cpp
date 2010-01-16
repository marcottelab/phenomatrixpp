/* 
 * File:   partialbayes.cpp
 * Author: jwoods
 * 
 * Created on August 1, 2009, 7:35 PM
 */

#include "partialbayes.h"

partialbayes::partialbayes(list<marshall>& mar, marshall& predmar, size_t kval, const set<gene_id_t>& test_set, const multimap<gene_id_t, phene_id_t>& cell_test_set, const string& identifier, const string& distance_function)
: knearest(mar, predmar, kval, test_set, cell_test_set, identifier, distance_function)
{
    // Set the method and column strings, which may vary depending on k
    if (kval == 0) {
        method_str = "all non-infinite-distance neighbors without k/n term";
        column_str = "pb-all";
    } else {
        ostringstream meth;
        meth << "k = " << kval << " nearest non-infinite distance neighbors without k/n term";
        method_str = meth.str();

        ostringstream col;
        col << "pb-knn" << kval;
        column_str = col.str();
    }
}

partialbayes::~partialbayes() {
}


// Predict for all phenotypes in the dest matrix.
void partialbayes::predict() {

    // Get the two items over which we will iterate
    const vector<phene_id_t>& phenes = dest->phene_vector();
    const vector<gene_id_t>&  genes  = dest->gene_vector();

    cout << "Beginning partial Naive Bayes calculations." << endl;
    // FORMULA FOR NAIVEBAYES:
    // P(gene is involved in disease X | all phenologs involving disease X) =
    //   1 - PRODUCT(for phenologs involving disease X)
    //              ( 1 - P(gene involved in disease X | phenolog is correct) * P(phenolog is correct) )
    // PARTIALBAYES leaves out the k/n term, ^^^^^ right above this ^^^^^
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
                    //double k_over_n = sg->second->common_items(to_phene,from_phene) / double(sg->second->num_genes_with_phene(from_phene));
                    //double k_over_m = sg->second->common_items(to_phene,from_phene) / double(sg->second->num_genes_with_phene(to_phene));
                    cumul *= dist;
                    // cumul *= (1.0 - k_over_n * (1.0 - dist));
                    //cumul *= (1.0 - k_over_m * (1.0-dist));
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

