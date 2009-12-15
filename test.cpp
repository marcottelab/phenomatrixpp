/* 
 * File:   test.cpp
 * Author: jwoods
 *
 * Created on June 19, 2009, 4:15 PM
 */

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

#include "mindist.h"
#include "avgmindist.h"
// #include "knearest.h"
#include "partialbayes.h"

namespace po = boost::program_options;


std::ostream& operator<<(std::ostream& out, const vector<string>& rhs) {

    vector<string>::const_iterator i = rhs.begin();
    out << *i; ++i;
    for ( ; i != rhs.end(); ++i)
        out << ' ' << *i;
    out << flush;
    return out;
}

// Convert comma-separated list to vector.
vector<string> make_source_species_vector(const string& commasep) {
    typedef tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(",");

    tokenizer tok(commasep, sep);
    vector<string> vec;
    for (tokenizer::iterator i = tok.begin(); i != tok.end(); ++i)
        vec.push_back(*i);

    
    return vec;
}


void save_predictions_for_crossvalidation(const oracle& o, const string& dir, const string& type, size_t num) {

    // Save to a specific directory (if a type is given).
    ostringstream predictions_directory;
    predictions_directory << dir << num;
    
    if (type == "row") {
        o.save_row_test_predictions(predictions_directory.str());
    } else if (type == "cell") {
        o.save_cell_test_predictions(predictions_directory.str());
    } else {
        cerr << "Error: Cross-validation type not recognized. Saving without test set information." << endl;
        o.save_predictions(dir);
    }
}


// Make sure they remembered to include a slash and stuff.
void normalize_directory(string& dir) {
    if (dir.size() == 0 || (dir.size() == 1 && dir == "/"))  dir = "./";
    else if (dir[dir.size()-1] != '/')                       dir += "/";
}


/*
 * MAIN: test genephene matrix class.
 */
int main(int argc, char** argv) {

    // Use Boost to handle command line arguments.
    po::options_description desc("Allowed options");
    string identifier, predict_genes_file, distance_measure, predict_root;
    desc.add_options()
            ("help,h", "produce help message")
            ("source-species,s", po::value<string>()->default_value("Hs,Mm,Dm,Ce,Sc,At"), "comma-separated list of source species")
            ("predict-species,S", po::value<string>()->default_value("Hs"), "species to predict")
            ("method,m", po::value<string>()->default_value("naivebayes"), "prediction method to use")
            ("k,k", po::value<size_t>(), "k for k-nearest neighbors")
            ("cross-validation,n", po::value<size_t>(), "n-fold cross-validation")
            ("cross-validation-type,t", po::value<string>()->default_value("row"), "type of cross-validation (row-based or cell-based")
            ("predict-phenotypes-file,p", po::value<string>(), "file containing list of phenotypes to predict (all in predict-species by default)")
            ("identifier,i", po::value<string>(&identifier)->default_value("no_id"), "identifier for a run")
            ("predict-genes-file,g", po::value<string>(&predict_genes_file)->default_value("predict_genes"), "file containing list of genes to go in the prediction matrix")
            ("distance-measure,d", po::value<string>(&distance_measure)->default_value("hypergeometric"), "distance function to use")
            ("predict-root,D", po::value<string>(&predict_root)->default_value("./"), "where to save predictions (./ or ../ recommended)")
    ;
    po::positional_options_description p;
    p.add("source-species", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    // Print the help message.
    if (vm.count("help")) {
        cerr << desc << endl;
        return 1;
    }

    // Make sure the directory for saving is acceptable.
    normalize_directory(predict_root);


    string predict_phenotypes_file;
    string predict_species;
    if (vm.count("predict-species"))
        predict_species = vm["predict-species"].as<string>();


    // If this is the default, fix it.
    if (vm.count("predict-phenotypes-file")) {
        predict_phenotypes_file = vm["predict-phenotypes-file"].as<string>();
        if (predict_phenotypes_file == "phenes.")
            predict_phenotypes_file += predict_species;
    } else {
        //predict_phenotypes_file = "phenes." + predict_species;
        // Changed for consistency with the Rails crossval app (since we use
        // predict_genes for destination, also makes it a little less confusing).
        predict_phenotypes_file = "predict_phenotypes";
    }

    cout << "vmsourcespecies is " << vm["source-species"].as<string>() << endl;
    vector<string> source_species = make_source_species_vector(vm["source-species"].as<string>());
    cout << "Source species are: " << source_species << endl;
    cout << "Predicting for: " << predict_species << endl;
    if (!vm.count("cross-validation"))
        cout << "Run identifier: " << identifier << endl;
    
    string method(vm["method"].as<string>());
    size_t kval = 0;
    if (vm.count("k")) {
        kval = vm["k"].as<size_t>();
    }
    cout << "Method: " << method;
    if (kval > 0) cout << "\tk = " << kval << endl;
    else          cout << "\tfor all non-infinite neighbors" << endl;

    if (vm.count("cross-validation")) {
        cout << vm["cross-validation"].as< size_t >() << "-fold cross-validation" << endl;
    }

    cout << "File containing list of phenotypes to predict: " << predict_phenotypes_file << endl;
    cout << "File containing list of genes to predict: " << predict_genes_file << endl;
    cout << endl;


    if (vm.count("cross-validation")) {
        size_t fold = vm["cross-validation"].as<size_t>();
        for (size_t i = 0; i < fold; ++i) {
            ostringstream testsetfile;
            testsetfile << "testset." << fold << "-" << i;
            string ident(testsetfile.str());
            
            // The 'marshall' object marshalls files into a genephene object, which it
            // creates via new. We'll let the oracle children handle allocations and
            // de-allocations; just create the marshalls from here.
            list<marshall> marshalls;
            for (vector<string>::iterator it = source_species.begin(); it != source_species.end(); ++it) {
                pair<string,string> src(*it, "genes_phenes." + *it);
                pair<string,string> dest;
                if (*it == predict_species) dest = pair<string,string>(*it, "");
                else                        dest = pair<string,string>(predict_species, "genes_phenes." + predict_species);

                // Create the marshall and add it.
                marshall m;
                m.distance_measure      = distance_measure;
                m.matrix_identifier     = "genes_phenes." + predict_species + *it;
                m.orthologs_filename    = "genes." + *it;
                m.source_species_info   = src;
                m.dest_species_info     = dest;
                m.construct();
  
                marshalls.push_back(m);
            }

            cout << "CROSS-VALIDATION: " << i << " of " << fold << endl;

            set<gene_id_t> row_test_set;
            multimap<gene_id_t, phene_id_t> cell_test_set;

            if (vm["cross-validation-type"].as<string>() == "row") {
                row_test_set = read_set<gene_id_t>(ident);
                cout << "Number of items in row test set: " << row_test_set.size() << endl;
                if (row_test_set.size() == 0) {
                    cerr << "Error: Empty row test set in file '" << ident << "'" << endl;
                    throw;
                }
            } else {
                cell_test_set = load_associations<gene_id_t,phene_id_t>(ident);
                cout << "Number of items in cell test set: " << cell_test_set.size() << endl;
                if (cell_test_set.size() == 0) {
                    cerr << "Error: Empty cell test set in file '" << ident << "'" << endl;
                    throw;
                }
            }

            // Create the marshall and add it.
            marshall predict;
            predict.distance_measure      = distance_measure;
            // predict.matrix_identifier     = "genes_phenes." + predict_species + *it;
            predict.orthologs_filename    = predict_genes_file;
            predict.phenotypes_filename   = predict_phenotypes_file;
            predict.dest_species_info     = make_pair<string,string>(predict_species, "");
            predict.construct(true);

            oracle* o = NULL;
            if (method == "naivebayes") {
                // Specify prediction matrix.
                // marshall predict(distance_measure, predict_species, predict_genes_file, predict_phenotypes_file, true);
                o = new knearest(marshalls, predict, kval, row_test_set, cell_test_set, ident, distance_measure);
            } else if (method == "partialbayes") {

                // marshall predict(distance_measure, predict_species, predict_genes_file, predict_phenotypes_file, true);
                o = new partialbayes(marshalls, predict, kval, row_test_set, cell_test_set, ident, distance_measure);
            }/* else if (method == "mindist") {
                o = new mindist(marshalls, predict, test_set, ident);
            } else if (method == "avgmindist") {
                o = new mindist(marshalls, predict, test_set, ident);
            }*/ else {
                cerr << "Error: Unrecognized method '" << method << "'." << endl;
                throw;
            }
            o->predict();
            save_predictions_for_crossvalidation(*o, predict_root + "predictions", vm["cross-validation-type"].as<string>(), i);
            delete o;
        }
    } else {

        // The 'marshall' object marshalls files into a genephene object, which it
        // creates via new. We'll let the oracle children handle allocations and
        // de-allocations; just create the marshalls from here.
        list<marshall> marshalls;
        for (vector<string>::iterator it = source_species.begin(); it != source_species.end(); ++it) {

            pair<string,string> dest(predict_species, "genes_phenes."+predict_species);
            pair<string,string> src(*it, "genes_phenes."+*it);
            if (*it == predict_species) dest.second = "";

            // Create the marshall and add it.
            marshall m;
            m.distance_measure      = distance_measure;
            m.matrix_identifier     = "genes_phenes." + predict_species + *it;
            m.orthologs_filename    = "genes." + *it;
            m.source_species_info   = src;
            m.dest_species_info     = dest;
            m.construct();

            marshalls.push_back(m);
        }

        marshall predict;
        predict.distance_measure         = distance_measure;
        predict.dest_species_info        = make_pair<string,string>("Hs", "");
        predict.orthologs_filename       = predict_genes_file;
        predict.phenotypes_filename      = predict_phenotypes_file;
        predict.construct(true);

        // empty row_test_set.
        set<gene_id_t> row_test_set;
        // empty cell_test_set
        multimap<gene_id_t,phene_id_t> cell_test_set;
        
        oracle* o = NULL;
        if (method == "naivebayes") {
            // Specify prediction matrix.
            o = new knearest(marshalls, predict, kval, row_test_set, cell_test_set, "no_id", distance_measure);
        } else if (method == "partialbayes") {
            o = new partialbayes(marshalls, predict, kval, row_test_set, cell_test_set, "no_id", distance_measure);
        }/* else if (method == "mindist") {
            o = new mindist(marshalls, predict, test_set, ident);
        } else if (method == "avgmindist") {
            o = new mindist(marshalls, predict, test_set, ident);
        }*/ else {
            cerr << "Error: Unrecognized method '" << method << "'." << endl;
            throw;
        }
        o->predict();
        o->save_predictions(predict_root + "predictions");

        delete o;
    }
    
    return 0;
}

