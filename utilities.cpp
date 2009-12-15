
#include "utilities.h"

size_t lines_in_file(const string& filename) {
    ifstream fin("wc-l_output.txt");

    size_t lines = 0;
    string file;

    fin >> lines;
    fin >> file;

    while (fin) {

        // Found it?
        if (file == filename) {
            fin.close();
            return lines;
        }

        // Not found -- keep looking.
        fin >> lines;
        fin >> file;

    }

    cerr << "Error: File '" << filename << "' not found in wc-l_output.txt. Please run 'wc -l * > wc-l_output.txt' before running this program." << endl;
    throw;
    fin.close();
    return 0;
}


// Given a species, get its phene set.
set<phene_id_t> species_phene_set(const species_id_t& species_id) {
    ostringstream filename_stream;
    filename_stream << "phenes." << species_id;
    string filename(filename_stream.str());
    if (!exists(filename)) {
        cerr << "Error: Species phene set file does not exist for species " << species_id << endl;
        throw;
    }

    cout << "Reading species phene set for '" << species_id << "'" << endl;

    ifstream fin(filename.c_str());

    phene_id_t phene;

    fin >> phene;
    set<phene_id_t> phene_set;
    set<phene_id_t>::iterator i = phene_set.insert(phene_set.begin(), phene);

    fin >> phene;
    while (fin) {
        // Add to set.
        i = phene_set.insert(i, phene);

        fin >> phene;
    }

    fin.close();

    return phene_set;
}


pair<size_t,size_t> num_genes_and_phenes_in_file(const string& filename) {
    ifstream fin(filename.c_str());

    set<gene_id_t> genes;
    set<phene_id_t> phenes;

    // ignore two lines
    string tmp;
    getline(fin,tmp);
    getline(fin,tmp);

    gene_id_t gene;
    phene_id_t phene;
    fin >> gene;
    fin >> phene;

    set<gene_id_t>::iterator i = genes.insert(gene).first;
    set<phene_id_t>::iterator j = phenes.insert(phene).first;

    fin >> gene;
    fin >> phene;

    while (fin) {
        i = genes.insert(i, gene);
        j = phenes.insert(j, phene);

        fin >> gene;
        fin >> phene;
    }

    fin.close();

    return std::pair<size_t,size_t>(genes.size(),phenes.size());
}


void read_gene_phene_line(ifstream& fin, gene_id_t& gene, phene_id_t& phene) {
    // Tokenize the line with spaces
    fin >> gene;
    fin >> phene;
}


// Given a string, count the number of tabs in it.
size_t count_tabs(const string& line) {
    size_t num = 0;
    for (size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '\t')
            ++num;
    }
    return num;
}


// Read the existing prediction file and return the number of tabs.
/*size_t read_prediction_file(const path& p, string& header_line1, string& header_line2, map<gene_id_t,string>& gene_line) {
    // If the file exists, read it into the map. We'll just add this as another column.
    if (exists(p)) {
        ifstream fin(p);
        getline(fin, header_line1);
        getline(fin, header_line2);

        gene_id_t g;
        fin >> g; // read the gene
        getline(fin, gene_line[g]); // read the rest of the line
        while (fin) {

            // keep reading until the end of the file.
            fin >> g;
            getline(fin, gene_line[g]);
        }
        fin.close();
    }

    // Pick a random gene and get the number of tabs
    return count_tabs(*(gene_line.begin()));
}*/


string make_phene_path(const string& dir, const phene_id_t& phene) {
    /*if (!exists(dir)) {
        cerr << "Error: Directory '" << dir << "' does not exist." << endl;
        return false;
    }*/

    ostringstream phene_path_stream;
    phene_path_stream << dir << "/" << phene;

    // create a path variable and let the calling function have it back as a ref.
    return phene_path_stream.str();
}


bool add_column_to_prediction_file(const string& dir, const phene_id_t& phene,
        const string& header,
        const string& column_header,
        unordered_map<gene_id_t,double> predicted)
{
    string phene_path = make_phene_path(dir, phene);
    if (!exists(phene_path)) {
        cerr << "Cannot add column to non-existent file." << endl;
        return false;
    }

    unordered_map<gene_id_t,string> gene_line;

    // First need to read everything in.
    ifstream fin(phene_path.c_str());
    string old_header, old_column_headers;
    getline(fin, old_header);
    getline(fin, old_column_headers);


    // Read in the key value (the gene id)
    gene_id_t gene;
    string line;

    // Count the number of columns for the first line.
    size_t cols = 1;
    
    fin >> gene;
    getline(fin, line);

    for (size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '\t')
            cols++;
    }

    while (fin) {
        gene_line[gene] = line;

        // priming read
        fin >> gene;
        getline(fin, line);
    }

    fin.close();

    // Now re-open the same file and dump it all back out.
    ofstream fout(phene_path.c_str());

    fout << old_header << " (" << cols << "); " << header << endl;
    fout << old_column_headers << '\t' << column_header << endl;
    
    for (unordered_map<gene_id_t,string>::const_iterator i = gene_line.begin(); i != gene_line.end(); ++i) {
        fout << i->first << '\t' << gene_line[i->first] << '\t';

        // Print - if not found in predicted hash
        unordered_map<gene_id_t,double>::iterator pred_it = predicted.find(i->first);
        if (pred_it != predicted.end()) {
            fout << pred_it->second << endl;
            // Get rid of it so it's not printed as a left-over.
            predicted.erase(pred_it);
        } else {
            fout << "-" << endl;
        }
    }

    // Print whatever is left in the map.
    for (unordered_map<gene_id_t,double>::const_iterator i = predicted.begin(); i != predicted.end(); ++i) {
        cerr << "Warning: Additional lines being added to file which were not originally in it: gene =" << i->first << endl;
        fout << i->first;
        // Print tabs to fill
        for (size_t t = 1; t < cols; ++t)
            fout << "\t-";
        fout << '\t' << i->second << endl;
    }


    fout.close();
    return true;
}


// Create a file for a phenotype which stores all of the predictions, one line
// per gene.
bool create_prediction_file(const string& dir, const phene_id_t& phene,
        const string& header,
        const string& column_header,
        const set<gene_id_t>& known,
        const unordered_map<gene_id_t,dist_t>& predicted)
{
    string phene_path = make_phene_path(dir,phene);
    if (exists(phene_path)) return false;

    // Open and create file
    ofstream fout(phene_path.c_str());

    // Print header and column headers
    fout << header << endl;
    fout << "gene\tknown\t" << column_header << endl;
    
    for (unordered_map<gene_id_t,dist_t>::const_iterator i = predicted.begin(); i != predicted.end(); ++i) {
        // Print gene \t known \t prediction
        fout << i->first << '\t';
        if (known.find(i->first) != known.end())    fout << "1\t";
        else                                        fout << "0\t";
        fout << i->second << endl;
    }

    fout.close();
    return true;
}

// Create a file for a phenotype which stores all of the predictions, one line
// per gene. This version doesn't print known values.
bool create_prediction_file(const string& dir, const phene_id_t& phene,
        const string& header,
        const string& column_header,
        const unordered_map<gene_id_t,dist_t>& predicted)
{
    string phene_path = make_phene_path(dir,phene);
    if (exists(phene_path)) return false;

    // Open and create file
    ofstream fout(phene_path.c_str());

    // Print header and column headers
    fout << header << endl;
    fout << "gene\t" << column_header << endl;

    for (unordered_map<gene_id_t,dist_t>::const_iterator i = predicted.begin(); i != predicted.end(); ++i)
        fout << i->first << '\t' << i->second << endl;

    fout.close();
    return true;
}

std::ostream& print_set(const unordered_set<size_t>& Set, std::ostream& out) {

    unordered_set<size_t>::const_iterator i = Set.begin();
    out << *i; ++i;
    
    for (; i != Set.end(); ++i) {
        out << ", " << *i;
    }
    out << flush;
    return out;
}


// Returns true if there are any digits in a text field.
bool has_numeric_component(const string& str) {
    for (size_t i = 0; i < str.size(); ++i) {
        if (str.at(i) >= '0' && str.at(i) <= '9') {
            return true;
        }
    }
    return false;
}

void ignore_set_header(istream& in) {
    // Figure out what line we're on.
    char c = in.get();
    char d = in.get();
    while ((c >= 'A' && c <= 'Z' && d >= 'a' && d <= 'z') || (c >= 'a' && d <= 'z')) {
        // Put them back so we can say what we ignored.
        in.putback(d);
        in.putback(c);

        // Ignore the header
        string header;
        getline(in, header);

        cerr << "Ignored header line '" << header << "'" << endl;

        if (has_numeric_component(header)) {
            cerr << "Error: Header line contains numeric component. This is like a gene or phenotype label, and ignore_header needs to be modified." << endl;
            throw;
        }

        // Get the next two characters.
        c = in.get();
        d = in.get();
    }
    // When the loop fails put them back.
    in.putback(d);
    in.putback(c);

    // Now we're done.
}

set<gene_id_t> load_test_set(const string& filename) {
    return read_set<gene_id_t>(filename, false);
}

