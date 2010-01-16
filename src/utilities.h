/* 
 * File:   utilities.h
 * Author: jwoods
 *
 * Created on July 16, 2009, 2:44 PM
 */

#ifndef _UTILITIES_H
#define	_UTILITIES_H

#include <map>
#include <set>
#include <string>
#include <limits>
#include <cmath>
using std::set;
using std::multimap;

#include "constants.h"
#include "distance.h"

size_t lines_in_file(const string&);
set<phene_id_t> species_phene_set(const species_id_t&);
pair<size_t,size_t> num_genes_and_phenes_in_file(const string&);
void read_gene_phene_line(ifstream&, gene_id_t&, phene_id_t&);
size_t count_tabs(const string& line);
string make_phene_path(const string& dir, const phene_id_t& phene);
//size_t read_prediction_file(const path& p, string& header_line1, string& header_line2, map<gene_id_t,string>& gene_line);
bool create_prediction_file(const string& dir, const phene_id_t& phene,
        const string& header,
        const string& column_header,
        const set<gene_id_t>& known,
        const unordered_map<gene_id_t,dist_t>& predicted);
bool create_prediction_file(const string& dir, const phene_id_t& phene,
        const string& header,
        const string& column_header,
        const unordered_map<gene_id_t,dist_t>& predicted);
bool add_column_to_prediction_file(const string& dir, const phene_id_t& phene,
        const string& header,
        const string& column_header,
        unordered_map<gene_id_t,double> predicted);

// Pipe a set through an output stream.
std::ostream& print_set(const unordered_set<size_t>& Set, std::ostream& out);

bool has_numeric_component(const string& str);

void ignore_set_header(istream& in);


// Find some sort of item in the counter and increment it, or add it to the
// counter.
template <typename T>
void increment_map_count(map<T,size_t>& counter, const T& item) {
    if (counter.find(item) == counter.end())
        counter[item] = 0;
    else
        counter[item]++;
}


template <typename T>
bool are_equal(T a, T b) {
    return std::abs<T>(a - b) < std::numeric_limits<T>::epsilon();
}


template <class S, class T>
multimap<S,T> read_associations(std::istream& in) {
    multimap<S,T> ret;

    S key;
    T val;
    in >> key;
    in >> val;
    while (in) {
        ret.insert(make_pair(key,val));
        //cerr << "Inserted " << key << "\t" << val << endl;
        //ret[key] = val;

        // Read the next line
        in >> key;
        in >> val;
    }
    return ret;
}

template <class S, class T>
multimap<S,T> load_associations(const string& filename) {
    ifstream f(filename.c_str());
    return read_associations<S,T>(f);
}

template <typename T>
void save_vector(const string& filename, const vector<T>& v, bool with_indeces = false) {
    ofstream out(filename.c_str());

    for (size_t idx = 0; idx < v.size(); ++idx) {
        if (with_indeces)
            out << idx << '\t';

        out << v[idx] << endl;
    }
    out.close();
}

template <typename T>
void save_vector_if_absent(const string& filename, const vector<T>& v, bool with_indeces = false) {
    if (!exists(filename)) save_vector(filename, v, with_indeces);
}
/*
template <class S, class T>
set<T> unique_values(const multimap<S,T>& test) {
    set<T> ret;
    for (multimap<S,T>::const_iterator it = test.begin(); it != test.end(); ++it) {
        ret.insert(it->second);
    }
    return ret;
}


template <typename S, typename T>
set<S> unique_keys(const multimap<S,T>& test) {
    set<S> ret;
    for (multimap<S,T>::const_iterator i = test.begin(); i != test.end(); ++i)
        ret.insert(i->first);
    return ret;
} */


// Uses find to determine if a value is in the set.
template <typename S, typename T>
bool in(const S& Set, const T& val) {
    return (Set.find(val) != Set.end());
}

// Read a set of items from a file.
template <typename T>
set<T> read_set(const string& filename, bool header = true) {
    set<T> ret;

    ifstream fin(filename.c_str());

    // ignore header
    string header_str;
    if (header)
        ignore_set_header(fin);

    T item;
    fin >> item;
    while (fin) {
        ret.insert(item);

        fin >> item;
    }

    fin.close();

    return ret;
}

// Same as above, but uses an open filestream (and keeps it open).
template <typename T>
set<T> read_set(istream& fin) {
    set<T> ret;

    T val;
    fin >> val;
    while (fin) {
        ret.insert(val);
        fin >> val;
    }
    
    return ret;
}

// Calls read_set.
set<gene_id_t> load_test_set(const string& filename);


#endif	/* _UTILITIES_H */

