/* 
 * File:   constants.h
 * Author: jwoods
 *
 * Created on July 15, 2009, 2:20 PM
 */

#ifndef _CONSTANTS_H
#define	_CONSTANTS_H

#include <list>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <limits>

#include <boost/progress.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

using std::map;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::ostringstream;
using std::ifstream;
using std::istream;
using std::ofstream;
using std::getline;
using std::vector;
using std::pair;
using std::set;
using std::list;
using std::multiset;
using std::multimap;
using boost::progress_display;
using boost::unordered_map;
using boost::unordered_set;
using boost::tokenizer;
using boost::lexical_cast;
using boost::numeric::ublas::compressed_matrix;
using boost::numeric::ublas::symmetric_matrix;
using boost::numeric::ublas::matrix;
using boost::filesystem::exists;
using boost::filesystem::path;
using boost::shared_ptr;
using boost::filesystem::is_directory;

#include "type_shield.h"

typedef unsigned int uint;
typedef double dist_t;
typedef unordered_map<uint,size_t>::iterator hash_iterator;
typedef string gene_id_t;
typedef string species_id_t;
typedef uint phene_id_t;
typedef type_shield<dist_t, phene_id_t> dist_phene_t;
typedef unordered_map<phene_id_t, species_id_t> phene_sp_map_t;
typedef std::numeric_limits<double> dbprec1;
typedef multiset<dist_phene_t> dist_phene_multiset_t;

const gene_id_t GENE_NONE = "";

#endif	/* _CONSTANTS_H */

