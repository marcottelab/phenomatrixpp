/* 
 * File:   knearest.h
 * Author: jwoods
 *
 * Created on July 23, 2009, 11:41 AM
 */

#ifndef _KNEAREST_H
#define	_KNEAREST_H

#include "oracle.h"

class knearest : public oracle {
    typedef vector<multiset<dist_phene_t> > sorted_t;
public:
    knearest(list<marshall>& mar, marshall& predmar, size_t kval, const set<gene_id_t>& test_set, const multimap<gene_id_t, phene_id_t>& cell_test_set,  const string& identifier, const string& distance_function);
    virtual ~knearest();

    virtual string identify_method() const {
        return method_str;
    }
    virtual string identify_column() const {
        return column_str;
    }

    virtual void predict();
protected:
    virtual bool save_sorted(const vector<phene_id_t>& j_to_phene, const string& filename, size_t k) const;
    virtual bool load_sorted(const vector<phene_id_t>& j_to_phene, const string& filename, size_t k);

    // Calculates sortings on the fly.
    sorted_t sorted;
    size_t   k_;
    // Column headers and method title depend on k in this child class.
    string method_str;
    string column_str;
};



#endif	/* _KNEAREST_H_ */

