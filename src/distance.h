/* 
 * File:   distance.h
 * Author: jwoods
 *
 * Created on November 16, 2009, 3:45 PM
 */

#ifndef _DISTANCE_H
#define	_DISTANCE_H

#include "euclidean.h"
#include "hypergeometric.h"
#include "utilities.h"

#include <iostream>
#include <string>
#include <map>
using std::map;


double (*switch_distance_function(const std::string& distance_measure))(size_t,size_t,size_t,size_t);

class Distance {
public:
    Distance(double (*fn_ptr)(size_t, size_t, size_t, size_t));
    ~Distance();

    double operator()(size_t m, size_t n, size_t k, size_t N);

    void test() const;
protected:
    double (*function_pointer)(size_t, size_t, size_t, size_t);
};


#endif	/* _DISTANCE_H */

