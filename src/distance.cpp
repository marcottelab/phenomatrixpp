#include "distance.h"

Distance::Distance(double (*fn_ptr)(size_t m, size_t n, size_t k, size_t N))
 : function_pointer(fn_ptr)
{
}

Distance::~Distance() {}

double Distance::operator()(size_t m, size_t n, size_t k, size_t N) {
    return (*function_pointer)(m, n, k, N);
}


// Returns a function pointer to a distance function based on a request made via
// a string.
double (*switch_distance_function(const std::string& distance_measure))(size_t,size_t,size_t,size_t) {
    map<std::string, double(*)(size_t,size_t,size_t,size_t)> choices;
    choices["hypergeometric"] = &hypergeometric;
    choices["euclidean"]      = &euclidean;
    choices["manhattan"]      = &manhattan;

    return choices[distance_measure];
}


void Distance::test() const {
    using std::cout;
    using std::flush;
    using std::endl;

    cout << "Testing distance function:" << endl;
    cout << " 1. Testing for 0 result..." << flush;
    if (are_equal<double>( (10, 10, 10, 100), 0.0 )) cout << "OK" << endl;
    else { cout << "FAIL" << endl; throw; }

    cout << " 1. Testing for 1 result..." << flush;
    if (are_equal<double>( (10, 10, 0, 100), 1.0 )) cout << "OK" << endl;
    else { cout << "FAIL" << endl; throw; }
}