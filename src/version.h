#ifndef _VERSION_H_
# define _VERSION_H_

# include <string>

const std::string VERSION  = "@VERSION@";
const std::string PLATFORM = "@ARCH@";

// Prints program version information.
void print_program_header() {
    cout << "Phenomatrix++ v" << VERSION << " on " << PLATFORM << endl;
    cout << "Copyright John O. Woods, 2009 - 2010" << endl;
    cout << "The Marcotte Lab, The University of Texas at Austin" << endl;
    cout << endl;
}

#endif

