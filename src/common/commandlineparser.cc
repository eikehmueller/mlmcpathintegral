#include "commandlineparser.hh"
/** @file commandlineparser.cc
 * @brief Implementation of commandlineparser.hh
 */

/* Construct new instance */
CommandLineParser::CommandLineParser(const int argc, char* argv[]) {
    for (int j=1;j<argc;++j) {
        bool malformatted = false;
        std::string sraw(argv[j]);
        
        malformatted = malformatted or (sraw.substr(0,2) != "--");
        int split_pos = sraw.find("=");
        malformatted = malformatted or (split_pos == std::string::npos);
        // Split command line parameter at the '=' sign
        std::string key = sraw.substr(2,split_pos-2);
        std::string value = sraw.substr(split_pos+1,sraw.size()-split_pos-1);
        if  (malformatted) {
            mpi_parallel::cout << "ERROR: malformatted command line option \'" << sraw << "\'" << std::endl;
            exit(-1);
        }
        dict[key] = value;
    }
    // Print out dictionary for debugging purposes
    bool debug = false;
    if (debug) {
        for (MapType::iterator it=dict.begin(); it!=dict.end(); ++it)
            std::cout << "\'" << it->first << "\' : \'" << it->second << "\'" << std::endl;
    }
}

/* Return option as string */
bool CommandLineParser::getopt_string(const std::string key, std::string& value) const {
    if (dict.count(key)==0) {
        return false;
    }
    value = dict[key];
    return true;
}

/* Return option as integer */
bool CommandLineParser::getopt_int(const std::string key, int& value) const {
    if (dict.count(key)==0) {
        return false;
    }
    try {
        value = std::stoi(dict[key]);
    } catch (std::exception& e) {
        mpi_parallel::cout << "WARNING: could not convert \'" << value << "\' to integer:" << std::endl;
        mpi_parallel::cout << e.what() << std::endl;
        return false;
    }
    return true;
}

/* Return option as unsigned long integer */
bool CommandLineParser::getopt_ulong(const std::string key, unsigned long& value) const {
    if (dict.count(key)==0) {
        return false;
    }
    try {
        value = std::stoul(dict[key]);
    } catch (std::exception& e) {
        mpi_parallel::cout << "WARNING: could not convert \'" << value << "\' to integer:" << std::endl;
        mpi_parallel::cout << e.what() << std::endl;
        return false;
    }
    return true;
}



/* Return option as double */
bool CommandLineParser::getopt_double(const std::string key, double& value) const {
    if (dict.count(key)==0) {
        return false;
    }
    try {
        value = std::stod(dict[key]);
    } catch (std::exception& e) {
        mpi_parallel::cout << "WARNING: could not convert \'" << value << "\' to double:" << std::endl;
        mpi_parallel::cout << e.what() << std::endl;
        return false;
    }
    return true;
}
