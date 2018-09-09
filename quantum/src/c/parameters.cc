/** @file parameters.cc
 * @brief Implementation of parameters.hh
 */
#include "parameters.hh"

/* Constructor for Parameters class */
Parameters::Parameters(const std::string section) :
  section_(section) {
  regcomp(&regexKeyword_,
          "^ *([a-zA-Z_]+): *(#.*)?$",REG_EXTENDED);
  regcomp(&regexComment_,
          "^ *(#.*)?$",REG_EXTENDED);
  std::string anyValue = "[+-]?[0-9]+|([+-]?[0-9]*\\.?[0-9]*([Ee][+-]?[0-9]+)?|(true|True|TRUE|false|False|FALSE)|[\'\"].*[\'\"])";
  std::string regexAny = "^ *([a-zA-Z0-9_]+) *= *(" + anyValue + "( *, *("+anyValue+"))*) *(#.*)?$";
  regcomp(&regexKeyValueAny_,regexAny.c_str(), REG_EXTENDED);
  std::string intValue = "[+-]?[0-9]+";
  std::string regexInt = "^ *([a-zA-Z0-9_]+) *= *(" + intValue + "( *, *(" + intValue+"))*) *(#.*)?$";
  regcomp(&regexKeyValueInt_,regexInt.c_str(), REG_EXTENDED);

  std::string doubleValue = "[+-]?[0-9]*\\.?[0-9]*([Ee][+-]?[0-9]+)?";
  std::string regexDouble = "^ *([a-zA-Z0-9_]+) *= *(" + doubleValue + "( *, *(" + doubleValue+"))*) *(#.*)?$";
  regcomp(&regexKeyValueDouble_, regexDouble.c_str(), REG_EXTENDED);
  std::string boolValue = "true|True|TRUE|false|False|FALSE";
  std::string regexBool = "^ *([a-zA-Z0-9_]+) *= *(" + boolValue + "( *, *(" + boolValue+"))*) *(#.*)?$";
  regcomp(&regexKeyValueBool_, regexBool.c_str(), REG_EXTENDED);
  std::string stringValue = "[\'\"].*[\'\"]";
  std::string regexString = "^ *([a-zA-Z0-9_]+) *= *(" + stringValue + "( *, *(" + stringValue+"))*) *(#.*)?$";
  regcomp(&regexKeyValueString_, regexString.c_str(), REG_EXTENDED);
}

/* Split a string by delimiter ',' and remove all whitespace */
std::vector<std::string> Parameters::split_string(const std::string s) {
  const char delim=',';
  std::string cur_str ="";
  std::vector<std::string> splitted_string;
  for (auto it=s.begin();it!=s.end();++it) {
    if ((*it != ' ') and (*it !=delim))
      cur_str += *it;
    if (*it == delim) {
      splitted_string.push_back(cur_str);
      cur_str = "";
    }
  }
  splitted_string.push_back(cur_str);
  return splitted_string;
}

/* Read parameters from file */
int Parameters::readFile(const std::string filename) {
  // reset stream
  std::ifstream is;
  is.open(filename.c_str());
  // read from stream until keyword is found
  bool lineValid;
  char line[256];
  bool parseSection = false;
  bool foundSection = false;
  // initialise found list to zero
  std::map<std::string,bool> found;
  typedef std::map<std::string,Datatype>::iterator Iterator;
  for (Iterator it = keywords.begin(); it!= keywords.end();++it) {
    found.insert(std::pair<std::string,bool>((*it).first,false));
  }
  while (is) {
    lineValid=false;
    is.getline(line,256);
    std::string sline(line);
    // Comment
    if (!regexec(&regexComment_,line,0,0,0)) {
      if (verbose > 0)
        std::cout << " >>> COMMENT: " << line << std::endl;
      lineValid = true;
    }
    // Section Keyword
    regmatch_t matchKeyword[3];
    if (!regexec(&regexKeyword_,line,3,matchKeyword,0)) {
      std::string section = sline.substr(matchKeyword[1].rm_so,
                                         matchKeyword[1].rm_eo);
      if (verbose > 0)
        std::cout << " >>> SECTION: \'" << section << "\'"<< std::endl;
      if (section == section_) {
        // Start parsing section
        parseSection = true;
        foundSection = true;
      } else {
        // Have reached a different section, stop parsing
        parseSection = false;
      }
      lineValid = true;
    }
    // Found a key-value pair
    regmatch_t matchKeyValueAny[4];
    if (!regexec(&regexKeyValueAny_,line,4,matchKeyValueAny,0) ) {
      // Check if this parameter is in the list of keywords
      std::string key = sline.substr(matchKeyValueAny[1].rm_so,
                                     matchKeyValueAny[1].rm_eo-
                                     matchKeyValueAny[1].rm_so);
      std::string value = sline.substr(matchKeyValueAny[2].rm_so,
                                       matchKeyValueAny[2].rm_eo-
                                       matchKeyValueAny[2].rm_so);
      std::vector<std::string> split_value = split_string(value);
      if (verbose > 0)
        std::cout << " >>> KEY: \'" << key << "\' VALUE: \'"
                  << value << "\'" << std::endl;
      if (parseSection) {
        if (keywords.count(key) > 0) {
          if (!found[key]) {
            // Integer
            if ( (!regexec(&regexKeyValueInt_,line,0,0,0)) &&
                 (keywords[key] == Integer) ) {
              for (auto it=split_value.begin();it!=split_value.end();++it) {
                contents[key].push_back(new IntParameter("0"));
                contents[key].back()->setValue(*it);
              }
              // Double
            } else if ( (!regexec(&regexKeyValueDouble_,line,0,0,0)) &&
                        (keywords[key] == Double) ) {
              for (auto it=split_value.begin();it!=split_value.end();++it) {
                contents[key].push_back(new DoubleParameter("0.0"));
                contents[key].back()->setValue(*it);
              }
              // Bool
            } else if ( (!regexec(&regexKeyValueBool_,line,0,0,0)) &&
                        (keywords[key] == Bool) ) {
              for (auto it=split_value.begin();it!=split_value.end();++it) {
                contents[key].push_back(new BoolParameter("false"));
                contents[key].back()->setValue(*it);
              }
              // String
            } else if ( (!regexec(&regexKeyValueString_,line,0,0,0)) &&
                        (keywords[key] == String) ) {
              for (auto it=split_value.begin();it!=split_value.end();++it) {
                contents[key].push_back(new StringParameter("blank"));
                contents[key].back()->setValue(*it);
              }
            } else {
              // Parameter has wrong type
              std::cerr << " ERROR reading parameters: Parameter: \'"
                        << key << " = " << value
                        << "\' in section \'" << section_ << "\' "
                        << " has wrong type." << std::endl;
              return 1;
            }
            found[key] = true;
          } else {
            std::cerr << " ERROR reading parameters: Duplicate parameter: \'"
                      << key << " = " << value
                      << "\' in section \'" << section_ << "\' " << std::endl;
            return 1;
          }
          // Unexpected key-value pair in section
        } else {
          std::cerr << " ERROR reading parameters: Invalid parameter: \'"
                    << key << " = " << value
                    << "\' in section \'" << section_ << "\'" << std::endl;
          return 1;
        }
      }
      lineValid = true;
    }
    if (!lineValid) {
      std::cerr << " ERROR reading parameters: Can not parse expression: \'"
                << line
                << "\' in section \'" << section_ << "\'" << std::endl;
      return 1;
    }
  }
  // ERROR handling
  // Check whether we found the section
  if (!foundSection) {
    std::cerr << " ERROR reading parameters: Can not find section: \'"
              << section_ << "\'" << std::endl;
    return 1;
  }
  // Check whether we have found all key-value pair of this section
  typedef std::map<std::string,bool>::iterator FIterator;
  for (FIterator it = found.begin(); it!= found.end();++it) {
    if (!(*it).second) {
      std::cerr << " ERROR reading parameters: parameter \'"
                << (*it).first << "\' in section \'" << section_
                << "\' has not been specified." << std::endl;
      return 1;
    }
  }
  is.close();
  return 0;
}

/* Register a key - value pair */
void Parameters::addKey(const std::string key, const Datatype datatype) {
  keywords.insert(std::pair<std::string,Datatype>(key,datatype));
  std::vector<Parameter*> empty;
  contents[key] = empty;
}

/* Write contents of parameters object to stream */
std::ostream& operator<<(std::ostream& output, const Parameters& p) {
  output << " -- " << p.section_ << " parameters --" << std::endl;
  std::map<std::string,std::vector<Parameter*>>::const_iterator it;
  std::map<std::string,Parameters::Datatype>::const_iterator pit;
  for (it=p.contents.begin(), pit=p.keywords.begin();
       it != p.contents.end(); ++it,++pit) {
    output << "    " << it->first << " = ";
    for (auto sit = (it->second).begin(); sit!=(it->second).end(); ++sit) {
      if (pit->second == Parameters::Integer) {
        output << (*sit)->getInt() << " ";
    } else if (pit->second == Parameters::Double) {
        output << std::scientific << (*sit)->getDouble();
    } else if (pit->second == Parameters::Bool) {
        if ((*sit)->getBool()) {
          output << "true";
        } else {
          output << "false";
        }
      } else if (pit->second == Parameters::String) {
        output << (*sit)->getString();
      }
      output << " ";
    }
    output << std::endl;
  }
  return output;
}
