#ifndef COMMANDLINEPARSER_HH
#define COMMANDLINEPARSER_HH COMMANDLINEPARSER_HH

#include "mpi/mpi_wrapper.hh"
#include <iostream>
#include <map>
#include <string>

/** @file commandlineparser.hh
 *
 * @brief Class for parsing command line options
 */

/** @class CommandLineParse
 *
 * @brief Simple class for parsing command line arguments
 *
 * Takes command line arguments of the form --OPTION=VALUE
 * and provides functionality for converting VALUE to different
 * formats.
 */

class CommandLineParser {
public:
  /** @brief
   * Construct new instance
   *
   * @param[in] argc Number of command line arguments
   * @param[in] argc Array with command line arguments
   */
  CommandLineParser(const int argc, char *argv[]);

  /** @brief Get value of option as string
   *
   * Returns false if option can not be found
   *
   * @param[in] key Name of option
   * @param[out] value Resulting value
   */
  bool getopt_string(const std::string key, std::string &value) const;

  /** @brief Get value of option as integer
   *
   * Returns false if option can not be found or can not be converted
   *
   * @param[in] key Name of option
   */
  bool getopt_int(const std::string key, int &value) const;

  /** @brief Get value of option as unsigned (long) integer
   *
   * Returns false if option can not be found or can not be converted
   *
   * @param[in] key Name of option
   */
  bool getopt_ulong(const std::string key, unsigned long &value) const;

  /** @brief Get value of option as double
   *
   * Returns false if option can not be found or can not be converted
   *
   * @param[in] key Name of option
   * @param[out] value Resulting value
   */
  bool getopt_double(const std::string key, double &result) const;

private:
  /** @brief Verify that a key is in the dictioary and abort if not
   *
   * @param[in] key Key to find
   */
  bool verify_key(const std::string key) const;

  /** @brief Type of map used for dictionary */
  typedef std::map<std::string, std::string> MapType;
  /** @brief Dictionary with key-value pairs */
  mutable MapType dict;
};

#endif // COMMANDLINEPARSER_HH
