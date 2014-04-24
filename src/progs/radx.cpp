/*    Copyright (C) 2014 University of Southern California and
 *                       Philip J Uren
 *                       Egor Dolzhenko
 *                       Andrew D Smith
 *
 *    Authors: Philip J Uren, Andrew D. Smith and Egor Dolzhenko
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>

// smithlab headers
#include "OptionParser.hpp"
#include "smithlab_os.hpp"

// headers for local code
#include "pipeline.hpp"

// using block
using std::string;
using std::vector;
using std::istringstream;
using std::cerr;
using std::endl;

int
main(int argc, const char **argv) {

try {
  string outfile, design_filename;
  string test_factor_name;
  bool VERBOSE = false;
  
  /****************** COMMAND LINE OPTIONS ********************/
  OptionParser opt_parse(strip_path(argv[0]), "Egor's program",
                           "<design-matrix> <data-matrix>");
  opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                    false, outfile);
  opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
  opt_parse.add_opt("design", 'd', "experiment design", true, design_filename);
  opt_parse.add_opt("factor", 'f', "a factor to test",
                    true, test_factor_name);
  vector<string> leftover_args;
  opt_parse.parse(argc, argv, leftover_args);
  if (argc == 1 || opt_parse.help_requested()) {
    cerr << opt_parse.help_message() << endl
         << opt_parse.about_message() << endl;
    return EXIT_SUCCESS;
  }
  if (opt_parse.about_requested()) {
    cerr << opt_parse.about_message() << endl;
    return EXIT_SUCCESS;
  }
  if (opt_parse.option_missing()) {
    cerr << opt_parse.option_missing_message() << endl;
    return EXIT_SUCCESS;
  }

  vector<string> input_fns(leftover_args);

  /****************** END COMMAND LINE OPTIONS *****************/
  
  std::ifstream design_file(design_filename.c_str());
  if (!design_file)
    throw SMITHLABException("could not open file: " + design_filename);
    
  //std::ifstream table_file(table_filename.c_str());
  //if (!table_file)
  //  throw SMITHLABException("could not open file: " + table_filename);
  
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  wand(design_file, input_fns, test_factor_name, out);
}
catch (const SMITHLABException &e) {
  cerr << e.what() << endl;
  return EXIT_FAILURE;
}
catch (std::bad_alloc &ba) {
  cerr << "ERROR: could not allocate memory" << endl;
  return EXIT_FAILURE;
}

  return EXIT_SUCCESS;
}
