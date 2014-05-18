/*    Copyright (C) 2014 University of Southern California and
 *                       Philip J. Uren, Egor Dolzhenko, Andrew D Smith
 *
 *    Authors: Philip J. Uren, Andrew D. Smith and Egor Dolzhenko
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

// STD headers. 
#include <set>
#include <iostream>
#include <algorithm>
#include <tr1/unordered_map>
#include <sstream>

// Smithlab headers.
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

// Local headers.
#include "Gene.hpp"
#include "Exon.hpp"
#include "design.hpp"
#include "pipeline.hpp"
#include "table_row.hpp"
#include "loglikratio_test.hpp"
#include "gsl_fitter.hpp"
#include "regression.hpp"

// bring the following into the default namespace
using std::istream;
using std::ifstream;
using std::ostream;
using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::make_pair;
using std::stringstream;
using std::istringstream;
using std::tr1::unordered_map;
using std::set;

/**
 * \brief TODO
 */
void
checkDesign (const Design &d, const vector<string> &sampleNames) {
  vector<string> designSampleNames(d.sample_names());

  // For output later.
  string designSamples_str = "", sampleNames_str = "";
  for (size_t i = 0; i < designSampleNames.size(); ++i) {
    if (i != 0) designSamples_str += ", ";
      designSamples_str += designSampleNames[i];
  }
  for (size_t i = 0; i < sampleNames.size(); ++i) {
    if (i != 0) sampleNames_str += ", ";
      sampleNames_str += sampleNames[i];
  }

  if (designSampleNames != sampleNames) {
    throw SMITHLABException("study design doesn't match input files. "
                            "Found these sample in input: " +\
                            sampleNames_str + " and these samples in design "
                            "matrix: " + designSamples_str);
  }
}

/**
 * \brief TODO
 * \param design_encoding     TODO
 * \param geneReadCount_fn    TODO
 * \param intronReadCount_fn  TODO
 * \param test_factor_name    TODO
 * \param out                 TODO
 * \param pThresh             output only those exons with a p-value less than
 *                            or equal to this value. Default is 0.01.
 * \param VERBOSE             TODO
 */
void
findDiffIncIntrons(istream &design_encoding, const string &geneReadCount_fn,
                   const string &intronReadCount_fn,
                   const string &test_factor_name,
                   ostream &out, const double pThresh, const bool VERBOSE) {

}

/**
 * \brief If exonReadCount_fns contains only a single entry, we assume it is
 *        in matrix format (see details below), otherwise we assume the input
 *        files are in BED format.
 *        TODO describe matrix format here too.
 * \param design_encoding   TODO
 * \param exonReadCount_fns TODO
 * \param test_factor_name  TODO
 * \param out               TODO
 * \param pThresh           output only those exons with a p-value less than
 *                          or equal to this value. Default is 0.01.
 * \param VERBOSE           TODO
 */
void
run(istream &design_encoding, vector<string> &exonReadCount_fns,
     string test_factor_name, ostream &out, const double pThresh,
     const bool VERBOSE) {
  // load the exon/gene counts and sample names
  unordered_map<string, Gene > genes;
  vector<string> sampleNames;
  readExons(exonReadCount_fns, genes, sampleNames, VERBOSE);
  
  // load the design matrix
  // TODO eliminate the need for the user to provide the 'base' factor
  Design full_design(design_encoding);
  vector<string> factor_names = full_design.factor_names();
  
  // make sure the design matrix sample names match the names of the read
  // count files.
  checkDesign(full_design, sampleNames);

  vector<string>::const_iterator test_factor_it = 
      std::find(factor_names.begin(), factor_names.end(), test_factor_name);
  
  // Checking that the provided test factor names exist.
  if (test_factor_it == factor_names.end())
    throw SMITHLABException(test_factor_name + " is not a part of the design" 
                            " specification.");
  
  // Factors are identified with their indexes to simplify naming.
  size_t test_factor = test_factor_it - factor_names.begin();
  
  Regression full_regression(full_design);
  Design null_design = full_design;
  null_design.remove_factor(test_factor);
  Regression null_regression(null_design);
  
  vector< std::pair<GenomicRegion, double> > results;
  for (unordered_map<string,Gene>::iterator gene_it = genes.begin();
       gene_it != genes.end(); gene_it++) {
    vector<size_t> geneReadCounts;
    gene_it->second.getReadcounts(sampleNames, geneReadCounts);
    cerr << "gene: " << gene_it->second.getName() << endl;
    for (Gene::const_iterator exon_it = gene_it->second.begin();
         exon_it != gene_it->second.end(); ++exon_it) {
      vector<size_t> exonReadCounts;
      exon_it->getReadcounts(sampleNames, exonReadCounts);

      full_regression.set_response(geneReadCounts, exonReadCounts);
      gsl_fitter(full_regression);

      null_regression.set_response(geneReadCounts, exonReadCounts);
      gsl_fitter(null_regression);

      // TODO skip tests that will obviously result in p-value of 1?
      double pval = loglikratio_test(null_regression.maximum_likelihood(),
                                     full_regression.maximum_likelihood());
      if (!std::isfinite(pval)) {
        std::cerr << "WARNING: log-likelihood ratio test failed for exon " 
                  << exon_it->get_name()
                  << "; setting p-value to 1.0" << endl;
        pval = 1.0;
      }

      GenomicRegion exonOut(*exon_it);
      exonOut.set_score(full_regression.log_fold_change(test_factor));
      results.push_back(std::make_pair<GenomicRegion, double>(exonOut, pval));
    }
  }

  std::sort(results.begin(), results.end());
  // TODO adjust p-values for multiple hypothesis testing
  for (size_t i = 0; i < results.size(); ++i) {
    if (results[i].second <= pThresh)
      out << results[i].first << "\t" << results[i].second << endl;
  }
}

/**
 * \brief Get a string representation of the exon names that are in the
 *        first std::set, but not in the second.
 */
string
getExonDifferenceString(const set<string> &exons1, const set<string> &exons2) {
  std::set<string> diff;
  std::set_difference(exons1.begin(), exons1.end(), exons2.begin(),
                      exons2.end(), std::inserter(diff, diff.end()));
  string dStr = "";
  for (std::set<string>::iterator it = diff.begin(); it != diff.end(); ++it) {
    if (it != diff.begin()) dStr += ",";
    dStr = dStr + (*it);
  }
  return dStr;
}

/**
 * \brief Load exon and gene read counts from a single read-count matrix file.
 *        The file is tab separated. The first row is a header and lists the
 *        sample names for each sample included in the file. Each subsequent
 *        row represents an exon. The first five columns in an exon row are:
 *        chromosome, start index, end index, exon name, and exon exon strand.
 *        The exon name must follow the format: exonName_exon_geneName,
 *        where exonName and geneName can be any strings as long as they
 *        (1) don't contain the tab character and (2) uniquely identify the
 *        exon in the file. Exons in a gene must not overlap each other. Each
 *        subsequent column i after the first five gives the read count for the
 *        exon in the ith sample (matching the order in the header row).
 * \param filenames    [in] filename to load the matrix from.
 * \param genes        [out] results (Gene objects) will be placed into this
 *                     map, where the keys (string) are gene names.
 * \param sampleNames  [out] the sample names from the matrix will be placed
 *                     into this vector. Anything that was in here before will
 *                     be cleared first.
 * \param VERBOSE      [in] print extra status messages if this is true.
 *                     Defaults to false if not set.
 */
void
readExons(const vector<string> filenames, unordered_map< string, Gene > &genes,
          vector<string> &sampleNames, const bool VERBOSE) {
  // note that the below two lines also populate the sampleNames vector
  vector<AugmentedGenomicRegion> regions;
  if (filenames.size() == 1)
    AugmentedGenomicRegion::readFromMatrixFile(filenames[0], regions,
                                               sampleNames, VERBOSE);
  else
    AugmentedGenomicRegion::readFromBEDFiles(filenames, regions,
                                             sampleNames, VERBOSE);

  for (size_t i = 0; i < regions.size(); ++i) {
    Exon e(regions[i]);
    // if we haven't seen this gene before, we make a new gene object for it
    unordered_map< string, Gene >::iterator loc = genes.find(e.getGeneName());
    if (loc == genes.end()) {
      genes.insert (make_pair<string, Gene>(e.getGeneName(), Gene(e)));
      loc = genes.find(e.getGeneName());
    }
    vector<size_t> readCounts;
    e.getReadcounts(sampleNames, readCounts);
    assert(readCounts.size() == sampleNames.size());
    for (size_t j = 0; j < readCounts.size(); j++)
      loc->second.addExonSampleCount(e, sampleNames[j], readCounts[j]);
  }
}
