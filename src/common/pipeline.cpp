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

void
wand(istream &design_encoding, vector<string> &exonReadCount_fns,
     string test_factor_name, ostream &out, const bool VERBOSE) {
  // load the exon/gene counts
  unordered_map<string, Gene > genes;
  if (exonReadCount_fns.size() == 1)
    readExons(exonReadCount_fns[0], genes, VERBOSE);
  else
    readExons(exonReadCount_fns, genes, VERBOSE);

  
  // load the design matrix
  // TODO eliminate the need for the user to provide the 'base' factor
  Design full_design(design_encoding);
  vector<string> factor_names = full_design.factor_names();
  
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
  
  for (unordered_map<string,Gene>::iterator gene_it = genes.begin();
       gene_it != genes.end(); gene_it++) {
    vector<size_t> geneReadCounts;
    gene_it->second.getReadcounts(exonReadCount_fns, geneReadCounts);
    cerr << "gene: " << gene_it->second.getName() << endl;
    for (Gene::const_iterator exon_it = gene_it->second.begin();
         exon_it != gene_it->second.end(); ++exon_it) {
      vector<size_t> exonReadcounts;
      exon_it->getReadcounts(exonReadCount_fns, exonReadcounts);

      full_regression.set_response(geneReadCounts, exonReadcounts);
      gsl_fitter(full_regression);

      null_regression.set_response(geneReadCounts, exonReadcounts);
      gsl_fitter(null_regression);

      // TODO skip tests that will obviously result in p-value of 1
      double pval = loglikratio_test(null_regression.maximum_likelihood(),
                                     full_regression.maximum_likelihood());

      // TODO add user-selected significance threshold
      // TODO collect up all p-values and adjust for multiple hypothesis testing
      if (pval < 0.01) {
        GenomicRegion exonOut(exon_it->getGenomicRegion());
        exonOut.set_score(full_regression.log_fold_change(test_factor));
        out << exonOut << "\t" << pval << endl;
      }
    }
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
 * \brief Load exon and gene read counts from a set of files. Each file
 *        represents a sample, and must be in BED format. Each entry in a file
 *        gives the genomic region associated with an exon, and must be named
 *        exonName_exon_geneName, where exonName and geneName can be any
 *        strings as long as they (1) don't contain the tab character and (2)
 *        uniquely identify the exon in the file. Exons in a gene must not
 *        overlap each other. The score field must give the number of reads
 *        mapping to that exon in the given sample.
 * \param filenames filenames to load sample data from
 * \param genes     results (Gene objects) will be placed into this map, where
 *                  the keys (string) are gene names.
 * \param VERBOSE   print extra status messages if this is true. Defaults to
 *                  false if not set.
 */
void
readExons(const vector<string> &filenames,
          unordered_map< string, Gene > &genes, const bool VERBOSE) {
  const string nameDelim = "_exon_";
  set<string> firstSampleExons;
  for (size_t i = 0; i < filenames.size(); ++i) {
    if (VERBOSE) cerr << "LOADING SAMPLE " << filenames[i] << " ... ";
    set<string> currentSampleExons;
    const string sampleName(filenames[i]);
    vector<GenomicRegion> sampleExons;
    ReadBEDFile(filenames[i], sampleExons);
    if (VERBOSE) cerr << "DONE. LOADED " << sampleExons.size() << " "
                      << "EXONS" << endl;
    for (size_t j = 0; j < sampleExons.size(); ++j) {
      const string fullExonName(sampleExons[j].get_name());
      string geneName(getGeneName(fullExonName));

      // if we haven't seen this gene before, we make a new gene object for it
      unordered_map< string, Gene >::iterator loc = genes.find(geneName);
      if (loc == genes.end()) {
        Exon e(sampleExons[j]);
        genes.insert (make_pair<string, Gene>(geneName, Gene(e)));
        loc = genes.find(geneName);
      }
      loc->second.addExonSampleCount(sampleExons[j], sampleName,
                                     sampleExons[j].get_score());

      // we remember a kind of manual hash of all the exons we see
      stringstream ss;
      ss << sampleExons[j].get_chrom() << sampleExons[j].get_start()
         << sampleExons[j].get_end() << sampleExons[j].get_name()
         << sampleExons[j].get_strand();
      currentSampleExons.insert(ss.str());
    }

    // finished loading all of the exons for this sample. If this is the first
    // sample, just remember what the exon keys were. Otherwise, we're going
    // check that the keys from this sample matched the ones we already saw in
    // the first sample.
    if (i == 0) std::swap(firstSampleExons, currentSampleExons);
    else {
      string m = getExonDifferenceString(firstSampleExons, currentSampleExons);
      string a = getExonDifferenceString(currentSampleExons, firstSampleExons);
      if (!m.empty())
        throw SMITHLABException("The following exons were present in " +\
                                filenames[0] + "but absent in " +\
                                filenames[i] + ": " + m);
      if (!a.empty())
        throw SMITHLABException("The following exons were present in " +\
                                filenames[i] + "but absent in " +\
                                filenames[0] + ": " + a);
    }
  }
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
 * \param filenames filename to load the matrix from.
 * \param genes     results (Gene objects) will be placed into this map, where
 *                  the keys (string) are gene names.
 * \param VERBOSE   print extra status messages if this is true. Defaults to
 *                  false if not set.
 */
void
readExons(const string filename, unordered_map< string, Gene > &genes,
          const bool VERBOSE) {
  vector<string> sampleNames;
  bool first = true;
  ifstream file(filename.c_str());
  if (!file.good()) throw SMITHLABException("Unable to open file :" + filename);
  string line;
  while (!file.eof()) {
    getline(file, line);
    line = smithlab::strip(line);
    if (line.empty()) continue;
    vector<string> parts(smithlab::split(line,"\t"));
    if (first) {
      first = false;
      std::swap(parts, sampleNames);
    } else {
      assert(parts.size() == sampleNames.size() + 5);
      const string joined (parts[0] + "\t" + parts[1] + "\t" +\
                           parts[2] + "\t" + parts[3] + "\t" + parts[4]);
      cerr << joined << endl;
      GenomicRegion r(joined);
      const string fullExonName(r.get_name());
      const string geneName(getGeneName(fullExonName));

      // if we haven't seen this gene before, we make a new gene object for it
      unordered_map< string, Gene >::iterator loc = genes.find(geneName);
      if (loc == genes.end()) {
        Exon e(r);
        genes.insert (make_pair<string, Gene>(geneName, Gene(e)));
        loc = genes.find(geneName);
      }
      for (size_t i = 0; i < sampleNames.size(); ++i) {
        int sampleCount(atoi(parts[i+5].c_str()));
        if (sampleCount < 0) {
          stringstream ss;
          ss << "line: " << line << " contains an invalid sample count: "
             << sampleCount;
          throw SMITHLABException(ss.str());
        }
        loc->second.addExonSampleCount(r, sampleNames[i],
                                       static_cast<size_t>(sampleCount));
      }
    }
  }
}
