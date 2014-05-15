/*    Copyright (C) 2014 University of Southern California and
 *                       Philip J Uren, Egor Dolzhenko, Andrew D Smith
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

// stl includes
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <tr1/unordered_map>

// Radix includes
#include "AugmentedGenomicRegion.hpp"

// Smithlab common cpp includes
#include "smithlab_os.hpp"

// using the following names
using std::string;
using std::vector;
using std::stringstream;
using std::tr1::unordered_map;
using std::ifstream;
using std::make_pair;


/******************************************************************************
 *        DEFINITION OF METHODS FOR THE AugmentedGenomicRegion CLASS
 ******************************************************************************/

/**
 * \brief Does this region partially overlap the given genomic region? This
 *        function doesn't consider an exact region match (i.e. same chrom,
 *        start, end, and strand) to be a 'partial overlap'.
 * \param r TODO
 * \return TODO
 */
bool
AugmentedGenomicRegion::partialOverlap(const GenomicRegion &r) const {
  if (r.get_chrom() != this->get_chrom()) return false;
  if (r.get_strand() != this->get_strand()) return false;
  // by now they must be one the same chromosome and strand... first see if
  // they match exactly
  if ((r.get_start() == this->get_start()) &&
      (r.get_end() == this->get_end())) return false;
  // not an exact match
  return (this->overlaps(r));
}

/**
 * \brief determine whether this augmented genomic region occupies the same
 *        genomic locus as r. Ignores names and scores of regions.
 * \param r the genomic region to compare to.
 * \return true if the chrom, start, end and strand of this region matches r
 */
bool
AugmentedGenomicRegion::sameRegion(const GenomicRegion &r) const {
  return ((r.get_chrom() == this->get_chrom()) &&
          (r.get_start() == this->get_start()) &&
          (r.get_end() == this->get_end()) &&
          (r.get_strand() == this->get_strand()));
}

/**
 * \brief add a read count for the given sample to this exon.
 * \param sampleName TODO
 * \param count      TODO
 * \throws SMITHLABException if a count for the given sample name already
 *         exists for this exon
 */
void
AugmentedGenomicRegion::addSampleReadCount(const string &sampleName,
                                           const size_t count) {
  if (sampleCounts.find(sampleName) != sampleCounts.end()) {
    stringstream ss;
    ss << "trying to add sample count for genomic region at " << (*this)
       << " in sample " << sampleName << ", but this exon already has a "
       << "count for that sample: " << sampleCounts.find(sampleName)->second;
    throw SMITHLABException(ss.str());
  }
  sampleCounts[sampleName] = count;
}

/**
 * \brief populate the given vector<size_t> with the read counts in this exon
 *        for the given sample names. The vector of resultant read counts will
 *        have the same order as the names that are presented in sNames
 * \param sNames   TODO
 * \param res      TODO
 */
void
AugmentedGenomicRegion::getReadcounts(const vector<string> &sNames,
                                      vector<size_t> &res) const {
  res.clear();
  for (size_t i = 0; i < sNames.size(); ++i) {
    if (this->sampleCounts.find(sNames[i]) == this->sampleCounts.end()) {
      stringstream ss;
      ss << "Lookup of read count for " << sNames[i] << " for region "
         << this->get_name() << " failed. " << "No read-count for that sample. "
         << "This region contains read counts for the following samples: ";
      vector<string> knownSamples;
      this->getSampleNames(knownSamples);
      if (knownSamples.size() == 0) ss << "<NONE>";
      for (size_t j = 0; j < knownSamples.size(); ++j) {
        if (j != 0) ss << ", ";
        ss << knownSamples[j];
      }
      throw SMITHLABException(ss.str());
    }
    res.push_back(this->sampleCounts.find(sNames[i])->second);
  }
}

/**
 * \brief get a vector of the sample names that this region has counts for.
 * \param names the results will be added to this vector. Note that anything
 *              already in here is cleared.
 */
void
AugmentedGenomicRegion::getSampleNames(vector<string> &names) const {
  names.clear();
  typedef unordered_map<string, size_t>::const_iterator mapItType;
  for (mapItType it = this->sampleCounts.begin();
       it != this->sampleCounts.end(); ++it) {
    names.push_back(it->first);
  }
}

/**
 * \brief TODO
 */
size_t
AugmentedGenomicRegion::getSampleCount() const {
  return this->sampleCounts.size();
}


/**
 * \brief             Load a collection of AugmentedGenomicRegions from a set
 *                    of BED format files. Regions that appear in multiple
 *                    files will have sample counts added for each file; they
 *                    won't be duplicated.
 * \param fns         The filenames of the BED format files to load the
 *                    AugmentedGenomicRegions and their counts from.
 * \param res         The loaded AugmentedGenomicRegions are added to this.
 *                    Anything in this vector beforehand is cleared.
 * \param sampleNames The sample names (defined as the names of the BED files,
 *                    minus their extension) will be added to this vector.
 *                    Anything already in here is cleared.
 * \param VERBOSE     If true, status messages for each file loaded (including
 *                    the number of regions loaded from it, and the number that
 *                    are newly observed in that file) will be printed to stderr
 */
void
AugmentedGenomicRegion::readFromBEDFiles(const vector<string> &fns,
                                         vector<AugmentedGenomicRegion> &res,
                                         vector<string> &sampleNames,
                                         const bool VERBOSE) {
  // although we're going to give the caller a vector at the end, we keep
  // an unordered map while we're loading the samples so when we see a new
  // region, we can know whether we already created an AugmentedGenomicRegion
  // for that region or not.
  unordered_map<string, AugmentedGenomicRegion> regions;
  sampleNames.clear();
  for (size_t i = 0; i < fns.size(); ++i) {
    const string sampleName(strip_path_and_suffix(fns[i]));
    sampleNames.push_back(sampleName);
    if (VERBOSE) std::cerr << "LOADING SAMPLE " << sampleName << " ... ";
    vector<GenomicRegion> regionsThisFile;
    ReadBEDFile(fns[i], regionsThisFile);
    size_t newRegions = 0;
    for (size_t j = 0; j < regionsThisFile.size(); ++j) {
      // regions are uniquely identified by their chrom, start, end, name,
      // and strand; obviously, score is not considered.
      stringstream ss;
      ss << regionsThisFile[j].get_chrom() << regionsThisFile[j].get_start()
         << regionsThisFile[j].get_end() << regionsThisFile[j].get_name()
         << regionsThisFile[j].get_strand();
      string key = ss.str();

      unordered_map<string, AugmentedGenomicRegion>::iterator loc;
      loc = regions.find(key);
      if (loc == regions.end()) {
        AugmentedGenomicRegion e(regionsThisFile[j]);
        regions.insert(make_pair<string, AugmentedGenomicRegion>(key, e));
        loc = regions.find(key);
        newRegions += 1;
      }
      loc->second.addSampleReadCount(sampleName,
                                     regionsThisFile[j].get_score());
    }

    if (VERBOSE) std::cerr << "DONE. LOADED " << regionsThisFile.size() << " "
                           << "REGIONS (" << newRegions
                           << " WERE NOT SEEN IN EARLIER FILES)" << std::endl;
  }

  // now package them up for the caller.
  res.clear();
  res.reserve(regions.size());
  typedef unordered_map<string, AugmentedGenomicRegion>::iterator mapIterType;
  for (mapIterType it = regions.begin(); it != regions.end(); ++it)
    res.push_back(it->second);
}

/**
 * \brief             Load a collection of AugmentedGenomicRegions from a
 *                    read-count matrix file. The file is tab separated. The
 *                    first row is a header and lists the sample names for
 *                    each sample included in the file. Each subsequent row
 *                    represents a region. The first five columns in a row are:
 *                    chromosome, start index, end index, exon name, and exon
 *                    strand. Each subsequent column i after the first five
 *                    gives the read count for the region in the ith sample
 *                    (matching the order in the header row).
 * \param fn          [in] filename to load the matrix from.
 * \param res         [out] resultant AugmentedGenomicRegions are placed into
 *                    this. Anything already in here is cleared.
 * \param sampleNames [out] the sample names from the matrix will be placed
 *                    into this vector. Anything that was in here before will
 *                    be cleared first.
 * \param VERBOSE     [in] print extra status messages if this is true.
 *                    Defaults to false if not set.
 */
void
AugmentedGenomicRegion::readFromMatrixFile(const string &fn,
                                           vector<AugmentedGenomicRegion> &res,
                                           vector<string> &sampleNames,
                                           const bool VERBOSE) {
  sampleNames.clear();
  res.clear();
  bool first = true;
  ifstream file(fn.c_str());
  if (!file.good()) throw SMITHLABException("Unable to open file :" + fn);
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
      if (parts.size() != sampleNames.size() + 5) {
        throw SMITHLABException("encountered malformed line when reading " +\
                                fn + ": " + line);
      }
      const string joined (parts[0] + "\t" + parts[1] + "\t" +\
                           parts[2] + "\t" + parts[3] + "\t" + "0" + "\t" +\
                           parts[4]);

      // extra parens below disambiguate variable defn. from func. defn.
      AugmentedGenomicRegion r((GenomicRegion(joined)));
      for (size_t i = 0; i < sampleNames.size(); ++i) {
        int sampleCount(atoi(parts[i+5].c_str()));
        if (sampleCount < 0) {
          stringstream ss;
          ss << "line: " << line << " contains an invalid sample count: "
             << sampleCount;
          throw SMITHLABException(ss.str());
        }
        r.addSampleReadCount(sampleNames[i], static_cast<size_t>(sampleCount));
      }
      res.push_back(r);
    }
  }
}

