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

// Radix includes
#include "AugmentedGenomicRegion.hpp"

// using the following names
using std::string;
using std::vector;
using std::stringstream;


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
      throw SMITHLABException("Lookup of read count for " + sNames[i] +\
                              " for region " + this->get_name() + " failed. "
                              "No read-count for that sample");
    }
    res.push_back(this->sampleCounts.find(sNames[i])->second);
  }
}

