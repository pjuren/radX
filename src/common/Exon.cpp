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

// stl includes
#include <string>
#include <vector>
#include <sstream>

// Radix includes
#include "Exon.hpp"

// smithlab common includes
#include "GenomicRegion.hpp"

// using the following names
using std::string;
using std::vector;
using std::stringstream;


/******************************************************************************
 *                          STATIC HELPER FUNCTIONS
 ******************************************************************************/

/**
 * \brief Count the number of times a string occurs in another string
 * \param needle      TODO
 * \param haystack    TODO
 * \param overlapping TODO
 * \TODO  this belongs in some library somewhere..
 */
size_t
countOccurrences(const string &needle, const string &haystack,
                 const bool overlapping) {
  int n = 0;
  size_t pos = 0;
  while ((pos = haystack.find(needle, pos)) != std::string::npos ) {
    n++;
    if (overlapping) pos += 1;
    else pos += needle.size();
  }
  return n;
}

/**
 * \brief split a full exon name of the form someGeneName_exon_someExonName
 *        into its pair of <someGeneName, someExonName>
 * \param fullName TODO
 * \throws SMITHLABException if the delimiter '_exon_' appears more than once
 *         or not at all in the full exon name.
 */
std::pair<string, string>
splitFullExonName(const string &fullName) {
  const string delim = "_exon_";
  const size_t numDelimOcc = countOccurrences(delim, fullName);
  if (numDelimOcc == 0) {
    stringstream ss;
    ss << "incorrectly formatted full exon name: " << fullName << "; "
       << "could not find the expected delimiter " << delim << " to separate "
       << "the gene name from the exon name";
    throw SMITHLABException(ss.str());
  }
  if (numDelimOcc > 1) {
    stringstream ss;
    ss << "incorrectly formatted full exon name: " << fullName << "; "
       << "found multiple occurrences of the delimiter " << delim << ", so "
       << "could not unambiguously separate the gene name from the exon name";
    throw SMITHLABException(ss.str());
  }
  size_t s_delim = fullName.find(delim);
  size_t e_delim = s_delim + delim.size();
  return std::make_pair<string,string>(fullName.substr(0, s_delim),
                                       fullName.substr(e_delim, string::npos));

}

/**
 * \brief Extract the gene name from a full exon name. The full exon name is
 *        expected to have the format someGeneName_exon_someExonName.
 * \param exonName TODO
 */
string
getGeneName(const std::string &exonName) {
  std::pair<string, string> splitName(splitFullExonName(exonName));
  return splitName.first;
}


/******************************************************************************
 *                 DEFINITION OF METHODS FOR THE EXON CLASS
 ******************************************************************************/

/**
 * \brief Does this exon partially overlap the given genomic region? This
 *        function doesn't consider an exact region match (i.e. same chrom,
 *        start, end, and strand) to be a 'partial overlap'.
 */
bool
Exon::partialOverlap(const GenomicRegion &r) const {
  if (r.get_chrom() != this->region.get_chrom()) return false;
  if (r.get_strand() != this->region.get_strand()) return false;
  // by now they must be one the same chromosome and strand... first see if
  // they match exactly
  if ((r.get_start() == this->region.get_start()) &&
      (r.get_end() == this->region.get_end())) return false;
  // not an exact match
  return (r.overlaps(this->region));
}


/**
 * \brief TODO
 */
bool
Exon::sameRegion(const GenomicRegion &r) const {
  return ((r.get_chrom() == this->region.get_chrom()) &&
          (r.get_start() == this->region.get_start()) &&
          (r.get_end() == this->region.get_end()) &&
          (r.get_strand() == this->region.get_strand()));
}

/**
 * \brief add a read count for the given sample to this exon.
 * \throws SMITHLABException if a count for the given sample name already
 *         exists for this exon
 */
void
Exon::addSampleReadCount(const string &sampleName, const size_t count) {
  if (sampleCounts.find(sampleName) != sampleCounts.end()) {
    stringstream ss;
    ss << "trying to add sample count for exon at " << this->region
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
Exon::getReadcounts(const vector<string> &sNames, vector<size_t> &res) const {
  res.clear();
  for (size_t i = 0; i < sNames.size(); ++i) {
    if (this->sampleCounts.find(sNames[i]) == this->sampleCounts.end()) {
      throw SMITHLABException("Lookup of read count for " + sNames[i] +\
                              " for exon " + this->getExonName() + " failed. "
                              "No read-count for that sample");
    }
    res.push_back(this->sampleCounts.find(sNames[i])->second);
  }
}

/**
 * \brief TODO
 */
GenomicRegion
Exon::getGenomicRegion() const { return this->region; }

/**
 * \brief get the name of this exon. The name of the exon is defined to be
 *        whatever the region name was when the Exon object was created.
 */
string
Exon::getExonName() const { return this->region.get_name(); }

/**
 * \brief TODO
 */
string
Exon::getChrom() const { return this->region.get_chrom(); }

/**
 * brief TODO
 */
char
Exon::getStrand() const { return this->region.get_strand(); }
