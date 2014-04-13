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

/**
 * \brief TODO
 */
std::string
getGeneName(const string &exonName) {
  vector<string> parts = smithlab::split(exonName, "_exon_");
  assert(parts.size() == 2);
  return parts[0];
}

/**
 * \brief Does this exon partially overlap the given genomic region? This
 *        function doesn't consider an exact region match (i.e. same chrom,
 *        start, end, and strand) to be a 'partial overlap'.
 */
bool
Exon::partialOverlap(const GenomicRegion &r) {
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
Exon::sameRegion(const GenomicRegion &r) {
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
 * \brief TODO
 */
GenomicRegion
Exon::getGenomicRegion() { return this->region; }

/**
 * \brief get the name of this exon. The name of the exon is defined to be
 *        whatever the region name was when the Exon object was created.
 */
string
Exon::getExonName() { return this->region.get_name(); }

/**
 * \brief TODO
 */
string
Exon::getChrom() { return this->region.get_chrom(); }

/**
 * brief TODO
 */
char
Exon::getStrand() { return this->region.get_strand(); }
