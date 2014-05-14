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

// stl imports
#include <sstream>
#include <vector>
#include <string>

// Radix includes
#include "Gene.hpp"

// Bring the following into the global namespace
using std::vector;
using std::string;
using std::stringstream;
using std::tr1::unordered_map;


/**
 * \brief adds a read count for the given sample name in the exon defined by
 *        the given genomic region to this gene. If the genomic region exactly
 *        coincides with an existing exon definition in this gene, then the
 *        new sample read count is assumed to go to that exon, otherwise if
 *        (1) the genomic region partially overlaps an existing exon defined
 *        in this gene an error is thrown (overlapping exons are not allowed),
 *        or (2) the exon genomic region is disjoint from all others in the
 *        gene, then a new exon is added and this sample read count is assigned
 *        to it.
 */
void
Gene::addExonSampleCount(GenomicRegion r, string sName, size_t count) {
  // first check that it doesn't partially overlap any other exons
  for (size_t i = 0; i < this->exons.size(); ++i) {
    if (this->exons[i].partialOverlap(r)) {
      stringstream ss;
      ss << "A given gene cannot have overlapping exons, since we define gene "
         << "read counts as the sum of exon read counts in that gene. Exon "
         << r << " was found tooverlap exon "
         << this->exons[i];
      throw SMITHLABException(ss.str());
    }
  }

  // does it match an existing exon in this gene? If so, just go ahead and
  // add the sample count to that exon. It can't possibly match more than one.
  bool foundMatch = false;
  for (size_t i = 0; i < this->exons.size(); ++i) {
    if (this->exons[i].sameRegion(r)) {
      this->exons[i].addSampleReadCount(sName, count);
      foundMatch = true;
    }
  }

  // if we didn't find a match, add a new exon
  if (!foundMatch) {
    Exon e(r);
    e.addSampleReadCount(sName, count);
    this->addExon(e);
  }

  // also, add the count for this exon to the total for this gene in the
  // given sample
  if (this->geneReadcounts.find(sName) == this->geneReadcounts.end()) {
    this->geneReadcounts.insert(std::make_pair<string, int>(sName, 0));
  }
  unordered_map<string, size_t>::iterator it = this->geneReadcounts.find(sName);
  it->second += count;
}


/**
 * \brief Add an exon to this gene. Exons added to this gene must be on the
 *        same chromosome as other exons in this gene, and must be on the
 *        same strand.
 * \throws SMITHLABException if the exon e is on a different chromosome or
 *         strand to the other exons that have already been added to this gene.
 */
void
Gene::addExon(Exon &e) {
  if (e.get_chrom() != this->chrom) {
    stringstream ss;
    ss << "trying to add exon to gene " << this->geneName << " on chromosome "
       << e.get_chrom() << ", but this gene already has exons on chrom "
       << this->chrom;
    throw SMITHLABException(ss.str());
  }
  if (e.get_chrom() != this->chrom) {
    stringstream ss;
    ss << "trying to add exon to gene " << this->geneName << " on strand "
       << e.get_strand() << ", but this gene already has exons on strand "
       << this->strand;
    throw SMITHLABException(ss.str());
  }
  this->exons.push_back(e);
}

/**
 * \brief populate the given vector<size_t> with the read counts in this gene
 *        for the given sample names. The vector of resultant read counts will
 *        have the same order as the names that are presented in sNames
 * \param sNames   TODO
 * \param res      TODO
 */
void
Gene::getReadcounts(const vector<string> &sNames, vector<size_t> &res) const {
  res.clear();
  for (size_t i = 0; i < sNames.size(); ++i) {
    if (this->geneReadcounts.find(sNames[i]) == this->geneReadcounts.end()) {
      throw SMITHLABException("Lookup of read count for " + sNames[i] +\
                              " for gene " + this->geneName + " failed. "
                              "No read-count for that sample");
    }
    res.push_back(this->geneReadcounts.find(sNames[i])->second);
  }
}


