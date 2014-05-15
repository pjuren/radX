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

#ifndef AUG_GEN_REG_HPP_
#define AUG_GEN_REG_HPP_

// smithlab common code includes
#include "GenomicRegion.hpp"

// stl includes
#include <vector>
#include <string>
#include <tr1/unordered_map>

/**
 * \brief An AugmentedGenomicRegion differs from a regular GenomicRegion by
 *        the inclusion of a set of read-counts that are indexed by the sample
 *        name they come from.
 */
class AugmentedGenomicRegion : public GenomicRegion {
public:
  /*--------------------------- CONSTRUCTORS --------------------------------*/
  /** construct an augmented genomic region from a regular genomic region */
  AugmentedGenomicRegion(const GenomicRegion& r) : GenomicRegion(r) {};

  /*---------------------------- INSPECTORS ---------------------------------*/
  bool partialOverlap(const GenomicRegion &r) const;
  bool sameRegion(const GenomicRegion &r) const;
  void getReadcounts(const std::vector<std::string> &sampleNames,
                     std::vector<size_t> &res) const;
  void getSampleNames(std::vector<std::string> &names) const;
  size_t getSampleCount() const;

  /*----------------------------- MUTATORS ----------------------------------*/
  void addSampleReadCount(const std::string &sampleName, const size_t count);

  /*----------------------- STATIC MEMBERS FOR IO ---------------------------*/
  static void readFromBEDFiles(const std::vector<std::string> &fns,
                               std::vector<AugmentedGenomicRegion> &res,
                               std::vector<std::string> &sampleNames,
                               const bool VERBOSE=false);
  static void readFromMatrixFile(const std::string &fn,
                                 std::vector<AugmentedGenomicRegion> &res,
                                 std::vector<std::string> &sampleNames,
                                 const bool VERBOSE=false);
private :
  /** map of read counts indexed by sample name for this region **/
  std::tr1::unordered_map<std::string, size_t> sampleCounts;
};

#endif // AUG_GEN_REG_HPP_
