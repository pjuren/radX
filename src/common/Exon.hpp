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

#ifndef EXON_HPP_
#define EXON_HPP_


// smithlab common code includes
#include "GenomicRegion.hpp"

// stl includes
#include <vector>
#include <string>
#include <tr1/unordered_map>

std::string getGeneName(const std::string &exonName);

class Exon {
public:
  Exon(const GenomicRegion& r) : region(r) {};
  bool partialOverlap(const GenomicRegion &r);
  bool sameRegion(const GenomicRegion &r);
  void addSampleReadCount(const std::string &sampleName, const size_t count);
  GenomicRegion getGenomicRegion();
  std::string getExonName();
  std::string getChrom();
  char getStrand();

private :
  GenomicRegion region;
  std::tr1::unordered_map<std::string, double> sampleCounts;

  Exon() {};
  void addSample(std::string sampleName, size_t count);


  //bool operator== (const Design &other_design) const;
  //bool operator!= (const Design &other_design) const;
};

#endif // EXON_HPP_
