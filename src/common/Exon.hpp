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

#ifndef EXON_HPP_
#define EXON_HPP_

// smithlab common code includes
#include "GenomicRegion.hpp"

// RADIX includes
#include "AugmentedGenomicRegion.hpp"

// stl includes
#include <vector>
#include <string>
#include <tr1/unordered_map>

// forward declare this
class Exon;

// helper functions
std::pair<std::string, std::string> splitFullExonName(const std::string &name);
void readBEDFile (const std::string &fn, std::vector<Exon> &res);

/**
 * \brief An Exon is just an AugmentedGenomicRegion with some additional
 *        constraints on naming, and functionality for manipulating the name
 */
class Exon : public AugmentedGenomicRegion {
public :
  Exon(const GenomicRegion r) : AugmentedGenomicRegion(r) {
    std::pair<std::string, std::string> spltNm(splitFullExonName(r.get_name()));
    this->geneName = spltNm.first;
    this->exonName = spltNm.second;
  }
  std::string getExonName() const { return this->exonName; }
  std::string getGeneName() const { return this->geneName; }
private :
  std::string geneName;
  std::string exonName;
};

#endif // EXON_HPP_
