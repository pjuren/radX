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

#ifndef GENE_HPP_
#define GENE_HPP_

// RADIX common includes
#include "Exon.hpp"

// Smithlab common includes
#include "GenomicRegion.hpp"

// stl includes
#include <vector>
#include <string>
#include <tr1/unordered_map>

class Gene {
public:
  Gene(Exon e) : geneName(getGeneName(e.getExonName())), chrom(e.getChrom()),
                 strand(e.getStrand()) {};
  void addExon(Exon &e);
  void addExonSampleCount(GenomicRegion r, std::string sName, size_t count);

private :
  std::string geneName;
  std::string chrom;
  char strand;
  std::vector<Exon> exons;
};

#endif // GENE_HPP_
