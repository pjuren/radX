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
  Gene(Exon e) : geneName(e.getGeneName()), chrom(e.get_chrom()),
                 strand(e.get_strand()) {};
  void addExon(Exon &e);
  void addExonSampleCount(GenomicRegion r, std::string sName, size_t count);
  void getReadcounts(const std::vector<std::string> &sampleNames,
                     std::vector<size_t> &res) const;
  std::string getName() const { return this->geneName; };
  typedef std::vector<Exon>::const_iterator  const_iterator;
  const_iterator begin() const { return exons.begin(); };
  const_iterator end() const { return exons.end(); };

private :
  std::string geneName;
  std::string chrom;
  char strand;
  std::vector<Exon> exons;
  std::tr1::unordered_map<std::string, size_t> geneReadcounts;
};

#endif // GENE_HPP_
