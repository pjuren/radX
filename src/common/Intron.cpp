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
#include "Intron.hpp"
#include "GenomicRegion.hpp"
#include "utils.hpp"

// using the following names
using std::string;
using std::vector;
using std::stringstream;

/******************************************************************************
 *                             HELPER FUNCTIONS
 ******************************************************************************/

/**
 * \brief split a full intron name of the form
 *        someGeneName_intron_someIntronName into its pair of
 *        <someGeneName, someIntronName>
 * \param fullName TODO
 * \return TODO
 * \throws SMITHLABException if the delimiter '_intron_' appears more than once
 *         or not at all in the full exon name.
 */
std::pair<string, string>
splitFullIntronName(const string &fullName) {
  const string delim = "_intron_";
  const size_t numDelimOcc = countOccurrences(delim, fullName);
  if (numDelimOcc == 0) {
    stringstream ss;
    ss << "incorrectly formatted full intron name: " << fullName << "; "
       << "could not find the expected delimiter " << delim << " to separate "
       << "the gene name from the intron name";
    throw SMITHLABException(ss.str());
  }
  if (numDelimOcc > 1) {
    stringstream ss;
    ss << "incorrectly formatted full intron name: " << fullName << "; "
       << "found multiple occurrences of the delimiter " << delim << ", so "
       << "could not unambiguously separate the gene name from the intron name";
    throw SMITHLABException(ss.str());
  }
  size_t s_delim = fullName.find(delim);
  size_t e_delim = s_delim + delim.size();
  return std::make_pair<string,string>(fullName.substr(0, s_delim),
                                       fullName.substr(e_delim, string::npos));
}

/**
 * \brief TODO
 * \param fn  TODO
 * \param res TODO
 * \return TODO
 */
void
readBEDFile (const string &fn, vector<Intron> &res) {
  res.clear();
  vector<GenomicRegion> tmp;
  ReadBEDFile(fn, tmp);
  for (size_t i = 0; i < tmp.size(); ++i) {
    res.push_back(Intron(tmp[i]));
  }
}

/******************************************************************************
 *                          EXON CLASS DEFINITION
 ******************************************************************************/

/**
 * \brief TODO
 */
Intron::Intron(const GenomicRegion& r) : AugmentedGenomicRegion(r) {
  std::pair<std::string, std::string> spltNm(splitFullIntronName(r.get_name()));
  this->geneName = spltNm.first;
  this->intronName = spltNm.second;
}

/**
 * \brief TODO
 */
Intron::Intron(const AugmentedGenomicRegion& r) : AugmentedGenomicRegion(r) {
  std::pair<std::string, std::string> spltNm(splitFullIntronName(r.get_name()));
  this->geneName = spltNm.first;
  this->intronName = spltNm.second;
}


