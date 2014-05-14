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

// radix includes
#include "utils.hpp"

// stl includes
#include <string>
#include <vector>

// using these names
using std::string;
using std::vector;

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

