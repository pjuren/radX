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

#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <string>
#include <sstream>

size_t countOccurrences(const std::string &needle, const std::string &haystack,
                        const bool overlapping=false);

template<typename T>
std::string join(const T &cntnr, const std::string &delim = ", ") {
  std::stringstream ss;
  for (typename T::const_iterator it = cntnr.begin(); it != cntnr.end(); ++it) {
    if (it != cntnr.begin()) ss << delim;
    ss << (*it);
  }
  return ss.str();
}

#endif
