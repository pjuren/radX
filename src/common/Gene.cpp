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
#include <tr1/unordered_map>
#include <tr1/unordered_set>

// Radix includes
#include "Gene.hpp"
#include "utils.hpp"

// Bring the following into the global namespace
using std::vector;
using std::string;
using std::stringstream;
using std::tr1::unordered_map;
using std::tr1::unordered_set;
using std::make_pair;

/******************************************************************************
 **                        STATIC HELPER FUNCTIONS                           **
 ******************************************************************************/

/**
 * \brief get the set of genes which appear in gm1, but not in gm2
 */
 static unordered_set<string>
 missingGenes(const unordered_map< string, Gene > &gm1,
              const unordered_map< string, Gene > &gm2) {
   typedef unordered_map< string, Gene >::const_iterator geneMapIt;
   unordered_set<string> geneNames_gm1, geneNames_gm2;

   for (geneMapIt it = gm1.begin(); it != gm1.end(); ++it)
     geneNames_gm1.insert(it->second.getName());
   for (geneMapIt it = gm2.begin(); it != gm2.end(); ++it)
        geneNames_gm2.insert(it->second.getName());

   unordered_set<string> genesMissing;
   std::set_difference(geneNames_gm1.begin(), geneNames_gm1.end(),
                       geneNames_gm2.begin(), geneNames_gm2.end(),
                       std::inserter(genesMissing, genesMissing.end()));
   return genesMissing;
}

 /**
  * \brief get the set of samples that appear in gm1, but not in gm2
  */
static unordered_set<string>
missingSamples(unordered_map< string, Gene > &gm1,
               unordered_map< string, Gene > &gm2) {
  typedef unordered_map< string, Gene >::const_iterator geneMapIt;
  std::tr1::unordered_set<string> smplNames_gm1, smplNames_gm2;

  for (geneMapIt it = gm1.begin(); it != gm1.end(); ++it) {
    vector<string> sNames (it->second.getSampleNames());
    for (size_t j = 0; j < sNames.size(); ++j) smplNames_gm1.insert(sNames[j]);
  }
  for (geneMapIt it = gm2.begin(); it != gm2.end(); ++it) {
    vector<string> sNames (it->second.getSampleNames());
    for (size_t j = 0; j < sNames.size(); ++j) smplNames_gm2.insert(sNames[j]);
  }

  unordered_set<string> samplesMissing;
  std::set_difference(smplNames_gm1.begin(), smplNames_gm1.end(),
                      smplNames_gm2.begin(), smplNames_gm2.end(),
                      std::inserter(samplesMissing, samplesMissing.end()));
  return samplesMissing;

}

/******************************************************************************
 **             FUNCTIONS FOR MANIPULATING COLLECTIONS OF GENES              **
 ******************************************************************************/

/**
 * \brief Load gene read-counts from a set of exon read counts in either a
 *        single read-count matrix file, or a set of BED format read-count
 *        files.
 *
 *        The matrix file is tab separated. The first row is a header and lists
 *        the sample names for each sample included in the file. Each subsequent
 *        row represents an exon. The first five columns in an exon row are:
 *        chromosome, start index, end index, exon name, and exon strand.
 *        The exon name must follow the format: exonName_exon_geneName,
 *        where exonName and geneName can be any strings as long as they
 *        (1) don't contain the tab character and (2) uniquely identify the
 *        exon in the file. Exons in a gene must not overlap each other. Each
 *        subsequent column i after the first five gives the read count for the
 *        exon in the ith sample (matching the order in the header row).
 *
 *        BED format files are standard 6 column BED: chrom, start, end, name
 *        score, and strand. The score field contains the read count. Names
 *        must be formatted in the same way as the matrix-format file (above).
 *        Sample names will be taken from the names of the files after stripping
 *        their path and extension/suffix (i.e. everything after the last
 *        period).
 * \param exon_fns     [in] filename(s) to load the matrix from.
 * \param genes        [out] results (Gene objects) will be placed into this
 *                     map, where the keys (string) are gene names.
 * \param sampleNames  [out] the sample names from the matrix will be placed
 *                     into this vector. Anything that was in here before will
 *                     be cleared first.
 * \param VERBOSE      [in] print extra status messages if this is true.
 *                     Defaults to false if not set.
 * \TODO allow caller to specify a suffix to strip from filenames when
 *       generating sample names.
 */
void
readGenesExonsOnly(const vector<string> &exon_fns,
                   unordered_map< string, Gene > &genes,
                   vector<string> &sampleNames, const bool VERBOSE) {
  // note that the below two lines also populate the sampleNames vector
  vector<AugmentedGenomicRegion> regions;
  if (exon_fns.size() == 1)
    AugmentedGenomicRegion::readFromMatrixFile(exon_fns[0], regions,
                                               sampleNames, VERBOSE);
  else
    AugmentedGenomicRegion::readFromBEDFiles(exon_fns, regions,
                                             sampleNames, VERBOSE);

  for (size_t i = 0; i < regions.size(); ++i) {
    Exon e(regions[i]);
    // if we haven't seen this gene before, we make a new gene object for it
    unordered_map< string, Gene >::iterator loc = genes.find(e.getGeneName());
    if (loc == genes.end()) {
      genes.insert (make_pair<string, Gene>(e.getGeneName(), Gene(e)));
      loc = genes.find(e.getGeneName());
    }
    vector<size_t> readCounts;
    e.getReadcounts(sampleNames, readCounts);
    assert(readCounts.size() == sampleNames.size());
    for (size_t j = 0; j < readCounts.size(); j++)
      loc->second.addExonSampleCount(e, sampleNames[j], readCounts[j]);
  }
}

/**
 * \brief Load gene read-counts from a set of intron read counts in either a
 *        single read-count matrix file, or a set of BED format read-count
 *        files.
 *
 *        The matrix file is tab separated. The first row is a header and lists
 *        the sample names for each sample included in the file. Each subsequent
 *        row represents an intron. The first five columns in an intron row are:
 *        chromosome, start index, end index, intron name, and intron strand.
 *        The intron name must follow the format: intronName_intron_geneName,
 *        where intronName and geneName can be any strings as long as they
 *        (1) don't contain the tab character and (2) uniquely identify the
 *        intron in the file. Introns in a gene must not overlap each other.
 *        Each subsequent column i after the first five gives the read count for
 *        the intron in the ith sample (matching the order in the header row).
 *
 *        BED format files are standard 6 column BED: chrom, start, end, name
 *        score, and strand. The score field contains the read count. Names
 *        must be formatted in the same way as the matrix-format file (above).
 *        Sample names will be taken from the names of the files after stripping
 *        their path and extension/suffix (i.e. everything after the last
 *        period).
 * \param exon_fns     [in] filename(s) to load the matrix from.
 * \param genes        [out] results (Gene objects) will be placed into this
 *                     map, where the keys (string) are gene names.
 * \param sampleNames  [out] the sample names from the matrix will be placed
 *                     into this vector. Anything that was in here before will
 *                     be cleared first.
 * \param VERBOSE      [in] print extra status messages if this is true.
 *                     Defaults to false if not set.
 * \TODO allow caller to specify a suffix to strip from filenames when
 *       generating sample names.
 */
void
readGenesIntronsOnly(const vector<string> &intron_fns,
                     unordered_map< string, Gene > &genes,
                     vector<string> &sampleNames, const bool VERBOSE) {
  // note that the below two lines also populate the sampleNames vector
  vector<AugmentedGenomicRegion> regions;
  if (intron_fns.size() == 1)
    AugmentedGenomicRegion::readFromMatrixFile(intron_fns[0], regions,
                                               sampleNames, VERBOSE);
  else
    AugmentedGenomicRegion::readFromBEDFiles(intron_fns, regions,
                                             sampleNames, VERBOSE);

  for (size_t i = 0; i < regions.size(); ++i) {
    Intron e(regions[i]);
    // if we haven't seen this gene before, we make a new gene object for it
    unordered_map< string, Gene >::iterator loc = genes.find(e.getGeneName());
    if (loc == genes.end()) {
      genes.insert (make_pair<string, Gene>(e.getGeneName(), Gene(e)));
      loc = genes.find(e.getGeneName());
    }
    vector<size_t> readCounts;
    e.getReadcounts(sampleNames, readCounts);
    assert(readCounts.size() == sampleNames.size());
    for (size_t j = 0; j < readCounts.size(); j++)
      loc->second.addIntronSampleCount(e, sampleNames[j], readCounts[j]);
  }
}

/**
 * \breif TODO
 * \param exon_filenames   TODO
 * \param intron_filenames TODO
 * \param genes            TODO
 * \param sampleNames      TODO
 * \param VERBOSE          TODO
 */
void
readGenes(const vector<string> &exon_filenames,
          const vector<string> &intron_filenames,
          unordered_map< string, Gene > &genes,
          vector<string> &sampleNames, const bool VERBOSE) {
  unordered_map< string, Gene > genesExons, genesIntrons;
  vector<string> exonSampleNames, intronSampleNames;
  readGenesExonsOnly(exon_filenames, genesExons, exonSampleNames, VERBOSE);
  readGenesIntronsOnly(intron_filenames, genesIntrons,
                       intronSampleNames, VERBOSE);

  // check that the introns and exons files describe the same samples, and
  // the same genes
  unordered_set<string> genesMissingInIntrons (missingGenes(genesExons,
                                                            genesIntrons));
  unordered_set<string> genesMissingInExons (missingGenes(genesIntrons,
                                                          genesExons));
  unordered_set<string> samplesMissingInIntrons (missingSamples(genesExons,
                                                                genesIntrons));
  unordered_set<string> samplesMissingInExons (missingSamples(genesIntrons,
                                                              genesExons));
  if (genesMissingInExons.size() != 0) {
    throw SMITHLABException("Loading gene counts for exons and introns "
                            "failed. Reason: the following genes had counts "
                            "present for introns, but not exons: "            +\
                            join(genesMissingInExons));
  }
  if (genesMissingInIntrons.size() != 0) {
    throw SMITHLABException("Loading gene counts for exons and introns "
                            "failed. Reason: the following genes had counts "
                            "present for exons, but not introns: "            +\
                            join(genesMissingInIntrons));
  }
  if (samplesMissingInExons.size() != 0) {
    throw SMITHLABException("Loading gene counts for exons and introns "
                            "failed. Reason: the following samples had  "
                            "counts present for introns, but not exons: "     +\
                            join(genesMissingInExons));
  }
  if (samplesMissingInIntrons.size() != 0) {
    throw SMITHLABException("Loading gene counts for exons and introns "
                            "failed. Reason: the following samples had  "
                            "counts present for exons, but not introns: "     +\
                            join(genesMissingInIntrons));
  }

  // now we know all of the genes are present for both exons and introns, it's
  // okay just to iterate over one of those sets of gene names
  typedef unordered_map< string, Gene >::const_iterator GeneIter;
  genes.clear();
  for (GeneIter eIt = genesExons.begin(); eIt != genesExons.end(); ++eIt) {
    GeneIter iIt = genesIntrons.find(eIt->first);
    Gene tmp(eIt->second);
    genes.insert(std::make_pair(eIt->first, iIt->second.merge(tmp)));
  }

  // finally, fix up the sample names; again, we know from above that they're
  // the same for both exons and introns, so it's okay to just use on of those
  sampleNames.swap(exonSampleNames);
}

/******************************************************************************
 **                            THE GENE CLASS                                **
 ******************************************************************************/

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
Gene::addExonSampleCount(const GenomicRegion &r, const string &sName,
                         const size_t count) {
  // first check that it doesn't partially overlap any other exons
  for (size_t i = 0; i < this->exons.size(); ++i) {
    if (this->exons[i].partialOverlap(r)) {
      stringstream ss;
      ss << "A given gene cannot have overlapping exons, since we define gene "
         << "read counts as the sum of exon read counts in that gene. Exon "
         << r << " was found to overlap exon "
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
Gene::addExon(const Exon &e) {
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
 * \brief adds a read count for the given sample name in the intron defined by
 *        the given genomic region to this gene. If the genomic region exactly
 *        coincides with an existing intron definition in this gene, then the
 *        new sample read count is assumed to go to that intron, otherwise if
 *        (1) the genomic region partially overlaps an existing intron defined
 *        in this gene an error is thrown (overlapping intron are not allowed),
 *        or (2) the intron genomic region is disjoint from all others in the
 *        gene, then a new intron is added and this sample read count is
 *        assigned to it.
 */
void
Gene::addIntronSampleCount(const GenomicRegion &r, const string &sName,
                           const size_t count) {
  // first check that it doesn't partially overlap any other introns
  for (size_t i = 0; i < this->introns.size(); ++i) {
    if (this->introns[i].partialOverlap(r)) {
      stringstream ss;
      ss << "A given gene cannot have overlapping introns, since we define gene"
         << "read counts as the sum of intron and exon read counts in that "
         << "gene. Intron " << r << " was found to overlap exon "
         << this->introns[i];
      throw SMITHLABException(ss.str());
    }
  }

  // does it match an existing exon in this gene? If so, just go ahead and
  // add the sample count to that exon. It can't possibly match more than one.
  bool foundMatch = false;
  for (size_t i = 0; i < this->introns.size(); ++i) {
    if (this->introns[i].sameRegion(r)) {
      this->introns[i].addSampleReadCount(sName, count);
      foundMatch = true;
    }
  }

  // if we didn't find a match, add a new intron
  if (!foundMatch) {
    Intron e(r);
    e.addSampleReadCount(sName, count);
    this->addIntron(e);
  }

  // also, add the count for this intron to the total for this gene in the
  // given sample
  if (this->geneReadcounts.find(sName) == this->geneReadcounts.end()) {
    this->geneReadcounts.insert(std::make_pair<string, int>(sName, 0));
  }
  unordered_map<string, size_t>::iterator it = this->geneReadcounts.find(sName);
  it->second += count;
}

/**
 * \brief Add an intron to this gene. Introns added to this gene must be on the
 *        same chromosome as other introns in this gene, and must be on the
 *        same strand.
 * \throws SMITHLABException if the intron e is on a different chromosome or
 *         strand to the other introns that have already been added to this
 *         gene.
 */
void
Gene::addIntron(const Intron &e) {
  if (e.get_chrom() != this->chrom) {
    stringstream ss;
    ss << "trying to add intron to gene " << this->geneName << " on chromosome "
       << e.get_chrom() << ", but this gene already has introns and/or exons "
       << "on chromosome" << this->chrom;
    throw SMITHLABException(ss.str());
  }
  if (e.get_chrom() != this->chrom) {
    stringstream ss;
    ss << "trying to add intron to gene " << this->geneName << " on strand "
       << e.get_strand() << ", but this gene already has introns and/or exons "
       << "on strand " << this->strand;
    throw SMITHLABException(ss.str());
  }
  this->introns.push_back(e);
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

/**
 * \brief TODO
 * \return TODO
 */
vector<string>
Gene::getSampleNames() const {
  vector<string> res;
  typedef unordered_map<string, size_t>::const_iterator sIter;
  for (sIter it = this->geneReadcounts.begin();
             it != this->geneReadcounts.begin(); ++it) {
    res.push_back(it->first);
  }
  return res;
}

/**
 * \brief Merge this Gene object with another and return a new Gene object
 *        representing the union of the two
 *
 * The following conditions must be met for a successful merge:
 *      (1) same gene name
 *      (2) no exons (introns) in this partially overlap exons (introns)
 *          in other; introns/exons that represent the same genomic region
 *          and have the same name are merged and added to the result as one
 *          exon/intron
 *      (3) all exons/introns on same strand and chromosome.
 */
Gene
Gene::merge(const Gene &other) const {
  Gene res(other);
  for (size_t i = 0; i < this->exons.size(); ++i) {
    vector<string> sampleNames;
    vector<size_t> readCounts;
    this->exons[i].getSampleNames(sampleNames);
    this->exons[i].getReadcounts(sampleNames, readCounts);
    assert(sampleNames.size() == readCounts.size());
    for (size_t j = 0; j < sampleNames.size(); j++) {
      res.addExonSampleCount(this->exons[i], sampleNames[j], readCounts[j]);
    }
  }
  for (size_t i = 0; i < this->introns.size(); ++i) {
    vector<string> sampleNames;
    vector<size_t> readCounts;
    this->introns[i].getSampleNames(sampleNames);
    this->introns[i].getReadcounts(sampleNames, readCounts);
    assert(sampleNames.size() == readCounts.size());
    for (size_t j = 0; j < sampleNames.size(); j++) {
      res.addExonSampleCount(this->introns[i], sampleNames[j], readCounts[j]);
    }
  }
  return res;
}


