#ifndef UTILS_HPP
#define UTILS_HPP
#include <stdint.h>
#include <string>
#include <vector>


namespace kProcessor::utils{
bool has_suffix(const std::string& s, const std::string& suffix);
// Taken from
// https://stackoverflow.com/questions/19189014/how-do-i-find-files-with-a-specific-extension-in-a-directory-that-is-provided-by
std::vector<std::string> GetFilesExt(const char *dir, const char *ext);
std::string last_part(std::string str, char c);
std::string first_part(std::string str, char c);
bool FileExists(std::string);

/**
 * Code adopted from MCcortex
 * Pick a cleaning threshold from kmer coverage histogram. Assumes low coverage
 * kmers are all due to error. Fits a poisson with a gamma distributed mean.
 * Then chooses a cleaning threshold such than FDR (uncleaned kmers) occur at a
 * rate of < the FDR paramater.
 *
 * Translated from Gil McVean's initial proposed method in R code
 *
 * @param kmer_covg Histogram of kmer counts at coverages 1,2,.. arrlen-1
 * @param arrlen    Length of array kmer_covg
 * @param alpha_est_ptr If not NULL, used to return estimate for alpha
 * @param beta_est_ptr  If not NULL, used to return estimate for beta
 * @return -1 if no cut-off satisfies FDR, otherwise returns coverage cutoff
 */
int cleaning_pick_kmer_threshold(const uint64_t *kmer_covg, size_t arrlen,
                                 double *alpha_est_ptr, double *beta_est_ptr,
                                 double *false_pos_ptr, double *false_neg_ptr);



}

#endif
