#include <stdint.h>
#include <string>
#include <gqf.h>
#include <sys/stat.h>
#include <vector>
#include <limits>
#include <math.h>
#include <algorithm>


namespace kProcessor::utils {
// void removeReadsWithN(std::string inputFilename,std::string outputFilename)
// {
//   SeqFileIn seqFileIn(inputFilename.c_str());
//   SeqFileOut seqFileOut(outputFilename.c_str());
//   CharString id;
//   CharString quals;
//   std::string readT;
//
//   while(!atEnd(seqFileIn))
//   {
//     bool hasN=false;
//     readRecord(id, readT,quals, seqFileIn);
//     for(uint64_t i=0;i<readT.size();i++)
//     {
//       if(readT[i]=='N')
//       {
//         hasN=true;
//         break;
//       }
//     }
//     if(!hasN)
//     {
//       writeRecord(seqFileOut,id,readT,quals);
//     }
//   }
//
//
// }
    bool has_suffix(const std::string &s, const std::string &suffix) {
        return (s.size() >= suffix.size()) && equal(suffix.rbegin(),
                                                    suffix.rend(), s.rbegin());
    }
// Taken from
// https://stackoverflow.com/questions/19189014/how-do-i-find-files-with-a-specific-extension-in-a-directory-that-is-provided-by
// TO BE REMOVED TODO V2
/*
std::vector<std::string> GetFilesExt(const char *dir, const char *ext)
  {
    DIR *folder = opendir(dir);
		if (!folder) {
			std::cerr << "Directory doesn't exist " << dir << std::endl;
			exit(1);
		}

		std::vector<std::string> ret;
		dirent *entry;
		while((entry = readdir(folder)) != NULL)
		{
			if(has_suffix(entry->d_name, ext))
			{
				std::string filename(entry->d_name);
				std::string dirname(dir);
				ret.push_back(std::string(dirname + filename));
			}
		}

		return ret;
}
*/

    std::string last_part(std::string str, char c) {
        uint64_t found = str.find_last_of(c);
        return str.substr(found + 1);
    }

    std::string first_part(std::string str, char c) {
        uint64_t found = str.find_first_of(c);
        return str.substr(0, found);
    }

// copied from mantis
// Taken from
// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
    bool FileExists(std::string Spath) {
        const char *path = Spath.c_str();
        struct stat fileStat;
        if (stat(path, &fileStat)) {
            return false;
        }
        if (!S_ISREG(fileStat.st_mode)) {
            return false;
        }
        return true;
    }


// Find cutoff by finding first coverage level where errors make up less than
// `fdr` of total coverage
// returns -1 if not found
    static inline int pick_cutoff_with_fdr_thresh(const double *e_covg,
                                                  const uint64_t *kmer_covg,
                                                  size_t arrlen, double fdr) {
        size_t i;
        for (i = 1; i < arrlen; i++) {
            // printf(" %zu: %f %zu test: %f < %f\n", i, e_covg[i], kmer_covg[i],
            //                                       e_covg[i] / kmer_covg[i], fdr);
            if (e_covg[i] / kmer_covg[i] <= fdr) {
                return i;
            }
        }
        return -1;
    }

    // Get highest cutoff where false-positives < false-negatives
    // i.e. proportion of real kmers we are removing is less than the
    //      proportion of bad kmers we are keeping
    // returns -1 if not found
    static inline int pick_cutoff_FP_lt_FN(const double *e_covg, double e_total,
                                           const uint64_t *kmer_covg, uint64_t d_total,
                                           size_t arrlen) {
        size_t i;
        // for(i = 1; i < arrlen; i++) { printf("  %zu", kmer_covg[i]); } printf("\n");
        // for(i = 1; i < arrlen; i++) { printf("  %f", e_covg[i]); } printf("\n");
        // printf(" e_total: %f d_total: %f\n", e_total, (double)d_total);
        double e_rem = e_total, d_rem = d_total;
        double e_sum = 0, d_sum = 0;
        for (i = 1; i < arrlen; i++) {
            e_sum += e_covg[i];
            d_sum += kmer_covg[i];
            e_rem -= e_covg[i];
            d_rem -= kmer_covg[i];
            // printf(" %zu: e_total: %f d_total: %f\n", i, e_rem, d_rem);
            if (1 - e_sum / d_sum > e_rem / d_rem) {
                return i;
            }
        }
        return -1;
    }

    static inline int pick_cutoff_loss_vs_error(const double *e_covg,
                                                double e_total,
                                                const uint64_t *kmer_covg,
                                                size_t arrlen) {
        size_t i;
        // for(i = 1; i < arrlen; i++) { printf("  %zu", kmer_covg[i]); } printf("\n");
        // for(i = 1; i < arrlen; i++) { printf("  %f", e_covg[i]); } printf("\n");
        // printf(" e_total: %f d_total: %f\n", e_total, (double)d_total);
        double e_rem = e_total;
        double e_sum = 0, d_sum = 0;
        for (i = 1; i < arrlen; i++) {
            e_sum += e_covg[i];
            d_sum += kmer_covg[i];
            e_rem -= e_covg[i];
            double lost_seq = (d_sum - e_sum);
            double rem_err = e_rem;
            // printf(" %zu: e_total: %f d_total: %f\n", i, e_rem, d_rem);
            if (lost_seq > rem_err) return i;
        }
        return -1;
    }

    static inline void cutoff_get_FP_FN(const double *e_covg, double e_total,
                                        const uint64_t *kmer_covg, uint64_t d_total,
                                        size_t cutoff,
                                        double *false_pos, double *false_neg) {
        size_t i;
        double e_rem = e_total, d_rem = d_total;
        double e_sum = 0, d_sum = 0;
        for (i = 1; i < cutoff; i++) {
            e_sum += e_covg[i];
            d_sum += kmer_covg[i];
            e_rem -= e_covg[i];
            d_rem -= kmer_covg[i];
        }
        *false_pos = 1 - e_sum / d_sum;
        *false_neg = e_rem / d_rem;
    }

    // Check if at least `frac_covg_kept` coverage is kept when using threshold
    static inline bool is_cutoff_good(const uint64_t *kmer_covg, size_t arrlen,
                                      size_t cutoff, double frac_covg_kept) {
        uint64_t kmers_below = 0, kmers_above = 0;
        size_t i;
        for (i = 0; i < cutoff; i++) kmers_below += kmer_covg[i] * i;
        for (i = cutoff; i < arrlen; i++) kmers_above += kmer_covg[i] * i;

        // At least 20% of kmers should be kept
        return !arrlen || // any cutoff is good if no kmers
               ((double) kmers_above / (kmers_below + kmers_above) >= frac_covg_kept);
    }


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
                                     double *false_pos_ptr, double *false_neg_ptr) {

        if (kmer_covg[0] == 0)
            throw std::logic_error("Shouldn't see any kmers with coverage zero");


        size_t i, min_a_est_idx = 0;
        double r1, r2, rr, min_a_est = std::numeric_limits<double>::max(), tmp;
        double aa, faa, a_est, b_est, c0;

        r1 = (double) kmer_covg[2] / kmer_covg[1];
        r2 = (double) kmer_covg[3] / kmer_covg[2];
        rr = r2 / r1;

        // printf("r1: %.2f r2: %.2f rr: %.2f\n", r1, r2, rr);

        // iterate aa = { 0.01, 0.02, ..., 1.99, 2.00 }
        // find aa value that minimises abs(faa-rr)
        for (i = 1; i <= 200; i++) {
            aa = i * 0.01;
            faa = tgamma(aa) * tgamma(aa + 2) / (2 * pow(tgamma(aa + 1), 2));
            tmp = fabs(faa - rr);
            if (tmp < min_a_est) {
                min_a_est = tmp;
                min_a_est_idx = i;
            }
        }

        // a_est, b_est are estimates for alpha, beta of gamma distribution
        a_est = min_a_est_idx * 0.01;
        b_est = tgamma(a_est + 1.0) / (r1 * tgamma(a_est)) - 1.0;
        b_est = std::max(b_est, 1.0); // Avoid beta values <1
        c0 = kmer_covg[1] * pow(b_est / (1 + b_est), -a_est);

        if (alpha_est_ptr) *alpha_est_ptr = a_est;
        if (beta_est_ptr) *beta_est_ptr = b_est;

        // printf("min_a_est_idx: %zu\n", min_a_est_idx);
        // printf("a_est: %f b_est %f c0: %f\n", a_est, b_est, c0);

        // keep coverage estimates on the stack - this should be ok
        double e_covg_tmp, e_covg[arrlen];
        double e_total = 0;
        uint64_t d_total = 0;

        // Calculate some values here for speed
        double log_b_est = log(b_est);
        double log_one_plus_b_est = log(1 + b_est);
        double lgamma_a_est = lgamma(a_est);

        // note: lfactorial(x) = lgamma(x+1)

        for (i = 1; i < arrlen; i++) {
            e_covg_tmp = a_est * log_b_est - lgamma_a_est - lgamma(i)
                         + lgamma(a_est + i - 1)
                         - (a_est + i - 1) * log_one_plus_b_est;
            e_covg[i] = exp(e_covg_tmp) * c0;
            e_total += e_covg[i];
            d_total += kmer_covg[i];
        }

        // for(i = 1; i < MIN2(arrlen,100); i++)
        //   printf("  %zu: %f %zu\n", i, e_covg[i], (size_t)kmer_covg[i]);

        int cutoff = -1;

        // Find cutoff by finding first coverage level where errors make up less than
        // 0.1% of total coverage
        cutoff = pick_cutoff_with_fdr_thresh(e_covg, kmer_covg, arrlen, 0.001);
        // printf("A cutoff: %i\n", cutoff);

        // Pick highest cutoff that keeps FP < FN
        if (cutoff < 0)
            cutoff = pick_cutoff_FP_lt_FN(e_covg, e_total, kmer_covg, d_total, arrlen);

        if (cutoff < 0)
            cutoff = pick_cutoff_loss_vs_error(e_covg, e_total, kmer_covg, arrlen);

        // printf("B cutoff: %i\n", cutoff);

        if (cutoff < 0) return -1;

        // printf("C cutoff: %i\n", cutoff);

        // Check cutoff keeps at least 20% of coverage
        // (WGS should be much higher, Exome sequencing needs low cutoff)
        if (!is_cutoff_good(kmer_covg, arrlen, cutoff, 0.2)) return -1;

        // printf("D cutoff: %i\n", cutoff);

        // Calculate FP,FN rates
        if (false_pos_ptr || false_neg_ptr) {
            double false_pos = 0, false_neg = 0;
            cutoff_get_FP_FN(e_covg, e_total, kmer_covg, d_total, cutoff,
                             &false_pos, &false_neg);
            // printf("  FP: %f, FN: %f\n", false_pos, false_neg);
            if (false_pos_ptr) *false_pos_ptr = false_pos;
            if (false_neg_ptr) *false_neg_ptr = false_neg;
        }

        // printf(" kmers_above : %zu / (%zu + %zu) = %f\n",
        //        kmers_above, kmers_below, kmers_above,
        //        (double)kmers_above/(kmers_below+kmers_above));

        // printf("cutoff: %i\n", cutoff);

        // printf(" cutoff: %zu fdr: %f fdr_limit: %f good: %i\n",
        //        cutoff, fdr, fdr_limit, (int)good_cutoff);

        return cutoff;

    }


}
