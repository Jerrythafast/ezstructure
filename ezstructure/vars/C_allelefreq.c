
#include "allelefreq.h"
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

void P_update_simple(const uint8_t* G, const double* zetabeta, const double* zetagamma, const double* xi, const double* beta, const double* gamma, double* var_beta, double* var_gamma, long N, long L, long K)
{
    uint8_t genotype;
    long idx, n, l, k;
    double theta_beta_sum, theta_gamma_sum;
    double *var_beta_tmp, *var_gamma_tmp;

    var_beta_tmp = (double*) malloc(K * sizeof(double));
    var_gamma_tmp = (double*) malloc(K * sizeof(double));

    // loop over loci
    for (l=0; l<L; l++) {

        for (k=0; k<K; k++) {

            var_beta_tmp[k] = 0.0;
            var_gamma_tmp[k] = 0.0;
        }

        // loop over samples
        for (n=0; n<N; n++) {

            genotype = G[n*L+l];

            // missing data do not contribute
            if (genotype!=3) {

                // compute xi*zeta_{beta,gamma}
                theta_beta_sum = 0.0;
                theta_gamma_sum = 0.0;
                for (k=0; k<K; k++) {
                    theta_beta_sum += xi[n*K+k] * zetabeta[l*K+k];
                    theta_gamma_sum += xi[n*K+k] * zetagamma[l*K+k];
                }

                // increment var_{beta,gamma}_tmp
                for (k=0; k<K; k++) {
                    var_beta_tmp[k] += (double) genotype * xi[n*K+k] / theta_beta_sum;
                    var_gamma_tmp[k] += (double) (2-genotype) * xi[n*K+k] / theta_gamma_sum;
                }
            }
        }

        // compute var_{beta,gamma}
        for (k=0; k<K; k++) {
            idx = l*K+k;
            var_beta[idx] = beta[idx] + zetabeta[idx] * var_beta_tmp[k];
            var_gamma[idx] = gamma[idx] + zetagamma[idx] * var_gamma_tmp[k];
        }
    }

    free( var_beta_tmp );
    free( var_gamma_tmp );
}
