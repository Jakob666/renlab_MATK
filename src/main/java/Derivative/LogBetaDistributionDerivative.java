package Derivative;

/**
 * logarithm Beta distribution function,
 *      log-Beta distribution = (w-1)*log(p) - (k-1)*log(1-p) - logΓ(w) - logΓ(k) + logΓ(w+k)
 * where p: methylation level or nonspecific enrichment ratio
 *       w: beta distribution scale parameter
 *       k: beta distribution scale parameter
 */
public class LogBetaDistributionDerivative {

    public static double logBetaSecondOrderDerivative(double w, double k, double level) {
        return (1 - w) / Math.pow(level, 2) - (k - 1) / Math.pow(1-level, 2);
    }
}
