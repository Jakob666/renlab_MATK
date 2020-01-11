package Derivative;

/**
 * logarithm Gamma distribution function,
 *      log-Gamma distribution = (k-1)*log(σ) - σ/θ - logΓ(k) - k*logθ
 * where σ: negative binomial distribution overdispersion
 *       k: gamma distribution shape
 *       θ: gamma distribution scale
 */
public class LogGammaDistributionDerivative {

    public static double logGammaSecondOrderDerivative(double shape, double overdispersion) {
        return (1 - shape) / Math.pow(overdispersion, 2);
    }
}
