package SeqDataModel;

import org.apache.commons.math3.special.Beta;

/**
 * Negative binomial distribution with parameter miu and overdispersion a, the density compute as
 *                      (a*miu/a*miu+1)^K (1/a*miu+1)^(-1/a)
 *      p(k|miu, a) = ----------------------------------------
 *                              Beta(k+1, 1/a) (k+1/a)
 */
public class NegativeBinomialDistribution {
    private double miu, overdispersion;

    /**
     * Constructor
     * @param miu expectation of negative binomial distribution
     * @param overdispersion overdispersion of negative binomial distribution
     */
    public NegativeBinomialDistribution(double miu, double overdispersion) {
        this.miu = miu;
        this.overdispersion = overdispersion;
    }

    /**
     * correspond density of a sample number
     * @param sampleNum sampling number
     * @return density
     */
    public double density(double sampleNum) {
        double numerator = Math.pow(this.miu * this.overdispersion/(this.miu * this.overdispersion+1), sampleNum) * Math.pow(1/(this.overdispersion*this.miu+1), 1/this.overdispersion);
        double beta = this.betaFunction(sampleNum+1, 1/this.overdispersion);
        double denominator = beta * (sampleNum+1/this.overdispersion);

        return numerator / denominator;
    }

    /**
     * correspond log density of a sample number, more robust than density method
     * @param sampleNum sampling number
     * @return log density
     */
    public double logDensity(double sampleNum) {
        double d = this.miu * this.overdispersion;
        double logNumerator = sampleNum * (Math.log(d) - Math.log(d + 1)) - Math.pow(this.overdispersion, -1) * Math.log(d+1);
        double logDenominator = 1 * this.logBetaFunction(sampleNum+1, Math.pow(this.overdispersion, -1)) + Math.log(sampleNum+Math.pow(this.overdispersion, -1));

        return logNumerator - logDenominator;
    }

    private double logBetaFunction(double a, double b) {
        return Beta.logBeta(a, b);
    }

    private double betaFunction(double a, double b) {
        return Math.exp(Beta.logBeta(a, b));
    }

}
