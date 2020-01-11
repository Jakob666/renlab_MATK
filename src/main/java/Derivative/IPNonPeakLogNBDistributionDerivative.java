package Derivative;

import org.apache.commons.math3.special.Gamma;

/**
 * IP sample non-peak region logarithm negative binomial distribution function,
 *      log-NB = x*log(Uσ) - y*log(Uσ+1) - 1/σ * log(1/(Uσ+1)) - logΓ(x+1) - logΓ(1/σ) + logΓ(x+1+1/σ) - log(x+1/σ)
 * where x: observed reads count,
 *       U: negative binomial distribution expectation, U=s*E*r;
 *          s is size factor, E is background expression, r is nonspecific enrichment
 *       σ: negative binomial distribution overdispersion
 */
public class IPNonPeakLogNBDistributionDerivative {
    private double backgroundExpression, nonspecificEnrichment, ipOverdispersion;

    public IPNonPeakLogNBDistributionDerivative(double backgroundExpression, double nonspecificEnrichment, double ipOverdispersion) {
        this.backgroundExpression = backgroundExpression;
        this.nonspecificEnrichment = nonspecificEnrichment;
        this.ipOverdispersion = ipOverdispersion;
    }

    /**
     * d2 f / d sigma 2
     * @param sizeFactor size factor
     * @param readsCount observed reads count
     * @return second-order derivative
     */
    public double logNBSecondOrderDerivativeOverdispersion(double sizeFactor, double readsCount) {
        // exp * r * sizeFactor
        double ers = this.backgroundExpression * this.nonspecificEnrichment * sizeFactor;
        // 1 + exp * (r+p) * sizeFactor * sigma
        double erso = 1 + ers * this.ipOverdispersion;
        double inverse = Math.pow(this.ipOverdispersion, -1);
        double n = readsCount + 1 + inverse;

        double item1 = Math.pow(ers, 2) / (this.ipOverdispersion * Math.pow(erso, 2));
        double item2 = 2 * ers / (Math.pow(this.ipOverdispersion, 2) * erso);
        double item3 = readsCount / Math.pow(this.ipOverdispersion, 2);
        double item4 = Math.pow(ers, 2) * readsCount / Math.pow(erso, 2);
        double item5 = 1 / (Math.pow(this.ipOverdispersion, 4) * Math.pow(inverse + readsCount, 2));
        double item6 = 2 / (Math.pow(this.ipOverdispersion, 3) * (inverse + readsCount));
        double item7 = 2 * Math.log(erso) / Math.pow(this.ipOverdispersion, 3);
        double item8 = 2 * Gamma.digamma(inverse) / Math.pow(this.ipOverdispersion, 3);
        double item9 = 2 * Gamma.digamma(n) / Math.pow(this.ipOverdispersion, 3);
        double item10 = Gamma.trigamma(inverse) / Math.pow(this.ipOverdispersion, 4);
        double item11 = Gamma.trigamma(n) / Math.pow(this.ipOverdispersion, 4);

        return item1 + item2 - item3 + item4 + item5 - item6 - item7 - item8 + item9 - item10 + item11;
    }

    /**
     * d2 f / d r 2
     * @param sizeFactor size factor
     * @param readsCount observed reads count
     * @return second-order derivative
     */
    public double logNBSecondOrderDerivativeNonspecificEnrich(double sizeFactor, double readsCount) {
        // exp * sizeFactor * sigma
        double eso = this.backgroundExpression * sizeFactor * this.ipOverdispersion;
        // exp * r * sizeFactor
        double erps = this.backgroundExpression * this.nonspecificEnrichment * sizeFactor;
        // 1 + exp * r * sizeFactor * sigma
        double erpso = 1 + erps * this.ipOverdispersion;

        return this.backgroundExpression * eso * sizeFactor / Math.pow(erpso, 2) -
               readsCount / Math.pow(this.nonspecificEnrichment, 2) +
               Math.pow(eso, 2) * readsCount / Math.pow(erpso, 2);
    }
}
