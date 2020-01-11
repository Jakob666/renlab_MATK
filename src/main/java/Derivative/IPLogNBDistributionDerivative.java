package Derivative;

import org.apache.commons.math3.special.Gamma;

/**
 * IP sample peak region logarithm negative binomial distribution function,
 *      log-NB = x*log(Uσ) - y*log(Uσ+1) - 1/σ * log(1/(Uσ+1)) - logΓ(x+1) - logΓ(1/σ) + logΓ(x+1+1/σ) - log(x+1/σ)
 * where x: observed reads count,
 *       U: negative binomial distribution expectation, U=s*E*(r+p)
 *          s is size factor, E is background expression, r is nonspecific enrichment, p is methylation level
 *       σ: negative binomial distribution overdispersion
 */
public class IPLogNBDistributionDerivative {
    private double backgroundExpression, nonspecificEnrichment, methylationLevel, ipOverdispersion;
    
    public IPLogNBDistributionDerivative(double backgroundExpression, double nonspecificEnrichment, 
                                         double methylationLevel, double ipOverdispersion) {
        this.backgroundExpression = backgroundExpression;
        this.nonspecificEnrichment = nonspecificEnrichment;
        this.methylationLevel = methylationLevel;
        this.ipOverdispersion = ipOverdispersion;
    }

    /**
     * d2 f / d sigma 2
     * @param sizeFactor size factor
     * @param readsCount observed reads count
     * @return second-order derivative
     */
    public double logNBSecondOrderDerivativeOverdispersion(double sizeFactor, double readsCount) {
        // exp * (r+p) * sizeFactor
        double erps = this.backgroundExpression * (this.nonspecificEnrichment + this.methylationLevel) * sizeFactor;
        // 1 + exp * (r+p) * sizeFactor * sigma
        double erpso = 1 + erps * this.ipOverdispersion;
        double inverse = Math.pow(this.ipOverdispersion, -1);
        double n = readsCount + 1 + inverse;

        double item1 = Math.pow(erps, 2) / (this.ipOverdispersion * Math.pow(erpso, 2));
        double item2 = 2 * erps / (Math.pow(this.ipOverdispersion, 2) * erpso);
        double item3 = readsCount / Math.pow(this.ipOverdispersion, 2);
        double item4 = Math.pow(erps, 2) * readsCount / Math.pow(erpso, 2);
        double item5 = 1 / (Math.pow(this.ipOverdispersion, 4) * Math.pow(inverse + readsCount, 2));
        double item6 = 2 / (Math.pow(this.ipOverdispersion, 3) * (inverse + readsCount));
        double item7 = 2 * Math.log(erpso) / Math.pow(this.ipOverdispersion, 3);
        double item8 = 2 * Gamma.digamma(inverse) / Math.pow(this.ipOverdispersion, 3);
        double item9 = 2 * Gamma.digamma(n) / Math.pow(this.ipOverdispersion, 3);
        double item10 = Gamma.trigamma(inverse) / Math.pow(this.ipOverdispersion, 4);
        double item11 = Gamma.trigamma(n) / Math.pow(this.ipOverdispersion, 4);

        return item1 + item2 - item3 + item4 + item5 - item6 - item7 - item8 + item9 - item10 + item11;
    }

    /**
     * d2 f / d sigma d r
     * @param sizeFactor size factor
     * @param readsCount observed reads count
     * @return second-order derivative
     */
    public double logNBOverdispersionNonspecificEnrich(double sizeFactor, double readsCount) {
        double erps = this.backgroundExpression * (this.nonspecificEnrichment + this.methylationLevel) * sizeFactor;
        double erpso = 1 + erps * this.ipOverdispersion;

        return this.backgroundExpression * erps * sizeFactor / Math.pow(erpso, 2) +
               this.backgroundExpression * (erpso-1) * sizeFactor * readsCount / Math.pow(erpso, 2) -
               this.backgroundExpression * sizeFactor * readsCount / erpso;
    }

    /**
     * d2 f / d sigma d p
     * @param sizeFactor size factor
     * @param readsCount observed reads count
     * @return second-order derivative
     */
    public double logNBOverdispersionMethylationLevel(double sizeFactor, double readsCount) {
        double erps = this.backgroundExpression * (this.nonspecificEnrichment + this.methylationLevel) * sizeFactor;
        double erpso = 1 + erps * this.ipOverdispersion;

        return this.backgroundExpression * erps * sizeFactor / Math.pow(erpso, 2) +
               this.backgroundExpression * (erpso-1) * sizeFactor * readsCount / Math.pow(erpso, 2) -
               this.backgroundExpression * sizeFactor * readsCount / erpso;
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
        // exp * (r+p) * sizeFactor
        double erps = this.backgroundExpression * (this.nonspecificEnrichment + this.methylationLevel) * sizeFactor;
        // 1 + exp * (r+p) * sizeFactor * sigma
        double erpso = 1 + erps * this.ipOverdispersion;

        return this.backgroundExpression * eso * sizeFactor / Math.pow(erpso, 2) -
               readsCount / Math.pow(this.methylationLevel + this.nonspecificEnrichment, 2) +
               Math.pow(eso, 2) * readsCount / Math.pow(erpso, 2);
    }

    /**
     * d2 f / d p 2
     * @param sizeFactor size factor
     * @param readsCount observed reads count
     * @return second-order derivative
     */
    public double logNBSecondOrderDerivativeMethylationLevel(double sizeFactor, double readsCount) {
        // exp * sizeFactor * sigma
        double eso = this.backgroundExpression * sizeFactor * this.ipOverdispersion;
        // exp * (r+p) * sizeFactor
        double erps = this.backgroundExpression * (this.nonspecificEnrichment + this.methylationLevel) * sizeFactor;
        // 1 + exp * (r+p) * sizeFactor * sigma
        double erpso = 1 + erps * this.ipOverdispersion;

        return this.backgroundExpression * eso * sizeFactor / Math.pow(erpso, 2) -
               readsCount / Math.pow(this.methylationLevel + this.nonspecificEnrichment, 2) +
               Math.pow(eso, 2) * readsCount / Math.pow(erpso, 2);
    }
}
