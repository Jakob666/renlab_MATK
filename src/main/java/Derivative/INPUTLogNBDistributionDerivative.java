package Derivative;

import org.apache.commons.math3.special.Gamma;

/**
 * INPUT sample logarithm negative binomial distribution function,
 *      log-NB = y*log(Vσ) - y*log(Vσ+1) - 1/σ * log(1/(Vσ+1)) - logΓ(y+1) - logΓ(1/σ) + logΓ(y+1+1/σ) - log(y+1/σ)
 * where y: observed reads count,
 *       V: negative binomial distribution expectation, V=s*E, s is size factor, E is background expression
 *       σ: negative binomial distribution overdispersion
 */
public class INPUTLogNBDistributionDerivative {
    private double backgroundExpression, inputOverdispersion;

    public INPUTLogNBDistributionDerivative(double backgroundExpression, double inputOverdispersion) {
        this.backgroundExpression = backgroundExpression;
        this.inputOverdispersion = inputOverdispersion;
    }

    /**
     * d2 f / d sigma 2
     * @param sizeFactor size factor
     * @param readsCount observed reads count
     * @return second-order derivative
     */
    public double logNBSecondOrderDerivativeOverdispersion(double sizeFactor, double readsCount) {
        // exp * sizeFactor
        double es = this.backgroundExpression * sizeFactor;
        // 1 + exp * r * sizeFactor * sigma
        double eso = 1 + es * this.inputOverdispersion;
        double inverse = Math.pow(this.inputOverdispersion, -1);
        double n = readsCount + 1 + inverse;

        double item1 = Math.pow(es, 2) / (this.inputOverdispersion * Math.pow(eso, 2));
        double item2 = 2 * es / (Math.pow(this.inputOverdispersion, 2) * eso);
        double item3 = readsCount / Math.pow(this.inputOverdispersion, 2);
        double item4 = Math.pow(es, 2) * readsCount / Math.pow(eso, 2);
        double item5 = 1 / (Math.pow(this.inputOverdispersion, 4) * Math.pow(inverse + readsCount, 2));
        double item6 = 2 / (Math.pow(this.inputOverdispersion, 3) * (inverse + readsCount));
        double item7 = 2 * Math.log(eso) / Math.pow(this.inputOverdispersion, 3);
        double item8 = 2 * Gamma.digamma(inverse) / Math.pow(this.inputOverdispersion, 3);
        double item9 = 2 * Gamma.digamma(n) / Math.pow(this.inputOverdispersion, 3);
        double item10 = Gamma.trigamma(inverse) / Math.pow(this.inputOverdispersion, 4);
        double item11 = Gamma.trigamma(n) / Math.pow(this.inputOverdispersion, 4);

        return item1 + item2 - item3 + item4 + item5 - item6 - item7 - item8 + item9 - item10 + item11;
    }
}
