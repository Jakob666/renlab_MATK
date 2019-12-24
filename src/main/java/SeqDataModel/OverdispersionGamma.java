package SeqDataModel;

import org.apache.commons.math3.distribution.GammaDistribution;

/**
 * Overdispersion of negative binomial distribution follow Gamma distribution
 */
public class OverdispersionGamma {
    private GammaDistribution gd;

    public OverdispersionGamma(double a, double b) {
        this.gd = new GammaDistribution(a, b);
    }

    public double density(double overdispersion) {
        return this.gd.density(overdispersion);
    }

    public double logDensity(double overdispersion) {
        return this.gd.logDensity(overdispersion);
    }

    public double sampling() {
        return this.gd.sample();
    }
}
