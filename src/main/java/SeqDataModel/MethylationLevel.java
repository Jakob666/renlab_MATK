package SeqDataModel;

import org.apache.commons.math3.distribution.BetaDistribution;

/**
 * Methylation level of a gene, follows Beta distribution
 */
public class MethylationLevel {
    private BetaDistribution bd;

    public MethylationLevel(double a, double b) {
        this.bd = new BetaDistribution(a, b);
    }

    public double density(double methylationLevel) {
        return this.bd.density(methylationLevel);
    }

    public double logDensity(double methylationLevel) {
        return this.bd.logDensity(methylationLevel);
    }

    public double sample() {
        return this.bd.sample();
    }
}
