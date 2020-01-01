package SeqDataModel;

import org.apache.commons.math3.distribution.BetaDistribution;

public class NonSpecificEnrichment {
    private BetaDistribution bd;

    public NonSpecificEnrichment(double a, double b) {
        this.bd = new BetaDistribution(a, b);
    }

    public double density(double nonSpecificEnrichmentRatio) {
        return this.bd.density(nonSpecificEnrichmentRatio);
    }

    public double logDensity(double nonSpecificEnrichmentRatio) {
        return this.bd.logDensity(nonSpecificEnrichmentRatio);
    }

    public double sample() {
        return this.bd.sample();
    }
}
