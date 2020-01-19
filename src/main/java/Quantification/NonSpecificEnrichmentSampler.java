package Quantification;

import SeqDataModel.NonSpecificEnrichment;
import org.apache.commons.math3.distribution.NormalDistribution;

@Deprecated
public class NonSpecificEnrichmentSampler extends MHSampling {
    private NonSpecificEnrichment nse;
    private NormalDistribution nd;

    public NonSpecificEnrichmentSampler(double w, double k) {
        this.nse = new NonSpecificEnrichment(w, k);
        this.nd = new NormalDistribution(0, 0.05);
    }

    public double randomInit() {
        return this.nse.sample();
    }

    public double randomSample(double curMethLevel) {
        return Math.min(1, Math.abs(curMethLevel + this.nd.sample()));
    }

    public double getDensity(double methLevel) {
        return this.nse.density(methLevel);
    }

    public double getLogDensity(double methLevel) {
        return this.nse.logDensity(methLevel);
    }
}
