package Quantification;

import SeqDataModel.OverdispersionGamma;
import org.apache.commons.math3.distribution.NormalDistribution;

public class OverdispersionSampler extends MHSampling {
    private OverdispersionGamma og;
    private NormalDistribution nd;

    public OverdispersionSampler(double shape, double scale) {
        this.og = new OverdispersionGamma(shape, scale);
        this.nd = new NormalDistribution(0, 0.05);
    }

    public double randomInit() {
        return this.og.sampling();
    }

    public double randomSample(double curMethLevel) {
        return Math.abs(curMethLevel + this.nd.sample());
    }

    public double getDensity(double overdispersion) {
        return this.og.density(overdispersion);
    }

    public double getLogDensity(double overdispersion) {
        return this.og.logDensity(overdispersion);
    }
}
