package Quantification;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

public class ExpansionEffectSampler extends MHSampling {
    private LogNormalDistribution lnd;
    private NormalDistribution nd;

    public ExpansionEffectSampler(double shape, double scale) {
        this.lnd = new LogNormalDistribution(scale, shape);
        this.nd = new NormalDistribution(0, 0.01);
    }

    public double randomSample(double curBackgroundExpression) {
        return Math.abs(curBackgroundExpression + this.nd.sample());
    }

    public double getDensity(double backgroundExpression) {
        return this.lnd.density(backgroundExpression);
    }

    public double getLogDensity(double backgroundExpression) {
        return this.lnd.logDensity(backgroundExpression);
    }
}
