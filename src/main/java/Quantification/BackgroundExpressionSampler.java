package Quantification;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * assume that the background expression of a gene follows log normal distribution as prior distribution
 */
public class BackgroundExpressionSampler extends MHSampling {
    private LogNormalDistribution lnd;
    private NormalDistribution nd;

    /**
     * Constructor
     * @param scale square root background expression mean
     * @param shape background expression std
     */
    public BackgroundExpressionSampler(double scale, double shape) {
        this.lnd = new LogNormalDistribution(scale, shape);
        this.nd = new NormalDistribution(0, 1);
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
