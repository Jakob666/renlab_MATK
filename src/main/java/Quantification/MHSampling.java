package Quantification;

import org.apache.commons.math3.distribution.UniformRealDistribution;

public class MHSampling {
    private UniformRealDistribution random;

    public MHSampling() {
        this.random = new UniformRealDistribution(0, 1);
    }

    /**
     * get new round sampling result via random u and receptance
     * @param curSamplingDensity posterior density of sampling result of new round
     * @param prevSamplingDensity posterior density of sampling result of last round
     * @param log true, if logarithm probability
     * @return sampling result
     */
    public boolean getSamplingRes(double curSamplingDensity, double prevSamplingDensity, boolean log) {
        double u, receptance;
        u = this.getRandomValue();
        receptance = this.getReceptance(curSamplingDensity, prevSamplingDensity, log);
        return u < receptance;
    }

    /**
     * random sample from uniform(0, 1)
     * @return random value
     */
    private double getRandomValue() {
        return this.random.sample();
    }

    /**
     * calculate MH sampling receptance
     * @param curSamplingDensity density of current sampling result
     * @param prevSamplingDensity density of previous sampling result
     * @param log true, if logarithm probability
     * @return receptance
     */
    private double getReceptance(double curSamplingDensity, double prevSamplingDensity, boolean log) {
        if (log) {
            return Math.min(1, Math.exp(curSamplingDensity - prevSamplingDensity));
        }
        return Math.min(1, curSamplingDensity/prevSamplingDensity);
    }

}
