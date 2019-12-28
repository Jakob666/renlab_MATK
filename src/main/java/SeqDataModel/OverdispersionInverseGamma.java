package SeqDataModel;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

/**
 * Inverse-Chi square distribution
 */
public class OverdispersionInverseChi {
    private double miu;
    private ChiSquaredDistribution csd;

    public OverdispersionInverseChi(double miu) {
        assert miu-2>0.00001;
        this.miu = miu;
        this.csd = new ChiSquaredDistribution(miu);
    }

    public double density(double x) {
        return this.csd.density(x) / Math.pow(x, this.miu) * Math.exp(x/2-1/(2*x));
    }

    public double logDensity(double x) {
        return this.csd.logDensity(x) - this.miu * Math.log(x) + x / 2 - 1 / (2*x);
    }

    public double sample() {
        return Math.pow(this.csd.sample(), -1);
    }
}
