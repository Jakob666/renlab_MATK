package Quantification;

import SeqDataModel.MethylationLevel;
import org.apache.commons.math3.distribution.NormalDistribution;

public class MethylationLevelSampler extends MHSampling {
    private MethylationLevel methLevel;
    private NormalDistribution nd;

    /**
     * Constructor
     * @param w parameter of Beta distribution
     * @param k parameter of Beta distribution
     */
    public MethylationLevelSampler(double w, double k) {
        this.methLevel = new MethylationLevel(w, k);
        this.nd = new NormalDistribution(0, 0.05);
    }

    public double randomInit() {
        return this.methLevel.sample();
    }

    public double randomSample(double curMethLevel) {
        return Math.min(1, Math.abs(curMethLevel + this.nd.sample()));
    }

    public double getDensity(double methLevel) {
        return this.methLevel.density(methLevel);
    }

    public double getLogDensity(double methLevel) {
        return this.methLevel.logDensity(methLevel);
    }
}
