package HessianMatrix;

import Derivative.LogBetaDistributionDerivative;

import java.math.BigDecimal;

/**
 * same methylation level model contains 6 parameters, the determinate contains 36 elements.
 */
public class SameMethylationModelDeterminate extends ModelDeterminate {
    private double betaScale1, betaScale2, methylationLevel;

    public SameMethylationModelDeterminate(double tretIPOverdispersion, double tretINPUTOverdispersion, double ctrlIPOverdispersion, double ctrlINPUTOverdispersion,
                                           double methylation, double tretBackgroundExpression, double ctrlBackgroundExpression,
                                           double tretNonPeakExpression, double ctrlNonPeakExpression, double nonspecificEnrich,
                                           double gammaShape, double betaScale1, double betaScale2) {
        super(tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion,
              methylation, methylation, tretBackgroundExpression, ctrlBackgroundExpression,
              tretNonPeakExpression, ctrlNonPeakExpression, nonspecificEnrich, gammaShape, betaScale1, betaScale2);
        this.betaScale1 = betaScale1;
        this.betaScale2 = betaScale2;
        this.methylationLevel = methylation;
    }

    /**
     * logarithm determinate
     * @return logarithm determinate
     */
    public BigDecimal modelDeterminate() {
        // d2 f / d tret_sigma_ip 2
        double secDerivTretIPOverdispersion = this.hm.secondDerivativeTreatmentIPOverdispersion();
        // d2 f / d tret_sigma_input 2
        double secDerivTretINPUTOverdispersion = this.hm.secondDerivativeTreatmentINPUTOverdispersion();
        // d2 f / d ctrl_sigma_ip 2
        double secDerivCtrlIPOverdispersion = this.hm.secondDerivativeControlIPOverdispersion();
        // d2 f / d ctrl_sigma_input 2
        double secDerivCtrlINPUTOverdispersion = this.hm.secondDerivativeControlINPUTOverdispersion();
        // d2 f / d r 2
        double secDerivNonspecificEnrich = this.hm.secondDerivativeNonspecificEnrichment();
        // d2 f / d p 2, here different with diff methylation level model
        double secDerivMethylationLevel = this.hm.secondDerivativeControlMethylationLevel() +
                                          this.hm.secondDerivativeTreatmentMethylationLevel() -
                                          LogBetaDistributionDerivative.logBetaSecondOrderDerivative(this.betaScale1, this.betaScale2, this.methylationLevel);

        double[] values = new double[] {secDerivTretIPOverdispersion, secDerivTretINPUTOverdispersion,
                                        secDerivCtrlIPOverdispersion, secDerivCtrlINPUTOverdispersion,
                                        secDerivMethylationLevel, secDerivNonspecificEnrich};

        BigDecimal res = new BigDecimal("1");
        for (int i=0; i<values.length; i++) {
            res = res.multiply(new BigDecimal(Double.toString(values[i])));
        }

        return res;
    }
}
