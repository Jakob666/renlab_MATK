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
        double secDerivTretIPOverdispersion = -1 * this.hm.secondDerivativeTreatmentIPOverdispersion();
        // d2 f / d tret_sigma_input 2
        double secDerivTretINPUTOverdispersion = -1 * this.hm.secondDerivativeTreatmentINPUTOverdispersion();
        // d2 f / d ctrl_sigma_ip 2
        double secDerivCtrlIPOverdispersion = -1 * this.hm.secondDerivativeControlIPOverdispersion();
        // d2 f / d ctrl_sigma_input 2
        double secDerivCtrlINPUTOverdispersion = -1 * this.hm.secondDerivativeControlINPUTOverdispersion();
        // d2 f / d r 2
        double secDerivNonspecificEnrich = -1 * this.hm.secondDerivativeNonspecificEnrichment();
        // d2 f / d p 2, here different with diff methylation level model
        double secDerivMethylationLevel = -1 * (this.hm.secondDerivativeControlMethylationLevel() +
                                          this.hm.secondDerivativeTreatmentMethylationLevel() -
                                          LogBetaDistributionDerivative.logBetaSecondOrderDerivative(this.betaScale1, this.betaScale2, this.methylationLevel));

        double[] values = new double[] {secDerivTretIPOverdispersion, secDerivCtrlIPOverdispersion,
                                        secDerivTretINPUTOverdispersion, secDerivCtrlINPUTOverdispersion,
                                        secDerivMethylationLevel, secDerivNonspecificEnrich};

        BigDecimal res = new BigDecimal("-1");
        for (double val: values) {
            res = res.multiply(new BigDecimal(val));
        }

        return res;
    }
}
