package HessianMatrix;

import java.math.BigDecimal;

public class DiffMethylationModelDeterminate extends ModelDeterminate{

    public DiffMethylationModelDeterminate(double tretIPOverdispersion, double tretINPUTOverdispersion, double ctrlIPOverdispersion, double ctrlINPUTOverdispersion,
                                           double tretMethylation, double ctrlMethylation, double tretBackgroundExpression, double ctrlBackgroundExpression,
                                           double tretNonPeakExpression, double ctrlNonPeakExpression, double nonspecificEnrich,
                                           double gammaShape, double betaScale1, double betaScale2) {
        super(tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion,
              tretMethylation, ctrlMethylation, tretBackgroundExpression, ctrlBackgroundExpression,
              tretNonPeakExpression, ctrlNonPeakExpression, nonspecificEnrich, gammaShape, betaScale1, betaScale2);
    }

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
        // d2 f / d pT 2
        double secDerivPTret = this.hm.secondDerivativeTreatmentMethylationLevel();
        // d2 f / d pC 2
        double secDerivPCtrl = this.hm.secondDerivativeControlMethylationLevel();

        double[] values = new double[] {secDerivTretIPOverdispersion, secDerivTretINPUTOverdispersion,
                                        secDerivCtrlIPOverdispersion, secDerivCtrlINPUTOverdispersion,
                                        secDerivPTret, secDerivPCtrl, secDerivNonspecificEnrich};

        BigDecimal res = new BigDecimal("1");
        for (int i=0; i<values.length; i++) {
            res = res.multiply(new BigDecimal(Double.toString(values[i])));
        }

        return res;
    }
}
