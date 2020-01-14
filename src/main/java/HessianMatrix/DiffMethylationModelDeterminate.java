package HessianMatrix;

import java.math.BigDecimal;
import java.util.Arrays;

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
        double secDerivTretIPOverdispersion = -1 * this.hm.secondDerivativeTreatmentIPOverdispersion();
        // d2 f / d tret_sigma_input 2
        double secDerivTretINPUTOverdispersion = -1 * this.hm.secondDerivativeTreatmentINPUTOverdispersion();
        // d2 f / d ctrl_sigma_ip 2
        double secDerivCtrlIPOverdispersion = -1 * this.hm.secondDerivativeControlIPOverdispersion();
        // d2 f / d ctrl_sigma_input 2
        double secDerivCtrlINPUTOverdispersion = -1 * this.hm.secondDerivativeControlINPUTOverdispersion();
        // d2 f / d r 2
        double secDerivNonspecificEnrich = -1 * this.hm.secondDerivativeNonspecificEnrichment();
        // d2 f / d pT 2
        double secDerivPTret = -1 * this.hm.secondDerivativeTreatmentMethylationLevel();
        // d2 f / d pC 2
        double secDerivPCtrl = -1 * this.hm.secondDerivativeControlMethylationLevel();

        double[] values = new double[] {secDerivTretIPOverdispersion, secDerivCtrlIPOverdispersion,
                                        secDerivTretINPUTOverdispersion, secDerivCtrlINPUTOverdispersion,
                                        secDerivPTret, secDerivPCtrl, secDerivNonspecificEnrich};

        BigDecimal res = new BigDecimal("1");
        for (double val: values) {
            res = res.multiply(new BigDecimal(Double.toString(val)));
        }

        return res;
    }
}
