package HessianMatrix;

import java.math.BigDecimal;

public abstract class ModelDeterminate {
    protected HessianMatrix hm;

    public ModelDeterminate(double tretIPOverdispersion, double tretINPUTOverdispersion, double ctrlIPOverdispersion, double ctrlINPUTOverdispersion,
                            double tretMethylation, double ctrlMethylation, double tretBackgroundExpression, double ctrlBackgroundExpression,
                            double tretNonPeakExpression, double ctrlNonPeakExpression, double nonspecificEnrich,
                            double gammaShape, double betaScale1, double betaScale2) {
        this.hm = new HessianMatrix(tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                tretMethylation, ctrlMethylation, tretBackgroundExpression, ctrlBackgroundExpression,
                tretNonPeakExpression, ctrlNonPeakExpression, nonspecificEnrich, gammaShape, betaScale1, betaScale2);
        this.hm.initializeDerivatives();
    }

    public void setReads(int[] tretIPReads, int[] tretINPUTReads, int[] tretIPNonPeak, int[] tretINPUTNonPeak,
                         int[] ctrlIPReads, int[] ctrlINPUTReads, int[] ctrlIPNonPeak, int[] ctrlINPUTNonPeak) {
        this.hm.setReads(tretIPReads, tretINPUTReads, tretIPNonPeak, tretINPUTNonPeak,
                ctrlIPReads, ctrlINPUTReads, ctrlIPNonPeak,ctrlINPUTNonPeak);
    }

    public void setSizeFactors(double[] tretSampleSizeFactor, double[] ctrlSampleSizeFactor) {
        this.hm.setSampleSizeFactor(tretSampleSizeFactor, ctrlSampleSizeFactor);
    }

    public abstract BigDecimal modelDeterminate();
}
