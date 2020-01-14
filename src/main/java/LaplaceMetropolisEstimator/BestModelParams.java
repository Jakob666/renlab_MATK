package LaplaceMetropolisEstimator;

import DifferentialMethylation.ModelSelection;

/**
 * parameters which lead to the maximum logarithm posterior probability
 */
public class BestModelParams {
    private double[] treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion,
                     treatmentMethylationLevel, controlMethylationLevel, nonspecificEnrichment;
    private double[] tretIPSizeFactors, tretINPUTSizeFactors, tretIPNonPeakSizeFactors, tretINPUTNonPeakSizeFactors,
                     ctrlIPSizeFactors, ctrlINPUTSizeFactors, ctrlIPNonPeakSizeFactors, ctrlINPUTNonPeakSizeFactors;
    private double maximumLogPosterior = -1 * Double.MAX_VALUE;
    private int samplingTime, burnIn, paramNum;
    private ModelSelection model;

    public BestModelParams(double[] treatmentIPOverdispersion, double[] treatmentINPUTOverdispersion,
                           double[] controlIPOverdispersion, double[] controlINPUTOverdispersion,
                           double[] treatmentMethylationLevel, double[] controlMethylationLevel, double[] nonspecificEnrichment,
                           ModelSelection model, int burnIn) {
        this.treatmentIPOverdispersion = treatmentIPOverdispersion;
        this.treatmentINPUTOverdispersion = treatmentINPUTOverdispersion;
        this.controlIPOverdispersion = controlIPOverdispersion;
        this.controlINPUTOverdispersion = controlINPUTOverdispersion;
        this.treatmentMethylationLevel = treatmentMethylationLevel;
        this.controlMethylationLevel = controlMethylationLevel;
        this.nonspecificEnrichment = nonspecificEnrichment;
        this.samplingTime = treatmentIPOverdispersion.length;
        this.burnIn = burnIn;
        if (controlMethylationLevel == null)
            this.paramNum = 6;      // same methylation level model parameter number
        else
            this.paramNum = 7;      // diff methylation level model parameter number
        this.model = model;
    }

    /**
     * set sample size factors
     */
    public void setSizeFactors(double[] tretIPSizeFactors, double[] tretINPUTSizeFactors, double[] tretIPNonPeakSizeFactors, double[] tretINPUTNonPeakSizeFactors,
                               double[] ctrlIPSizeFactors, double[] ctrlINPUTSizeFactors, double[] ctrlIPNonPeakSizeFactors, double[] ctrlINPUTNonPeakSizeFactors) {
        this.tretIPSizeFactors = tretIPSizeFactors;
        this.tretINPUTSizeFactors = tretINPUTSizeFactors;
        this.tretIPNonPeakSizeFactors = tretIPNonPeakSizeFactors;
        this.tretINPUTNonPeakSizeFactors = tretINPUTNonPeakSizeFactors;
        this.ctrlIPSizeFactors = ctrlIPSizeFactors;
        this.ctrlINPUTSizeFactors = ctrlINPUTSizeFactors;
        this.ctrlIPNonPeakSizeFactors = ctrlIPNonPeakSizeFactors;
        this.ctrlINPUTNonPeakSizeFactors = ctrlINPUTNonPeakSizeFactors;
    }

    private void calcParamsProbability() {
        int keepTime = this.samplingTime - this.burnIn;

        double logPosteriorProba, logLikelihood, logPrior;
        double[] tretIPExpectations, tretIPNonPeakExpect, tretINPUTExpectations, tretINPUTNonPeakExpect,
                 ctrlIPExpectations, ctrlIPNonPeakExpect, ctrlINPUTExpectations, ctrlINPUTNonPeakExpect;
        double[][] tretReads, ctrlReads;
        for (int t=0; t<keepTime; t++) {
            double tempTretMeth, tempCtrlMeth, tempNonspecificEnrich,
                   tempTretIPOverdispersion, tempTretINPUTOverdispersion,
                   tempCtrlIPOverdispersion, tempCtrlINPUTOverdispersion;

            tempTretMeth = this.treatmentMethylationLevel[t];
            if (this.paramNum == 7)
                tempCtrlMeth = this.controlMethylationLevel[t];
            else
                tempCtrlMeth = 0;
            tempNonspecificEnrich = this.nonspecificEnrichment[t];
            tempTretIPOverdispersion = this.treatmentIPOverdispersion[t];
            tempTretINPUTOverdispersion = this.treatmentINPUTOverdispersion[t];
            tempCtrlIPOverdispersion = this.controlIPOverdispersion[t];
            tempCtrlINPUTOverdispersion = this.controlINPUTOverdispersion[t];

            tretReads = this.model.renewReadsExpectationViaNonspecificEnrich(true, new double[]{tempTretMeth}, new double[]{tempNonspecificEnrich},
                                                                             this.tretIPSizeFactors, this.tretINPUTSizeFactors, this.tretIPNonPeakSizeFactors, this.tretINPUTNonPeakSizeFactors);
            tretIPExpectations = tretReads[0];
            tretIPNonPeakExpect = tretReads[1];
            tretINPUTExpectations = tretReads[2];
            tretINPUTNonPeakExpect = tretReads[3];
            if (this.paramNum == 7)
                ctrlReads = this.model.renewReadsExpectationViaNonspecificEnrich(false, new double[]{tempCtrlMeth}, new double[]{tempNonspecificEnrich},
                                                                                 this.ctrlIPSizeFactors, this.ctrlINPUTSizeFactors, this.ctrlIPNonPeakSizeFactors, this.ctrlINPUTNonPeakSizeFactors);
            else
                ctrlReads = this.model.renewReadsExpectationViaNonspecificEnrich(false, new double[]{tempTretMeth}, new double[]{tempNonspecificEnrich},
                                                                                 this.ctrlIPSizeFactors, this.ctrlINPUTSizeFactors, this.ctrlIPNonPeakSizeFactors, this.ctrlINPUTNonPeakSizeFactors);
            ctrlIPExpectations = ctrlReads[0];
            ctrlIPNonPeakExpect = ctrlReads[1];
            ctrlINPUTExpectations = ctrlReads[2];
            ctrlINPUTNonPeakExpect = ctrlReads[3];

            logLikelihood = this.model.logLikelihood(tretIPExpectations, tretINPUTExpectations, ctrlIPExpectations, ctrlINPUTExpectations,
                                                     tretIPNonPeakExpect, tretINPUTNonPeakExpect, ctrlIPNonPeakExpect, ctrlINPUTNonPeakExpect,
                                                     tempTretIPOverdispersion, tempTretINPUTOverdispersion, tempCtrlIPOverdispersion, tempCtrlINPUTOverdispersion);
            logPrior = this.model.logPriority(tempTretIPOverdispersion, tempTretINPUTOverdispersion,
                                              tempCtrlIPOverdispersion, tempCtrlINPUTOverdispersion,
                                              tempTretMeth, tempCtrlMeth, tempNonspecificEnrich, 0,0,0,0);

            logPosteriorProba = logLikelihood + logPrior;
            if (logPosteriorProba - this.maximumLogPosterior > 0.00001)
                this.maximumLogPosterior = logPosteriorProba;
        }
    }

    public double getMaximumPosterior() {
        this.calcParamsProbability();
        return Math.exp(this.maximumLogPosterior);
    }
}
