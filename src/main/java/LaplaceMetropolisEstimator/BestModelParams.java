package LaplaceMetropolisEstimator;

import DifferentialMethylation.ModelSelection;
import Quantification.BackgroundExpressionSampler;
import Quantification.MethylationLevelSampler;

/**
 * parameters which lead to the maximum logarithm posterior probability
 */
public class BestModelParams {
    private double[] treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion,
                     treatmentMethylationLevel, controlMethylationLevel, treatmentExpression, controlExpression;
    private double[] tretIPSizeFactors, tretINPUTSizeFactors, ctrlIPSizeFactors, ctrlINPUTSizeFactors;
    private double maximumLogPosterior = -1 * Double.MAX_VALUE, expansionEffect;
    private int samplingTime, burnIn, paramNum, peakIdx;
    private ModelSelection model;

    /**
     * Constructor
     * @param treatmentIPOverdispersion treatment IP overdispersion, shape 1 × samplingTime
     * @param treatmentINPUTOverdispersion treatment INPUT overdispersion, shape 1 × samplingTime
     * @param controlIPOverdispersion control IP overdispersion, shape 1 × samplingTime
     * @param controlINPUTOverdispersion control INPUT overdispersion, shape 1 × samplingTime
     * @param treatmentMethylationLevel treatment IP methylation level, shape 1 × samplingTime
     * @param controlMethylationLevel control IP methylation level, shape 1 × samplingTime
     * @param treatmentExpression treatment IP methylation level, shape 1 × samplingTime
     * @param controlExpression control IP methylation level, shape 1 × samplingTime
     * @param model model
     * @param expansionEffect quantified antigen expansion effect
     * @param burnIn burn-in time
     */
    public BestModelParams(double[] treatmentIPOverdispersion, double[] treatmentINPUTOverdispersion,
                           double[] controlIPOverdispersion, double[] controlINPUTOverdispersion,
                           double[] treatmentMethylationLevel, double[] controlMethylationLevel,
                           double[] treatmentExpression, double[] controlExpression,
                           ModelSelection model, double expansionEffect, int burnIn, int peakIdx) {
        this.treatmentIPOverdispersion = treatmentIPOverdispersion;
        this.treatmentINPUTOverdispersion = treatmentINPUTOverdispersion;
        this.controlIPOverdispersion = controlIPOverdispersion;
        this.controlINPUTOverdispersion = controlINPUTOverdispersion;
        this.treatmentMethylationLevel = treatmentMethylationLevel;
        this.controlMethylationLevel = controlMethylationLevel;
        this.treatmentExpression = treatmentExpression;
        this.controlExpression = controlExpression;
        this.samplingTime = treatmentIPOverdispersion.length;
        this.expansionEffect = expansionEffect;
        this.burnIn = burnIn;
        this.peakIdx = peakIdx;
        if (controlMethylationLevel == null)
            this.paramNum = 7;      // same methylation level model parameter number
        else
            this.paramNum = 8;      // diff methylation level model parameter number
        this.model = model;
    }

    /**
     * set sample size factors
     */
    public void setSizeFactors(double[] tretIPSizeFactors, double[] tretINPUTSizeFactors, double[] ctrlIPSizeFactors, double[] ctrlINPUTSizeFactors) {
        this.tretIPSizeFactors = tretIPSizeFactors;
        this.tretINPUTSizeFactors = tretINPUTSizeFactors;
        this.ctrlIPSizeFactors = ctrlIPSizeFactors;
        this.ctrlINPUTSizeFactors = ctrlINPUTSizeFactors;
    }

    /**
     * calculate parameter posterior probability for each peak. log-posterior = log-likelihood + log-prior
     */
    private void calcParamsProbability() {
        int keepTime = this.samplingTime - this.burnIn,
            tretIndividualNumber = this.model.getTretIndividualNumber(), ctrlIndividualNumber = this.model.getCtrlIndividualNumber();
        int[] tretIPReads = new int[tretIndividualNumber], tretINPUTReads = new int[tretIndividualNumber],
              ctrlIPReads = new int[ctrlIndividualNumber], ctrlINPUTReads = new int[ctrlIndividualNumber];
        int[][] modelTretIPReads = this.model.getTreatmentIPReads(), modelTretINPUTReads = this.model.getTreatmentINPUTReads(),
                modelCtrlIPReads = this.model.getControlIPReads(), modelCtrlINPUTReads = this.model.getControlINPUTReads();

        for (int sampleIdx=0; sampleIdx<tretIndividualNumber; sampleIdx++) {
            tretIPReads[sampleIdx] = modelTretIPReads[sampleIdx][this.peakIdx];
            tretINPUTReads[sampleIdx] = modelTretINPUTReads[sampleIdx][this.peakIdx];
        }
        for (int sampleIdx=0; sampleIdx<ctrlIndividualNumber; sampleIdx++) {
            ctrlIPReads[sampleIdx] = modelCtrlIPReads[sampleIdx][this.peakIdx];
            ctrlINPUTReads[sampleIdx] = modelCtrlINPUTReads[sampleIdx][this.peakIdx];
        }
        BackgroundExpressionSampler tretExpSampler = this.model.getTretBackgroundExpressionSampler(this.peakIdx),
                                    ctrlExpSampler = this.model.getCtrlBackgroundExpressionSampler(this.peakIdx);
        MethylationLevelSampler tretMethSampler = this.model.getTretMethylationLevelSampler(),
                                ctrlMethSampler = this.model.getCtrlMethylationLevelSampler();

        double logPosteriorProba, logLikelihood, logPrior;
        double[] tretIPExpectations, tretINPUTExpectations, ctrlIPExpectations, ctrlINPUTExpectations;
        for (int t=0; t<keepTime; t++) {
            double tempTretMeth, tempCtrlMeth, tempTretExp, tempCtrlExp,
                   tempTretIPOverdispersion, tempTretINPUTOverdispersion,
                   tempCtrlIPOverdispersion, tempCtrlINPUTOverdispersion;

            tempTretMeth = this.treatmentMethylationLevel[t];
            if (this.paramNum == 8)
                tempCtrlMeth = this.controlMethylationLevel[t];
            else
                tempCtrlMeth = 0;
            tempTretExp = this.treatmentExpression[t];
            tempCtrlExp = this.controlExpression[t];
            tempTretIPOverdispersion = this.treatmentIPOverdispersion[t];
            tempTretINPUTOverdispersion = this.treatmentINPUTOverdispersion[t];
            tempCtrlIPOverdispersion = this.controlIPOverdispersion[t];
            tempCtrlINPUTOverdispersion = this.controlINPUTOverdispersion[t];

            tretIPExpectations = this.model.renewReadsExpectation(this.tretIPSizeFactors, this.expansionEffect, tempTretMeth, tempTretExp, true);
            tretINPUTExpectations = this.model.renewReadsExpectation(this.tretINPUTSizeFactors, this.expansionEffect, tempTretMeth, tempTretExp, true);
            if (this.paramNum == 8) {
                ctrlIPExpectations = this.model.renewReadsExpectation(this.ctrlIPSizeFactors, this.expansionEffect, tempCtrlMeth, tempCtrlExp, false);
                ctrlINPUTExpectations = this.model.renewReadsExpectation(this.ctrlINPUTSizeFactors, this.expansionEffect, tempCtrlMeth, tempCtrlExp, false);
            } else {
                ctrlIPExpectations = this.model.renewReadsExpectation(this.ctrlIPSizeFactors, this.expansionEffect, tempTretMeth, tempCtrlExp, false);
                ctrlINPUTExpectations = this.model.renewReadsExpectation(this.ctrlINPUTSizeFactors, this.expansionEffect, tempTretMeth, tempCtrlExp, false);
            }

            logLikelihood = this.model.singlePeakLogLikelihood(tretIPReads, tretINPUTReads, ctrlIPReads, ctrlINPUTReads,
                                                               tretIPExpectations, tretINPUTExpectations, ctrlIPExpectations, ctrlINPUTExpectations,
                                                               tempTretIPOverdispersion, tempTretINPUTOverdispersion, tempCtrlIPOverdispersion, tempCtrlINPUTOverdispersion);

            logPrior = tretExpSampler.getLogDensity(tempTretExp) + ctrlExpSampler.getLogDensity(tempCtrlExp) +
                       tretMethSampler.getLogDensity(tempTretMeth) + ctrlMethSampler.getLogDensity(tempCtrlMeth);

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
