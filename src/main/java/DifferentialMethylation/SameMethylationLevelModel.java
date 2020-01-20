package DifferentialMethylation;

import Quantification.*;

import java.util.concurrent.CountDownLatch;

/**
 * treatment and control group has same methylation level. For a single gene i, the model parameter vector is consisted of 5 components.
 *  theta ~ (treatmentINPUTOverdispersion_i, treatmentIPOverdispersion_i,
 *           controlINPUTOverdispersion_i, controlIPOverdispersion_i,
 *           methylationLevel_i)
 *
 * for gene i in an individual j, let D_j = {X_treatment_ip_i_j, X_treatment_input_i_j,
 *                                                    X_control_ip_i_j, X_control_input_i_j}
 * denote as the reads count in IP and INPUT of treatment and control group, respectively.
 *
 * suppose there exists totally m individuals, for gene i, the dataset can represent as D = {D_1, D_2, ..., D_m}
 * the posterior probability p(D|theta) approximates to
 *
 */
public class SameMethylationLevelModel extends ModelSelection {

    public SameMethylationLevelModel(MethylationLevelSampler tretMethylationLevelSampler, MethylationLevelSampler ctrlMethylationLevelSampler,
                                     OverdispersionSampler tretIPOverdispersionSampler, OverdispersionSampler tretINPUTOverdispersionSampler,
                                     OverdispersionSampler ctrlIPOverdispersionSampler, OverdispersionSampler ctrlINPUTOverdispersionSampler,
                                     BackgroundExpressionSampler[] treatmentBackgroundExpressionSampler, BackgroundExpressionSampler[] controlBackgroundExpressionSampler,
                                     ExpansionEffectSampler expansionEffectSampler, int samplingTime, int burnIn, int threadNumber, String tmpDir) {
        super(tretMethylationLevelSampler, ctrlMethylationLevelSampler,
              tretIPOverdispersionSampler, tretINPUTOverdispersionSampler, ctrlIPOverdispersionSampler, ctrlINPUTOverdispersionSampler,
              treatmentBackgroundExpressionSampler, controlBackgroundExpressionSampler, expansionEffectSampler, samplingTime, burnIn, threadNumber, tmpDir);
    }

    /**
     * MH sampling for treatment group methylation level, in these model we assume treatment and control group
     * have same methylation levels
     * @param time iteration time
     */
    @Override
    protected double treatmentMethylationSampling(int time, double logPrevParamsProba) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion, curExpansionEffect;
        // shape 1 × peakNumber
        double[] tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp, newTretMethLevel,
                 tretIPSizeFactors, tretINPUTSizeFactors, ctrlIPSizeFactors, ctrlINPUTSizeFactors;
        double logCurParamsProba, logLikeProba, logPrior;

        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion()[time];
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion()[time];
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion()[time];
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion()[time];
        curExpansionEffect = this.specificExpansionEffect? this.quantifiedExpansionEffect: this.parameters.getExpansionEffect()[time];
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        newTretMethLevel = new double[this.peakNumber];

        tretIPSizeFactors = this.parameters.getTretIPSizeFactors();
        tretINPUTSizeFactors = this.parameters.getTretINPUTSizeFactors();
        ctrlIPSizeFactors = this.parameters.getCtrlIPSizeFactors();
        ctrlINPUTSizeFactors = this.parameters.getCtrlINPUTSizeFactors();

        CountDownLatch countDownLatch = new CountDownLatch(this.peakNumber);
        MethLevelSampling taskMaker = ((prevValue, curValue, idx) -> () -> {
            try {
                int[] tretIPReadsCount = new int[this.tretIndividualNumber], tretINPUTReadsCount = new int[this.tretIndividualNumber],
                      ctrlIPReadsCount = new int[this.ctrlIndividualNumber], ctrlINPUTReadsCount = new int[this.ctrlIndividualNumber];
                double[] tretIPExpectation = new double[this.tretIndividualNumber], tretINPUTExpectation = new double[this.tretIndividualNumber],
                         ctrlIPExpectation = new double[this.ctrlIndividualNumber], ctrlINPUTExpectation = new double[this.ctrlIndividualNumber],
                         newTretIPExpectation, newTretINPUTExpectation, newCtrlIPExpectation, newCtrlINPUTExpectation;
                double tretExp = tretBkgExp[idx];
                double ctrlExp = ctrlBkgExp[idx];
                double prevLogLikelihood, curLogLikelihood, prevLogPrior, curLogPrior, prevLogPosterior, curLogPosterior;
                boolean samplingRes;

                for (int sampleIdx=0; sampleIdx<this.tretIndividualNumber; sampleIdx++) {
                    tretIPReadsCount[sampleIdx] = this.treatmentIPReads[sampleIdx][idx];
                    tretINPUTReadsCount[sampleIdx] = this.treatmentINPUTReads[sampleIdx][idx];
                    tretIPExpectation[sampleIdx] = this.treatmentIPExpectations[sampleIdx][idx];
                    tretINPUTExpectation[sampleIdx] = this.treatmentINPUTExpectations[sampleIdx][idx];
                }

                for (int sampleIdx=0; sampleIdx<this.ctrlIndividualNumber; sampleIdx++) {
                    ctrlIPReadsCount[sampleIdx] = this.controlIPReads[sampleIdx][idx];
                    ctrlINPUTReadsCount[sampleIdx] = this.controlINPUTReads[sampleIdx][idx];
                    ctrlIPExpectation[sampleIdx] = this.controlIPExpectations[sampleIdx][idx];
                    ctrlINPUTExpectation[sampleIdx] = this.controlINPUTExpectations[sampleIdx][idx];
                }

                newTretIPExpectation = this.renewReadsExpectation(tretIPSizeFactors, curExpansionEffect, curValue, tretExp, true);
                newTretINPUTExpectation = this.renewReadsExpectation(tretINPUTSizeFactors, 1, 1, tretExp, true);
                newCtrlIPExpectation = this.renewReadsExpectation(ctrlIPSizeFactors, curExpansionEffect, curValue, ctrlExp, false);
                newCtrlINPUTExpectation = this.renewReadsExpectation(ctrlINPUTSizeFactors, 1, 1, ctrlExp, false);

                prevLogLikelihood = this.singlePeakLogLikelihood(tretIPReadsCount, tretINPUTReadsCount, ctrlIPReadsCount, ctrlINPUTReadsCount,
                                                                 tretIPExpectation, tretINPUTExpectation, ctrlIPExpectation, ctrlINPUTExpectation,
                                                                 tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion);
                prevLogPrior = this.tretMethylationLevelSampler.getLogDensity(prevValue);
                prevLogPosterior = prevLogLikelihood + prevLogPrior;

                curLogLikelihood = this.singlePeakLogLikelihood(tretIPReadsCount, tretINPUTReadsCount, ctrlIPReadsCount, ctrlINPUTReadsCount,
                                                                newTretIPExpectation, newTretINPUTExpectation, newCtrlIPExpectation, newCtrlINPUTExpectation,
                                                                tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion);
                curLogPrior = this.tretMethylationLevelSampler.getLogDensity(curValue);
                curLogPosterior = curLogLikelihood + curLogPrior;

                samplingRes = this.mhSampling.getSamplingRes(curLogPosterior, prevLogPosterior, true);
                this.lock.lock();
                if (samplingRes) {
                    newTretMethLevel[idx] = curValue;
                    for (int i=0; i<this.tretIndividualNumber; i++) {
                        this.treatmentIPExpectations[i][idx] = newTretIPExpectation[i];
                        this.treatmentINPUTExpectations[i][idx] = newTretINPUTExpectation[i];
                    }
                    for (int i=0; i<this.ctrlIndividualNumber; i++) {
                        this.controlIPExpectations[i][idx] = newCtrlIPExpectation[i];
                        this.controlINPUTExpectations[i][idx] = newCtrlINPUTExpectation[i];
                    }
                }
                else
                    newTretMethLevel[idx] = prevValue;
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                countDownLatch.countDown();
                this.lock.unlock();
            }
        });

        double methValue, newMethValue;
        Runnable runnable;
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            methValue = tretMethLevel[peakIdx];
            newMethValue = this.tretMethylationLevelSampler.randomSample(methValue);
            runnable = taskMaker.createTask(methValue, newMethValue, peakIdx);
            this.executorService.submit(runnable);
        }

        try{
            countDownLatch.await();
        } catch (InterruptedException ie) {
            ie.printStackTrace();
            this.executorService.shutdown();
        }

        logLikeProba = this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion,
                                          ctrlIPOverdispersion, ctrlINPUTOverdispersion);
        logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                    ctrlIPOverdispersion, ctrlINPUTOverdispersion, curExpansionEffect,
                                    newTretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp);
        logCurParamsProba = logLikeProba + logPrior;

        this.parameters.setTretMethylation(newTretMethLevel);

        return logCurParamsProba;
    }

    /**
     * MH sampling for control group methylation level, in these model we assume treatment and control group
     * have same methylation levels
     * @param time iteration time
     */
    @Override
    protected double controlMethylationSampling(int time, double logPrevParamsProba) {
        this.parameters.setCtrlMethylation(this.parameters.getTretMethylation());

        return logPrevParamsProba;
    }

    /**
     * parameters logarithm prior probability
     *      log-gamma(tret_ip_overdispersion) + log-gamma(tret_input_overdispersion) +
     *      log-gamma(ctrl_ip_overdispersion) + log-gamma(ctrl_input_overdispersion) +
     *      log-beta(methylation_level)
     * tret_methylation_level = ctrl_methylation_level
     * @param treatmentIPOverdispersion treatment group IP reads count overdispersion of a gene
     * @param treatmentINPUTOverdispersion treatment group INPUT reads count overdispersion of a gene
     * @param controlIPOverdispersion control group IP reads count overdispersion of a gene
     * @param controlINPUTOverdispersion control group INPUT reads count overdispersion of a gene
     * @param treatmentMethylationLevel methylation level of a treatment group gene
     * @param controlMethylationLevel methylation level of a control group gene
     * @param treatmentBackgroundExpression background expression of a gene in treatment group
     * @param controlBackgroundExpression background expression of a gene in control group
     * @return logarithm prior probability
     */
    @Override
    public double logPriority(double treatmentIPOverdispersion, double treatmentINPUTOverdispersion,
                              double controlIPOverdispersion, double controlINPUTOverdispersion,
                              double expansionEffect, double[] treatmentMethylationLevel, double[] controlMethylationLevel,
                              double[] treatmentBackgroundExpression, double[] controlBackgroundExpression) {
        double proba = 0;
        proba += this.tretIPOverdispersionSampler.getLogDensity(treatmentIPOverdispersion);
        proba += this.tretINPUTOverdispersionSampler.getLogDensity(treatmentINPUTOverdispersion);
        proba += this.ctrlIPOverdispersionSampler.getLogDensity(controlIPOverdispersion);
        proba += this.ctrlINPUTOverdispersionSampler.getLogDensity(controlINPUTOverdispersion);
        proba += this.expansionEffectSampler.getLogDensity(expansionEffect);

        double tretMeth, tretExp, ctrlExp;
        BackgroundExpressionSampler tretSampler, ctrlSampler;
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            tretMeth = treatmentMethylationLevel[peakIdx];
            tretExp = treatmentBackgroundExpression[peakIdx];
            ctrlExp = controlBackgroundExpression[peakIdx];
            tretSampler = this.treatmentBackgroundExpressionSamplers[peakIdx];
            ctrlSampler = this.controlBackgroundExpressionSamplers[peakIdx];

            proba += this.tretMethylationLevelSampler.getLogDensity(tretMeth);
            proba += tretSampler.getLogDensity(tretExp);
            proba += ctrlSampler.getLogDensity(ctrlExp);
        }

        return proba;
    }
}
