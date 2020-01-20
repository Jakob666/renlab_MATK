package DifferentialMethylation;

import ProbabilityCalculation.ProbabilityCalculator;
import Quantification.*;

import java.io.*;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.locks.ReentrantLock;

public abstract class ModelSelection {
    protected int tretIndividualNumber, ctrlIndividualNumber, peakNumber, samplingTime, burnIn, threadNumber;
    protected double quantifiedExpansionEffect;
    protected boolean specificExpansionEffect = false;
    private String tmpDir;
    protected Parameters parameters;
    protected OverdispersionSampler tretIPOverdispersionSampler, tretINPUTOverdispersionSampler, ctrlIPOverdispersionSampler, ctrlINPUTOverdispersionSampler;
    protected MethylationLevelSampler tretMethylationLevelSampler, ctrlMethylationLevelSampler;
    protected BackgroundExpressionSampler[] treatmentBackgroundExpressionSamplers, controlBackgroundExpressionSamplers;
    protected ExpansionEffectSampler expansionEffectSampler;
    // shape individualNumber × peakNumber
    protected int[][] treatmentIPReads, treatmentINPUTReads, controlIPReads, controlINPUTReads;
    protected double[][] treatmentIPExpectations, treatmentINPUTExpectations, controlIPExpectations, controlINPUTExpectations;
    protected MHSampling mhSampling = new MHSampling();
    protected ExecutorService executorService;
    protected ReentrantLock lock;
    private DecimalFormat df = new DecimalFormat("0.000000");

    /**
     * Construct with default prior probabilities of same methylation level model and variant methylation level model
     * @param tretMethylationLevelSampler methylation level sampler follows Beta distribution
     * @param ctrlMethylationLevelSampler methylation level sampler follows Beta distribution
     * @param tretIPOverdispersionSampler IP reads count dispersion sampler follows Gamma distribution
     * @param tretINPUTOverdispersionSampler INPUT reads count dispersion sampler follows Gamma distribution
     * @param ctrlIPOverdispersionSampler IP reads count dispersion sampler follows Gamma distribution
     * @param ctrlINPUTOverdispersionSampler INPUT reads count dispersion sampler follows Gamma distribution
     * @param treatmentBackgroundExpressionSampler treatment background expression sampler follows log-normal distribution
     * @param controlBackgroundExpressionSampler control background expression sampler follows log-normal distribution
     */
    public ModelSelection(MethylationLevelSampler tretMethylationLevelSampler, MethylationLevelSampler ctrlMethylationLevelSampler,
                          OverdispersionSampler tretIPOverdispersionSampler, OverdispersionSampler tretINPUTOverdispersionSampler,
                          OverdispersionSampler ctrlIPOverdispersionSampler, OverdispersionSampler ctrlINPUTOverdispersionSampler,
                          BackgroundExpressionSampler[] treatmentBackgroundExpressionSampler, BackgroundExpressionSampler[] controlBackgroundExpressionSampler,
                          ExpansionEffectSampler expansionEffectSampler, int samplingTime, int burnIn, int threadNumber, String tmpDir) {
        this.tretMethylationLevelSampler = tretMethylationLevelSampler;
        this.ctrlMethylationLevelSampler = ctrlMethylationLevelSampler;
        this.tretIPOverdispersionSampler = tretIPOverdispersionSampler;
        this.tretINPUTOverdispersionSampler = tretINPUTOverdispersionSampler;
        this.ctrlIPOverdispersionSampler = ctrlIPOverdispersionSampler;
        this.ctrlINPUTOverdispersionSampler = ctrlINPUTOverdispersionSampler;
        this.controlBackgroundExpressionSamplers = controlBackgroundExpressionSampler;
        this.treatmentBackgroundExpressionSamplers = treatmentBackgroundExpressionSampler;
        this.expansionEffectSampler = expansionEffectSampler;
        this.threadNumber = threadNumber;
        this.tmpDir = tmpDir;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.executorService = Executors.newFixedThreadPool(this.threadNumber);
        this.lock = new ReentrantLock();
    }

    /**
     * set sample size factors
     */
    public void setSizeFactors(double[] tretIPSizeFactors, double[] tretINPUTSizeFactors,
                               double[] ctrlIPSizeFactors, double[] ctrlINPUTSizeFactors) {
        this.parameters.setTretIPSizeFactors(tretIPSizeFactors);
        this.parameters.setTretINPUTSizeFactors(tretINPUTSizeFactors);
        this.parameters.setCtrlIPSizeFactors(ctrlIPSizeFactors);
        this.parameters.setCtrlINPUTSizeFactors(ctrlINPUTSizeFactors);
    }

    /**
     * set sampling result list for a single gene
     * @param treatmentIPOverdispersion treatment IP overdispersion sampling list, shape 1 × samplingTime
     * @param treatmentINPUTOverdispersion treatment INPUT overdispersion sampling list, shape 1 × samplingTime
     * @param controlIPOverdispersion control IP overdispersion sampling list, shape 1 × samplingTime
     * @param controlINPUTOverdispersion control INPUT overdispersion sampling list, shape 1 × samplingTime
     * @param treatmentBkgExp treatment background expression sampling list, shape 1 × peakNumber
     * @param controlBkgExp control background expression sampling list, shape 1 × peakNumber
     * @param treatmentMethylationLevel treatment methylation level sampling list, shape 1 × peakNumber
     * @param controlMethylationLevel control methylation level sampling list, shape 1 × peakNumber
     */
    public void setSamplingList(double[] treatmentIPOverdispersion, double[] treatmentINPUTOverdispersion, double[] controlIPOverdispersion, double[] controlINPUTOverdispersion,
                                double[] expansionEffect, double[] treatmentBkgExp, double[] controlBkgExp, double[] treatmentMethylationLevel, double[] controlMethylationLevel) {
        this.parameters = new Parameters();
        this.parameters.setTretINPUTOverdispersion(treatmentINPUTOverdispersion);
        this.parameters.setTretIPOverdispersion(treatmentIPOverdispersion);
        this.parameters.setCtrlINPUTOverdispersion(controlINPUTOverdispersion);
        this.parameters.setCtrlIPOverdispersion(controlIPOverdispersion);
        this.parameters.setExpansionEffect(expansionEffect);
        this.parameters.setTretMethylation(treatmentMethylationLevel);
        this.parameters.setCtrlMethylation(controlMethylationLevel);
        this.parameters.setTretBkgExp(treatmentBkgExp);
        this.parameters.setCtrlBkgExp(controlBkgExp);
    }

    /**
     * treatment and control group reads count of a single gene
     * @param treatmentIPReads treatment IP gene reads count, shape individualNumber × peakNumber
     * @param treatmentINPUTReads treatment INPUT gene reads count, shape individualNumber × peakNumber
     * @param controlIPReads control IP gene reads count, shape individualNumber × peakNumber
     * @param controlINPUTReads control INPUT gene reads count, shape individualNumber × peakNumber
     */
    public void setTreatmentControlGeneReads(int[][] treatmentIPReads, int[][] treatmentINPUTReads,
                                             int[][] controlIPReads, int[][] controlINPUTReads) {
        this.tretIndividualNumber = treatmentIPReads.length;
        this.ctrlIndividualNumber = controlIPReads.length;
        this.peakNumber = treatmentIPReads[0].length;
        this.treatmentIPReads = treatmentIPReads;
        this.treatmentINPUTReads = treatmentINPUTReads;
        this.controlIPReads = controlIPReads;
        this.controlINPUTReads = controlINPUTReads;

        double[] tretIPSizeFactors = this.parameters.getTretIPSizeFactors();
        double[] tretINPUTSizeFactors = this.parameters.getTretINPUTSizeFactors();
        double[] ctrlIPSizeFactors = this.parameters.getCtrlIPSizeFactors();
        double[] ctrlINPUTSizeFactors = this.parameters.getCtrlINPUTSizeFactors();
        double[] tretBkgExpressions = this.parameters.getTretBkgExp();
        double[] ctrlBkgExpressions = this.parameters.getCtrlBkgExp();
        double[] tretMethylationLevels = this.parameters.getTretMethylation();
        double[] ctrlMethylationLevels = this.parameters.getCtrlMethylation();
        double expansion = this.parameters.getExpansionEffect()[0];

        // initial reads expectations
        this.treatmentIPExpectations = this.calcReadsExpectation(tretBkgExpressions, tretMethylationLevels, expansion, tretIPSizeFactors, false, true);
        this.treatmentINPUTExpectations = this.calcReadsExpectation(tretBkgExpressions, tretMethylationLevels, expansion, tretINPUTSizeFactors, true, true);
        this.controlIPExpectations = this.calcReadsExpectation(ctrlBkgExpressions, ctrlMethylationLevels, expansion, ctrlIPSizeFactors, false, false);
        this.controlINPUTExpectations = this.calcReadsExpectation(ctrlBkgExpressions, ctrlMethylationLevels, expansion, ctrlINPUTSizeFactors, true, false);
    }

    /**
     * calculate reads expectations for each peak in each sample
     * @param expressions expression values, shape 1 × peakNumber
     * @param methylationLevels methylation values, shape 1 × peakNumber
     * @param expansion antigen expansion effect
     * @param sampleSizeFactor sample size factors, shape 1 × individualNumber
     * @param input true, if calculate INPUT data reads expectation; otherwise, false
     * @return reads expectations, shape individualNumber × peakNumber
     */
    private double[][] calcReadsExpectation(double[] expressions, double[] methylationLevels, double expansion, double[] sampleSizeFactor, boolean input, boolean treatment) {
        double factor, expression, methylationLevel, readsExpectation;
        int individualNumber = treatment? this.tretIndividualNumber: this.ctrlIndividualNumber;
        double[][] expectations = new double[individualNumber][this.peakNumber];

        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            expression = expressions[peakIdx];
            methylationLevel = methylationLevels[peakIdx];
            for (int sampleIdx=0; sampleIdx<individualNumber; sampleIdx++) {
                factor = sampleSizeFactor[sampleIdx];
                if (input)
                    readsExpectation = this.calcExpectation(factor, expression, methylationLevel, expansion);
                else
                    readsExpectation = this.calcExpectation(factor, expression, 1, 1);
                expectations[sampleIdx][peakIdx] = readsExpectation;
            }
        }

        return expectations;
    }

    /**
     * calculate IP and INPUT reads count expectation
     * @param sizeFactor size factor
     * @param expression background expression
     * @param methylationLevel methylation level
     * @param expansion expansion effect of antigen
     * @return reads count expectation
     */
    private double calcExpectation(double sizeFactor, double expression, double methylationLevel, double expansion) {
        return sizeFactor * expansion * expression * methylationLevel;
    }

    /**
     * set quantified expansion effect
     * @param quantifiedExpansionEffect quantified expansion effect
     */
    public void setQuantifiedExpansionEffect(double quantifiedExpansionEffect) {
        this.quantifiedExpansionEffect = quantifiedExpansionEffect;
        this.specificExpansionEffect = true;
    }

    /**
     * an sampling iteration
     * @param time iteration time
     */
    public double iterate(int time, double logPrevParamsProba) {
        logPrevParamsProba = this.treatmentINPUTOverdispersionSampling(time, logPrevParamsProba);
        logPrevParamsProba = this.treatmentIPOverdispersionSampling(time, logPrevParamsProba);
        logPrevParamsProba = this.controlINPUTOverdispersionSampling(time, logPrevParamsProba);
        logPrevParamsProba = this.controlIPOverdispersionSampling(time, logPrevParamsProba);
        if (this.specificExpansionEffect)
            logPrevParamsProba = this.expansionEffectSampling(time, logPrevParamsProba);
        logPrevParamsProba = this.treatmentBkgExpSampling(time, logPrevParamsProba);
        logPrevParamsProba = this.controlBkgExpSampling(time, logPrevParamsProba);
        logPrevParamsProba = this.treatmentMethylationSampling(time, logPrevParamsProba);
        logPrevParamsProba = this.controlMethylationSampling(time, logPrevParamsProba);

        return logPrevParamsProba;
    }

    /**
     * MH sampling for treatment INPUT overdispersion
     * @param time iteration time
     * @param logPrevParamsProba log-likehood + log-prior
     */
    private double treatmentINPUTOverdispersionSampling(int time, double logPrevParamsProba) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion, curExpansionEffect,
               newTretINPUTOverdispersion;
        // shape 1 × peakNumber
        double[] tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion()[time-1];
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion()[time-1];
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion()[time-1];
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion()[time-1];
        curExpansionEffect = this.specificExpansionEffect? this.quantifiedExpansionEffect: this.parameters.getExpansionEffect()[time-1];
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        newTretINPUTOverdispersion = this.tretINPUTOverdispersionSampler.randomSample(tretINPUTOverdispersion);

        logLikeProba = this.logLikelihood(tretIPOverdispersion, newTretINPUTOverdispersion,
                                          ctrlIPOverdispersion, ctrlINPUTOverdispersion);
        logPrior = this.logPriority(tretIPOverdispersion, newTretINPUTOverdispersion,
                                    ctrlIPOverdispersion, ctrlINPUTOverdispersion, curExpansionEffect,
                                    tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp);
        logCurParamsProba = logLikeProba + logPrior;

        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        double[] tretINPUTOverdispersions = this.parameters.getTretINPUTOverdispersion();
        if (samplingRes) {
            tretINPUTOverdispersions[time] = newTretINPUTOverdispersion;
            this.parameters.setTretINPUTOverdispersion(tretINPUTOverdispersions);
            return logCurParamsProba;
        }
        else {
            tretINPUTOverdispersions[time] = tretINPUTOverdispersion;
            this.parameters.setTretINPUTOverdispersion(tretINPUTOverdispersions);
        }

        return logPrevParamsProba;
    }

    /**
     * MH sampling for treatment IP overdispersion
     * @param time iteration time
     */
    private double treatmentIPOverdispersionSampling(int time, double logPrevParamsProba) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion, curExpansionEffect,
                newTretIPOverdispersion;
        // shape 1 × peakNumber
        double[] tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion()[time];
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion()[time-1];
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion()[time-1];
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion()[time-1];
        curExpansionEffect = this.specificExpansionEffect? this.quantifiedExpansionEffect: this.parameters.getExpansionEffect()[time-1];
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        newTretIPOverdispersion = this.tretIPOverdispersionSampler.randomSample(tretIPOverdispersion);

        logLikeProba = this.logLikelihood(newTretIPOverdispersion, tretINPUTOverdispersion,
                                          ctrlIPOverdispersion, ctrlINPUTOverdispersion);
        logPrior = this.logPriority(newTretIPOverdispersion, tretINPUTOverdispersion,
                                    ctrlIPOverdispersion, ctrlINPUTOverdispersion, curExpansionEffect,
                                    tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp);
        logCurParamsProba = logLikeProba + logPrior;

        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        double[] tretIPOverdispersions = this.parameters.getTretIPOverdispersion();
        if (samplingRes) {
            tretIPOverdispersions[time] = newTretIPOverdispersion;
            this.parameters.setTretIPOverdispersion(tretIPOverdispersions);
            return logCurParamsProba;
        } else {
            tretIPOverdispersions[time] = tretIPOverdispersion;
            this.parameters.setTretIPOverdispersion(tretIPOverdispersions);
        }

        return logPrevParamsProba;
    }

    /**
     * MH sampling for control INPUT overdispersion
     * @param time iteration time
     */
    private double controlINPUTOverdispersionSampling(int time, double logPrevParamsProba){
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion, curExpansionEffect,
                newCtrlINPUTOverdispersion;
        // shape 1 × peakNumber
        double[] tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion()[time];
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion()[time];
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion()[time-1];
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion()[time-1];
        curExpansionEffect = this.specificExpansionEffect? this.quantifiedExpansionEffect: this.parameters.getExpansionEffect()[time-1];
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        newCtrlINPUTOverdispersion = this.ctrlINPUTOverdispersionSampler.randomSample(ctrlINPUTOverdispersion);
        logLikeProba = this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion,
                                          ctrlIPOverdispersion, newCtrlINPUTOverdispersion);
        logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                    ctrlIPOverdispersion, newCtrlINPUTOverdispersion, curExpansionEffect,
                                    tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp);
        logCurParamsProba = logLikeProba + logPrior;

        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        double[] ctrlINPUTOverdispersions = this.parameters.getCtrlINPUTOverdispersion();
        if (samplingRes) {
            ctrlINPUTOverdispersions[time] = newCtrlINPUTOverdispersion;
            this.parameters.setCtrlINPUTOverdispersion(ctrlINPUTOverdispersions);
            return logCurParamsProba;
        } else {
            ctrlINPUTOverdispersions[time] = ctrlINPUTOverdispersion;
            this.parameters.setCtrlINPUTOverdispersion(ctrlINPUTOverdispersions);
        }

        return logPrevParamsProba;
    }

    /**
     * MH sampling for control IP overdispersion
     * @param time iteration time
     */
    private double controlIPOverdispersionSampling(int time, double logPrevParamsProba){
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion, curExpansionEffect,
                newCtrlIPOverdispersion;
        // shape 1 × peakNumber
        double[] tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion()[time];
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion()[time];
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion()[time];
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion()[time-1];
        curExpansionEffect = this.specificExpansionEffect? this.quantifiedExpansionEffect: this.parameters.getExpansionEffect()[time-1];
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        newCtrlIPOverdispersion = this.ctrlIPOverdispersionSampler.randomSample(ctrlIPOverdispersion);
        logLikeProba = this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion,
                                          newCtrlIPOverdispersion, ctrlINPUTOverdispersion);
        logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                    newCtrlIPOverdispersion, ctrlINPUTOverdispersion, curExpansionEffect,
                                    tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp);
        logCurParamsProba = logLikeProba + logPrior;

        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        double[] ctrlIPOverdispersions = this.parameters.getCtrlIPOverdispersion();
        if (samplingRes) {
            ctrlIPOverdispersions[time] = newCtrlIPOverdispersion;
            this.parameters.setCtrlIPOverdispersion(ctrlIPOverdispersions);
            return logCurParamsProba;
        } else {
            ctrlIPOverdispersions[time] = ctrlIPOverdispersion;
            this.parameters.setCtrlIPOverdispersion(ctrlIPOverdispersions);
        }

        return logPrevParamsProba;
    }

    /**
     * MH sampling for control IP overdispersion
     * @param time iteration time
     */
    private double expansionEffectSampling(int time, double logPrevParamsProba){
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion, curExpansionEffect,
                newExpansionEffect;
        // shape 1 × peakNumber
        double[] tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion()[time];
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion()[time];
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion()[time];
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion()[time];
        curExpansionEffect = this.parameters.getExpansionEffect()[time-1];
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        newExpansionEffect = this.ctrlIPOverdispersionSampler.randomSample(curExpansionEffect);

        logLikeProba = this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion,
                                          ctrlIPOverdispersion, ctrlINPUTOverdispersion);
        logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                    ctrlIPOverdispersion, ctrlINPUTOverdispersion, newExpansionEffect,
                                    tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp);
        logCurParamsProba = logLikeProba + logPrior;

        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        double[] expansionEffects = this.parameters.getExpansionEffect();
        if (samplingRes) {
            expansionEffects[time] = newExpansionEffect;
            this.parameters.setExpansionEffect(expansionEffects);
            return logCurParamsProba;
        } else {
            expansionEffects[time] = curExpansionEffect;
            this.parameters.setCtrlIPOverdispersion(expansionEffects);
        }

        return logPrevParamsProba;
    }

    /**
     * MH sampling for treatment group methylation level
     * @param time iteration time
     */
    protected abstract double treatmentMethylationSampling(int time, double logPrevParamsProba);

    /**
     * MH sampling for control group methylation level
     * @param time iteration time
     */
    protected abstract double controlMethylationSampling(int time, double logPrevParamsProba);

    /**
     * MH sampling for treatment group background expression
     * @param time iteration time
     */
    private double treatmentBkgExpSampling(int time, double logPrevParamsProba) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion, curExpansionEffect;
        // shape 1 × peakNumber
        double[] tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp, newTretBkgExp,
                 tretIPSizeFactors, tretINPUTSizeFactors;
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
        newTretBkgExp = new double[this.peakNumber];

        tretIPSizeFactors = this.parameters.getTretIPSizeFactors();
        tretINPUTSizeFactors = this.parameters.getTretINPUTSizeFactors();

        CountDownLatch countDownLatch = new CountDownLatch(this.peakNumber);
        BkgExpSamplingTask taskMaker = ((prevValue, curValue, bkgSampler, idx) -> () -> {
            try {
                int[] tretIPReadsCount = new int[this.tretIndividualNumber], tretINPUTReadsCount = new int[this.tretIndividualNumber],
                        ctrlIPReadsCount = new int[this.ctrlIndividualNumber], ctrlINPUTReadsCount = new int[this.ctrlIndividualNumber];
                double[] tretIPExpectation = new double[this.tretIndividualNumber], tretINPUTExpectation = new double[this.tretIndividualNumber],
                        ctrlIPExpectation = new double[this.ctrlIndividualNumber], ctrlINPUTExpectation = new double[this.ctrlIndividualNumber],
                        newTretIPExpectation, newTretINPUTExpectation;
                double tretMethValue = tretMethLevel[idx];
                double prevLogLikelihood, curLogLikelihood, prevLogPrior, curLogPrior, prevLogPosterior, curLogPosterior;
                boolean samplingRes;

                for (int sampleIdx=0; sampleIdx<this.tretIndividualNumber; sampleIdx++) {
                    tretIPReadsCount[sampleIdx] = this.treatmentIPReads[sampleIdx][idx];
                    tretINPUTReadsCount[sampleIdx] = this.treatmentINPUTReads[sampleIdx][idx];
                    tretIPExpectation[sampleIdx] = this.treatmentIPExpectations[sampleIdx][idx];
                    tretINPUTExpectation[sampleIdx] = this.treatmentINPUTExpectations[sampleIdx][idx];
                }

                for (int sampleIdx=0; sampleIdx<this.ctrlIndividualNumber; sampleIdx++) {
                    ctrlIPReadsCount[sampleIdx] = this.treatmentIPReads[sampleIdx][idx];
                    ctrlINPUTReadsCount[sampleIdx] = this.treatmentINPUTReads[sampleIdx][idx];
                    ctrlIPExpectation[sampleIdx] = this.controlIPExpectations[sampleIdx][idx];
                    ctrlINPUTExpectation[sampleIdx] = this.controlINPUTExpectations[sampleIdx][idx];
                }

                newTretIPExpectation = this.renewReadsExpectation(tretIPSizeFactors, curExpansionEffect, tretMethValue, curValue, true);
                newTretINPUTExpectation = this.renewReadsExpectation(tretINPUTSizeFactors, 1, 1, curValue, true);

                prevLogLikelihood = this.singlePeakLogLikelihood(tretIPReadsCount, tretINPUTReadsCount, ctrlIPReadsCount, ctrlINPUTReadsCount,
                                                                 tretIPExpectation, tretINPUTExpectation, ctrlIPExpectation, ctrlINPUTExpectation,
                                                                 tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion);
                prevLogPrior = bkgSampler.getLogDensity(prevValue);
                prevLogPosterior = prevLogLikelihood + prevLogPrior;

                curLogLikelihood = this.singlePeakLogLikelihood(tretIPReadsCount, tretINPUTReadsCount, ctrlIPReadsCount, ctrlINPUTReadsCount,
                                                                newTretIPExpectation, newTretINPUTExpectation, ctrlIPExpectation, ctrlINPUTExpectation,
                                                                tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion);
                curLogPrior = bkgSampler.getLogDensity(curValue);
                curLogPosterior = curLogLikelihood + curLogPrior;

                samplingRes = this.mhSampling.getSamplingRes(curLogPosterior, prevLogPosterior, true);
                this.lock.lock();
                if (samplingRes) {
                    newTretBkgExp[idx] = curValue;
                    for (int i=0; i<this.tretIndividualNumber; i++) {
                        this.treatmentIPExpectations[i][idx] = newTretIPExpectation[i];
                        this.treatmentINPUTExpectations[i][idx] = newTretINPUTExpectation[i];
                    }
                }
                else
                    newTretBkgExp[idx] = prevValue;
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                countDownLatch.countDown();
                this.lock.unlock();
            }
        });

        double expression, newExpression;
        BackgroundExpressionSampler sampler;
        Runnable runnable;
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            expression = tretBkgExp[peakIdx];
            sampler = this.treatmentBackgroundExpressionSamplers[peakIdx];
            newExpression = sampler.randomSample(expression);
            runnable = taskMaker.createSamplingTask(expression, newExpression, sampler, peakIdx);
            this.executorService.submit(runnable);
        }
        // wait all peaks renew its background expression
        try {
            countDownLatch.await();
        } catch (InterruptedException ie) {
            ie.printStackTrace();
            this.executorService.shutdown();
        }

        logLikeProba = this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion,
                                          ctrlIPOverdispersion, ctrlINPUTOverdispersion);
        logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                    ctrlIPOverdispersion, ctrlINPUTOverdispersion, curExpansionEffect,
                                    tretMethLevel, ctrlMethLevel, newTretBkgExp, ctrlBkgExp);
        logCurParamsProba = logLikeProba + logPrior;

        this.parameters.setTretBkgExp(newTretBkgExp);

        return logCurParamsProba;
    }

    /**
     * MH sampling for control group background expression
     * @param time iteration time
     */
    private double controlBkgExpSampling(int time, double logPrevParamsProba) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion, curExpansionEffect;
        // shape 1 × peakNumber
        double[] tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp, newCtrlBkgExp,
                 ctrlIPSizeFactors, ctrlINPUTSizeFactors;
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
        newCtrlBkgExp = new double[this.peakNumber];

        ctrlIPSizeFactors = this.parameters.getCtrlIPSizeFactors();
        ctrlINPUTSizeFactors = this.parameters.getCtrlINPUTSizeFactors();

        CountDownLatch countDownLatch = new CountDownLatch(this.peakNumber);
        BkgExpSamplingTask taskMaker = ((prevValue, curValue, bkgSampler, idx) -> () -> {
            try {
                int[] tretIPReadsCount = new int[this.tretIndividualNumber], tretINPUTReadsCount = new int[this.tretIndividualNumber],
                        ctrlIPReadsCount = new int[this.ctrlIndividualNumber], ctrlINPUTReadsCount = new int[this.ctrlIndividualNumber];
                double[] tretIPExpectation = new double[this.tretIndividualNumber], tretINPUTExpectation = new double[this.tretIndividualNumber],
                        ctrlIPExpectation = new double[this.ctrlIndividualNumber], ctrlINPUTExpectation = new double[this.ctrlIndividualNumber],
                        newCtrlIPExpectation, newCtrlINPUTExpectation;
                double ctrlMethValue = ctrlMethLevel[idx];
                double prevLogLikelihood, curLogLikelihood, prevLogPrior, curLogPrior, prevLogPosterior, curLogPosterior;
                boolean samplingRes;

                for (int sampleIdx=0; sampleIdx<this.tretIndividualNumber; sampleIdx++) {
                    tretIPReadsCount[sampleIdx] = this.treatmentIPReads[sampleIdx][idx];
                    tretINPUTReadsCount[sampleIdx] = this.treatmentINPUTReads[sampleIdx][idx];
                    tretIPExpectation[sampleIdx] = this.treatmentIPExpectations[sampleIdx][idx];
                    tretINPUTExpectation[sampleIdx] = this.treatmentINPUTExpectations[sampleIdx][idx];
                }

                for (int sampleIdx=0; sampleIdx<this.ctrlIndividualNumber; sampleIdx++) {
                    ctrlIPReadsCount[sampleIdx] = this.treatmentIPReads[sampleIdx][idx];
                    ctrlINPUTReadsCount[sampleIdx] = this.treatmentINPUTReads[sampleIdx][idx];
                    ctrlIPExpectation[sampleIdx] = this.controlIPExpectations[sampleIdx][idx];
                    ctrlINPUTExpectation[sampleIdx] = this.controlINPUTExpectations[sampleIdx][idx];
                }

                newCtrlIPExpectation = this.renewReadsExpectation(ctrlIPSizeFactors, curExpansionEffect, ctrlMethValue, curValue, false);
                newCtrlINPUTExpectation = this.renewReadsExpectation(ctrlINPUTSizeFactors, 1, 1, curValue, false);

                prevLogLikelihood = this.singlePeakLogLikelihood(tretIPReadsCount, tretINPUTReadsCount, ctrlIPReadsCount, ctrlINPUTReadsCount,
                                                                 tretIPExpectation, tretINPUTExpectation, ctrlIPExpectation, ctrlINPUTExpectation,
                                                                 tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion);
                prevLogPrior = bkgSampler.getLogDensity(prevValue);
                prevLogPosterior = prevLogLikelihood + prevLogPrior;

                curLogLikelihood = this.singlePeakLogLikelihood(tretIPReadsCount, tretINPUTReadsCount, ctrlIPReadsCount, ctrlINPUTReadsCount,
                                                                tretIPExpectation, tretINPUTExpectation, newCtrlIPExpectation, newCtrlINPUTExpectation,
                                                                tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion);
                curLogPrior = bkgSampler.getLogDensity(curValue);
                curLogPosterior = curLogLikelihood + curLogPrior;

                samplingRes = this.mhSampling.getSamplingRes(curLogPosterior, prevLogPosterior, true);
                this.lock.lock();
                if (samplingRes) {
                    newCtrlBkgExp[idx] = curValue;
                    for (int i=0; i<this.ctrlIndividualNumber; i++) {
                        this.controlIPExpectations[i][idx] = newCtrlIPExpectation[i];
                        this.controlINPUTExpectations[i][idx] = newCtrlINPUTExpectation[i];
                    }
                }
                else
                    newCtrlBkgExp[idx] = prevValue;
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                countDownLatch.countDown();
                this.lock.unlock();
            }
        });

        double expression, newExpression;
        BackgroundExpressionSampler sampler;
        Runnable runnable;
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            expression = ctrlBkgExp[peakIdx];
            sampler = this.treatmentBackgroundExpressionSamplers[peakIdx];
            newExpression = sampler.randomSample(expression);
            runnable = taskMaker.createSamplingTask(expression, newExpression, sampler, peakIdx);
            this.executorService.submit(runnable);
        }
        // wait all peaks renew its background expression
        try {
            countDownLatch.await();
        } catch (InterruptedException ie) {
            ie.printStackTrace();
            this.executorService.shutdown();
        }

        logLikeProba = this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion,
                                          ctrlIPOverdispersion, ctrlINPUTOverdispersion);
        logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                    ctrlIPOverdispersion, ctrlINPUTOverdispersion, curExpansionEffect,
                                    tretMethLevel, ctrlMethLevel, tretBkgExp, newCtrlBkgExp);
        logCurParamsProba = logLikeProba + logPrior;

        this.parameters.setCtrlBkgExp(newCtrlBkgExp);

        return logCurParamsProba;
    }

    /**
     * calculate reads expectations
     * @param sampleSizeFactor samples' size factors, shape 1 × individualNumber
     * @param expansion antigen expansion effect
     * @param methLevel methylation level
     * @param bkgExp background expression
     * @return reads expectations, shape 1 × individualNumber
     */
    public double[] renewReadsExpectation(double[] sampleSizeFactor, double expansion, double methLevel,double bkgExp, boolean treatment) {
        double readsExpectation, factor;
        int individualNumber = treatment? this.tretIndividualNumber: this.ctrlIndividualNumber;
        double[] readsExpectations = new double[individualNumber];
        for (int sampleIdx=0; sampleIdx<individualNumber; sampleIdx++) {
            factor = sampleSizeFactor[sampleIdx];
            readsExpectation = this.calcReadsExpectation(factor, bkgExp, methLevel, expansion);
            readsExpectations[sampleIdx] = readsExpectation;
        }

        return readsExpectations;
    }

    /**
     * calculate IP and INPUT reads count expectation
     * @param sizeFactor size factor
     * @param expression background expression
     * @param methylationLevel methylation level
     * @param expansion expansion effect of antigen
     * @return reads count expectation
     */
    private double calcReadsExpectation(double sizeFactor, double expression, double methylationLevel, double expansion) {
        return sizeFactor * expansion * expression * methylationLevel;
    }

    /**
     * calculate log-likelihood  p(D|params) with given parameters,
     * where D denotes observed data (X_tret_ip_i_j, miu_tret_input_i_j, X_ctrl_ip_i_j, miu_ctrl_input_i_j)
     * params represents parameters (miu_tret_ip_i_j, miu_tret_input_i_j, miu_ctrl_ip_i_j, miu_ctrl_input_i_j, sigma_ip_i, sigma_input_i)
     *
     *      sum(log NB(X_tret_ip_i_j| miu_tret_ip_i_j, sigma_ip_i)) + sum(log NB(X_tret_input_i_j| miu_tret_input_i_j, sigma_input_i))
     *    + sum(log NB(X_ctrl_ip_i_j| miu_ctrl_ip_i_j, sigma_ip_i)) + sum(log NB(X_ctrl_input_i_j| miu_ctrl_input_i_j, sigma_input_i))
     * where NB denotes as negative binomial distribution. X, miu and sigma represent reads count, reads count expectation and overdispersion
     * of a gene, respectively
     * @param treatmentIPOverdispersion treatment group IP reads count overdispersion of a gene
     * @param treatmentINPUTOverdispersion treatment group INPUT reads count overdispersion of a gene
     * @param controlIPOverdispersion control group IP reads count overdispersion of a gene
     * @param controlINPUTOverdispersion control group INPUT reads count overdispersion of a gene
     */
    protected double logLikelihood(double treatmentIPOverdispersion, double treatmentINPUTOverdispersion,
                                   double controlIPOverdispersion, double controlINPUTOverdispersion) {
        double proba = 0;
        // shape 1 × individualNumber
        int[] tretIPReadsCount = new int[this.tretIndividualNumber], tretINPUTReadsCount = new int[this.tretIndividualNumber],
              ctrlIPReadsCount = new int[this.ctrlIndividualNumber], ctrlINPUTReadsCount = new int[this.ctrlIndividualNumber];
        double[] tretIPExpectation = new double[this.tretIndividualNumber], tretINPUTExpectation = new double[this.tretIndividualNumber],
                 ctrlIPExpectation = new double[this.ctrlIndividualNumber], ctrlINPUTExpectation = new double[this.ctrlIndividualNumber];

        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            for (int sampleIdx=0; sampleIdx<this.tretIndividualNumber; sampleIdx++) {
                tretIPReadsCount[sampleIdx] = this.treatmentIPReads[sampleIdx][peakIdx];
                tretINPUTReadsCount[sampleIdx] = this.treatmentINPUTReads[sampleIdx][peakIdx];
                tretIPExpectation[sampleIdx] = this.treatmentIPExpectations[sampleIdx][peakIdx];
                tretINPUTExpectation[sampleIdx] = this.treatmentINPUTExpectations[sampleIdx][peakIdx];
            }

            for (int sampleIdx=0; sampleIdx<this.ctrlIndividualNumber; sampleIdx++) {
                ctrlIPReadsCount[sampleIdx] = this.controlIPReads[sampleIdx][peakIdx];
                ctrlINPUTReadsCount[sampleIdx] = this.controlINPUTReads[sampleIdx][peakIdx];
                ctrlIPExpectation[sampleIdx] = this.controlIPExpectations[sampleIdx][peakIdx];
                ctrlINPUTExpectation[sampleIdx] = this.controlINPUTExpectations[sampleIdx][peakIdx];
            }
            proba += this.singlePeakLogLikelihood(tretIPReadsCount, tretINPUTReadsCount, ctrlIPReadsCount, ctrlINPUTReadsCount,
                                                  tretIPExpectation, tretINPUTExpectation, ctrlIPExpectation, ctrlINPUTExpectation,
                                                  treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion);
        }

        return proba;
    }

    /**
     * calculate log-likelihood  p(D|params) with given parameters and reads count,
     * where D denotes observed data (X_tret_ip_i_j, miu_tret_input_i_j, X_ctrl_ip_i_j, miu_ctrl_input_i_j)
     * params represents parameters (miu_tret_ip_i_j, miu_tret_input_i_j, miu_ctrl_ip_i_j, miu_ctrl_input_i_j, sigma_ip_i, sigma_input_i)
     *
     *      sum(log NB(X_tret_ip_i_j| miu_tret_ip_i_j, sigma_ip_i)) + sum(log NB(X_tret_input_i_j| miu_tret_input_i_j, sigma_input_i))
     *    + sum(log NB(X_ctrl_ip_i_j| miu_ctrl_ip_i_j, sigma_ip_i)) + sum(log NB(X_ctrl_input_i_j| miu_ctrl_input_i_j, sigma_input_i))
     * where NB denotes as negative binomial distribution. X, miu and sigma represent reads count, reads count expectation and overdispersion
     * of a gene, respectively
     * @param treatmentIPOverdispersion treatment group IP reads count overdispersion of a gene
     * @param treatmentINPUTOverdispersion treatment group INPUT reads count overdispersion of a gene
     * @param controlIPOverdispersion control group IP reads count overdispersion of a gene
     * @param controlINPUTOverdispersion control group INPUT reads count overdispersion of a gene
     */
    public double logLikelihood(double[][] treatmentIPExpectations, double[][] treatmentINPUTExpectations,
                                double[][] controlIPExpectations, double[][] controlINPUTExpectations,
                                double treatmentIPOverdispersion, double treatmentINPUTOverdispersion,
                                double controlIPOverdispersion, double controlINPUTOverdispersion) {
        double proba = 0;
        // shape 1 × individualNumber
        int[] tretIPReadsCount = new int[this.tretIndividualNumber], tretINPUTReadsCount = new int[this.tretIndividualNumber],
              ctrlIPReadsCount = new int[this.ctrlIndividualNumber], ctrlINPUTReadsCount = new int[this.ctrlIndividualNumber];
        double[] tretIPExpectation = new double[this.tretIndividualNumber], tretINPUTExpectation = new double[this.tretIndividualNumber],
                 ctrlIPExpectation = new double[this.ctrlIndividualNumber], ctrlINPUTExpectation = new double[this.ctrlIndividualNumber];
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            for (int sampleIdx=0; sampleIdx<this.tretIndividualNumber; sampleIdx++) {
                tretIPReadsCount[sampleIdx] = this.treatmentIPReads[sampleIdx][peakIdx];
                tretINPUTReadsCount[sampleIdx] = this.treatmentINPUTReads[sampleIdx][peakIdx];
                tretIPExpectation[sampleIdx] = treatmentIPExpectations[sampleIdx][peakIdx];
                tretINPUTExpectation[sampleIdx] = treatmentINPUTExpectations[sampleIdx][peakIdx];
            }

            for (int sampleIdx=0; sampleIdx<this.ctrlIndividualNumber; sampleIdx++) {
                ctrlIPReadsCount[sampleIdx] = this.controlIPReads[sampleIdx][peakIdx];
                ctrlINPUTReadsCount[sampleIdx] = this.controlINPUTReads[sampleIdx][peakIdx];
                ctrlIPExpectation[sampleIdx] = controlIPExpectations[sampleIdx][peakIdx];
                ctrlINPUTExpectation[sampleIdx] = controlINPUTExpectations[sampleIdx][peakIdx];
            }
            proba += this.singlePeakLogLikelihood(tretIPReadsCount, tretINPUTReadsCount, ctrlIPReadsCount, ctrlINPUTReadsCount,
                                                  tretIPExpectation, tretINPUTExpectation, ctrlIPExpectation, ctrlINPUTExpectation,
                                                  treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion);
        }

        return proba;
    }

    /**
     * log likelihood for a single peak region
     * @param tretIPReadsCount treatment IP reads count, shape 1 × individualNumber
     * @param tretINPUTReadsCount treatment INPUT reads count, shape 1 × individualNumber
     * @param ctrlIPReadsCount control IP reads count, shape 1 × individualNumber
     * @param ctrlINPUTReadsCount control INPUT reads count, shape 1 × individualNumber
     * @param tretIPExpectation treatment IP reads expectation, shape 1 × individualNumber
     * @param tretINPUTExpectation treatment INPUT reads expectation, shape 1 × individualNumber
     * @param ctrlIPExpectation control IP reads expectation, shape 1 × individualNumber
     * @param ctrlINPUTExpectation control INPUT reads expectation, shape 1 × individualNumber
     * @param treatmentIPOverdispersion treatment IP overdispersion
     * @param treatmentINPUTOverdispersion treatment INPUT overdispersion
     * @param controlIPOverdispersion control IP overdispersion
     * @param controlINPUTOverdispersion control INPUT overdispersion
     * @return log likelihood for a single peak
     */
    public double singlePeakLogLikelihood(int[] tretIPReadsCount, int[] tretINPUTReadsCount, int[] ctrlIPReadsCount, int[] ctrlINPUTReadsCount,
                                          double[] tretIPExpectation, double[] tretINPUTExpectation, double[] ctrlIPExpectation, double[] ctrlINPUTExpectation,
                                          double treatmentIPOverdispersion, double treatmentINPUTOverdispersion,  double controlIPOverdispersion, double controlINPUTOverdispersion) {
        double proba = 0;
        proba += ProbabilityCalculator.logNegativeProbability(tretIPReadsCount, tretIPExpectation, treatmentIPOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(tretINPUTReadsCount, tretINPUTExpectation, treatmentINPUTOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(ctrlIPReadsCount, ctrlIPExpectation, controlINPUTOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(ctrlINPUTReadsCount, ctrlINPUTExpectation, controlIPOverdispersion);

        return proba;
    }

    /**
     * parameters logarithm prior probability
     * @param treatmentIPOverdispersion treatment group IP reads count overdispersion of a gene
     * @param treatmentINPUTOverdispersion treatment group INPUT reads count overdispersion of a gene
     * @param controlIPOverdispersion control group IP reads count overdispersion of a gene
     * @param controlINPUTOverdispersion control group INPUT reads count overdispersion of a gene
     * @param treatmentMethylationLevel methylation level of a treatment group gene, shape 1 × peakNumber
     * @param controlMethylationLevel methylation level of a control group gene, shape 1 × peakNumber
     * @param treatmentBackgroundExpression background expression of a gene in treatment group, shape 1 × peakNumber
     * @param controlBackgroundExpression background expression of a gene in control group, shape 1 × peakNumber
     * @return logarithm prior probability
     */
    public abstract double logPriority(double treatmentIPOverdispersion, double treatmentINPUTOverdispersion,
                                       double controlIPOverdispersion, double controlINPUTOverdispersion,
                                       double expansionEffect, double[] treatmentMethylationLevel, double[] controlMethylationLevel,
                                       double[] treatmentBackgroundExpression, double[] controlBackgroundExpression);

    /**
     * calculate log-posterior = log-likelihood + log-prior at the end of an iteration
     * @return log-posterior
     */
    public double paramsLogFullConditionalProbability(int time) {
        // shape 1 × peakNumber
        double[] tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp;
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion, expansionEffect;

        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();

        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion()[time];
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion()[time];
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion()[time];
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion()[time];
        expansionEffect = this.parameters.getExpansionEffect()[time];

        double logLike= this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion);
        double prior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        expansionEffect, tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp);
        return logLike + prior;
    }

    /**
     * renew sampling files for each peak
     */
    public void renewRecordFiles() {
        // shape 1 × peakNumber
        double[] tretMethylationLevels = this.parameters.getTretMethylation();
        double[] ctrlMethylationLevels = this.parameters.getCtrlMethylation();
        double[] tretExpressions = this.parameters.getTretBkgExp();
        double[] ctrlExpressions = this.parameters.getCtrlBkgExp();

        CountDownLatch countDownLatch = new CountDownLatch(this.peakNumber);
        FileRecordTask taskMaker = ((tretMeth, ctrlMeth, tretExp, ctrlExp, idx) -> () -> {
            BufferedWriter bfw = null;
            File targetFile = new File(this.tmpDir, idx+".txt");
            try {
                bfw = new BufferedWriter(new FileWriter(targetFile, true));
                String newLine = String.join("\t", new String[] {this.df.format(tretMeth), this.df.format(ctrlMeth),
                                                                           this.df.format(tretExp), this.df.format(ctrlExp)});
                bfw.write(newLine);
                bfw.newLine();
            } catch (IOException ie) {
                ie.printStackTrace();
            } finally {
                countDownLatch.countDown();
                if (bfw != null) {
                    try {
                        bfw.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        });

        double tretMethLevel, ctrlMethLevel, tretExpression, ctrlExpression;
        Runnable runnable;
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            tretMethLevel = tretMethylationLevels[peakIdx];
            ctrlMethLevel = ctrlMethylationLevels[peakIdx];
            tretExpression = tretExpressions[peakIdx];
            ctrlExpression = ctrlExpressions[peakIdx];

            runnable = taskMaker.createTask(tretMethLevel, ctrlMethLevel, tretExpression, ctrlExpression, peakIdx);
            this.executorService.submit(runnable);
        }

        try {
            countDownLatch.await();
        } catch (InterruptedException ie) {
            ie.printStackTrace();
            this.executorService.shutdown();
        }
    }

    /**
     * quantify methylation level and expression for each peak
     */
    public void quantify() {
        double[] tretIPOverdispersions = this.parameters.getTretIPOverdispersion();
        double[] tretINPUTOverdispersions = this.parameters.getTretINPUTOverdispersion();
        double[] ctrlIPOverdispersions = this.parameters.getCtrlIPOverdispersion();
        double[] ctrlINPUTOverdispersions = this.parameters.getCtrlINPUTOverdispersion();

        double quantifiedTretIPOverdispersion = this.calcMedian(tretIPOverdispersions);
        double quantifiedTretINPUTOverdispersion = this.calcMedian(tretINPUTOverdispersions);
        double quantifiedCtrlIPOverdispersion = this.calcMedian(ctrlIPOverdispersions);
        double quantifiedCtrlINPUTOverdispersion = this.calcMedian(ctrlINPUTOverdispersions);
        this.parameters.setOverdispersions(quantifiedTretIPOverdispersion, quantifiedTretINPUTOverdispersion,
                                           quantifiedCtrlIPOverdispersion, quantifiedCtrlINPUTOverdispersion);
        CountDownLatch countDownLatch = new CountDownLatch(this.peakNumber);

        double[] quantifiedTretMethLevels = new double[this.peakNumber];
        double[] quantifiedCtrlMethLevels = new double[this.peakNumber];
        double[] quantifiedTretExpressions = new double[this.peakNumber];
        double[] quantifiedCtrlExpressions = new double[this.peakNumber];
        QuantificationTask taskMaker = (idx -> () -> {
            File targetFile = new File(this.tmpDir, idx+".txt");
            BufferedReader bfr = null;
            double[] tretMethValues, ctrlMethValues, tretExpValues, ctrlExpValues;

            try {
                bfr = new BufferedReader(new InputStreamReader(new FileInputStream(targetFile)));
                tretMethValues = new double[this.samplingTime];
                ctrlMethValues = new double[this.samplingTime];
                tretExpValues = new double[this.samplingTime];
                ctrlExpValues = new double[this.samplingTime];

                int lineNum = 0;
                String line = "";
                String[] info;
                double tretMeth, ctrlMeth, tretExp, ctrlExp;
                while (line!=null) {
                    line = bfr.readLine();
                    if (line!=null) {
                        info = line.split("\t");
                        tretMeth = Double.valueOf(info[0]);
                        ctrlMeth = Double.valueOf(info[1]);
                        tretExp = Double.valueOf(info[2]);
                        ctrlExp = Double.valueOf(info[3]);

                        tretMethValues[lineNum] = tretMeth;
                        ctrlMethValues[lineNum] = ctrlMeth;
                        tretExpValues[lineNum] = tretExp;
                        ctrlExpValues[lineNum] = ctrlExp;
                    }
                }

                double tretMethylationValue = this.calcMedian(tretMethValues);
                double ctrlMethylationValue = this.calcMedian(ctrlMethValues);
                double tretExpressionValue = this.calcMedian(tretExpValues);
                double ctrlExpressionValue = this.calcMedian(ctrlExpValues);
                this.lock.lock();
                quantifiedTretMethLevels[idx] = tretMethylationValue;
                quantifiedCtrlMethLevels[idx] = ctrlMethylationValue;
                quantifiedTretExpressions[idx] = tretExpressionValue;
                quantifiedCtrlExpressions[idx] = ctrlExpressionValue;
            } catch (IOException ie) {
                ie.printStackTrace();
            } finally {
                countDownLatch.countDown();
                this.lock.unlock();
                if (bfr!=null) {
                    try {
                        bfr.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                tretMethValues = null;
                ctrlMethValues = null;
                tretExpValues = null;
                ctrlExpValues = null;
            }
        });

        Runnable runnable;
        for (int peakIdx=0; peakIdx<this.peakNumber; peakIdx++) {
            runnable = taskMaker.createTask(peakIdx);
            this.executorService.submit(runnable);
        }

        try {
            countDownLatch.await();
        } catch (InterruptedException ie) {
            ie.printStackTrace();
        } finally {
            this.executorService.shutdown();
        }
        this.parameters.setQuantifiedTretMethLevel(quantifiedTretMethLevels);
        this.parameters.setQuantifiedCtrlMethLevel(quantifiedCtrlMethLevels);
        this.parameters.setQuantifiedTretExpression(quantifiedTretExpressions);
        this.parameters.setQuantifiedCtrlExpression(quantifiedCtrlExpressions);
    }

    /**
     * calculate median values
     * @param values array
     * @return median
     */
    public double calcMedian(double[] values) {
        if (values.length == 1)
            return values[0];
        values = Arrays.stream(values).skip(this.burnIn).toArray();
        Arrays.sort(values);
        int medianIdx = values.length / 2;
        if (values.length % 2 == 0)
            return (values[medianIdx-1] + values[medianIdx]) * 0.5;
        else
            return values[medianIdx];
    }

    public int[][] getTreatmentIPReads() {
        return this.treatmentIPReads;
    }

    public int[][] getTreatmentINPUTReads() {
        return this.treatmentINPUTReads;
    }

    public int[][] getControlIPReads() {
        return this.controlIPReads;
    }

    public int[][] getControlINPUTReads() {
        return this.controlINPUTReads;
    }

    public int getTretIndividualNumber() {
        return this.tretIndividualNumber;
    }

    public int getCtrlIndividualNumber() {
        return this.ctrlIndividualNumber;
    }

    public BackgroundExpressionSampler getTretBackgroundExpressionSampler(int peakIdx) {
        return this.treatmentBackgroundExpressionSamplers[peakIdx];
    }

    public BackgroundExpressionSampler getCtrlBackgroundExpressionSampler(int peakIdx) {
        return this.controlBackgroundExpressionSamplers[peakIdx];
    }

    public MethylationLevelSampler getTretMethylationLevelSampler() {
        return this.tretMethylationLevelSampler;
    }

    public MethylationLevelSampler getCtrlMethylationLevelSampler() {
        return this.ctrlMethylationLevelSampler;
    }
}
