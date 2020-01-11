package DifferentialMethylation;

import ProbabilityCalculation.ProbabilityCalculator;
import Quantification.*;
import SeqDataModel.ReadsExpectation;

import java.util.ArrayList;

public abstract class ModelSelection {
    private int individualNumber;
    protected Parameters parameters;
    protected OverdispersionSampler tretIPOverdispersionSampler, tretINPUTOverdispersionSampler, ctrlIPOverdispersionSampler, ctrlINPUTOverdispersionSampler,
              tretIPOverdispersionPseudoSampler, tretINPUTOverdispersionPseudoSampler, ctrlIPOverdispersionPseudoSampler, ctrlINPUTOverdispersionPseudoSampler;
    protected NonSpecificEnrichmentSampler nonspecificEnrichSampler, nonspecificEnrichPseudoSampler;
    protected MethylationLevelSampler tretMethylationLevelSampler, ctrlMethylationLevelSampler, tretMethylationLevelPseudoSampler, ctrlMethylationLevelPseudoSampler;
    protected BackgroundExpressionSampler treatmentBackgroundExpressionSampler, controlBackgroundExpressionSampler, treatmentBackgroundExpressionPseudoSampler, controlBackgroundExpressionPseudoSampler,
                                          treatmentNonPeakExpressionSampler, controlNonPeakExpressionSampler, treatmentNonPeakExpressionPseudoSampler, controlNonPeakExpressionPseudoSampler;
    // shape 1 × samplingTime
    protected double[] treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion,
                       treatmentBkgExp, treatmentNonPeakBkgExp, controlBkgExp, controlNonPeakBkgExp, treatmentMethylationLevel, controlMethylationLevel, nonspecificEnrichment;
    // shape 1 × individualNumber
    protected int[] treatmentIPReads, treatmentINPUTReads, controlIPReads, controlINPUTReads,
                    treatmentIPNonPeakReads, treatmentINPUTNonPeakReads, controlIPNonPeakReads, controlINPUTNonPeakReads;
    protected double[] treatmentIPExpectations, treatmentINPUTExpectations, controlIPExpectations, controlINPUTExpectations,
                       treatmentIPNonPeakExpect, treatmentINPUTNonPeakExpect, controlIPNonPeakExpect, controlINPUTNonPeakExpect;
    protected MHSampling mhSampling = new MHSampling();
    protected ArrayList<Double> tretIPOverdispersionList, tretINPUTOverdispersionList, ctrlIPOverdispersionList, ctrlINPUTOverdispersionList,
                                tretMethyLationList, ctrlMethylationList, nonspecificEnrichList,
                                tretBkgExpList, ctrlBkgExpList, tretNonPeakExpList, ctrlNonPeakExpList;

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
                          NonSpecificEnrichmentSampler nonspecificEnrichSampler,
                          OverdispersionSampler tretIPOverdispersionSampler, OverdispersionSampler tretINPUTOverdispersionSampler,
                          OverdispersionSampler ctrlIPOverdispersionSampler, OverdispersionSampler ctrlINPUTOverdispersionSampler,
                          BackgroundExpressionSampler treatmentBackgroundExpressionSampler, BackgroundExpressionSampler treatmentNonPeakExpressionSampler,
                          BackgroundExpressionSampler controlBackgroundExpressionSampler, BackgroundExpressionSampler controlNonPeakExpressionSampler) {
        this.tretMethylationLevelSampler = tretMethylationLevelSampler;
        this.ctrlMethylationLevelSampler = ctrlMethylationLevelSampler;
        this.nonspecificEnrichSampler = nonspecificEnrichSampler;
        this.tretIPOverdispersionSampler = tretIPOverdispersionSampler;
        this.tretINPUTOverdispersionSampler = tretINPUTOverdispersionSampler;
        this.ctrlIPOverdispersionSampler = ctrlIPOverdispersionSampler;
        this.ctrlINPUTOverdispersionSampler = ctrlINPUTOverdispersionSampler;
        this.controlBackgroundExpressionSampler = controlBackgroundExpressionSampler;
        this.treatmentBackgroundExpressionSampler = treatmentBackgroundExpressionSampler;
        this.treatmentNonPeakExpressionSampler = treatmentNonPeakExpressionSampler;
        this.controlNonPeakExpressionSampler = controlNonPeakExpressionSampler;
    }

    /**
     * treatment and control group reads count of a single gene
     * @param treatmentIPReads treatment IP gene reads count, shape 1 × individualNumber
     * @param treatmentINPUTReads treatment INPUT gene reads count, shape 1 × individualNumber
     * @param controlIPReads control IP gene reads count, shape 1 × individualNumber
     * @param controlINPUTReads control INPUT gene reads count, shape 1 × individualNumber
     * @param controlIPExpectations treatment IP gene reads count expectation, shape 1 × individualNumber
     * @param controlINPUTExpectations treatment INPUT gene reads count expectation, shape 1 × individualNumber
     * @param treatmentIPExpectations control IP gene reads count expectation, shape 1 × individualNumber
     * @param treatmentINPUTExpectations control INPUT gene reads count expectation, shape 1 × individualNumber
     */
    public void setTreatmentControlGeneReads(int[] treatmentIPReads, int[] treatmentINPUTReads,
                                             int[] controlIPReads, int[] controlINPUTReads,
                                             int[] treatmentIPNonPeakReads, int[] treatmentINPUTNonPeakReads,
                                             int[] controlIPNonPeakReads, int[] controlINPUTNonPeakReads,
                                             double[] treatmentIPExpectations, double[] treatmentINPUTExpectations,
                                             double[] controlIPExpectations, double[] controlINPUTExpectations,
                                             double[] treatmentIPNonPeakExpect, double[] treatmentINPUTNonPeakExpect,
                                             double[] controlIPNonPeakExpect, double[] controlINPUTNonPeakExpect) {
        this.individualNumber = treatmentIPReads.length;

        this.treatmentIPReads = new int[this.individualNumber];
        System.arraycopy(treatmentIPReads, 0, this.treatmentIPReads, 0, this.individualNumber);
        this.treatmentINPUTReads = new int[this.individualNumber];
        System.arraycopy(treatmentINPUTReads, 0, this.treatmentINPUTReads, 0, this.individualNumber);
        this.controlIPReads = new int[this.individualNumber];
        System.arraycopy(controlIPReads, 0, this.controlIPReads, 0, this.individualNumber);
        this.controlINPUTReads = new int[this.individualNumber];
        System.arraycopy(controlINPUTReads, 0, this.controlINPUTReads, 0, this.individualNumber);
        this.treatmentIPNonPeakReads = new int[this.individualNumber];
        System.arraycopy(treatmentIPNonPeakReads, 0, this.treatmentIPNonPeakReads, 0, this.individualNumber);
        this.treatmentINPUTNonPeakReads = new int[this.individualNumber];
        System.arraycopy(treatmentINPUTNonPeakReads, 0, this.treatmentINPUTNonPeakReads, 0, this.individualNumber);
        this.controlIPNonPeakReads = new int[this.individualNumber];
        System.arraycopy(controlIPNonPeakReads, 0, this.controlIPNonPeakReads, 0, this.individualNumber);
        this.controlINPUTNonPeakReads = new int[this.individualNumber];
        System.arraycopy(controlINPUTNonPeakReads, 0, this.controlINPUTNonPeakReads, 0, this.individualNumber);

        this.treatmentIPExpectations = new double[this.individualNumber];
        System.arraycopy(treatmentIPExpectations, 0, this.treatmentIPExpectations, 0, this.individualNumber);
        this.treatmentINPUTExpectations = new double[this.individualNumber];
        System.arraycopy(treatmentINPUTExpectations, 0, this.treatmentINPUTExpectations, 0, this.individualNumber);
        this.controlIPExpectations = new double[this.individualNumber];
        System.arraycopy(controlIPExpectations, 0, this.controlIPExpectations, 0, this.individualNumber);
        this.controlINPUTExpectations = new double[this.individualNumber];
        System.arraycopy(controlINPUTExpectations, 0, this.controlINPUTExpectations, 0, this.individualNumber);
        this.treatmentIPNonPeakExpect = new double[this.individualNumber];
        System.arraycopy(treatmentIPNonPeakExpect, 0, this.treatmentIPNonPeakExpect, 0, this.individualNumber);
        this.treatmentINPUTNonPeakExpect = new double[this.individualNumber];
        System.arraycopy(treatmentINPUTNonPeakExpect, 0, this.treatmentINPUTNonPeakExpect, 0, this.individualNumber);
        this.controlIPNonPeakExpect = new double[this.individualNumber];
        System.arraycopy(controlIPNonPeakExpect, 0, this.controlIPNonPeakExpect, 0, this.individualNumber);
        this.controlINPUTNonPeakExpect = new double[this.individualNumber];
        System.arraycopy(controlINPUTNonPeakExpect, 0, this.controlINPUTNonPeakExpect, 0, this.individualNumber);
    }

    /**
     * set sampling result list for a single gene
     * @param treatmentIPOverdispersion treatment IP overdispersion sampling list, shape 1 × samplingNumber
     * @param treatmentINPUTOverdispersion treatment INPUT overdispersion sampling list, shape 1 × samplingNumber
     * @param controlIPOverdispersion control IP overdispersion sampling list, shape 1 × samplingNumber
     * @param controlINPUTOverdispersion control INPUT overdispersion sampling list, shape 1 × samplingNumber
     * @param treatmentBkgExp treatment background expression sampling list, shape 1 × samplingNumber
     * @param controlBkgExp control background expression sampling list, shape 1 × samplingNumber
     * @param treatmentMethylationLevel treatment methylation level sampling list, shape 1 × samplingNumber
     * @param controlMethylationLevel control methylation level sampling list, shape 1 × samplingNumber
     */
    public void setSamplingList(double[]treatmentIPOverdispersion, double[] treatmentINPUTOverdispersion, double[] controlIPOverdispersion, double[] controlINPUTOverdispersion,
                                double[] treatmentBkgExp, double[] controlBkgExp, double[] treatmentNonPeakBkgExp, double[] controlNonPeakBkgExp,
                                double[] treatmentMethylationLevel, double[] controlMethylationLevel, double[] nonspecificEnrichment) {
        this.treatmentIPOverdispersion = treatmentIPOverdispersion;
        this.treatmentINPUTOverdispersion = treatmentINPUTOverdispersion;
        this.controlIPOverdispersion = controlIPOverdispersion;
        this.controlINPUTOverdispersion = controlINPUTOverdispersion;
        this.treatmentBkgExp = treatmentBkgExp;
        this.treatmentNonPeakBkgExp = treatmentNonPeakBkgExp;
        this.controlBkgExp = controlBkgExp;
        this.controlNonPeakBkgExp = controlNonPeakBkgExp;
        this.treatmentMethylationLevel = treatmentMethylationLevel;
        this.controlMethylationLevel = controlMethylationLevel;
        this.nonspecificEnrichment = nonspecificEnrichment;

        this.parameters = new Parameters();
        this.parameters.setTretINPUTOverdispersion(treatmentINPUTOverdispersion[0]);
        this.parameters.setTretIPOverdispersion(treatmentIPOverdispersion[0]);
        this.parameters.setCtrlINPUTOverdispersion(controlINPUTOverdispersion[0]);
        this.parameters.setCtrlIPOverdispersion(controlIPOverdispersion[0]);
        this.parameters.setTretMethylation(treatmentMethylationLevel[0]);
        this.parameters.setCtrlMethylation(controlMethylationLevel[0]);
        this.parameters.setTretBkgExp(treatmentBkgExp[0]);
        this.parameters.setTretNonPeakBkgExp(treatmentNonPeakBkgExp[0]);
        this.parameters.setCtrlBkgExp(controlBkgExp[0]);
        this.parameters.setCtrlNonPeakBkgExp(controlNonPeakBkgExp[0]);
        this.parameters.setNonspecificEnrichment(nonspecificEnrichment[0]);
    }

    /**
     * an sampling iteration
     * @param time iteration time
     */
    public double iterate(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior) {
        logPrevParamsProba = this.treatmentINPUTOverdispersionSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        logPrevParamsProba = this.treatmentIPOverdispersionSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        logPrevParamsProba = this.controlINPUTOverdispersionSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        logPrevParamsProba = this.controlIPOverdispersionSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        // skip expression sampling, saving time
        // logPrevParamsProba = this.treatmentBkgExpSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        // logPrevParamsProba = this.controlBkgExpSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        // logPrevParamsProba = this.treatmentNonPeakExpSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        // logPrevParamsProba = this.controlNonPeakExpSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        logPrevParamsProba = this.treatmentMethylationSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        logPrevParamsProba = this.controlMethylationSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        logPrevParamsProba = this.nonspecificEnrichSampling(time, logPrevParamsProba, curModel, usePseudoPrior);

        return logPrevParamsProba;
    }

    /**
     * MH sampling for treatment INPUT overdispersion
     * @param time iteration time
     * @param logPrevParamsProba log-likehood + log-prior
     * @param curModel true, if M=j; otherwise false
     */
    private double treatmentINPUTOverdispersionSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, nonspecificEnrich,
               tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp, newTretINPUTOverdispersion;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        nonspecificEnrich = this.parameters.getNonspecificEnrichment();
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        tretNonPeakExp = this.parameters.getTretNonPeakBkgExp();
        ctrlNonPeakExp = this.parameters.getCtrlNonPeakBkgExp();
        newTretINPUTOverdispersion = this.tretINPUTOverdispersionSampler.randomSample(tretINPUTOverdispersion);
        if (curModel) {
            logLikeProba = this.logLikelihood(tretIPOverdispersion, newTretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, newTretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, nonspecificEnrich,
                                        tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp);
            logCurParamsProba = logLikeProba + logPrior;
        } else {
            if (!usePseudoPrior) {
                logPrevParamsProba = this.tretINPUTOverdispersionSampler.getLogDensity(tretINPUTOverdispersion);
                logCurParamsProba = this.tretINPUTOverdispersionSampler.getLogDensity(newTretINPUTOverdispersion);
            } else {
                logPrevParamsProba = this.tretINPUTOverdispersionPseudoSampler.getLogDensity(tretINPUTOverdispersion);
                logCurParamsProba = this.tretINPUTOverdispersionPseudoSampler.getLogDensity(newTretINPUTOverdispersion);
            }
        }
        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        if (samplingRes) {
            this.parameters.setTretINPUTOverdispersion(newTretINPUTOverdispersion);
            this.treatmentINPUTOverdispersion[time] = newTretINPUTOverdispersion;
            return logCurParamsProba;
        }
        else
            this.treatmentINPUTOverdispersion[time] = tretINPUTOverdispersion;

        return logPrevParamsProba;
    }

    /**
     * MH sampling for treatment IP overdispersion
     * @param time iteration time
     */
    private double treatmentIPOverdispersionSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, nonspecificEnrich,
               tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp, newTretIPOverdispersion;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        nonspecificEnrich = this.parameters.getNonspecificEnrichment();
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        tretNonPeakExp = this.parameters.getTretNonPeakBkgExp();
        ctrlNonPeakExp = this.parameters.getCtrlNonPeakBkgExp();
        newTretIPOverdispersion = this.tretIPOverdispersionSampler.randomSample(tretIPOverdispersion);
        if (curModel) {
            logLikeProba = this.logLikelihood(newTretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(newTretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, nonspecificEnrich,
                                        tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp);
            logCurParamsProba = logLikeProba + logPrior;
        } else {
            if (!usePseudoPrior) {
                logPrevParamsProba = this.tretIPOverdispersionSampler.getLogDensity(tretIPOverdispersion);
                logCurParamsProba = this.tretIPOverdispersionSampler.getLogDensity(newTretIPOverdispersion);
            } else {
                logPrevParamsProba = this.tretIPOverdispersionPseudoSampler.getLogDensity(tretIPOverdispersion);
                logCurParamsProba = this.tretIPOverdispersionPseudoSampler.getLogDensity(newTretIPOverdispersion);
            }

        }
        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        if (samplingRes) {
            this.treatmentIPOverdispersion[time] = newTretIPOverdispersion;
            this.parameters.setTretIPOverdispersion(newTretIPOverdispersion);
            return logCurParamsProba;
        } else
            this.treatmentIPOverdispersion[time] = tretIPOverdispersion;
        return logPrevParamsProba;
    }

    /**
     * MH sampling for control INPUT overdispersion
     * @param time iteration time
     */
    private double controlINPUTOverdispersionSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior){
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, nonspecificEnrich,
               tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp, newCtrlINPUTOverdispersion;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        nonspecificEnrich = this.parameters.getNonspecificEnrichment();
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        tretNonPeakExp = this.parameters.getTretNonPeakBkgExp();
        ctrlNonPeakExp = this.parameters.getCtrlNonPeakBkgExp();
        newCtrlINPUTOverdispersion = this.ctrlINPUTOverdispersionSampler.randomSample(ctrlINPUTOverdispersion);
        if (curModel) {
            logLikeProba = this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, newCtrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, newCtrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, nonspecificEnrich,
                                        tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp);
            logCurParamsProba = logLikeProba + logPrior;

        } else {
            if (!usePseudoPrior) {
                logPrevParamsProba = this.ctrlINPUTOverdispersionSampler.getLogDensity(ctrlINPUTOverdispersion);
                logCurParamsProba = this.ctrlINPUTOverdispersionSampler.getLogDensity(newCtrlINPUTOverdispersion);
            } else {
                logPrevParamsProba = this.ctrlINPUTOverdispersionPseudoSampler.getLogDensity(ctrlINPUTOverdispersion);
                logCurParamsProba = this.ctrlINPUTOverdispersionPseudoSampler.getLogDensity(newCtrlINPUTOverdispersion);
            }
        }
        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        if (samplingRes) {
            this.controlINPUTOverdispersion[time] = newCtrlINPUTOverdispersion;
            this.parameters.setCtrlINPUTOverdispersion(newCtrlINPUTOverdispersion);
            return logCurParamsProba;
        } else
            this.controlINPUTOverdispersion[time] = ctrlINPUTOverdispersion;

        return logPrevParamsProba;
    }

    /**
     * MH sampling for control IP overdispersion
     * @param time iteration time
     */
    private double controlIPOverdispersionSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior){
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, nonspecificEnrich,
               tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp, newCtrlIPOverdispersion;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        nonspecificEnrich = this.parameters.getNonspecificEnrichment();
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        tretNonPeakExp = this.parameters.getTretNonPeakBkgExp();
        ctrlNonPeakExp = this.parameters.getCtrlNonPeakBkgExp();
        newCtrlIPOverdispersion = this.ctrlIPOverdispersionSampler.randomSample(ctrlIPOverdispersion);
        if (curModel) {
            logLikeProba = this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion,
                                              newCtrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        newCtrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, nonspecificEnrich,
                                        tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp);
            logCurParamsProba = logLikeProba + logPrior;
        } else {
            if (!usePseudoPrior) {
                logPrevParamsProba = this.ctrlIPOverdispersionSampler.getLogDensity(ctrlIPOverdispersion);
                logCurParamsProba = this.ctrlIPOverdispersionSampler.getLogDensity(newCtrlIPOverdispersion);
            } else {
                logPrevParamsProba = this.ctrlIPOverdispersionPseudoSampler.getLogDensity(ctrlIPOverdispersion);
                logCurParamsProba = this.ctrlIPOverdispersionPseudoSampler.getLogDensity(newCtrlIPOverdispersion);
            }
        }
        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        if (samplingRes) {
            this.controlIPOverdispersion[time] = newCtrlIPOverdispersion;
            this.parameters.setCtrlIPOverdispersion(newCtrlIPOverdispersion);
            return logCurParamsProba;
        } else
            this.controlIPOverdispersion[time] = ctrlIPOverdispersion;

        return logPrevParamsProba;
    }

    /**
     * MH sampling for treatment group methylation level
     * @param time iteration time
     */
    protected abstract double treatmentMethylationSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior);

    /**
     * MH sampling for control group methylation level
     * @param time iteration time
     */
    protected abstract double controlMethylationSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior);

    /**
     * MH sampling for treatment group background expression
     * @param time iteration time
     */
    private double treatmentBkgExpSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, nonspecificEnrich,
               tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp, newTretBkgExp;
        boolean samplingRes;
        double[][] expectationReads;
        double[] newIPReadsExpectation, newINPUTReadsExpectation;
        double logCurParamsProba, logLikeProba, logPrior;

        nonspecificEnrich = this.parameters.getNonspecificEnrichment();
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        tretNonPeakExp = this.parameters.getTretNonPeakBkgExp();
        ctrlNonPeakExp = this.parameters.getCtrlNonPeakBkgExp();
        newTretBkgExp = this.treatmentBackgroundExpressionSampler.randomSample(tretBkgExp);
        // renew IP and INPUT reads expectations via new sampling background expression IP and INPUT reads count expectations, shape 1 × individualNumber
        expectationReads = this.renewReadsExpectationViaBkgExp(true, true, new double[]{tretMethLevel}, new double[]{nonspecificEnrich}, new double[]{newTretBkgExp});
        newIPReadsExpectation = expectationReads[0];
        newINPUTReadsExpectation = expectationReads[1];
        if (curModel) {
            logLikeProba = this.logLikelihood(newIPReadsExpectation, newINPUTReadsExpectation,
                                              this.controlIPExpectations, this.controlINPUTExpectations,
                                              this.treatmentIPNonPeakExpect, this.treatmentINPUTNonPeakExpect,
                                              this.controlIPNonPeakExpect, this.controlINPUTNonPeakExpect,
                                              tretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, nonspecificEnrich,
                                        newTretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp);
            logCurParamsProba = logLikeProba + logPrior;

        } else {
            if (!usePseudoPrior) {
                logPrevParamsProba = this.treatmentBackgroundExpressionSampler.getLogDensity(tretBkgExp);
                logCurParamsProba = this.treatmentBackgroundExpressionSampler.getLogDensity(newTretBkgExp);
            } else {
                logPrevParamsProba = this.treatmentBackgroundExpressionPseudoSampler.getLogDensity(tretBkgExp);
                logCurParamsProba = this.treatmentBackgroundExpressionPseudoSampler.getLogDensity(newTretBkgExp);
            }
        }
        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        if (samplingRes) {
            this.parameters.setTretBkgExp(newTretBkgExp);
            this.treatmentBkgExp[time] = newTretBkgExp;
            this.treatmentIPExpectations = newIPReadsExpectation;
            this.treatmentINPUTExpectations = newINPUTReadsExpectation;
            return logCurParamsProba;
        } else
            this.treatmentBkgExp[time] = tretBkgExp;

        return logPrevParamsProba;
    }

    /**
     * MH sampling for control group background expression
     * @param time iteration time
     */
    private double controlBkgExpSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, nonspecificEnrich,
               tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp, newCtrlBkgExp;
        boolean samplingRes;
        double[][] expectationReads;
        double[] newIPReadsExpectation, newINPUTReadsExpectation;
        double logCurParamsProba, logLikeProba, logPrior;

        nonspecificEnrich = this.parameters.getNonspecificEnrichment();
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        tretNonPeakExp = this.parameters.getTretNonPeakBkgExp();
        ctrlNonPeakExp = this.parameters.getCtrlNonPeakBkgExp();
        newCtrlBkgExp = this.controlBackgroundExpressionSampler.randomSample(ctrlBkgExp);
        // renew IP and INPUT reads expectations via new sampling background expression IP and INPUT reads count expectations, shape 1 × individualNumber
        expectationReads = this.renewReadsExpectationViaBkgExp(false, true, new double[]{ctrlMethLevel}, new double[] {nonspecificEnrich}, new double[]{newCtrlBkgExp});
        newIPReadsExpectation = expectationReads[0];
        newINPUTReadsExpectation = expectationReads[1];
        if (curModel) {
            logLikeProba = this.logLikelihood(this.treatmentIPExpectations, this.treatmentINPUTExpectations,
                                              newIPReadsExpectation, newINPUTReadsExpectation,
                                              this.treatmentIPNonPeakExpect, this.treatmentINPUTNonPeakExpect,
                                              this.controlIPNonPeakExpect, this.controlINPUTNonPeakExpect,
                                              tretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, nonspecificEnrich,
                                        tretBkgExp, newCtrlBkgExp, tretNonPeakExp, ctrlNonPeakExp);
            logCurParamsProba = logLikeProba + logPrior;
        } else {
            if (!usePseudoPrior) {
                logPrevParamsProba = this.controlBackgroundExpressionSampler.getLogDensity(ctrlBkgExp);
                logCurParamsProba = this.controlBackgroundExpressionSampler.getLogDensity(newCtrlBkgExp);
            } else {
                logPrevParamsProba = this.controlBackgroundExpressionPseudoSampler.getLogDensity(ctrlBkgExp);
                logCurParamsProba = this.controlBackgroundExpressionPseudoSampler.getLogDensity(newCtrlBkgExp);
            }
        }
        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        if (samplingRes) {
            this.parameters.setCtrlBkgExp(newCtrlBkgExp);
            this.controlBkgExp[time] = newCtrlBkgExp;
            this.controlIPExpectations = newIPReadsExpectation;
            this.controlINPUTExpectations = newINPUTReadsExpectation;
            return logCurParamsProba;
        } else
            this.controlBkgExp[time] = ctrlBkgExp;

        return logPrevParamsProba;
    }

    /**
     * MH sampling for treatment group non-peak region expression
     * @param time iteration time
     */
    private double treatmentNonPeakExpSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, nonspecificEnrich,
               tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp, newTretNonPeakExp;
        boolean samplingRes;
        double[][] expectationReads;
        double[] newIPReadsExpectation, newINPUTReadsExpectation;
        double logCurParamsProba, logLikeProba, logPrior;

        nonspecificEnrich = this.parameters.getNonspecificEnrichment();
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        tretNonPeakExp = this.parameters.getTretNonPeakBkgExp();
        ctrlNonPeakExp = this.parameters.getCtrlNonPeakBkgExp();
        newTretNonPeakExp = this.treatmentBackgroundExpressionSampler.randomSample(tretNonPeakExp);
        // renew IP and INPUT reads expectations via new sampling background expression IP and INPUT reads count expectations, shape 1 × individualNumber
        expectationReads = this.renewReadsExpectationViaBkgExp(true, false, new double[]{tretMethLevel}, new double[]{nonspecificEnrich}, new double[]{newTretNonPeakExp});
        newIPReadsExpectation = expectationReads[0];
        newINPUTReadsExpectation = expectationReads[1];
        if (curModel) {
            logLikeProba = this.logLikelihood(this.treatmentIPExpectations, this.treatmentINPUTExpectations,
                                              this.controlIPExpectations, this.controlINPUTExpectations,
                                              newIPReadsExpectation, newINPUTReadsExpectation,
                                              this.controlIPNonPeakExpect, this.controlINPUTNonPeakExpect,
                                              tretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, nonspecificEnrich,
                                        tretBkgExp, ctrlBkgExp, newTretNonPeakExp, ctrlNonPeakExp);
            logCurParamsProba = logLikeProba + logPrior;

        } else {
            if (!usePseudoPrior) {
                logPrevParamsProba = this.treatmentNonPeakExpressionSampler.getLogDensity(tretNonPeakExp);
                logCurParamsProba = this.treatmentNonPeakExpressionSampler.getLogDensity(newTretNonPeakExp);
            } else {
                logPrevParamsProba = this.treatmentNonPeakExpressionPseudoSampler.getLogDensity(tretNonPeakExp);
                logCurParamsProba = this.treatmentNonPeakExpressionPseudoSampler.getLogDensity(newTretNonPeakExp);
            }
        }
        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        if (samplingRes) {
            this.parameters.setTretNonPeakBkgExp(newTretNonPeakExp);
            this.treatmentNonPeakBkgExp[time] = newTretNonPeakExp;
            this.treatmentIPNonPeakExpect = newIPReadsExpectation;
            this.treatmentINPUTNonPeakExpect = newINPUTReadsExpectation;
            return logCurParamsProba;
        } else
            this.treatmentNonPeakBkgExp[time] = tretNonPeakExp;

        return logPrevParamsProba;
    }

    /**
     * MH sampling for treatment group non-peak region expression
     * @param time iteration time
     */
    private double controlNonPeakExpSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, nonspecificEnrich,
               tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp, newCtrlNonPeakExp;
        boolean samplingRes;
        double[][] expectationReads;
        double[] newIPReadsExpectation, newINPUTReadsExpectation;
        double logCurParamsProba, logLikeProba, logPrior;

        nonspecificEnrich = this.parameters.getNonspecificEnrichment();
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        tretNonPeakExp = this.parameters.getTretNonPeakBkgExp();
        ctrlNonPeakExp = this.parameters.getCtrlNonPeakBkgExp();
        newCtrlNonPeakExp = this.treatmentBackgroundExpressionSampler.randomSample(ctrlNonPeakExp);
        // renew IP and INPUT reads expectations via new sampling background expression IP and INPUT reads count expectations, shape 1 × individualNumber
        expectationReads = this.renewReadsExpectationViaBkgExp(true, false, new double[]{tretMethLevel}, new double[]{nonspecificEnrich}, new double[]{newCtrlNonPeakExp});
        newIPReadsExpectation = expectationReads[0];
        newINPUTReadsExpectation = expectationReads[1];
        if (curModel) {
            logLikeProba = this.logLikelihood(this.treatmentIPExpectations, this.treatmentINPUTExpectations,
                                              this.controlIPExpectations, this.controlINPUTExpectations,
                                              this.treatmentIPNonPeakExpect, this.treatmentINPUTNonPeakExpect,
                                              newIPReadsExpectation, newINPUTReadsExpectation,
                                              tretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, nonspecificEnrich,
                                        tretBkgExp, ctrlBkgExp, tretNonPeakExp, newCtrlNonPeakExp);
            logCurParamsProba = logLikeProba + logPrior;

        } else {
            if (!usePseudoPrior) {
                logPrevParamsProba = this.controlNonPeakExpressionSampler.getLogDensity(ctrlNonPeakExp);
                logCurParamsProba = this.controlNonPeakExpressionSampler.getLogDensity(newCtrlNonPeakExp);
            } else {
                logPrevParamsProba = this.controlNonPeakExpressionPseudoSampler.getLogDensity(ctrlNonPeakExp);
                logCurParamsProba = this.controlNonPeakExpressionPseudoSampler.getLogDensity(newCtrlNonPeakExp);
            }
        }
        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        if (samplingRes) {
            this.parameters.setCtrlNonPeakBkgExp(newCtrlNonPeakExp);
            this.controlNonPeakBkgExp[time] = newCtrlNonPeakExp;
            this.treatmentIPNonPeakExpect = newIPReadsExpectation;
            this.treatmentINPUTNonPeakExpect = newINPUTReadsExpectation;
            return logCurParamsProba;
        } else
            this.controlNonPeakBkgExp[time] = ctrlNonPeakExp;

        return logPrevParamsProba;
    }

    /**
     * MH sampling for treatment and control group nonspecific enrichment ratio
     * @param time iteration time
     */
    protected double nonspecificEnrichSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, nonspecificEnrich,
               tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp, newNonspecificEnrich;
        double[][] newExpectations;
        double[] newTretIPReadsExpectation, newTretIPNonPeakExpectation, newCtrlIPReadsExpectation, newCtrlIPNonPeakExpectation;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        nonspecificEnrich = this.parameters.getNonspecificEnrichment();
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        tretNonPeakExp = this.parameters.getTretNonPeakBkgExp();
        ctrlNonPeakExp = this.parameters.getCtrlNonPeakBkgExp();
        newNonspecificEnrich = this.nonspecificEnrichSampler.randomSample(nonspecificEnrich);
        // renew IP and INPUT reads expectations via new sampling nonspecific enrichment ratio, shape individualNumber × geneNumber
        newExpectations = this.renewReadsExpectationViaNonspecificEnrich(true, new double[]{tretMethLevel}, new double[]{newNonspecificEnrich});
        newTretIPReadsExpectation = newExpectations[0];
        newTretIPNonPeakExpectation = newExpectations[1];

        newExpectations = this.renewReadsExpectationViaNonspecificEnrich(false, new double[]{ctrlMethLevel}, new double[]{newNonspecificEnrich});
        newCtrlIPReadsExpectation = newExpectations[0];
        newCtrlIPNonPeakExpectation = newExpectations[1];

        if (curModel) {
            logLikeProba = this.logLikelihood(newTretIPReadsExpectation, this.treatmentINPUTExpectations,
                                              newCtrlIPReadsExpectation, this.controlINPUTExpectations,
                                              newTretIPNonPeakExpectation, this.treatmentINPUTNonPeakExpect,
                                              newCtrlIPNonPeakExpectation, this.controlINPUTNonPeakExpect,
                                              tretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, newNonspecificEnrich,
                                        tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp);
            logCurParamsProba = logLikeProba + logPrior;
        } else {
            if (!usePseudoPrior) {
                logPrevParamsProba = this.nonspecificEnrichSampler.getLogDensity(nonspecificEnrich);
                logCurParamsProba = this.nonspecificEnrichSampler.getLogDensity(newNonspecificEnrich);
            } else {
                logPrevParamsProba = this.nonspecificEnrichPseudoSampler.getLogDensity(nonspecificEnrich);
                logCurParamsProba = this.nonspecificEnrichPseudoSampler.getLogDensity(newNonspecificEnrich);
            }
        }
        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        if (samplingRes) {
            this.parameters.setNonspecificEnrichment(newNonspecificEnrich);
            this.nonspecificEnrichment[time] = newNonspecificEnrich;
            this.treatmentIPExpectations = newTretIPReadsExpectation;
            this.treatmentIPNonPeakExpect = newTretIPNonPeakExpectation;
            this.controlIPExpectations = newCtrlIPReadsExpectation;
            this.controlIPNonPeakExpect = newCtrlIPNonPeakExpectation;
            return logCurParamsProba;
        } else
            this.nonspecificEnrichment[time] = nonspecificEnrich;

        return logPrevParamsProba;
    }

    /**
     * calculate reads expectations via new nonspecific enrichment ratio and methylation level
     * @param treatment true, if calculate for treatment group; otherwise, for control group
     * @param methLevel methylation level 1×1
     * @return reads count expectations, shape 1×individualNumber
     */
    protected double[][] renewReadsExpectationViaNonspecificEnrich(boolean treatment, double[] methLevel, double[] nonspecificEnrich) {
        double[][] newIPReadsExpectation, newIPNonPeakExpectation, newINPUTReadsExpectation, newINPUTNonPeakExpectation;
        // shape individualNumber × geneNumber
        int[][] ipReads, inputReads, ipNonpeak, inputNonpeak;
        if (treatment) {
            ipReads = this.dataForBkgExpCalculation(this.treatmentIPReads);
            inputReads = this.dataForBkgExpCalculation(this.treatmentINPUTReads);
            ipNonpeak = this.dataForBkgExpCalculation(this.treatmentIPNonPeakReads);
            inputNonpeak = this.dataForBkgExpCalculation(this.treatmentINPUTNonPeakReads);
        } else {
            ipReads = this.dataForBkgExpCalculation(this.controlIPReads);
            inputReads = this.dataForBkgExpCalculation(this.controlINPUTReads);
            ipNonpeak = this.dataForBkgExpCalculation(this.controlIPNonPeakReads);
            inputNonpeak = this.dataForBkgExpCalculation(this.controlINPUTNonPeakReads);
        }

        ReadsExpectation re = new ReadsExpectation(ipReads, inputReads, ipNonpeak, inputNonpeak, methLevel, nonspecificEnrich);
        // shape individualNumber × 1
        newIPReadsExpectation = re.getIPReadsExepectation();
        newIPNonPeakExpectation = re.getIPNonPeakExepectation();
        newINPUTReadsExpectation = re.getINPUTReadsExpectation();
        newINPUTNonPeakExpectation = re.getINPUTNonPeakExpectation();

        double[] newExpectations = new double[this.individualNumber], newNonPeakExpect = new double[this.individualNumber],
                 newINPUTExpectations = new double[this.individualNumber], newINPUTNonPeakExpect = new double[this.individualNumber];
        for (int i=0; i<this.individualNumber; i++) {
            newExpectations[i] = newIPReadsExpectation[i][0];
            newNonPeakExpect[i] = newIPNonPeakExpectation[i][0];
            newINPUTExpectations[i] = newINPUTReadsExpectation[i][0];
            newINPUTNonPeakExpect[i] = newINPUTNonPeakExpectation[i][0];
        }
        ipReads = null;
        inputReads = null;
        ipNonpeak = null;
        inputNonpeak = null;
        newIPReadsExpectation = null;
        re = null;

        return new double[][] {newExpectations, newNonPeakExpect, newINPUTExpectations, newINPUTNonPeakExpect};
    }

    /**
     * calculate reads expectations via new methylation level
     * @param treatment true, if calculate for treatment group; otherwise, for control group
     * @param methLevel methylation level 1×1
     * @return reads count expectations, shape 1×individualNumber
     */
    protected double[] renewReadsExpectationViaMethLevel(boolean treatment, double[] methLevel, double[] nonspecificEnrich) {
        double[][] newIPReadsExpectation;
        // shape individualNumber × geneNumber
        int[][] ipReads, inputReads, ipNonpeak, inputNonpeak;
        if (treatment) {
            ipReads = this.dataForBkgExpCalculation(this.treatmentIPReads);
            inputReads = this.dataForBkgExpCalculation(this.treatmentINPUTReads);
            ipNonpeak = this.dataForBkgExpCalculation(this.treatmentIPNonPeakReads);
            inputNonpeak = this.dataForBkgExpCalculation(this.treatmentINPUTNonPeakReads);
        } else {
            ipReads = this.dataForBkgExpCalculation(this.controlIPReads);
            inputReads = this.dataForBkgExpCalculation(this.controlINPUTReads);
            ipNonpeak = this.dataForBkgExpCalculation(this.controlIPNonPeakReads);
            inputNonpeak = this.dataForBkgExpCalculation(this.controlINPUTNonPeakReads);
        }

        ReadsExpectation re = new ReadsExpectation(ipReads, inputReads, ipNonpeak, inputNonpeak, methLevel, nonspecificEnrich);
        // shape individualNumber × 1
        newIPReadsExpectation = re.getIPReadsExepectation();

        double[] newExpectations = new double[this.individualNumber];
        for (int i=0; i<this.individualNumber; i++) {
            newExpectations[i] = newIPReadsExpectation[i][0];
        }
        ipReads = null;
        inputReads = null;
        newIPReadsExpectation = null;
        re = null;

        return newExpectations;
    }

    /**
     * calculate reads expectations via methylation level
     * @param treatment true, if calculate for treatment group; otherwise, for control group
     * @param bkgExp background expression 1×1
     * @return reads count expectations, shape 1×individualNumber
     */
    private double[][] renewReadsExpectationViaBkgExp(boolean treatment, boolean inPeak, double[] methLevel, double[] nonspecificEnrich,double[] bkgExp) {

        double[][] newIPReadsExpectation, newINPUTReadsExpectation;
        int[][] ipReads, inputReads, ipNonpeak, inputNonpeak;
        if (treatment) {
            ipReads = this.dataForBkgExpCalculation(this.treatmentIPReads);
            inputReads = this.dataForBkgExpCalculation(this.treatmentINPUTReads);
            ipNonpeak = this.dataForBkgExpCalculation(this.treatmentIPNonPeakReads);
            inputNonpeak = this.dataForBkgExpCalculation(this.treatmentINPUTNonPeakReads);
        } else {
            ipReads = this.dataForBkgExpCalculation(this.controlIPReads);
            inputReads = this.dataForBkgExpCalculation(this.controlINPUTReads);
            ipNonpeak = this.dataForBkgExpCalculation(this.controlIPNonPeakReads);
            inputNonpeak = this.dataForBkgExpCalculation(this.controlINPUTNonPeakReads);
        }

        ReadsExpectation re = new ReadsExpectation(ipReads, inputReads, ipNonpeak, inputNonpeak, methLevel, nonspecificEnrich, bkgExp, inPeak);
        // shape individualNumber × 1
        if (inPeak) {
            newIPReadsExpectation = re.getIPReadsExepectation();
            newINPUTReadsExpectation = re.getINPUTReadsExpectation();
        } else {
            newIPReadsExpectation = re.getIPNonPeakExepectation();
            newINPUTReadsExpectation = re.getINPUTNonPeakExpectation();
        }


        double[] ipExpectations = new double[this.individualNumber], inputExpectations = new double[this.individualNumber];
        for (int i=0; i<this.individualNumber; i++) {
            ipExpectations[i] = newIPReadsExpectation[i][0];
            inputExpectations[i] = newINPUTReadsExpectation[i][0];
        }
        ipReads = null;
        inputReads = null;
        newIPReadsExpectation = null;
        newINPUTReadsExpectation = null;
        re = null;

        return new double[][] {ipExpectations, inputExpectations};
    }

    /**
     * data preparation for reads expectation calculation
     * @param readsCount shape 1 × individualNumber, will be transformed to individualNumber × 1
     * @return data for expectation
     */
    private int[][] dataForBkgExpCalculation(int[] readsCount) {
        int[][] data = new int[readsCount.length][1];
        for (int i=0; i<readsCount.length; i++) {
            data[i][0] = readsCount[i];
        }

        return data;
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
        double proba = 0, logTretIPProba, logTretINPUTProba, logCtrlIPProba, logCtrlINPUTProba,
               logTretIPNonPeakProba, logTretINPUTNonPeakProba, logCtrlIPNonPeakProba, logCtrlINPUTNonPeakProba;
        logTretIPProba = ProbabilityCalculator.logNegativeProbability(this.treatmentIPReads, this.treatmentIPExpectations, treatmentIPOverdispersion);
        logTretINPUTProba = ProbabilityCalculator.logNegativeProbability(this.treatmentINPUTReads, this.treatmentINPUTExpectations, treatmentINPUTOverdispersion);
        logCtrlINPUTProba = ProbabilityCalculator.logNegativeProbability(this.controlINPUTReads, this.controlINPUTExpectations, controlINPUTOverdispersion);
        logCtrlIPProba = ProbabilityCalculator.logNegativeProbability(this.controlIPReads, this.controlIPExpectations, controlIPOverdispersion);
        logTretIPNonPeakProba = ProbabilityCalculator.logNegativeProbability(this.treatmentIPNonPeakReads, this.treatmentIPNonPeakExpect, treatmentIPOverdispersion);
        logTretINPUTNonPeakProba = ProbabilityCalculator.logNegativeProbability(this.treatmentINPUTNonPeakReads, this.treatmentINPUTNonPeakExpect, treatmentINPUTOverdispersion);
        logCtrlIPNonPeakProba = ProbabilityCalculator.logNegativeProbability(this.controlIPNonPeakReads, this.controlIPNonPeakExpect, controlIPOverdispersion);
        logCtrlINPUTNonPeakProba = ProbabilityCalculator.logNegativeProbability(this.controlINPUTNonPeakReads, this.controlINPUTNonPeakExpect, controlINPUTOverdispersion);

        proba += logTretIPProba + logTretINPUTProba + logCtrlIPProba + logCtrlINPUTProba +
                 logTretIPNonPeakProba + logTretINPUTNonPeakProba + logCtrlIPNonPeakProba + logCtrlINPUTNonPeakProba;

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
    protected double logLikelihood(double[] treatmentIPExpectations, double[] treatmentINPUTExpectations,
                                   double[] controlIPExpectations, double[] controlINPUTExpectations,
                                   double[] treatmentIPNonPeakExpect, double[] treatmentINPUTNonPeakExpect,
                                   double[] controlIPNonPeakExpect, double[] controlINPUTNonPeakExpect,
                                   double treatmentIPOverdispersion, double treatmentINPUTOverdispersion,
                                   double controlIPOverdispersion, double controlINPUTOverdispersion) {
        double proba = 0;
        proba += ProbabilityCalculator.logNegativeProbability(this.treatmentIPReads, treatmentIPExpectations, treatmentIPOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(this.treatmentINPUTReads, treatmentINPUTExpectations, treatmentINPUTOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(this.controlINPUTReads, controlINPUTExpectations, controlINPUTOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(this.controlIPReads, controlIPExpectations, controlIPOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(this.treatmentIPNonPeakReads, treatmentIPNonPeakExpect, treatmentIPOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(this.treatmentINPUTNonPeakReads, treatmentINPUTNonPeakExpect, treatmentINPUTOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(this.controlIPNonPeakReads, controlIPNonPeakExpect, controlIPOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(this.controlINPUTNonPeakReads, controlINPUTNonPeakExpect, controlINPUTOverdispersion);

        return proba;
    }

    /**
     * parameters logarithm prior probability
     *      log-gamma(tret_ip_overdispersion) + log-gamma(tret_input_overdispersion) +
     *      log-gamma(ctrl_ip_overdispersion) + log-gamma(ctrl_input_overdispersion) +
     *      log-beta(tret_methylation_level) + log-beta(ctrl_methylation_level) +
     *      log-normal(tret_bkg_exp) + log-normal(ctrl_bkg_exp)
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
    protected abstract double logPriority(double treatmentIPOverdispersion, double treatmentINPUTOverdispersion,
                                 double controlIPOverdispersion, double controlINPUTOverdispersion,
                                 double treatmentMethylationLevel, double controlMethylationLevel,
                                 double nonspecificEnrich,
                                 double treatmentBackgroundExpression, double controlBackgroundExpression,
                                 double treatmentNonPeakExpression, double controlNonPeakExpression);

    /**
     * model parameters prior probability in pseudo distribution
     * @return pseudo prior probability
     */
    protected abstract double paramsPseudoDistributionProbability(boolean usePseudoPrior);

    /**
     * calculate log-posterior = log-likelihood + log-prior at the end of an iteration
     * @return log-posterior
     */
    public double paramsLogFullConditionalProbability() {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, nonspecificEnrich,
               tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp;

        nonspecificEnrich = this.parameters.getNonspecificEnrichment();
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        tretNonPeakExp = this.parameters.getTretNonPeakBkgExp();
        ctrlNonPeakExp = this.parameters.getCtrlNonPeakBkgExp();

        double logLike= this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion);
        double prior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, nonspecificEnrich,
                                        tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp);
        return logLike + prior;
    }

    /**
     * record current iteration sampling parameters
     */
    public void recordCurrentIterationParams() {
        if (this.tretIPOverdispersionList == null)
            this.tretIPOverdispersionList = new ArrayList<>();
        this.tretIPOverdispersionList.add(this.parameters.getTretIPOverdispersion());
        if (this.tretINPUTOverdispersionList == null)
            this.tretINPUTOverdispersionList = new ArrayList<>();
        this.tretINPUTOverdispersionList.add(this.parameters.getTretINPUTOverdispersion());

        if (this.ctrlIPOverdispersionList == null)
            this.ctrlIPOverdispersionList = new ArrayList<>();
        this.ctrlIPOverdispersionList.add(this.parameters.getCtrlIPOverdispersion());
        if (this.ctrlINPUTOverdispersionList == null)
            this.ctrlINPUTOverdispersionList = new ArrayList<>();
        this.ctrlINPUTOverdispersionList.add(this.parameters.getCtrlINPUTOverdispersion());

        if (this.tretMethyLationList == null)
            this.tretMethyLationList = new ArrayList<>();
        this.tretMethyLationList.add(this.parameters.getTretMethylation());
        if (this.ctrlMethylationList == null)
            this.ctrlMethylationList = new ArrayList<>();
        this.ctrlMethylationList.add(this.parameters.getCtrlMethylation());

        if (this.nonspecificEnrichList == null)
            this.nonspecificEnrichList = new ArrayList<>();
        this.nonspecificEnrichList.add(this.parameters.getNonspecificEnrichment());

        if (this.tretBkgExpList == null)
            this.tretBkgExpList = new ArrayList<>();
        this.tretBkgExpList.add(this.parameters.getTretBkgExp());
        if (this.ctrlBkgExpList == null)
            this.ctrlBkgExpList = new ArrayList<>();
        this.ctrlBkgExpList.add(this.parameters.getCtrlBkgExp());

        if (this.tretNonPeakExpList == null)
            this.tretNonPeakExpList = new ArrayList<>();
        this.tretNonPeakExpList.add(this.parameters.getTretNonPeakBkgExp());
        if (this.ctrlNonPeakExpList == null)
            this.ctrlNonPeakExpList = new ArrayList<>();
        this.ctrlNonPeakExpList.add(this.parameters.getCtrlNonPeakBkgExp());
    }

    /**
     * use the first-time sampling result to regress pseudo prior for each parameters
     */
    protected abstract void setModelParameterPseudoPrior();

    /**
     * cleaning up sampling result to save memory
     */
    protected void cleanUpSamplingList() {
        if (this.tretBkgExpList == null)
            return;
        this.tretMethyLationList.clear();
        this.ctrlMethylationList.clear();
        this.tretBkgExpList.clear();
        this.ctrlBkgExpList.clear();
        this.tretNonPeakExpList.clear();
        this.ctrlNonPeakExpList.clear();
        this.tretIPOverdispersionList.clear();
        this.tretINPUTOverdispersionList.clear();
        this.ctrlIPOverdispersionList.clear();
        this.ctrlINPUTOverdispersionList.clear();
        this.nonspecificEnrichList.clear();
    }

    /**
     * release redundant instance for saving memory
     */
    protected void helpGC() {
        this.cleanUpSamplingList();
        this.tretBkgExpList = null;
        this.ctrlBkgExpList = null;
        this.tretNonPeakExpList = null;
        this.ctrlNonPeakExpList = null;
        this.tretMethyLationList = null;
        this.ctrlMethylationList = null;
        this.nonspecificEnrichList = null;
        this.tretIPOverdispersionList = null;
        this.tretINPUTOverdispersionList = null;
        this.ctrlIPOverdispersionList = null;
        this.ctrlINPUTOverdispersionList = null;
        this.treatmentIPOverdispersion = null;
        this.treatmentINPUTOverdispersion = null;
        this.controlIPOverdispersion = null;
        this.controlINPUTOverdispersion = null;
        this.treatmentMethylationLevel = null;
        this.controlMethylationLevel = null;
        this.treatmentBkgExp = null;
        this.controlBkgExp = null;
    }
}
