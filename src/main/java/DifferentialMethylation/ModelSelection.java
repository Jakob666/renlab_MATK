package DifferentialMethylation;

import ProbabilityCalculation.ProbabilityCalculator;
import Quantification.BackgroundExpressionSampler;
import Quantification.MHSampling;
import Quantification.MethylationLevelSampler;
import Quantification.OverdispersionSampler;
import SeqDataModel.ReadsExpectation;

import java.util.ArrayList;

public abstract class ModelSelection {
    private int individualNumber;
    protected Parameters parameters;
    protected OverdispersionSampler tretIPOverdispersionSampler, tretINPUTOverdispersionSampler, ctrlIPOverdispersionSampler, ctrlINPUTOverdispersionSampler,
              tretIPOverdispersionPseudoSampler, tretINPUTOverdispersionPseudoSampler, ctrlIPOverdispersionPseudoSampler, ctrlINPUTOverdispersionPseudoSampler;
    protected MethylationLevelSampler tretMethylationLevelSampler, ctrlMethylationLevelSampler, tretMethylationLevelPseudoSampler, ctrlMethylationLevelPseudoSampler;
    protected BackgroundExpressionSampler treatmentBackgroundExpressionSampler, controlBackgroundExpressionSampler, treatmentBackgroundExpressionPseudoSampler, controlBackgroundExpressionPseudoSampler;
    // shape 1 × samplingTime
    protected double[] treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion,
                       treatmentBkgExp, controlBkgExp, treatmentMethylationLevel, controlMethylationLevel;
    // shape 1 × individualNumber
    protected int[] treatmentIPReads, treatmentINPUTReads, controlIPReads, controlINPUTReads;
    protected double[] treatmentIPExpectations, treatmentINPUTExpectations, controlIPExpectations, controlINPUTExpectations;
    protected MHSampling mhSampling = new MHSampling();
    protected ArrayList<Double> tretIPOverdispersionList, tretINPUTOverdispersionList, ctrlIPOverdispersionList, ctrlINPUTOverdispersionList,
                                tretMethyLationList, ctrlMethylationList, tretBkgExpList, ctrlBkgExpList;

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
                          BackgroundExpressionSampler treatmentBackgroundExpressionSampler,
                          BackgroundExpressionSampler controlBackgroundExpressionSampler) {
        this.tretMethylationLevelSampler = tretMethylationLevelSampler;
        this.ctrlMethylationLevelSampler = ctrlMethylationLevelSampler;
        this.tretIPOverdispersionSampler = tretIPOverdispersionSampler;
        this.tretINPUTOverdispersionSampler = tretINPUTOverdispersionSampler;
        this.ctrlIPOverdispersionSampler = ctrlIPOverdispersionSampler;
        this.ctrlINPUTOverdispersionSampler = ctrlINPUTOverdispersionSampler;
        this.controlBackgroundExpressionSampler = controlBackgroundExpressionSampler;
        this.treatmentBackgroundExpressionSampler = treatmentBackgroundExpressionSampler;
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
                                             double[] treatmentIPExpectations, double[] treatmentINPUTExpectations,
                                             double[] controlIPExpectations, double[] controlINPUTExpectations) {
        this.individualNumber = treatmentIPReads.length;

        this.treatmentIPReads = new int[this.individualNumber];
        System.arraycopy(treatmentIPReads, 0, this.treatmentIPReads, 0, this.individualNumber);
        this.treatmentINPUTReads = new int[this.individualNumber];
        System.arraycopy(treatmentINPUTReads, 0, this.treatmentINPUTReads, 0, this.individualNumber);
        this.controlIPReads = new int[this.individualNumber];
        System.arraycopy(controlIPReads, 0, this.controlIPReads, 0, this.individualNumber);
        this.controlINPUTReads = new int[this.individualNumber];
        System.arraycopy(controlINPUTReads, 0, this.controlINPUTReads, 0, this.individualNumber);

        this.treatmentIPExpectations = new double[this.individualNumber];
        System.arraycopy(treatmentIPExpectations, 0, this.treatmentIPExpectations, 0, this.individualNumber);
        this.treatmentINPUTExpectations = new double[this.individualNumber];
        System.arraycopy(treatmentINPUTExpectations, 0, this.treatmentINPUTExpectations, 0, this.individualNumber);
        this.controlIPExpectations = new double[this.individualNumber];
        System.arraycopy(controlIPExpectations, 0, this.controlIPExpectations, 0, this.individualNumber);
        this.controlINPUTExpectations = new double[this.individualNumber];
        System.arraycopy(controlINPUTExpectations, 0, this.controlINPUTExpectations, 0, this.individualNumber);
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
                                double[] treatmentBkgExp, double[] controlBkgExp, double[] treatmentMethylationLevel, double[] controlMethylationLevel) {
        this.treatmentIPOverdispersion = treatmentIPOverdispersion;
        this.treatmentINPUTOverdispersion = treatmentINPUTOverdispersion;
        this.controlIPOverdispersion = controlIPOverdispersion;
        this.controlINPUTOverdispersion = controlINPUTOverdispersion;
        this.treatmentBkgExp = treatmentBkgExp;
        this.controlBkgExp = controlBkgExp;
        this.treatmentMethylationLevel = treatmentMethylationLevel;
        this.controlMethylationLevel = controlMethylationLevel;

        this.parameters = new Parameters();
        this.parameters.setTretINPUTOverdispersion(treatmentINPUTOverdispersion[0]);
        this.parameters.setTretIPOverdispersion(treatmentIPOverdispersion[0]);
        this.parameters.setCtrlINPUTOverdispersion(controlINPUTOverdispersion[0]);
        this.parameters.setCtrlIPOverdispersion(controlIPOverdispersion[0]);
        this.parameters.setTretMethylation(treatmentMethylationLevel[0]);
        this.parameters.setCtrlMethylation(controlMethylationLevel[0]);
        this.parameters.setTretBkgExp(treatmentBkgExp[0]);
        this.parameters.setCtrlBkgExp(controlBkgExp[0]);
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
        logPrevParamsProba = this.treatmentBkgExpSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        logPrevParamsProba = this.controlBkgExpSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        logPrevParamsProba = this.treatmentMethylationSampling(time, logPrevParamsProba, curModel, usePseudoPrior);
        logPrevParamsProba = this.controlMethylationSampling(time, logPrevParamsProba, curModel, usePseudoPrior);

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
               tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp, newTretINPUTOverdispersion;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        newTretINPUTOverdispersion = this.tretINPUTOverdispersionSampler.randomSample(tretINPUTOverdispersion);
        if (curModel) {
            logLikeProba = this.logLikelihood(tretIPOverdispersion, newTretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, newTretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp);
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
               tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp, newTretIPOverdispersion;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        newTretIPOverdispersion = this.tretIPOverdispersionSampler.randomSample(tretIPOverdispersion);
        if (curModel) {
            logLikeProba = this.logLikelihood(newTretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(newTretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp);
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
               tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp, newCtrlINPUTOverdispersion;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        newCtrlINPUTOverdispersion = this.ctrlINPUTOverdispersionSampler.randomSample(ctrlINPUTOverdispersion);
        if (curModel) {
            logLikeProba = this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, newCtrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, newCtrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp);
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
               tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp, newCtrlIPOverdispersion;
        boolean samplingRes;
        double logCurParamsProba, logLikeProba, logPrior;

        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        newCtrlIPOverdispersion = this.ctrlIPOverdispersionSampler.randomSample(ctrlIPOverdispersion);
        if (curModel) {
            logLikeProba = this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion,
                                              newCtrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        newCtrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp);
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
               tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp, newTretBkgExp;
        boolean samplingRes;
        double[][] expectationReads;
        double[] newIPReadsExpectation, newINPUTReadsExpectation;
        double logCurParamsProba, logLikeProba, logPrior;

        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        newTretBkgExp = this.treatmentBackgroundExpressionSampler.randomSample(tretBkgExp);
        // renew IP and INPUT reads expectations via new sampling background expression IP and INPUT reads count expectations, shape 1 × individualNumber
        expectationReads = this.renewReadsExpectationViaBkgExp(true, new double[]{tretMethLevel}, new double[]{newTretBkgExp});
        newIPReadsExpectation = expectationReads[0];
        newINPUTReadsExpectation = expectationReads[1];
        if (curModel) {
            logLikeProba = this.logLikelihood(newIPReadsExpectation, newINPUTReadsExpectation,
                                              this.controlIPExpectations, this.controlINPUTExpectations,
                                              tretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, newTretBkgExp, ctrlBkgExp);
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
               tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp, newCtrlBkgExp;
        boolean samplingRes;
        double[][] expectationReads;
        double[] newIPReadsExpectation, newINPUTReadsExpectation;
        double logCurParamsProba, logLikeProba, logPrior;

        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        newCtrlBkgExp = this.controlBackgroundExpressionSampler.randomSample(ctrlBkgExp);
        // renew IP and INPUT reads expectations via new sampling background expression IP and INPUT reads count expectations, shape 1 × individualNumber
        expectationReads = this.renewReadsExpectationViaBkgExp(false, new double[]{ctrlMethLevel}, new double[]{newCtrlBkgExp});
        newIPReadsExpectation = expectationReads[0];
        newINPUTReadsExpectation = expectationReads[1];
        if (curModel) {
            logLikeProba = this.logLikelihood(this.treatmentIPExpectations, this.treatmentINPUTExpectations,
                                              newIPReadsExpectation, newINPUTReadsExpectation,
                                              tretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, tretBkgExp, newCtrlBkgExp);
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
     * calculate reads expectations via methylation level
     * @param treatment true, if calculate for treatment group; otherwise, for control group
     * @param methLevel methylation level 1×1
     * @return reads count expectations, shape 1×individualNumber
     */
    protected double[] renewReadsExpectationViaMethLevel(boolean treatment, double[] methLevel) {
        double[][] newIPReadsExpectation;
        // shape individualNumber × geneNumber
        int[][] ipReads, inputReads;
        if (treatment) {
            ipReads = this.dataForBkgExpCalculation(this.treatmentIPReads);
            inputReads = this.dataForBkgExpCalculation(this.treatmentINPUTReads);
        } else {
            ipReads = this.dataForBkgExpCalculation(this.controlIPReads);
            inputReads = this.dataForBkgExpCalculation(this.controlINPUTReads);
        }
        ReadsExpectation re = new ReadsExpectation(ipReads, inputReads, methLevel);
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
    private double[][] renewReadsExpectationViaBkgExp(boolean treatment, double[] methLevel, double[] bkgExp) {

        double[][] newIPReadsExpectation, newINPUTReadsExpectation;
        int[][] ipReads, inputReads;
        if (treatment) {
            ipReads = this.dataForBkgExpCalculation(this.treatmentIPReads);
            inputReads = this.dataForBkgExpCalculation(this.treatmentINPUTReads);
        } else {
            ipReads = this.dataForBkgExpCalculation(this.controlIPReads);
            inputReads = this.dataForBkgExpCalculation(this.controlINPUTReads);
        }
        ReadsExpectation re = new ReadsExpectation(ipReads, inputReads, methLevel, bkgExp);
        // shape individualNumber × 1
        newIPReadsExpectation = re.getIPReadsExepectation();
        newINPUTReadsExpectation = re.getINPUTReadsExpectation();

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
        double proba = 0, logTretIPProba, logTretINPUTProba, logCtrlIPProba, logCtrlINPUTProba;
        logTretIPProba = ProbabilityCalculator.logNegativeProbability(this.treatmentIPReads, this.treatmentIPExpectations, treatmentIPOverdispersion);
        logTretINPUTProba = ProbabilityCalculator.logNegativeProbability(this.treatmentINPUTReads, this.treatmentINPUTExpectations, treatmentINPUTOverdispersion);
        logCtrlINPUTProba = ProbabilityCalculator.logNegativeProbability(this.controlINPUTReads, this.controlINPUTExpectations, controlINPUTOverdispersion);
        logCtrlIPProba = ProbabilityCalculator.logNegativeProbability(this.controlIPReads, this.controlIPExpectations, controlIPOverdispersion);
        proba += logTretIPProba + logTretINPUTProba + logCtrlIPProba + logCtrlINPUTProba;

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
                                 double treatmentIPOverdispersion, double treatmentINPUTOverdispersion,
                                 double controlIPOverdispersion, double controlINPUTOverdispersion) {
        double proba = 0;
        proba += ProbabilityCalculator.logNegativeProbability(this.treatmentIPReads, treatmentIPExpectations, treatmentIPOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(this.treatmentINPUTReads, treatmentINPUTExpectations, treatmentINPUTOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(this.controlINPUTReads, controlINPUTExpectations, controlINPUTOverdispersion);
        proba += ProbabilityCalculator.logNegativeProbability(this.controlIPReads, controlIPExpectations, controlIPOverdispersion);

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
                                 double treatmentBackgroundExpression, double controlBackgroundExpression);

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
               tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp;

        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();

        double logLike= this.logLikelihood(tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion);
        double prior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp);
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

        if (this.tretBkgExpList == null)
            this.tretBkgExpList = new ArrayList<>();
        this.tretBkgExpList.add(this.parameters.getTretBkgExp());
        if (this.ctrlBkgExpList == null)
            this.ctrlBkgExpList = new ArrayList<>();
        this.ctrlBkgExpList.add(this.parameters.getCtrlBkgExp());
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
        this.tretIPOverdispersionList.clear();
        this.tretINPUTOverdispersionList.clear();
        this.ctrlIPOverdispersionList.clear();
        this.ctrlINPUTOverdispersionList.clear();
    }

    /**
     * release redundant instance for saving memory
     */
    protected void helpGC() {
        this.cleanUpSamplingList();
        this.tretBkgExpList = null;
        this.ctrlBkgExpList = null;
        this.tretMethyLationList = null;
        this.ctrlMethylationList = null;
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
