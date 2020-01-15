package DifferentialMethylation;

import DifferentialMethylation.PseudoPrior.PseudoBetaDistribution;
import DifferentialMethylation.PseudoPrior.PseudoInverseGammaDistribution;
import DifferentialMethylation.PseudoPrior.PseudoLogNormalDistribution;
import Quantification.BackgroundExpressionSampler;
import Quantification.MethylationLevelSampler;
import Quantification.NonSpecificEnrichmentSampler;
import Quantification.OverdispersionSampler;

/**
 * treatment and control group has different methylation level. For a single gene, the model parameter vector is consisted of 6 components.
 *  theta ~ (treatmentINPUTOverdispersion, treatmentIPOverdispersion,
 *           controlINPUTOverdispersion, controlIPOverdispersion,
 *           treatmentMethylationLevel, controlMethylationLevel)
 *
 * for a single gene, let D = {X_treatment_ip, X_treatment_input, X_control_ip, X_control_input} denote as the reads count
 * in IP and INPUT of treatment and control group, respectively.
 */
public class DiffMethylationLevelModel extends ModelSelection {

    public DiffMethylationLevelModel(MethylationLevelSampler tretMethylationLevelSampler, MethylationLevelSampler ctrlMethylationLevelSampler,
                                     NonSpecificEnrichmentSampler nonspecificEnrichSampler,
                                     OverdispersionSampler tretIPOverdispersionSampler, OverdispersionSampler tretINPUTOverdispersionSampler,
                                     OverdispersionSampler ctrlIPOverdispersionSampler, OverdispersionSampler ctrlINPUTOverdispersionSampler,
                                     BackgroundExpressionSampler treatmentBackgroundExpressionSampler, BackgroundExpressionSampler treatmentNonPeakExpressionSampler,
                                     BackgroundExpressionSampler controlBackgroundExpressionSampler, BackgroundExpressionSampler controlNonPeakExpressionSampler) {
        super(tretMethylationLevelSampler, ctrlMethylationLevelSampler, nonspecificEnrichSampler,
              tretIPOverdispersionSampler, tretINPUTOverdispersionSampler, ctrlIPOverdispersionSampler, ctrlINPUTOverdispersionSampler,
              treatmentBackgroundExpressionSampler, treatmentNonPeakExpressionSampler, controlBackgroundExpressionSampler, controlNonPeakExpressionSampler);
    }

    /**
     * MH sampling for treatment group methylation level
     * @param time iteration time
     */
    @Override
    protected double treatmentMethylationSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, nonspecificEnrich,
               tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp, newTretMethLevel;
        boolean samplingRes;
        double[] newIPReadsExpectation, tretIPSizeFactors;
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
        tretIPSizeFactors = this.parameters.getTretIPSizeFactors();
        newTretMethLevel = this.tretMethylationLevelSampler.randomSample(tretMethLevel);
        // renew IP and INPUT reads expectations via new sampling methylation level, shape individualNumber × geneNumber
        newIPReadsExpectation = this.renewReadsExpectationViaMethLevel(true, new double[]{newTretMethLevel}, new double[]{nonspecificEnrich}, tretIPSizeFactors);
        if (curModel) {
            logLikeProba = this.logLikelihood(newIPReadsExpectation, this.treatmentINPUTExpectations,
                                              this.controlIPExpectations, this.controlINPUTExpectations,
                                              this.treatmentIPNonPeakExpect, this.treatmentINPUTNonPeakExpect,
                                              this.controlIPNonPeakExpect, this.controlINPUTNonPeakExpect,
                                              tretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        newTretMethLevel, ctrlMethLevel, nonspecificEnrich,
                                        tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp);
            logCurParamsProba = logLikeProba + logPrior;

        } else {
            if (!usePseudoPrior) {
                logPrevParamsProba = this.tretMethylationLevelSampler.getLogDensity(tretMethLevel);
                logCurParamsProba = this.tretMethylationLevelSampler.getLogDensity(newTretMethLevel);
            } else {
                logPrevParamsProba = this.tretMethylationLevelPseudoSampler.getLogDensity(tretMethLevel);
                logCurParamsProba = this.tretMethylationLevelPseudoSampler.getLogDensity(newTretMethLevel);
            }
        }
        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        if (samplingRes) {
            this.parameters.setTretMethylation(newTretMethLevel);
            this.treatmentMethylationLevel[time] = newTretMethLevel;
            this.treatmentIPExpectations = newIPReadsExpectation;
            return logCurParamsProba;
        } else
            this.treatmentMethylationLevel[time] = tretMethLevel;

        return logPrevParamsProba;
    }

    /**
     * MH sampling for control group methylation level
     * @param time iteration time
     */
    @Override
    protected double controlMethylationSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, nonspecificEnrich,
               tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp, newCtrlMethLevel;
        boolean samplingRes;
        double[] newIPReadsExpectation, ctrlIPSizeFactors;
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
        ctrlIPSizeFactors = this.parameters.getCtrlIPSizeFactors();
        newCtrlMethLevel = this.ctrlMethylationLevelSampler.randomSample(ctrlMethLevel);
        // renew IP and INPUT reads expectations via new sampling methylation level, shape individualNumber × geneNumber
        newIPReadsExpectation = this.renewReadsExpectationViaMethLevel(false, new double[]{newCtrlMethLevel}, new double[]{nonspecificEnrich}, ctrlIPSizeFactors);
        if (curModel) {
            logLikeProba = this.logLikelihood(this.treatmentIPExpectations, this.treatmentINPUTExpectations,
                                              newIPReadsExpectation, this.controlINPUTExpectations,
                                              this.treatmentIPNonPeakExpect, this.treatmentINPUTNonPeakExpect,
                                              this.controlIPNonPeakExpect, this.controlINPUTNonPeakExpect,
                                              tretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        tretMethLevel, newCtrlMethLevel, nonspecificEnrich,
                                        tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp);
            logCurParamsProba = logLikeProba + logPrior;
        } else {
            if (!usePseudoPrior) {
                logPrevParamsProba = this.ctrlMethylationLevelSampler.getLogDensity(ctrlMethLevel);
                logCurParamsProba = this.ctrlMethylationLevelSampler.getLogDensity(newCtrlMethLevel);
            } else {
                logPrevParamsProba = this.ctrlMethylationLevelPseudoSampler.getLogDensity(ctrlMethLevel);
                logCurParamsProba = this.ctrlMethylationLevelPseudoSampler.getLogDensity(newCtrlMethLevel);
            }
        }
        // if samplingRes is true, means the new sampling value was accepted
        samplingRes = this.mhSampling.getSamplingRes(logCurParamsProba, logPrevParamsProba, true);
        if (samplingRes) {
            this.parameters.setCtrlMethylation(newCtrlMethLevel);
            this.controlMethylationLevel[time] = newCtrlMethLevel;
            this.controlIPExpectations = newIPReadsExpectation;
            return logCurParamsProba;
        } else
            this.controlMethylationLevel[time] = ctrlMethLevel;

        return logPrevParamsProba;
    }

    /**
     * model parameters prior probability in pseudo distribution
     * @return pseudo prior probability
     */
    @Override
    public double paramsPseudoDistributionProbability(boolean usePseudoPrior) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
                tretMethLevel, ctrlMethLevel, nonspecificEnrich;// , tretBkgExp, ctrlBkgExp, tretNonPeakExp, ctrlNonPeakExp;
        nonspecificEnrich = this.parameters.getNonspecificEnrichment();
        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        // tretBkgExp = this.parameters.getTretBkgExp();
        // ctrlBkgExp = this.parameters.getCtrlBkgExp();
        // tretNonPeakExp = this.parameters.getTretNonPeakBkgExp();
        // ctrlNonPeakExp = this.parameters.getCtrlNonPeakBkgExp();

        double proba, tretIPPseudo, tretINPUTPseudo, ctrlIPPseudo, ctrlINPUTPseudo, tretMethPseudo, ctrlMethPseudo,
               nonspecificPseudo; // , ctrlNonspecificPseudo, tretBkgExpPseudo, ctrlBkgExpPseudo, tretNonPeakExpPseudo, ctrlNonPeakExpPseudo;
        if (!usePseudoPrior) {
            nonspecificPseudo = this.nonspecificEnrichSampler.getLogDensity(nonspecificEnrich);
            tretMethPseudo = this.tretMethylationLevelSampler.getLogDensity(tretMethLevel);
            ctrlMethPseudo = this.ctrlMethylationLevelSampler.getLogDensity(ctrlMethLevel);
            tretIPPseudo = this.tretIPOverdispersionSampler.getLogDensity(tretIPOverdispersion);
            tretINPUTPseudo = this.tretINPUTOverdispersionSampler.getLogDensity(tretINPUTOverdispersion);
            ctrlIPPseudo = this.ctrlIPOverdispersionSampler.getLogDensity(ctrlIPOverdispersion);
            ctrlINPUTPseudo = this.ctrlINPUTOverdispersionSampler.getLogDensity(ctrlINPUTOverdispersion);
            // tretBkgExpPseudo = this.treatmentBackgroundExpressionSampler.getLogDensity(tretBkgExp);
            // ctrlBkgExpPseudo = this.controlBackgroundExpressionSampler.getLogDensity(ctrlBkgExp);
            // tretNonPeakExpPseudo = this.treatmentNonPeakExpressionSampler.getLogDensity(tretNonPeakExp);
            // ctrlNonPeakExpPseudo = this.treatmentNonPeakExpressionSampler.getLogDensity(ctrlNonPeakExp);
        } else {
            nonspecificPseudo = this.nonspecificEnrichPseudoSampler.getLogDensity(nonspecificEnrich);
            tretMethPseudo = this.tretMethylationLevelPseudoSampler.getLogDensity(tretMethLevel);
            ctrlMethPseudo = this.ctrlMethylationLevelPseudoSampler.getLogDensity(ctrlMethLevel);
            tretIPPseudo = this.tretIPOverdispersionPseudoSampler.getLogDensity(tretIPOverdispersion);
            tretINPUTPseudo = this.tretINPUTOverdispersionPseudoSampler.getLogDensity(tretINPUTOverdispersion);
            ctrlIPPseudo = this.ctrlIPOverdispersionPseudoSampler.getLogDensity(ctrlIPOverdispersion);
            ctrlINPUTPseudo = this.ctrlINPUTOverdispersionPseudoSampler.getLogDensity(ctrlINPUTOverdispersion);
            // tretBkgExpPseudo = this.treatmentBackgroundExpressionPseudoSampler.getLogDensity(tretBkgExp);
            // ctrlBkgExpPseudo = this.controlBackgroundExpressionPseudoSampler.getLogDensity(ctrlBkgExp);
            // tretNonPeakExpPseudo = this.treatmentNonPeakExpressionSampler.getLogDensity(tretNonPeakExp);
            // ctrlNonPeakExpPseudo = this.treatmentNonPeakExpressionSampler.getLogDensity(ctrlNonPeakExp);
        }

        proba = tretIPPseudo + tretINPUTPseudo + ctrlIPPseudo + ctrlINPUTPseudo + tretMethPseudo + ctrlMethPseudo +
                nonspecificPseudo;//  + ctrlNonspecificPseudo + tretBkgExpPseudo + ctrlBkgExpPseudo + tretNonPeakExpPseudo + ctrlNonPeakExpPseudo;

        return proba;
    }

    /**
     * parameters logarithm prior probability
     *      log-gamma(tret_ip_overdispersion) + log-gamma(tret_input_overdispersion) +
     *      log-gamma(ctrl_ip_overdispersion) + log-gamma(ctrl_input_overdispersion) +
     *      log-beta(tret_methylation_level) + log-beta(ctrl_methylation_level)
     * tret_methylation_level != ctrl_methylation_level
     * @param treatmentIPOverdispersion treatment group IP reads count overdispersion of a gene
     * @param treatmentINPUTOverdispersion treatment group INPUT reads count overdispersion of a gene
     * @param controlIPOverdispersion control group IP reads count overdispersion of a gene
     * @param controlINPUTOverdispersion control group INPUT reads count overdispersion of a gene
     * @param treatmentMethylationLevel methylation level of a treatment group gene
     * @param controlMethylationLevel methylation level of a control group gene
     * @return logarithm prior probability
     */
    @Override
    public double logPriority(double treatmentIPOverdispersion, double treatmentINPUTOverdispersion,
                              double controlIPOverdispersion, double controlINPUTOverdispersion,
                              double treatmentMethylationLevel, double controlMethylationLevel,
                              double nonspecificEnrich,
                              double treatmentBackgroundExpression, double controlBackgroundExpression,
                              double treatmentNonPeakExpression, double controlNonPeakExpression) {
        double tretIPOverdispersionProba = this.tretIPOverdispersionSampler.getLogDensity(treatmentIPOverdispersion);
        double tretINPUTOverdispersionProba = this.tretINPUTOverdispersionSampler.getLogDensity(treatmentINPUTOverdispersion);
        double ctrlIPOverdispersionProba = this.ctrlIPOverdispersionSampler.getLogDensity(controlIPOverdispersion);
        double ctrlINPUTOverdispersionProba= this.ctrlINPUTOverdispersionSampler.getLogDensity(controlINPUTOverdispersion);
        double tretMethProba = this.tretMethylationLevelSampler.getLogDensity(treatmentMethylationLevel);
        double ctrlMethProba = this.ctrlMethylationLevelSampler.getLogDensity(controlMethylationLevel);
        // double tretBkgExpProba = this.treatmentBackgroundExpressionSampler.getLogDensity(treatmentBackgroundExpression);
        // double ctrlBkgExpProba = this.controlBackgroundExpressionSampler.getLogDensity(controlBackgroundExpression);
        double nonspecificEnrichProba = this.nonspecificEnrichSampler.getLogDensity(nonspecificEnrich);
        // double tretNonPeakExpProba = this.treatmentNonPeakExpressionSampler.getLogDensity(treatmentNonPeakExpression);
        // double ctrlNonPeakExpProba = this.controlNonPeakExpressionSampler.getLogDensity(controlNonPeakExpression);

        double proba = tretINPUTOverdispersionProba + tretIPOverdispersionProba + ctrlIPOverdispersionProba + ctrlINPUTOverdispersionProba
                       + tretMethProba + ctrlMethProba + nonspecificEnrichProba;
                       // + ctrlNonspecificEnrichProba+ tretBkgExpProba + ctrlBkgExpProba  + tretNonPeakExpProba + ctrlNonPeakExpProba;

        return proba;
    }

    /**
     * use the first-time sampling result to regress pseudo prior for each parameters
     * Deprecated!
     */
    @Override
    protected void setModelParameterPseudoPrior() {
        if (tretMethyLationList == null) {
            this.cleanUpSamplingList();
            this.tretIPOverdispersionPseudoSampler = this.tretIPOverdispersionSampler;
            this.tretINPUTOverdispersionPseudoSampler = this.tretINPUTOverdispersionSampler;
            this.ctrlIPOverdispersionPseudoSampler = this.ctrlIPOverdispersionSampler;
            this.ctrlINPUTOverdispersionPseudoSampler = this.ctrlINPUTOverdispersionSampler;
            this.tretMethylationLevelPseudoSampler = this.tretMethylationLevelSampler;
            this.ctrlMethylationLevelPseudoSampler = this.ctrlMethylationLevelSampler;
            this.nonspecificEnrichPseudoSampler = this.nonspecificEnrichSampler;
            this.treatmentBackgroundExpressionPseudoSampler = this.treatmentBackgroundExpressionSampler;
            this.controlBackgroundExpressionPseudoSampler = this.controlBackgroundExpressionSampler;
            this.treatmentNonPeakExpressionPseudoSampler = this.treatmentNonPeakExpressionSampler;
            this.controlNonPeakExpressionPseudoSampler = this.controlNonPeakExpressionSampler;
        } else {
            double[] tretMethylationParams = PseudoBetaDistribution.estimate(this.tretMethyLationList);
            if (tretMethylationParams[0] < 0.00001 | tretMethylationParams[1] < 0.00001)
                this.tretMethylationLevelPseudoSampler = this.tretMethylationLevelSampler;
            else
                this.tretMethylationLevelPseudoSampler = new MethylationLevelSampler(tretMethylationParams[0], tretMethylationParams[1]);
            double[] ctrlMethylationParams = PseudoBetaDistribution.estimate(this.ctrlMethylationList);
            if (ctrlMethylationParams[0] < 0.00001 | ctrlMethylationParams[1] < 0.00001)
                this.ctrlMethylationLevelPseudoSampler = this.ctrlMethylationLevelSampler;
            else
                this.ctrlMethylationLevelPseudoSampler = new MethylationLevelSampler(ctrlMethylationParams[0], ctrlMethylationParams[1]);

            double[] nonspecificEnrichParams = PseudoBetaDistribution.estimate(this.nonspecificEnrichList);
            if (nonspecificEnrichParams[0] < 0.00001 | nonspecificEnrichParams[1] < 0.00001)
                this.nonspecificEnrichPseudoSampler = this.nonspecificEnrichSampler;
            else
                this.nonspecificEnrichPseudoSampler = new NonSpecificEnrichmentSampler(nonspecificEnrichParams[0], nonspecificEnrichParams[1]);

            double[] tretIPOverdispersionParams = PseudoInverseGammaDistribution.estimate(this.tretIPOverdispersionList);
            if (tretIPOverdispersionParams == null)
                this.tretIPOverdispersionPseudoSampler = this.tretIPOverdispersionSampler;
            else
                this.tretIPOverdispersionPseudoSampler = new OverdispersionSampler(tretIPOverdispersionParams[0], tretIPOverdispersionParams[1]);
            double[] tretINPUTOverdispersionParams = PseudoInverseGammaDistribution.estimate(this.tretINPUTOverdispersionList);
            if (tretINPUTOverdispersionParams == null)
                this.tretINPUTOverdispersionPseudoSampler = this.tretINPUTOverdispersionSampler;
            else
                this.tretINPUTOverdispersionPseudoSampler = new OverdispersionSampler(tretINPUTOverdispersionParams[0], tretINPUTOverdispersionParams[1]);

            double[] ctrlIPOverdispersionParams = PseudoInverseGammaDistribution.estimate(this.ctrlIPOverdispersionList);
            if (ctrlIPOverdispersionParams == null)
                this.ctrlIPOverdispersionPseudoSampler = this.ctrlIPOverdispersionSampler;
            else
                this.ctrlIPOverdispersionPseudoSampler = new OverdispersionSampler(ctrlIPOverdispersionParams[0], ctrlIPOverdispersionParams[1]);
            double[] ctrlINPUTOverdispersionParams = PseudoInverseGammaDistribution.estimate(this.ctrlINPUTOverdispersionList);
            if (ctrlINPUTOverdispersionParams == null)
                this.ctrlINPUTOverdispersionPseudoSampler = this.ctrlINPUTOverdispersionSampler;
            else
                this.ctrlINPUTOverdispersionPseudoSampler = new OverdispersionSampler(ctrlINPUTOverdispersionParams[0], ctrlINPUTOverdispersionParams[1]);

            double[] tretBkgExpParams = PseudoLogNormalDistribution.estimate(this.tretBkgExpList);
            if (tretBkgExpParams[0] < 0.00001 | tretBkgExpParams[1] < 0.00001)
                this.treatmentBackgroundExpressionPseudoSampler = this.treatmentBackgroundExpressionSampler;
            else
                this.treatmentBackgroundExpressionPseudoSampler = new BackgroundExpressionSampler(tretBkgExpParams[0], tretBkgExpParams[1]);
            double[] ctrlBkgExpParams = PseudoLogNormalDistribution.estimate(this.ctrlBkgExpList);
            if (ctrlBkgExpParams[0] < 0.00001 | ctrlBkgExpParams[1] < 0.00001)
                this.controlBackgroundExpressionPseudoSampler = this.controlBackgroundExpressionSampler;
            else
                this.controlBackgroundExpressionPseudoSampler = new BackgroundExpressionSampler(ctrlBkgExpParams[0], ctrlBkgExpParams[1]);

            double[] tretNonPeakExpParams = PseudoLogNormalDistribution.estimate(this.tretNonPeakExpList);
            if (tretNonPeakExpParams[0] < 0.00001 | tretNonPeakExpParams[1] < 0.00001)
                this.treatmentNonPeakExpressionPseudoSampler = this.treatmentNonPeakExpressionSampler;
            else
                this.treatmentNonPeakExpressionPseudoSampler = new BackgroundExpressionSampler(tretNonPeakExpParams[0], tretNonPeakExpParams[1]);
            double[] ctrlNonPeakExpParams = PseudoLogNormalDistribution.estimate(this.ctrlNonPeakExpList);
            if (ctrlNonPeakExpParams[0] < 0.00001 | ctrlNonPeakExpParams[1] < 0.00001)
                this.controlNonPeakExpressionPseudoSampler = this.controlNonPeakExpressionSampler;
            else
                this.controlNonPeakExpressionPseudoSampler = new BackgroundExpressionSampler(ctrlNonPeakExpParams[0], ctrlNonPeakExpParams[1]);

            this.cleanUpSamplingList();
        }
    }
}
