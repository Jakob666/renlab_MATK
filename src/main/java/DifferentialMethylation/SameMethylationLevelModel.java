package DifferentialMethylation;

import DifferentialMethylation.PseudoPrior.PseudoBetaDistribution;
import DifferentialMethylation.PseudoPrior.PseudoInverseGammaDistribution;
import DifferentialMethylation.PseudoPrior.PseudoLogNormalDistribution;
import Quantification.BackgroundExpressionSampler;
import Quantification.MethylationLevelSampler;
import Quantification.OverdispersionSampler;

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
                                     BackgroundExpressionSampler treatmentBackgroundExpressionSampler,
                                     BackgroundExpressionSampler controlBackgroundExpressionSampler) {
        super(tretMethylationLevelSampler, ctrlMethylationLevelSampler, tretIPOverdispersionSampler, tretINPUTOverdispersionSampler,
              ctrlIPOverdispersionSampler, ctrlINPUTOverdispersionSampler, treatmentBackgroundExpressionSampler, controlBackgroundExpressionSampler);
    }

    /**
     * MH sampling for treatment group methylation level, in these model we assume treatment and control group
     * have same methylation levels
     * @param time iteration time
     */
    @Override
    protected double treatmentMethylationSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, ctrlMethLevel, tretBkgExp, ctrlBkgExp, newTretMethLevel;
        boolean samplingRes;
        double[] newTretIPReadsExpectation, newCtrlIPReadsExpectation;
        double logCurParamsProba, logLikeProba, logPrior;

        tretMethLevel = this.parameters.getTretMethylation();
        ctrlMethLevel = this.parameters.getCtrlMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        assert Math.abs(tretMethLevel-ctrlMethLevel)<0.00001;
        newTretMethLevel = this.tretMethylationLevelSampler.randomSample(tretMethLevel);
        // renew IP and INPUT reads expectations via new sampling methylation level, shape individualNumber × geneNumber
        newTretIPReadsExpectation = this.renewReadsExpectationViaMethLevel(true, new double[]{newTretMethLevel});
        newCtrlIPReadsExpectation = this.renewReadsExpectationViaMethLevel(false, new double[]{newTretMethLevel});
        if (curModel) {
            logLikeProba = this.logLikelihood(newTretIPReadsExpectation, this.treatmentINPUTExpectations,
                                              newCtrlIPReadsExpectation, this.controlINPUTExpectations,
                                              tretIPOverdispersion, tretINPUTOverdispersion,
                                              ctrlIPOverdispersion, ctrlINPUTOverdispersion);
            logPrior = this.logPriority(tretIPOverdispersion, tretINPUTOverdispersion,
                                        ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                                        newTretMethLevel, newTretMethLevel, tretBkgExp, ctrlBkgExp);
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
            this.treatmentMethylationLevel[time] = newTretMethLevel;
            this.treatmentIPExpectations = newTretIPReadsExpectation;
            this.controlIPExpectations = newCtrlIPReadsExpectation;
            this.parameters.setTretMethylation(newTretMethLevel);
            return logCurParamsProba;
        } else
            this.treatmentMethylationLevel[time] = tretMethLevel;

        return logPrevParamsProba;
    }

    /**
     * MH sampling for control group methylation level, in these model we assume treatment and control group
     * have same methylation levels
     * @param time iteration time
     */
    @Override
    protected double controlMethylationSampling(int time, double logPrevParamsProba, boolean curModel, boolean usePseudoPrior) {
        this.controlMethylationLevel[time] = this.treatmentMethylationLevel[time];
        this.parameters.setCtrlMethylation(this.parameters.getTretMethylation());

        return logPrevParamsProba;
    }

    /**
     * model parameters prior probability in pseudo distribution
     * @return pseudo prior probability
     */
    @Override
    public double paramsPseudoDistributionProbability(boolean usePseudoPrior) {
        double tretINPUTOverdispersion, tretIPOverdispersion, ctrlINPUTOverdispersion, ctrlIPOverdispersion,
               tretMethLevel, tretBkgExp, ctrlBkgExp;
        tretMethLevel = this.parameters.getTretMethylation();
        tretINPUTOverdispersion = this.parameters.getTretINPUTOverdispersion();
        tretIPOverdispersion = this.parameters.getTretIPOverdispersion();
        ctrlINPUTOverdispersion = this.parameters.getCtrlINPUTOverdispersion();
        ctrlIPOverdispersion = this.parameters.getCtrlIPOverdispersion();
        tretBkgExp = this.parameters.getTretBkgExp();
        ctrlBkgExp = this.parameters.getCtrlBkgExp();
        double proba, tretIPPseudo, tretINPUTPseudo, ctrlIPPseudo, ctrlINPUTPseudo, tretMethPseudo, tretBkgExpPseudo, ctrlBkgExpPseudo;

        if (!usePseudoPrior) {
            tretIPPseudo = this.tretIPOverdispersionSampler.getLogDensity(tretIPOverdispersion);
            tretINPUTPseudo = this.tretINPUTOverdispersionSampler.getLogDensity(tretINPUTOverdispersion);
            ctrlIPPseudo = this.ctrlIPOverdispersionSampler.getLogDensity(ctrlIPOverdispersion);
            ctrlINPUTPseudo = this.ctrlINPUTOverdispersionSampler.getLogDensity(ctrlINPUTOverdispersion);
            tretBkgExpPseudo = this.treatmentBackgroundExpressionSampler.getLogDensity(tretBkgExp);
            ctrlBkgExpPseudo = this.controlBackgroundExpressionSampler.getLogDensity(ctrlBkgExp);
            tretMethPseudo = this.tretMethylationLevelSampler.getLogDensity(tretMethLevel);
        } else {
            tretIPPseudo = this.tretIPOverdispersionPseudoSampler.getLogDensity(tretIPOverdispersion);
            tretINPUTPseudo = this.tretINPUTOverdispersionPseudoSampler.getLogDensity(tretINPUTOverdispersion);
            ctrlIPPseudo = this.ctrlIPOverdispersionPseudoSampler.getLogDensity(ctrlIPOverdispersion);
            ctrlINPUTPseudo = this.ctrlINPUTOverdispersionPseudoSampler.getLogDensity(ctrlINPUTOverdispersion);
            tretBkgExpPseudo = this.treatmentBackgroundExpressionPseudoSampler.getLogDensity(tretBkgExp);
            ctrlBkgExpPseudo = this.controlBackgroundExpressionPseudoSampler.getLogDensity(ctrlBkgExp);
            tretMethPseudo = this.tretMethylationLevelPseudoSampler.getLogDensity(tretMethLevel);
        }

        proba = tretIPPseudo + tretINPUTPseudo + ctrlIPPseudo + ctrlINPUTPseudo + tretMethPseudo + tretBkgExpPseudo + ctrlBkgExpPseudo;

        return proba;
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
    protected double logPriority(double treatmentIPOverdispersion, double treatmentINPUTOverdispersion,
                                 double controlIPOverdispersion, double controlINPUTOverdispersion,
                                 double treatmentMethylationLevel, double controlMethylationLevel,
                                 double treatmentBackgroundExpression, double controlBackgroundExpression) {
        double proba = 0;
        double tretIPOverdispersionProba = this.tretIPOverdispersionSampler.getLogDensity(treatmentIPOverdispersion);
        double tretINPUTOverdispersionProba = this.tretINPUTOverdispersionSampler.getLogDensity(treatmentINPUTOverdispersion);
        double ctrlIPOverdispersionProba = this.ctrlIPOverdispersionSampler.getLogDensity(controlIPOverdispersion);
        double ctrlINPUTOverdispersionProba= this.ctrlINPUTOverdispersionSampler.getLogDensity(controlINPUTOverdispersion);
        double tretMethProba = this.tretMethylationLevelSampler.getLogDensity(treatmentMethylationLevel);
        double tretBkgExpProba = this.treatmentBackgroundExpressionSampler.getLogDensity(treatmentBackgroundExpression);
        double ctrlBkgExpProba = this.controlBackgroundExpressionSampler.getLogDensity(controlBackgroundExpression);

        proba += tretINPUTOverdispersionProba + tretIPOverdispersionProba + ctrlIPOverdispersionProba + ctrlINPUTOverdispersionProba
                + tretMethProba + tretBkgExpProba + ctrlBkgExpProba;
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
            this.treatmentBackgroundExpressionPseudoSampler = this.treatmentBackgroundExpressionSampler;
            this.controlBackgroundExpressionPseudoSampler = this.controlBackgroundExpressionSampler;
        } else {
            double[] methylationParams = PseudoBetaDistribution.estimate(this.tretMethyLationList);
            if (methylationParams[0] < 0.00001 | methylationParams[1] < 0.00001) {
                this.tretMethylationLevelPseudoSampler = this.tretMethylationLevelSampler;
                this.ctrlMethylationLevelPseudoSampler = this.ctrlMethylationLevelSampler;
            } else {
                this.tretMethylationLevelPseudoSampler = new MethylationLevelSampler(methylationParams[0], methylationParams[1]);
                this.ctrlMethylationLevelPseudoSampler = new MethylationLevelSampler(methylationParams[0], methylationParams[1]);
            }
            // TODO: 测试
//            double[] tretIPOverdispersionParams = PseudoGammaDistribution.estimate(this.tretIPOverdispersionList);
//            if (tretIPOverdispersionParams[0] < 0.00001 | tretIPOverdispersionParams[1] < 0.00001)
//                this.tretIPOverdispersionPseudoSampler = this.tretIPOverdispersionSampler;
//            else
//                this.tretIPOverdispersionPseudoSampler = new OverdispersionSampler(tretIPOverdispersionParams[0], tretIPOverdispersionParams[1]);
//            double[] tretINPUTOverdispersionParams = PseudoGammaDistribution.estimate(this.tretINPUTOverdispersionList);
//            if (tretINPUTOverdispersionParams[0] < 0.00001 | tretINPUTOverdispersionParams[1] < 0.00001)
//                this.tretINPUTOverdispersionPseudoSampler = this.tretINPUTOverdispersionSampler;
//            else
//                this.tretINPUTOverdispersionPseudoSampler = new OverdispersionSampler(tretINPUTOverdispersionParams[0], tretINPUTOverdispersionParams[1]);
//
//            double[] ctrlIPOverdispersionParams = PseudoGammaDistribution.estimate(this.ctrlIPOverdispersionList);
//            if (ctrlIPOverdispersionParams[0] < 0.00001 | ctrlIPOverdispersionParams[1] < 0.00001)
//                this.ctrlIPOverdispersionPseudoSampler = this.ctrlIPOverdispersionSampler;
//            else
//                this.ctrlIPOverdispersionPseudoSampler = new OverdispersionSampler(ctrlIPOverdispersionParams[0], ctrlIPOverdispersionParams[1]);
//            double[] ctrlINPUTOverdispersionParams = PseudoGammaDistribution.estimate(this.ctrlINPUTOverdispersionList);
//            if (ctrlINPUTOverdispersionParams[0] < 0.00001 | ctrlINPUTOverdispersionParams[1] < 0.00001)
//                this.ctrlINPUTOverdispersionPseudoSampler = this.ctrlINPUTOverdispersionSampler;
//            else
//                this.ctrlINPUTOverdispersionPseudoSampler = new OverdispersionSampler(ctrlINPUTOverdispersionParams[0], ctrlINPUTOverdispersionParams[1]);
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

            this.cleanUpSamplingList();
        }
    }
}
