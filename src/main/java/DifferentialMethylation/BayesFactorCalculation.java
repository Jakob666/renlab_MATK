package DifferentialMethylation;

import HessianMatrix.DiffMethylationModelDeterminate;
import HessianMatrix.ModelDeterminate;
import HessianMatrix.SameMethylationModelDeterminate;

import java.math.BigDecimal;
import java.math.RoundingMode;

public class BayesFactorCalculation {
    private int model0ParamNumber, model1ParamNumber;
    private ModelDeterminate sameModelDeterminate, diffModelDeterminate;
    private double model0Likelihood, model0Prior, model1Likelihood, model1Prior;
    private double tretBackgroundExpression, ctrlBackgroundExpression, tretNonPeakExpression, ctrlNonPeakExpression,
                   gammaShape, betaScale1, betaScale2;
    private double[] model0Quantification, model1Quantification;

    public BayesFactorCalculation(double model0LogLikelihood, double model1LogLikelihood,
                                  double model0LogPrior, double model1LogPrior,
                                  double[] sameMethModelQuantification, double[] diffMethModelQuantification,
                                  double tretBackgroundExpression, double ctrlBackgroundExpression,
                                  double tretNonPeakExpression, double ctrlNonPeakExpression,
                                  double gammaShape, double betaScale1, double betaScale2) {
        this.model0Likelihood = model0LogLikelihood;
        this.model0Prior = model0LogPrior;
        this.model1Likelihood = model1LogLikelihood;
        this.model1Prior = model1LogPrior;
        this.model0Quantification = sameMethModelQuantification;
        this.model1Quantification = diffMethModelQuantification;
        this.model0ParamNumber = sameMethModelQuantification.length - 1;
        this.model1ParamNumber = diffMethModelQuantification.length;

        this.tretBackgroundExpression = tretBackgroundExpression;
        this.ctrlBackgroundExpression = ctrlBackgroundExpression;
        this.tretNonPeakExpression = tretNonPeakExpression;
        this.ctrlNonPeakExpression = ctrlNonPeakExpression;
        this.gammaShape = gammaShape;
        this.betaScale1 = betaScale1;
        this.betaScale2 = betaScale2;

        this.createModelDeterminate(true);
        this.createModelDeterminate(false);
    }

    public void createModelDeterminate(boolean sameMethLevel) {
        double tretMeth, ctrlMeth, nonspecificEnrich,
                tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion;
        if (sameMethLevel) {
            tretMeth = this.model0Quantification[0];
            nonspecificEnrich = this.model0Quantification[2];
            tretIPOverdispersion = this.model0Quantification[3];
            tretINPUTOverdispersion = this.model0Quantification[4];
            ctrlIPOverdispersion = this.model0Quantification[5];
            ctrlINPUTOverdispersion = this.model0Quantification[6];
            this.sameModelDeterminate = new SameMethylationModelDeterminate(tretIPOverdispersion, tretINPUTOverdispersion,
                                                                            ctrlIPOverdispersion, ctrlINPUTOverdispersion, tretMeth,
                                                                            this.tretBackgroundExpression, this.ctrlBackgroundExpression,
                                                                            this.tretNonPeakExpression, this.ctrlNonPeakExpression, nonspecificEnrich,
                                                                            this.gammaShape, this.betaScale1, this.betaScale2);
        } else {
            tretMeth = this.model1Quantification[0];
            ctrlMeth = this.model1Quantification[1];
            nonspecificEnrich = this.model1Quantification[2];
            tretIPOverdispersion = this.model1Quantification[3];
            tretINPUTOverdispersion = this.model1Quantification[4];
            ctrlIPOverdispersion = this.model1Quantification[5];
            ctrlINPUTOverdispersion = this.model1Quantification[6];
            this.diffModelDeterminate = new DiffMethylationModelDeterminate(tretIPOverdispersion, tretINPUTOverdispersion,
                                                                            ctrlIPOverdispersion, ctrlINPUTOverdispersion, tretMeth, ctrlMeth,
                                                                            this.tretBackgroundExpression, this.ctrlBackgroundExpression,
                                                                            this.tretNonPeakExpression, this.ctrlNonPeakExpression, nonspecificEnrich,
                                                                            this.gammaShape, this.betaScale1, this.betaScale2);
        }
    }

    public void setReads(int[] tretIPReads, int[] tretINPUTReads, int[] tretIPNonPeak, int[] tretINPUTNonPeak,
                         int[] ctrlIPReads, int[] ctrlINPUTReads, int[] ctrlIPNonPeak, int[] ctrlINPUTNonPeak) {
        this.sameModelDeterminate.setReads(tretIPReads, tretINPUTReads, tretIPNonPeak, tretINPUTNonPeak,
                                           ctrlIPReads, ctrlINPUTReads, ctrlIPNonPeak, ctrlINPUTNonPeak);
        this.diffModelDeterminate.setReads(tretIPReads, tretINPUTReads, tretIPNonPeak, tretINPUTNonPeak,
                                           ctrlIPReads, ctrlINPUTReads, ctrlIPNonPeak, ctrlINPUTNonPeak);
    }

    public void setSizeFactors(double[] tretSampleSizeFactor, double[] ctrlSampleSizeFactor) {
        this.sameModelDeterminate.setSizeFactors(tretSampleSizeFactor, ctrlSampleSizeFactor);
        this.diffModelDeterminate.setSizeFactors(tretSampleSizeFactor, ctrlSampleSizeFactor);
    }

    /**
     * logBF = model1LogLikelihood - model0LogLikelihood + model1LogPrior - model0LogPrior + (param1-param0)*log(2PI) + log|-H0| - log|-H1|
     * @return Bayes factor
     */
    public double calcBayesFactor() {
        BigDecimal model0Determinate = this.sameModelDeterminate.modelDeterminate();
        BigDecimal model1Determinate = this.diffModelDeterminate.modelDeterminate();

        double likelihood = Math.exp(this.model1Likelihood - this.model0Likelihood);
        double prior = Math.exp(this.model1Prior - this.model0Prior);
        double penalty = Math.pow(2*Math.PI, (this.model1ParamNumber - this.model0ParamNumber) * 0.5);
        BigDecimal determinate = model0Determinate.divide(model1Determinate, 60, RoundingMode.HALF_UP);

        double bayesFactor = likelihood * prior * penalty * Math.pow(Math.abs(determinate.doubleValue()), 0.5);

        return bayesFactor;
    }
}
