package HessianMatrix;

import Derivative.*;

public class HessianMatrix {
    // model quantification result
    private double tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                   tretMethylation, ctrlMethylation, tretBackgroundExpression, ctrlBackgroundExpression,
                   tretNonPeakExpression, ctrlNonPeakExpression, nonspecificEnrich;
    // prior distribution parameters
    private double gammaShape, betaScale1, betaScale2;
    // shape 1 × individualNumber
    private double[] tretSampleSizeFactor, ctrlSampleSizeFactor;
    // individual number
    private int tretIndividualNumber, ctrlIndividualNumber;
    // reads count
    private int[] tretIPReads, tretINPUTReads, tretIPNonPeak, tretINPUTNonPeak,
                  ctrlIPReads, ctrlINPUTReads, ctrlIPNonPeak, ctrlINPUTNonPeak;
    private IPLogNBDistributionDerivative tretIPLogNB, ctrlIPLogNB;
    private IPNonPeakLogNBDistributionDerivative tretIPNonPeakLogNB, ctrlIPNonPeakLogNB;
    private INPUTLogNBDistributionDerivative tretINPUTLogNB, ctrlINPUTLogNB, tretINPUTNonPeakLogNB, ctrlINPUTNonPeakLogNB;

    public HessianMatrix(double tretIPOverdispersion, double tretINPUTOverdispersion, double ctrlIPOverdispersion, double ctrlINPUTOverdispersion,
                         double tretMethylation, double ctrlMethylation, double tretBackgroundExpression, double ctrlBackgroundExpression,
                         double tretNonPeakExpression, double ctrlNonPeakExpression, double nonspecificEnrich,
                         double gammaShape, double betaScale1, double betaScale2) {
        this.tretIPOverdispersion = tretIPOverdispersion;
        this.tretINPUTOverdispersion = tretINPUTOverdispersion;
        this.ctrlIPOverdispersion = ctrlIPOverdispersion;
        this.ctrlINPUTOverdispersion = ctrlINPUTOverdispersion;
        this.tretMethylation = tretMethylation;
        this.ctrlMethylation = ctrlMethylation;
        this.tretBackgroundExpression = tretBackgroundExpression;
        this.ctrlBackgroundExpression = ctrlBackgroundExpression;
        this.tretNonPeakExpression = tretNonPeakExpression;
        this.ctrlNonPeakExpression = ctrlNonPeakExpression;
        this.nonspecificEnrich = nonspecificEnrich;
        this.gammaShape = gammaShape;
        this.betaScale1 = betaScale1;
        this.betaScale2 = betaScale2;
    }

    public void setReads(int[] tretIPReads, int[] tretINPUTReads, int[] tretIPNonPeak, int[] tretINPUTNonPeak,
                         int[] ctrlIPReads, int[] ctrlINPUTReads, int[] ctrlIPNonPeak, int[] ctrlINPUTNonPeak) {
        this.tretIPReads = tretIPReads;
        this.tretINPUTReads = tretINPUTReads;
        this.tretIPNonPeak = tretIPNonPeak;
        this.tretINPUTNonPeak = tretINPUTNonPeak;
        this.ctrlIPReads = ctrlIPReads;
        this.ctrlINPUTReads = ctrlINPUTReads;
        this.ctrlIPNonPeak = ctrlIPNonPeak;
        this.ctrlINPUTNonPeak = ctrlINPUTNonPeak;

        this.tretIndividualNumber = tretIPReads.length;
        this.ctrlIndividualNumber = ctrlIPReads.length;
    }

    public void setSampleSizeFactor(double[] tretSampleSizeFactor, double[] ctrlSampleSizeFactor) {
        assert this.tretIndividualNumber == tretSampleSizeFactor.length;
        assert this.ctrlIndividualNumber == ctrlSampleSizeFactor.length;
        this.tretSampleSizeFactor = tretSampleSizeFactor;
        this.ctrlSampleSizeFactor = ctrlSampleSizeFactor;
    }

    public void initializeDerivatives() {
        this.tretIPLogNB = new IPLogNBDistributionDerivative(this.tretBackgroundExpression, this.nonspecificEnrich, this.tretMethylation, this.tretIPOverdispersion);
        this.ctrlIPLogNB = new IPLogNBDistributionDerivative(this.ctrlBackgroundExpression, this.nonspecificEnrich, this.ctrlMethylation, this.ctrlIPOverdispersion);
        this.tretIPNonPeakLogNB = new IPNonPeakLogNBDistributionDerivative(this.tretNonPeakExpression, this.nonspecificEnrich, this.tretIPOverdispersion);
        this.ctrlIPNonPeakLogNB = new IPNonPeakLogNBDistributionDerivative(this.ctrlNonPeakExpression, this.nonspecificEnrich, this.ctrlIPOverdispersion);
        this.tretINPUTLogNB = new INPUTLogNBDistributionDerivative(this.tretBackgroundExpression, this.tretINPUTOverdispersion);
        this.ctrlINPUTLogNB = new INPUTLogNBDistributionDerivative(this.ctrlBackgroundExpression, this.ctrlINPUTOverdispersion);
        this.tretINPUTNonPeakLogNB = new INPUTLogNBDistributionDerivative(this.tretNonPeakExpression, this.tretINPUTOverdispersion);
        this.ctrlINPUTNonPeakLogNB = new INPUTLogNBDistributionDerivative(this.ctrlNonPeakExpression, this.ctrlINPUTOverdispersion);
    }

    /**
     * ∂^2 f / ∂ sigma^2, where f=log[likelihood * prior], sigma denote as treatment group IP overdispersion
     * @return second order derivative
     */
    public double secondDerivativeTreatmentIPOverdispersion() {
        double sizeFactor, peakReadsCount, nonPeakReadsCount, result = 0;
        for (int i=0; i<this.tretIndividualNumber; i++) {
            sizeFactor = this.tretSampleSizeFactor[i];
            peakReadsCount = this.tretIPReads[i];
            nonPeakReadsCount = this.tretIPNonPeak[i];
            // IP peak region:  ∂^2 logNB / ∂ sigma^2
            result += this.tretIPLogNB.logNBSecondOrderDerivativeOverdispersion(sizeFactor, peakReadsCount);
            // IP non-peak region:  ∂^2 logNB / ∂ sigma^2
            result += this.tretIPNonPeakLogNB.logNBSecondOrderDerivativeOverdispersion(sizeFactor, nonPeakReadsCount);
        }
        // ∂^2 log-Gamma / ∂ sigma^2
        result += LogGammaDistributionDerivative.logGammaSecondOrderDerivative(this.gammaShape, this.tretIPOverdispersion);

        return result;
    }

    /**
     * ∂^2 f / ∂ sigma^2, where f=log[likelihood * prior], sigma denote as control group IP overdispersion
     * @return second order derivative
     */
    public double secondDerivativeControlIPOverdispersion() {
        double sizeFactor, peakReadsCount, nonPeakReadsCount, result = 0;
        for (int i=0; i<this.ctrlIndividualNumber; i++) {
            sizeFactor = this.ctrlSampleSizeFactor[i];
            peakReadsCount = this.ctrlIPReads[i];
            nonPeakReadsCount = this.ctrlIPNonPeak[i];
            // IP peak region:  ∂^2 logNB / ∂ sigma^2
            result += this.ctrlIPLogNB.logNBSecondOrderDerivativeOverdispersion(sizeFactor, peakReadsCount);
            // IP non-peak region:  ∂^2 logNB / ∂ sigma^2
            result += this.ctrlIPNonPeakLogNB.logNBSecondOrderDerivativeOverdispersion(sizeFactor, nonPeakReadsCount);
        }
        // ∂^2 log-Gamma / ∂ sigma^2
        result += LogGammaDistributionDerivative.logGammaSecondOrderDerivative(this.gammaShape, this.ctrlIPOverdispersion);

        return result;
    }

    /**
     * ∂^2 f / ∂ sigma^2, where f=log[likelihood * prior], sigma denote as treatment group INPUT overdispersion
     * @return second order derivative
     */
    public double secondDerivativeTreatmentINPUTOverdispersion() {
        double sizeFactor, peakReadsCount, nonPeakReadsCount, result = 0;
        for (int i=0; i<this.tretIndividualNumber; i++) {
            sizeFactor = this.tretSampleSizeFactor[i];
            peakReadsCount = this.tretINPUTReads[i];
            nonPeakReadsCount = this.tretINPUTNonPeak[i];
            // INPUT peak region:  ∂^2 logNB / ∂ sigma^2
            result += this.tretINPUTLogNB.logNBSecondOrderDerivativeOverdispersion(sizeFactor, peakReadsCount);
            // INPUT non-peak region:  ∂^2 logNB / ∂ sigma^2
            result += this.tretINPUTNonPeakLogNB.logNBSecondOrderDerivativeOverdispersion(sizeFactor, nonPeakReadsCount);
        }
        // ∂^2 log-Gamma / ∂ sigma^2
        double secondDerivGamma = LogGammaDistributionDerivative.logGammaSecondOrderDerivative(this.gammaShape, this.tretINPUTOverdispersion);
        result += secondDerivGamma;

        return result;
    }

    /**
     * ∂^2 f / ∂ sigma^2, where f=log[likelihood * prior], sigma denote as control group INPUT overdispersion
     * @return second order derivative
     */
    public double secondDerivativeControlINPUTOverdispersion() {
        double sizeFactor, peakReadsCount, nonPeakReadsCount, result = 0;
        for (int i=0; i<this.ctrlIndividualNumber; i++) {
            sizeFactor = this.ctrlSampleSizeFactor[i];
            peakReadsCount = this.ctrlINPUTReads[i];
            nonPeakReadsCount = this.ctrlINPUTNonPeak[i];
            // INPUT peak region:  ∂^2 logNB / ∂ sigma^2
            result += this.ctrlINPUTLogNB.logNBSecondOrderDerivativeOverdispersion(sizeFactor, peakReadsCount);
            // INPUT non-peak region:  ∂^2 logNB / ∂ sigma^2
            result += this.ctrlINPUTNonPeakLogNB.logNBSecondOrderDerivativeOverdispersion(sizeFactor, nonPeakReadsCount);
        }
        // ∂^2 log-Gamma / ∂ sigma^2
        result += LogGammaDistributionDerivative.logGammaSecondOrderDerivative(this.gammaShape, this.ctrlINPUTOverdispersion);

        return result;
    }

    /**
     * ∂^2 f / ∂ pT^2, where f=log[likelihood * prior], pT denote as treatment group methylation level
     * @return second order derivative
     */
    public double secondDerivativeTreatmentMethylationLevel() {
        double sizeFactor, peakReadsCount, result = 0;
        for (int i=0; i<this.tretIndividualNumber; i++) {
            sizeFactor = this.tretSampleSizeFactor[i];
            peakReadsCount = this.tretIPReads[i];
            // peak region:  ∂^2 logNB / ∂ pT^2
            result += this.tretIPLogNB.logNBSecondOrderDerivativeMethylationLevel(sizeFactor, peakReadsCount);
        }
        // ∂^2 log-Beta / ∂ pT^2
        result += LogBetaDistributionDerivative.logBetaSecondOrderDerivative(this.betaScale1, this.betaScale2, this.tretMethylation);

        return result;
    }

    /**
     * ∂^2 f / ∂ pC^2, where f=log[likelihood * prior], pC denote as control group methylation level
     * @return second order derivative
     */
    public double secondDerivativeControlMethylationLevel() {
        double sizeFactor, peakReadsCount, result = 0;
        for (int i=0; i<this.ctrlIndividualNumber; i++) {
            sizeFactor = this.ctrlSampleSizeFactor[i];
            peakReadsCount = this.ctrlIPReads[i];
            // peak region:  ∂^2 logNB / ∂ pC^2
            result += this.ctrlIPLogNB.logNBSecondOrderDerivativeMethylationLevel(sizeFactor, peakReadsCount);
        }
        // ∂^2 log-Beta / ∂ pC^2
        result += LogBetaDistributionDerivative.logBetaSecondOrderDerivative(this.betaScale1, this.betaScale2, this.ctrlMethylation);

        return result;
    }

    /**
     * ∂^2 f / ∂ r^2, where f=log[likelihood * prior], r denote as treatment/control group nonspecific enrichment ratio
     * @return second order derivative
     */
    public double secondDerivativeNonspecificEnrichment() {
        double sizeFactor, peakReadsCount, nonPeakReadsCount, result = 0;
        for (int i=0; i<this.tretIndividualNumber; i++) {
            sizeFactor = this.tretSampleSizeFactor[i];
            peakReadsCount = this.tretIPReads[i];
            nonPeakReadsCount = this.tretIPNonPeak[i];
            // IP peak region ∂^2 logNB / ∂ r^2
            result += this.tretIPLogNB.logNBSecondOrderDerivativeNonspecificEnrich(sizeFactor, peakReadsCount);
            // IP non-peak region ∂^2 logNB / ∂ r^2
            result += this.tretIPNonPeakLogNB.logNBSecondOrderDerivativeNonspecificEnrich(sizeFactor, nonPeakReadsCount);
        }

        for (int i=0; i<this.ctrlIndividualNumber; i++) {
            sizeFactor = this.ctrlSampleSizeFactor[i];
            peakReadsCount = this.ctrlIPReads[i];
            nonPeakReadsCount = this.ctrlIPNonPeak[i];
            // IP peak region ∂^2 logNB / ∂ r^2
            result += this.ctrlIPLogNB.logNBSecondOrderDerivativeNonspecificEnrich(sizeFactor, peakReadsCount);
            // IP non-peak region ∂^2 logNB / ∂ r^2
            result += this.ctrlIPNonPeakLogNB.logNBSecondOrderDerivativeNonspecificEnrich(sizeFactor, nonPeakReadsCount);
        }
        // ∂^2 log-Beta / ∂ r^2
        result += LogBetaDistributionDerivative.logBetaSecondOrderDerivative(this.betaScale1, this.betaScale2, this.nonspecificEnrich);

        return result;
    }
}
