package LaplaceMetropolisEstimator;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Arrays;

/**
 *  covariance matrix V
 *      V = average {(θj - θ_bar)^T * (θj - θ_bar)} , where θ_bar = average(θj) and j={1,2,3,...,N}
 *  reference to http://hedibert.org/wp-content/uploads/2013/12/bayesianmodelcriticism.pdf?tdsourcetag=s_pctim_aiomsg
 */
public class CovarianceMatrix {
    private RealMatrix covarianceMatrix;
    private double[] treatmentIPOverdispersion, treatmentINPUTOverdispersion, controlIPOverdispersion, controlINPUTOverdispersion,
                     treatmentMethylationLevel, controlMethylationLevel, nonspecificEnrichment;
    private double[] averageParams;
    private int samplingTime, burnIn, paramNum;

    public CovarianceMatrix(double[] treatmentIPOverdispersion, double[] treatmentINPUTOverdispersion,
                            double[] controlIPOverdispersion, double[] controlINPUTOverdispersion,
                            double[] treatmentMethylationLevel, double[] controlMethylationLevel, double[] nonspecificEnrichment,
                            int burnIn) {
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
    }

    public double[][] getCovarianceMatrix() {
        this.calcAverageParams();
        this.calcCovarianceMatrix();

        return this.covarianceMatrix.getData();
    }

    /**
     * calculate θ_bar, shape 1 × paramNum
     */
    private void calcAverageParams() {
        double averageTretMeth = this.calcAverage(this.treatmentMethylationLevel);
        double averageCtrlMeth;
        if (this.paramNum == 7)
            averageCtrlMeth = this.calcAverage(this.controlMethylationLevel);
        else
            averageCtrlMeth = Double.MIN_VALUE;
        double averageNonspecificEnrich = this.calcAverage(this.nonspecificEnrichment);
        double averageTretIPOverdispersion = this.calcAverage(this.treatmentIPOverdispersion);
        double averageCtrlIPOverdispersion = this.calcAverage(this.controlIPOverdispersion);
        double averageTretINPUTOverdispersion = this.calcAverage(this.treatmentINPUTOverdispersion);
        double averageCtrlINPUTOverdispersion = this.calcAverage(this.controlINPUTOverdispersion);

        if (this.paramNum == 7)
            this.averageParams = new double[] {averageTretMeth, averageCtrlMeth, averageNonspecificEnrich,
                                               averageTretIPOverdispersion, averageTretINPUTOverdispersion,
                                               averageCtrlIPOverdispersion, averageCtrlINPUTOverdispersion};
        else
            this.averageParams = new double[] {averageTretMeth, 0, averageNonspecificEnrich,
                                               averageTretIPOverdispersion, averageTretINPUTOverdispersion,
                                               averageCtrlIPOverdispersion, averageCtrlINPUTOverdispersion};
    }

    /**
     * calculate average of sampling numbers
     * @param samplingValues sampling numbers
     * @return average
     */
    private double calcAverage(double[] samplingValues) {
        return Arrays.stream(samplingValues).skip(this.burnIn).average().getAsDouble();
    }

    private void calcCovarianceMatrix() {
        int keepTime = this.samplingTime - this.burnIn;
        double[][] matrix = new double[keepTime][this.paramNum];   // shape keepNum × paramNum
        for (int t=0; t<this.samplingTime; t++) {
            if (t<this.burnIn)
                continue;
            calcMatrix(matrix, t-this.burnIn);
        }
        Array2DRowRealMatrix samples = new Array2DRowRealMatrix(matrix);
        RealMatrix transSamples = samples.transpose();
        this.covarianceMatrix = transSamples.multiply(samples).scalarMultiply(Math.pow(keepTime, -1));
    }

    private void calcMatrix(double[][] matrix, int time) {
        double tempTretMeth, tempCtrlMeth, tempNonspecificEnrich,
               tempTretIPOverdispersion, tempTretINPUTOverdispersion,
               tempCtrlIPOverdispersion, tempCtrlINPUTOverdispersion;

        tempTretMeth = this.treatmentMethylationLevel[time] - this.averageParams[0];
        if (this.paramNum == 7)
            tempCtrlMeth = this.controlMethylationLevel[time] - this.averageParams[1];
        else
            tempCtrlMeth = 0;
        tempNonspecificEnrich = this.nonspecificEnrichment[time] - this.averageParams[2];
        tempTretIPOverdispersion = this.treatmentIPOverdispersion[time] - this.averageParams[3];
        tempTretINPUTOverdispersion = this.treatmentINPUTOverdispersion[time] - this.averageParams[4];
        tempCtrlIPOverdispersion = this.controlIPOverdispersion[time] - this.averageParams[5];
        tempCtrlINPUTOverdispersion = this.controlINPUTOverdispersion[time] - this.averageParams[6];

        if (this.paramNum == 7)
            matrix[time] = new double[] {tempTretMeth, tempCtrlMeth, tempNonspecificEnrich,
                                         tempTretIPOverdispersion, tempTretINPUTOverdispersion,
                                         tempCtrlIPOverdispersion, tempCtrlINPUTOverdispersion};
        else
            matrix[time] = new double[] {tempTretMeth, tempNonspecificEnrich,
                                         tempTretIPOverdispersion, tempTretINPUTOverdispersion,
                                         tempCtrlIPOverdispersion, tempCtrlINPUTOverdispersion};
    }
}
