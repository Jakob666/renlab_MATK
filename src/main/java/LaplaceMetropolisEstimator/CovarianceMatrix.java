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
                     treatmentMethylationLevel, controlMethylationLevel, treatmentExpression, controlExpression;
    private double[] averageParams;
    private int samplingTime, burnIn, paramNum;

    /**
     * Constructor
     * @param treatmentIPOverdispersion treatment IP overdispersion sampling, shape 1 × samplingTime
     * @param treatmentINPUTOverdispersion treatment INPUT overdispersion sampling, shape 1 × samplingTime
     * @param controlIPOverdispersion control IP overdispersion sampling, shape 1 × samplingTime
     * @param controlINPUTOverdispersion control INPUT overdispersion sampling, shape 1 × samplingTime
     * @param treatmentMethylationLevel treatment IP methylation level sampling, shape 1 × samplingTime
     * @param controlMethylationLevel control IP methylation level sampling, shape 1 × samplingTime
     * @param treatmentExpression treatment background expression sampling, shape 1 × samplingTime
     * @param controlExpression control background expression sampling, shape 1 × samplingTime
     * @param burnIn burn-in time
     */
    public CovarianceMatrix(double[] treatmentIPOverdispersion, double[] treatmentINPUTOverdispersion,
                            double[] controlIPOverdispersion, double[] controlINPUTOverdispersion,
                            double[] treatmentMethylationLevel, double[] controlMethylationLevel,
                            double[] treatmentExpression, double[] controlExpression,
                            int burnIn) {
        this.treatmentIPOverdispersion = treatmentIPOverdispersion;
        this.treatmentINPUTOverdispersion = treatmentINPUTOverdispersion;
        this.controlIPOverdispersion = controlIPOverdispersion;
        this.controlINPUTOverdispersion = controlINPUTOverdispersion;
        this.treatmentMethylationLevel = treatmentMethylationLevel;
        this.controlMethylationLevel = controlMethylationLevel;
        this.treatmentExpression = treatmentExpression;
        this.controlExpression = controlExpression;
        this.samplingTime = treatmentIPOverdispersion.length;
        this.burnIn = burnIn;
        if (controlMethylationLevel == null)
            this.paramNum = 7;      // same methylation level model parameter number
        else
            this.paramNum = 8;      // diff methylation level model parameter number
    }

    /**
     * get model covariance matrix, shape paramNum × paramNum
     * @return covariance matrix
     */
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
        if (this.paramNum == 8)
            averageCtrlMeth = this.calcAverage(this.controlMethylationLevel);
        else
            averageCtrlMeth = Double.MIN_VALUE;
        double averageTretExp = this.calcAverage(this.treatmentExpression);
        double averageCtrlExp = this.calcAverage(this.controlExpression);
        double averageTretIPOverdispersion = this.calcAverage(this.treatmentIPOverdispersion);
        double averageCtrlIPOverdispersion = this.calcAverage(this.controlIPOverdispersion);
        double averageTretINPUTOverdispersion = this.calcAverage(this.treatmentINPUTOverdispersion);
        double averageCtrlINPUTOverdispersion = this.calcAverage(this.controlINPUTOverdispersion);

        if (this.paramNum == 7)
            this.averageParams = new double[] {averageTretMeth, averageCtrlMeth, averageTretExp, averageCtrlExp,
                                               averageTretIPOverdispersion, averageTretINPUTOverdispersion,
                                               averageCtrlIPOverdispersion, averageCtrlINPUTOverdispersion};
        else
            this.averageParams = new double[] {averageTretMeth, 0, averageTretExp, averageCtrlExp,
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

    /**
     * calculate elements in covariance matrix, shape paramNum × paramNum
     */
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

    /**
     * calculate matrix element
     * @param matrix matrix
     * @param time iteration time
     */
    private void calcMatrix(double[][] matrix, int time) {
        double tempTretMeth, tempCtrlMeth, tempTretExp, tempCtrlExp,
               tempTretIPOverdispersion, tempTretINPUTOverdispersion,
               tempCtrlIPOverdispersion, tempCtrlINPUTOverdispersion;

        tempTretMeth = this.treatmentMethylationLevel[time] - this.averageParams[0];
        if (this.paramNum == 7)
            tempCtrlMeth = this.controlMethylationLevel[time] - this.averageParams[1];
        else
            tempCtrlMeth = 0;
        tempTretExp = this.treatmentExpression[time] - this.averageParams[2];
        tempCtrlExp = this.controlExpression[time] - this.averageParams[3];
        tempTretIPOverdispersion = this.treatmentIPOverdispersion[time] - this.averageParams[3];
        tempTretINPUTOverdispersion = this.treatmentINPUTOverdispersion[time] - this.averageParams[4];
        tempCtrlIPOverdispersion = this.controlIPOverdispersion[time] - this.averageParams[5];
        tempCtrlINPUTOverdispersion = this.controlINPUTOverdispersion[time] - this.averageParams[6];

        if (this.paramNum == 8)
            matrix[time] = new double[] {tempTretMeth, tempCtrlMeth, tempTretExp, tempCtrlExp,
                                         tempTretIPOverdispersion, tempTretINPUTOverdispersion,
                                         tempCtrlIPOverdispersion, tempCtrlINPUTOverdispersion};
        else
            matrix[time] = new double[] {tempTretMeth, tempTretExp, tempCtrlExp,
                                         tempTretIPOverdispersion, tempTretINPUTOverdispersion,
                                         tempCtrlIPOverdispersion, tempCtrlINPUTOverdispersion};
    }
}
