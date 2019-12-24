package Quantification;

import ProbabilityCalculation.ProbabilityCalculator;
import SeqDataModel.BackgroundExpression;
import SeqDataModel.ReadsExpectation;

import java.util.Arrays;

/**
 * Quantify the methylation level with component-wise MH sampling
 * reference https://theclevermachine.wordpress.com/2012/11/04/mcmc-multivariate-distributions-block-wise-component-wise-updates/
 */
public class MethylationQuantification {
    private MethylationLevelSampler methLevelSampler = null;
    private OverdispersionSampler ipOverdispersionSampler = null, inputOverdispersionSampler = null;
    private BackgroundExpressionSampler[] geneBackgroundExpressionSamplers = null;
    private int individualNumber, geneNumber, samplingTime, burnIn;
    private int[][] ipReads, inputReads;

    private double a, b, c, d, w, k;
    private double[] geneMethylationLevel, geneIPOverdispersion, geneINPUTOverdispersion, geneBackgroundExpression;
    private double[][] ipReadsExpectation, inputReadsExpectation, methLevel, ipOverdispersions, inputOverdispersions,
                      backgroundExpression;


    /**
     * Constructor
     * @param ipReads genes' IP reads of each individual, shape individualNumber × geneNumber
     * @param inputReads genes' INPUT reads of each individual, shape individualNumber × geneNumber
     * @param ipShape IP data overdispersion Gamma distribution shape parameter
     * @param ipScale IP data overdispersion Gamma distribution scale parameter
     * @param inputShape INPUT data overdispersion Gamma distribution shape parameter
     * @param inputScale INPUT data overdispersion Gamma distribution scale parameter
     * @param methLevelShape1 methylation level Beta distribution shape parameter
     * @param methLevelShape2 methylation level Beta distribution shape parameter
     * @param samplingTime MH sampling time
     * @param burnIn MH burn-in time
     */
    public MethylationQuantification(int[][] ipReads, int[][] inputReads,
                                     double ipShape, double ipScale, double inputShape, double inputScale,
                                     double methLevelShape1, double methLevelShape2, int samplingTime, int burnIn) {
        assert samplingTime > burnIn;
        // assert arrays contains same sample and gene number
        assert inputReads.length == ipReads.length && inputReads[0].length == ipReads[0].length;
        this.individualNumber = inputReads.length;
        this.geneNumber = inputReads[0].length;
        this.inputReads = inputReads;
        this.ipReads = ipReads;
        this.a = ipShape;
        this.b = ipScale;
        this.c = inputShape;
        this.d = inputScale;
        this.w = methLevelShape1;
        this.k = methLevelShape2;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
    }

    public void sampling() {
        this.initialize();
        this.generatePosteriorDistribution();
        this.quantify();
    }

    public double[] getGeneMethylationLevel() {
        return this.geneMethylationLevel;
    }

    public double[] getGeneIPOverdispersion() {
        return this.geneIPOverdispersion;
    }

    public double[] getGeneINPUTOverdispersion() {
        return this.geneINPUTOverdispersion;
    }

    public double[] getGeneBackgroundExpression() {
        return this.geneBackgroundExpression;
    }

    /**
     * MH sampling process
     */
    private void generatePosteriorDistribution() {
        int time = 0;
        while (time < this.samplingTime-1) {
            time += 1;
            // sampling methylation level
            this.methylationLevelSampling(time);
            // sampling INPUT data overdispersion
            this.inputOverdispersionSampling(time);
            // sampling IP data overdispersion
            this.ipOverdispersionSampling(time);
            // sampling background expression
            this.geneBackgroundExpressionSampling(time);
        }
    }

    /**
     * initialize sampling parameters with given parameters
     */
    public void initialize() {
        this.methLevelSampler = new MethylationLevelSampler(this.w, this.k);
        this.inputOverdispersionSampler = new OverdispersionSampler(this.c, this.d);
        this.ipOverdispersionSampler = new OverdispersionSampler(this.a, this.b);
        this.backgroundExpression = new double[this.samplingTime][this.geneNumber];
        BackgroundExpression be = new BackgroundExpression(this.ipReads, this.inputReads);
        // background expression expectation and standard deviation, shape 1 × geneNumber
        double[] backgroundExpressionMean = be.geneBackgroundExp();
        double[] backgroundExpressionStd = be.geneExpressionStd();

        // for each gene, initialize a background expression sampler.
        this.geneBackgroundExpressionSamplers = new BackgroundExpressionSampler[this.geneNumber];
        BackgroundExpressionSampler sampler;
        for (int i=0; i<this.geneNumber; i++) {
            // shape and scale parameters of log-normal distribution
            double scale = Math.log(backgroundExpressionMean[i]);
            double shape = Math.log(backgroundExpressionStd[i]);
            sampler = new BackgroundExpressionSampler(scale, shape);
            // random init background expression from corresponding distribution
            this.backgroundExpression[0][i] = backgroundExpressionMean[i];
            this.geneBackgroundExpressionSamplers[i] = sampler;
        }

        // methylation level, IP overdispersion, INPUT overdispersion of each gene, shape samplingTime × geneNumber
        this.methLevel = new double[this.samplingTime][this.geneNumber];
        this.ipOverdispersions = new double[this.samplingTime][this.geneNumber];
        this.inputOverdispersions = new double[this.samplingTime][this.geneNumber];
        for (int i=0; i<this.geneNumber; i++) {
//            this.methLevel[0][i] = this.methLevelSampler.randomInit();
            this.methLevel[0][i] = 0.1;
//            this.ipOverdispersions[0][i] = this.ipOverdispersionSampler.randomInit();
            this.ipOverdispersions[0][i] = 1;
//            this.inputOverdispersions[0][i] = this.inputOverdispersionSampler.randomInit();
            this.inputOverdispersions[0][i] = 1;
        }

        // initialize genes' IP and INPUT reads expectation of each individual with methylation level,
        // shape individualNumber × geneNumber. INPUT reads expectations will not be changed any more,
        // but IP reads change with different methylation level
        ReadsExpectation re = new ReadsExpectation(this.ipReads, this.inputReads, this.methLevel[0]);
        this.ipReadsExpectation = re.getIPReadsExepectation();
        this.inputReadsExpectation = re.getINPUTReadsExpectation();
        re = null;
    }

    /**
     * new sampling iteration of Methylation level for each gene
     */
    private void methylationLevelSampling(int time) {

        double curMethLevel, curPosteriorProba, prevMethLevel, prevPosteriorProba;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation, prevIPExpectation, prevINPUTExpectation;
        double ipOverdispersion, inputOverdispersion;
        boolean samplingRes;

        // genes' new methylation level
        double[] newMethLevel = new double[this.geneNumber];
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            newMethLevel[geneIdx] = this.methLevelSampler.randomSample(this.methLevel[time-1][geneIdx]);
        }

        // use the new methylation level to calculate IP and INPUT reads count expectations
        // shape individualNumber × geneNumber
        double[][] newIPReadsExpectation;
        ReadsExpectation re = new ReadsExpectation(this.ipReads, this.inputReads, newMethLevel);
        newIPReadsExpectation = re.getIPReadsExepectation();
        re = null;

        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            prevMethLevel = this.methLevel[time-1][geneIdx];
            curMethLevel = newMethLevel[geneIdx];

            // get IP and INPUT gene reads count expectation of each individual, shape 1 × individualNumber
            readsData = this.getReadsDataForGene(this.ipReads, this.inputReads, this.ipReadsExpectation, this.inputReadsExpectation, geneIdx);
            prevIPExpectation = readsData[2];
            prevINPUTExpectation = readsData[3];

            readsData = this.getReadsDataForGene(this.ipReads, this.inputReads, newIPReadsExpectation, this.inputReadsExpectation, geneIdx);
            ipCount = readsData[0];
            inputCount = readsData[1];
            ipExpectation = readsData[2];
            inputExpectation = readsData[3];
            ipOverdispersion = this.ipOverdispersions[time-1][geneIdx];
            inputOverdispersion = this.inputOverdispersions[time-1][geneIdx];

            BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[geneIdx];
            double curGeneBkgExp = this.backgroundExpression[time-1][geneIdx];
            prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, prevIPExpectation, prevINPUTExpectation,
                    ipOverdispersion, inputOverdispersion, prevMethLevel, sampler, curGeneBkgExp);
            curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                    ipOverdispersion, inputOverdispersion, curMethLevel, sampler, curGeneBkgExp);

            samplingRes = this.methLevelSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
            // if samplingRes is true, means the new sampling value was accepted, also need to renew reads count expectation
            if (samplingRes) {
                this.methLevel[time][geneIdx] = curMethLevel;
                this.renewReadsExpectation(ipExpectation, geneIdx);
            } else { // otherwise, the methylation level of new iteration is same as the previous iteration
                this.methLevel[time][geneIdx] = prevMethLevel;
            }
        }
        readsData = null;
        newIPReadsExpectation = null;
        ipCount = null;
        inputCount = null;
        ipExpectation = null;
        inputExpectation = null;
    }

    /**
     * new sampling iteration of INPUT overdispersion for genes
     */
    private void inputOverdispersionSampling(int time) {

        double curMethLevel, curINPUTOverdispersion, curPosteriorProba, prevINPUTOverdispersion, prevPosteriorProba, curIPOverdispersion;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation;
        boolean samplingRes;
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            curMethLevel = this.methLevel[time][geneIdx];
            readsData = this.getReadsDataForGene(this.ipReads, this.inputReads, this.ipReadsExpectation, this.inputReadsExpectation, geneIdx);
            ipCount = readsData[0];
            inputCount = readsData[1];
            ipExpectation = readsData[2];
            inputExpectation = readsData[3];
            prevINPUTOverdispersion = this.inputOverdispersions[time-1][geneIdx];
            curINPUTOverdispersion = this.inputOverdispersionSampler.randomSample(prevINPUTOverdispersion);
            curIPOverdispersion = this.ipOverdispersions[time-1][geneIdx];

            BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[geneIdx];
            double curGeneBkgExp = this.backgroundExpression[time-1][geneIdx];
            prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                    curIPOverdispersion, prevINPUTOverdispersion, curMethLevel, sampler, curGeneBkgExp);
            curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                    curIPOverdispersion, curINPUTOverdispersion, curMethLevel, sampler, curGeneBkgExp);

            samplingRes = this.inputOverdispersionSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
            // if samplingRes is true, means the new sampling value was accepted
            if (samplingRes) {
                this.inputOverdispersions[time][geneIdx] = curINPUTOverdispersion;
            } else { // otherwise, rejected
                this.inputOverdispersions[time][geneIdx] = prevINPUTOverdispersion;
            }
        }
        readsData = null;
        ipCount = null;
        inputCount = null;
        ipExpectation = null;
        inputExpectation = null;
    }

    /**
     * new sampling iteration of IP overdispersion for genes
     */
    private void ipOverdispersionSampling(int time) {

        double curMethLevel, curIPOverdispersion, curPosteriorProba, prevIPOverdispersion, prevPosteriorProba, curINPUTOverdispersion;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation;
        boolean samplingRes;
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            curMethLevel = this.methLevel[time][geneIdx];
            readsData = this.getReadsDataForGene(this.ipReads, this.inputReads, this.ipReadsExpectation, this.inputReadsExpectation, geneIdx);
            ipCount = readsData[0];
            inputCount = readsData[1];
            ipExpectation = readsData[2];
            inputExpectation = readsData[3];
            prevIPOverdispersion = this.ipOverdispersions[time-1][geneIdx];
            curIPOverdispersion = this.ipOverdispersionSampler.randomSample(prevIPOverdispersion);
            curINPUTOverdispersion = this.inputOverdispersions[time][geneIdx];

            BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[geneIdx];
            double curGeneBkgExp = this.backgroundExpression[time-1][geneIdx];
            prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                    prevIPOverdispersion, curINPUTOverdispersion, curMethLevel, sampler, curGeneBkgExp);
            curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                    curIPOverdispersion, curINPUTOverdispersion, curMethLevel, sampler, curGeneBkgExp);

            samplingRes = this.ipOverdispersionSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
            // if samplingRes is true, means the new sampling value was accepted
            if (samplingRes) {
                this.ipOverdispersions[time][geneIdx] = curIPOverdispersion;
            } else { // otherwise, rejected
                this.ipOverdispersions[time][geneIdx] = prevIPOverdispersion;
            }
        }
        readsData = null;
        ipCount = null;
        inputCount = null;
        ipExpectation = null;
        inputExpectation = null;
    }

    /**
     * new sampling iteration of background expression for genes
     */
    private void geneBackgroundExpressionSampling(int time) {
        double curBackgroundExp, curPosteriorProba, prevBackgroundExp, prevPosteriorProba;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation, prevIPExpectation, prevINPUTExpectation;
        double ipOverdispersion, inputOverdispersion;
        boolean samplingRes;
        BackgroundExpressionSampler sampler;

        // genes' new background expression
        double[] newBkgExp = new double[this.geneNumber];
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            sampler = this.geneBackgroundExpressionSamplers[geneIdx];
            newBkgExp[geneIdx] = sampler.randomSample(this.backgroundExpression[time-1][geneIdx]);
        }

        // use the new background expression value to calculate IP and INPUT reads count expectations
        // shape individualNumber × geneNumber
        double[] curMethLevel = this.methLevel[time];
        double[][] newIPReadsExpectation, newINPUTReadsExpectation;
        ReadsExpectation re = new ReadsExpectation(this.ipReads, this.inputReads, curMethLevel, newBkgExp);
        newIPReadsExpectation = re.getIPReadsExepectation();
        newINPUTReadsExpectation = re.getINPUTReadsExpectation();
        re = null;

        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            // background expression sampler for corresponding gene
            sampler = this.geneBackgroundExpressionSamplers[geneIdx];
            prevBackgroundExp = this.backgroundExpression[time-1][geneIdx];
            curBackgroundExp = newBkgExp[geneIdx];

            // get IP and INPUT gene reads count expectation of each individual, shape 1 × individualNumber
            readsData = this.getReadsDataForGene(this.ipReads, this.inputReads, this.ipReadsExpectation, this.inputReadsExpectation, geneIdx);
            prevIPExpectation = readsData[2];
            prevINPUTExpectation = readsData[3];

            readsData = this.getReadsDataForGene(this.ipReads, this.inputReads, newIPReadsExpectation, newINPUTReadsExpectation, geneIdx);
            ipCount = readsData[0];
            inputCount = readsData[1];
            ipExpectation = readsData[2];
            inputExpectation = readsData[3];
            ipOverdispersion = this.ipOverdispersions[time][geneIdx];
            inputOverdispersion = this.inputOverdispersions[time][geneIdx];

            double meth = this.methLevel[time][geneIdx];
            prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, prevIPExpectation, prevINPUTExpectation,
                    ipOverdispersion, inputOverdispersion, meth, sampler, prevBackgroundExp);
            curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                    ipOverdispersion, inputOverdispersion, meth, sampler, curBackgroundExp);

            samplingRes = sampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
            // if samplingRes is true, means the new sampling value was accepted, also need to renew reads count expectation
            if (samplingRes) {
                this.backgroundExpression[time][geneIdx] = curBackgroundExp;
                this.renewReadsExpectation(ipExpectation, geneIdx);
            } else { // otherwise, the methylation level of new iteration is same as the previous iteration
                this.backgroundExpression[time][geneIdx] = prevBackgroundExp;
            }
        }
        readsData = null;
        newIPReadsExpectation = null;
        newINPUTReadsExpectation = null;
        ipCount = null;
        inputCount = null;
        ipExpectation = null;
        inputExpectation = null;
    }

    /**
     * get gene reads count of each individual
     */
    private double[][] getReadsDataForGene(int[][] ipReadsCount, int[][] inputReadsCount,
                                           double[][] ipCountExpectation, double[][] inputCountExpectation, int geneIdx) {
        double[] ipCount = new double[this.individualNumber];
        double[] inputCount = new double[this.individualNumber];
        double[] ipExpectation = new double[this.individualNumber];
        double[] inputExpectation = new double[this.individualNumber];
        for (int j=0; j<this.individualNumber; j++) {
            ipCount[j] = ipReadsCount[j][geneIdx];
            inputCount[j] = inputReadsCount[j][geneIdx];
            ipExpectation[j] = ipCountExpectation[j][geneIdx];
            inputExpectation[j] = inputCountExpectation[j][geneIdx];
        }

        return new double[][] {ipCount, inputCount, ipExpectation, inputExpectation};
    }

    /**
     * renew IP reads expectation if a new methylation level value is accepted
     */
    private void renewReadsExpectation(double[] ipReadsExpectation, int geneIdx) {
        for (int i=0; i<this.individualNumber; i++) {
            this.ipReadsExpectation[i][geneIdx] = ipReadsExpectation[i];
        }
    }

    /**
     * posterior distribution computation
     * @param ipReadsCount reads count of each individual cover on a particular gene in IP data, shape 1 × individualNumber
     * @param inputReadsCount reads count of each individual cover on a particular gene in INPUT data, shape 1 × individualNumber
     * @param ipExpectation reads count expectation of each individual cover on a particular gene in IP data, shape 1 × individualNumber
     * @param inputExpectation reads count expectation of each individual cover on a particular gene in INPUT data, shape 1 × individualNumber
     * @param ipOverdispersion reads count overdispersion of IP data
     * @param inputOverdispersion reads count overdispersion of INPUT data
     * @param methLevel methylation level
     * @param sampler background expression sampler of a single gene
     * @param backgroundExp background expression value
     * @return posterior distribution
     */
    private double logPosteriorProbability(double[] ipReadsCount, double[] inputReadsCount,
                                           double[] ipExpectation, double[] inputExpectation,
                                           double ipOverdispersion, double inputOverdispersion, double methLevel,
                                           BackgroundExpressionSampler sampler, double backgroundExp) {
        double methLevelProba = this.methLevelSampler.getLogDensity(methLevel);
        double ipOverdispersionProba = this.ipOverdispersionSampler.getLogDensity(ipOverdispersion);
        double inputOverdispersionProba = this.inputOverdispersionSampler.getLogDensity(inputOverdispersion);
        double backgroundExpressionProba = sampler.getLogDensity(backgroundExp);

        double ipIndividualProba = ProbabilityCalculator.logNegativeProbability(ipReadsCount, ipExpectation, ipOverdispersion);
        double inputIndividualProba = ProbabilityCalculator.logNegativeProbability(inputReadsCount, inputExpectation, inputOverdispersion);
        return methLevelProba + ipOverdispersionProba + inputOverdispersionProba + ipIndividualProba + inputIndividualProba + backgroundExpressionProba;
    }

    /**
     * use the median of methylation level sampling result as methylation quantification result, shape 1 × geneNumber
     */
    private void quantify() {
        int keepNum = this.samplingTime - this.burnIn;
        double[][] methLevelRemainValues = new double[keepNum][this.geneNumber];
        double[][] ipOverdispersionRemainValues = new double[keepNum][this.geneNumber];
        double[][] inputOverdispersionRemainValues = new double[keepNum][this.geneNumber];
        System.arraycopy(this.methLevel, this.burnIn-1, methLevelRemainValues, 0, keepNum);
        System.arraycopy(this.ipOverdispersions, this.burnIn-1, ipOverdispersionRemainValues, 0, keepNum);
        System.arraycopy(this.inputOverdispersions, this.burnIn-1, inputOverdispersionRemainValues, 0, keepNum);

        double[][] geneBackgroundExpressionRemainValue = null;

        this.geneBackgroundExpression = new double[this.geneNumber];
        geneBackgroundExpressionRemainValue = new double[keepNum][this.geneNumber];
        System.arraycopy(this.backgroundExpression, this.burnIn-1, geneBackgroundExpressionRemainValue, 0, keepNum);

        this.geneMethylationLevel = new double[this.geneNumber];
        this.geneIPOverdispersion = new double[this.geneNumber];
        this.geneINPUTOverdispersion = new double[this.geneNumber];

        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            this.geneMethylationLevel[geneIdx] = this.distributionMedian(methLevelRemainValues, geneIdx);
            this.geneIPOverdispersion[geneIdx] = this.distributionMedian(ipOverdispersionRemainValues, geneIdx);
            this.geneINPUTOverdispersion[geneIdx] = this.distributionMedian(inputOverdispersionRemainValues, geneIdx);
            this.geneBackgroundExpression[geneIdx] = this.distributionMedian(geneBackgroundExpressionRemainValue, geneIdx);
        }
    }

    private double distributionMedian(double[][] samplingResult, int geneIdx) {
        int keepNum = samplingResult.length;
        double[] methylationLevels = new double[keepNum];
        for (int i=0; i<keepNum; i++) {
            methylationLevels[i] = samplingResult[i][geneIdx];
        }

        Arrays.sort(methylationLevels);

        int medianIdx = keepNum / 2;
        if (keepNum % 2 == 0) {
            return (methylationLevels[medianIdx] + methylationLevels[medianIdx+1]) / 2;
        } else {
            return methylationLevels[medianIdx];
        }
    }
}
