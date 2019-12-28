package Quantification;

import DifferentialMethylation.CreateTask;
import ProbabilityCalculation.ProbabilityCalculator;
import SeqDataModel.BackgroundExpression;
import SeqDataModel.ReadsExpectation;

import java.util.Arrays;
import java.util.concurrent.*;
import java.util.stream.Collectors;

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
    private double[] geneBackgroundExpression, backgroundExpressionMean;
    private double[][] ipReadsExpectation, inputReadsExpectation;
    private ConcurrentHashMap<Integer, String> testRecord;

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

    /**
     * API function
     */
    public void runExecution() {
        this.initialize();
        this.generatePosteriorDistribution();

    }

    public double[] getGeneBackgroundExpression() {
        return this.geneBackgroundExpression;
    }

    /**
     * initialize sampling parameters with given parameters
     */
    private void initialize() {
        this.methLevelSampler = new MethylationLevelSampler(this.w, this.k);
        this.inputOverdispersionSampler = new OverdispersionSampler(this.c, this.d);//
        this.ipOverdispersionSampler = new OverdispersionSampler(this.a, this.b);//
        BackgroundExpression be = new BackgroundExpression(this.ipReads, this.inputReads);
        // background expression expectation and standard deviation, shape 1 × geneNumber
        this.backgroundExpressionMean = be.geneBackgroundExp();
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
            this.geneBackgroundExpressionSamplers[i] = sampler;
        }
        // initialize genes' IP and INPUT reads expectation of each individual with methylation level,
        // shape individualNumber × geneNumber. INPUT reads expectations will not be changed any more,
        // but IP reads change with different methylation level
        double[] initMethLevel = new double[this.geneNumber];
        for (int i=0; i<this.geneNumber; i++) {
            initMethLevel[i] = 0.1;
        }
        ReadsExpectation re = new ReadsExpectation(this.ipReads, this.inputReads, initMethLevel);
        this.ipReadsExpectation = re.getIPReadsExepectation();
        this.inputReadsExpectation = re.getINPUTReadsExpectation();
        re = null;
        initMethLevel = null;
    }

    /**
     * MH sampling process
     */
    private void generatePosteriorDistribution() {
        this.testRecord = new ConcurrentHashMap<>();
        CountDownLatch countDownLatch = new CountDownLatch(this.geneNumber);
        ExecutorService threadPool = Executors.newFixedThreadPool(5);
        QuantifyTask task = (sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation, geneIdx) -> {
            return () -> {
                try {
                    double[] methylationLevel = new double[this.samplingTime];
                    methylationLevel[0] = 0.1;
                    double[] ipOverdispersion = new double[this.samplingTime];
                    ipOverdispersion[0] = 1;
                    double[] inputOverdispersion = new double[this.samplingTime];
                    inputOverdispersion[0] = 1;
                    double[] backgroundExpression = new double[this.samplingTime];
                    backgroundExpression[0] = this.backgroundExpressionMean[geneIdx];

                    int time = 0;
                    while (time < this.samplingTime-1) {
                        time += 1;
                        // sampling methylation level
                        this.methylationLevelSampling(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation,
                                methylationLevel, ipOverdispersion, inputOverdispersion, backgroundExpression, time, geneIdx);
                        // sampling INPUT data overdispersion
                        this.inputOverdispersionSampling(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation,
                                methylationLevel, ipOverdispersion, inputOverdispersion, backgroundExpression, time, geneIdx);
                        // sampling IP data overdispersion
                        this.ipOverdispersionSampling(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation,
                                methylationLevel, ipOverdispersion, inputOverdispersion, backgroundExpression, time, geneIdx);
                        // sampling background expression
                        this.geneBackgroundExpressionSampling(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation,
                                methylationLevel, ipOverdispersion, inputOverdispersion, backgroundExpression, time, geneIdx);
                    }

                    double[] quantificationResult = this.quantify(methylationLevel, ipOverdispersion, inputOverdispersion, backgroundExpression);
                    String[] res = new String[quantificationResult.length];
                    for (int i=0; i<quantificationResult.length; i++) {
                        res[i] = String.valueOf(quantificationResult[i]);
                    }
                    String record = String.join("\t", res);
                    System.out.println(record);
                    this.testRecord.put(geneIdx, record);
                } catch (Exception e) {
                    e.printStackTrace();
                } finally {
                    countDownLatch.countDown();
                }
            };
        };

        double[] initMethLevel = new double[] {0.1};
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            int[][] sampleIPReads = new int[this.individualNumber][1];
            int[][] sampleINPUTReads = new int[this.individualNumber][1];
            for (int individualIdx=0; individualIdx<this.individualNumber; individualIdx++) {
                sampleIPReads[individualIdx][0] = this.ipReads[individualIdx][geneIdx];
                sampleINPUTReads[individualIdx][0] = this.inputReads[individualIdx][geneIdx];
            }
            double[][] sampleIPExpectation, sampleINPUTExpectation;

            ReadsExpectation re = new ReadsExpectation(sampleIPReads, sampleINPUTReads, initMethLevel);
            sampleIPExpectation = re.getIPReadsExepectation();
            sampleINPUTExpectation = re.getINPUTReadsExpectation();
            re = null;

            Runnable runnable = task.getTask(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation, geneIdx);
            threadPool.submit(runnable);
        }

        try {
            countDownLatch.await();
            if (!threadPool.awaitTermination(1000, TimeUnit.MILLISECONDS))
                threadPool.shutdown();
        } catch (InterruptedException ie) {
            ie.printStackTrace();
        } finally {
            threadPool.shutdownNow();
        }
    }

    /**
     * new sampling iteration of Methylation level for each gene
     */
    private void methylationLevelSampling(int[][] sampleIPReads, int[][] sampleINPUTReads,
                                          double[][] sampleIPExpectation, double[][] sampleINPUTExpectation,
                                          double[] methylationLevels, double[] ipOverdispersions,
                                          double[] inputOverdispersions, double[] backgroundExpressions, int time, int geneIdx) {

        double curMethLevel, curPosteriorProba, prevMethLevel, prevPosteriorProba;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation, prevIPExpectation, prevINPUTExpectation;
        double ipOverdispersion, inputOverdispersion;
        boolean samplingRes;

        // gene's new methylation level
        curMethLevel =  this.methLevelSampler.randomSample(methylationLevels[time-1]);

        // use the new methylation level to calculate IP and INPUT reads count expectations
        // shape individualNumber × geneNumber
        double[][] newIPReadsExpectation;
        ReadsExpectation re = new ReadsExpectation(sampleIPReads, sampleINPUTReads, new double[] {curMethLevel});
        newIPReadsExpectation = re.getIPReadsExepectation();
        re = null;

        prevMethLevel = methylationLevels[time-1];

        // get IP and INPUT gene reads count expectation of each individual, shape 1 × individualNumber
        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation);
        prevIPExpectation = readsData[2];
        prevINPUTExpectation = readsData[3];

        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, newIPReadsExpectation, sampleINPUTExpectation);
        ipCount = readsData[0];
        inputCount = readsData[1];
        ipExpectation = readsData[2];
        inputExpectation = readsData[3];
        ipOverdispersion = ipOverdispersions[time-1];
        inputOverdispersion = inputOverdispersions[time-1];

        BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[geneIdx];
        double curGeneBkgExp = backgroundExpressions[time-1];
        prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, prevIPExpectation, prevINPUTExpectation,
                ipOverdispersion, inputOverdispersion, prevMethLevel, sampler, curGeneBkgExp);
        curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                ipOverdispersion, inputOverdispersion, curMethLevel, sampler, curGeneBkgExp);

        samplingRes = this.methLevelSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
        // if samplingRes is true, means the new sampling value was accepted, also need to renew reads count expectation
        if (samplingRes) {
            methylationLevels[time] = curMethLevel;
            this.renewReadsExpectation(sampleIPExpectation, ipExpectation);
        } else { // otherwise, the methylation level of new iteration is same as the previous iteration
            methylationLevels[time] = prevMethLevel;
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
    private void inputOverdispersionSampling(int[][] sampleIPReads, int[][] sampleINPUTReads,
                                             double[][] sampleIPExpectation, double[][] sampleINPUTExpectation,
                                             double[] methylationLevels, double[] ipOverdispersions,
                                             double[] inputOverdispersions, double[] backgroundExpressions, int time, int geneIdx) {

        double curMethLevel, curINPUTOverdispersion, curPosteriorProba, prevINPUTOverdispersion, prevPosteriorProba, curIPOverdispersion;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation;
        boolean samplingRes;

        curMethLevel = methylationLevels[time];
        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation);
        ipCount = readsData[0];
        inputCount = readsData[1];
        ipExpectation = readsData[2];
        inputExpectation = readsData[3];
        prevINPUTOverdispersion = inputOverdispersions[time-1];
        curINPUTOverdispersion = this.inputOverdispersionSampler.randomSample(prevINPUTOverdispersion);
        curIPOverdispersion = ipOverdispersions[time-1];

        BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[geneIdx];
        double curGeneBkgExp = backgroundExpressions[time-1];
        prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                curIPOverdispersion, prevINPUTOverdispersion, curMethLevel, sampler, curGeneBkgExp);
        curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                curIPOverdispersion, curINPUTOverdispersion, curMethLevel, sampler, curGeneBkgExp);

        samplingRes = this.inputOverdispersionSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
        // if samplingRes is true, means the new sampling value was accepted
        if (samplingRes) {
            inputOverdispersions[time] = curINPUTOverdispersion;
        } else { // otherwise, rejected
            inputOverdispersions[time] = prevINPUTOverdispersion;
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
    private void ipOverdispersionSampling(int[][] sampleIPReads, int[][] sampleINPUTReads,
                                          double[][] sampleIPExpectation, double[][] sampleINPUTExpectation,
                                          double[] methylationLevels, double[] ipOverdispersions,
                                          double[] inputOverdispersions, double[] backgroundExpressions, int time, int geneIdx) {

        double curMethLevel, curIPOverdispersion, curPosteriorProba, prevIPOverdispersion, prevPosteriorProba, curINPUTOverdispersion;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation;
        boolean samplingRes;

        curMethLevel = methylationLevels[time];
        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation);
        ipCount = readsData[0];
        inputCount = readsData[1];
        ipExpectation = readsData[2];
        inputExpectation = readsData[3];
        prevIPOverdispersion = ipOverdispersions[time-1];
        curIPOverdispersion = this.ipOverdispersionSampler.randomSample(prevIPOverdispersion);
        curINPUTOverdispersion = inputOverdispersions[time];

        BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[geneIdx];
        double curGeneBkgExp = backgroundExpressions[time-1];
        prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                prevIPOverdispersion, curINPUTOverdispersion, curMethLevel, sampler, curGeneBkgExp);
        curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                curIPOverdispersion, curINPUTOverdispersion, curMethLevel, sampler, curGeneBkgExp);

        samplingRes = this.ipOverdispersionSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
        // if samplingRes is true, means the new sampling value was accepted
        if (samplingRes) {
            ipOverdispersions[time] = curIPOverdispersion;
        } else { // otherwise, rejected
            ipOverdispersions[time] = prevIPOverdispersion;
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
    private void geneBackgroundExpressionSampling(int[][] sampleIPReads, int[][] sampleINPUTReads,
                                                  double[][] sampleIPExpectation, double[][] sampleINPUTExpectation,
                                                  double[] methylationLevels, double[] ipOverdispersions,
                                                  double[] inputOverdispersions, double[] backgroundExpressions, int time, int geneIdx) {
        double curBackgroundExp, curPosteriorProba, prevBackgroundExp, prevPosteriorProba;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation, prevIPExpectation, prevINPUTExpectation;
        double ipOverdispersion, inputOverdispersion;
        boolean samplingRes;
        BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[geneIdx];

        // gene's new background expression
        curBackgroundExp = sampler.randomSample(backgroundExpressions[time-1]);

        // use the new background expression value to calculate IP and INPUT reads count expectations
        // shape individualNumber × geneNumber
        double[] curMethLevel = new double[] {methylationLevels[time]};
        double[][] newIPReadsExpectation, newINPUTReadsExpectation;
        ReadsExpectation re = new ReadsExpectation(sampleIPReads, sampleINPUTReads, curMethLevel, new double[] {curBackgroundExp});
        newIPReadsExpectation = re.getIPReadsExepectation();
        newINPUTReadsExpectation = re.getINPUTReadsExpectation();
        re = null;

        // background expression sampler for corresponding gene
        sampler = this.geneBackgroundExpressionSamplers[geneIdx];
        prevBackgroundExp = backgroundExpressions[time-1];

        // get IP and INPUT gene reads count expectation of each individual, shape 1 × individualNumber
        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation);
        prevIPExpectation = readsData[2];
        prevINPUTExpectation = readsData[3];

        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, newIPReadsExpectation, newINPUTReadsExpectation);
        ipCount = readsData[0];
        inputCount = readsData[1];
        ipExpectation = readsData[2];
        inputExpectation = readsData[3];
        ipOverdispersion = ipOverdispersions[time];
        inputOverdispersion = inputOverdispersions[time];

        double meth = methylationLevels[time];
        prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, prevIPExpectation, prevINPUTExpectation,
                ipOverdispersion, inputOverdispersion, meth, sampler, prevBackgroundExp);
        curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                ipOverdispersion, inputOverdispersion, meth, sampler, curBackgroundExp);

        samplingRes = sampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
        // if samplingRes is true, means the new sampling value was accepted, also need to renew reads count expectation
        if (samplingRes) {
            backgroundExpressions[time] = curBackgroundExp;
            this.renewReadsExpectation(sampleIPExpectation, ipExpectation);
        } else { // otherwise, the methylation level of new iteration is same as the previous iteration
            backgroundExpressions[time] = prevBackgroundExp;
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
                                           double[][] ipCountExpectation, double[][] inputCountExpectation) {
        double[] ipCount = new double[this.individualNumber];
        double[] inputCount = new double[this.individualNumber];
        double[] ipExpectation = new double[this.individualNumber];
        double[] inputExpectation = new double[this.individualNumber];
        for (int j=0; j<this.individualNumber; j++) {
            ipCount[j] = ipReadsCount[j][0];
            inputCount[j] = inputReadsCount[j][0];
            ipExpectation[j] = ipCountExpectation[j][0];
            inputExpectation[j] = inputCountExpectation[j][0];
        }

        return new double[][] {ipCount, inputCount, ipExpectation, inputExpectation};
    }

    /**
     * renew IP reads expectation if a new methylation level value is accepted
     */
    private void renewReadsExpectation(double[][] sampleIPExpectation, double[] ipReadsExpectation) {
        for (int i=0; i<this.individualNumber; i++) {
            sampleIPExpectation[i][0] = ipReadsExpectation[i];
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
    private double[] quantify(double[] methylationLevels, double[] ipOverdispersions,
                          double[] inputOverdispersions, double[] backgroundExpressions) {
        int keepNum = this.samplingTime - this.burnIn;
        double[] methLevelRemainValues = new double[keepNum];
        double[] ipOverdispersionRemainValues = new double[keepNum];
        double[] inputOverdispersionRemainValues = new double[keepNum];
        System.arraycopy(methylationLevels, this.burnIn-1, methLevelRemainValues, 0, keepNum);
        System.arraycopy(ipOverdispersions, this.burnIn-1, ipOverdispersionRemainValues, 0, keepNum);
        System.arraycopy(inputOverdispersions, this.burnIn-1, inputOverdispersionRemainValues, 0, keepNum);

        double[] geneBackgroundExpressionRemainValue = null;

        this.geneBackgroundExpression = new double[this.geneNumber];
        geneBackgroundExpressionRemainValue = new double[keepNum];
        System.arraycopy(backgroundExpressions, this.burnIn-1, geneBackgroundExpressionRemainValue, 0, keepNum);

        double methylationValue, ipOverdispersionValue, inputOverdispersionValue, backgroundExpressionValue;
        methylationValue = this.distributionMedian(methLevelRemainValues);
        ipOverdispersionValue = this.distributionMedian(ipOverdispersionRemainValues);
        inputOverdispersionValue = this.distributionMedian(inputOverdispersionRemainValues);
        backgroundExpressionValue = this.distributionMedian(geneBackgroundExpressionRemainValue);

        return new double[] {methylationValue, ipOverdispersionValue, inputOverdispersionValue, backgroundExpressionValue};
    }

    /**
     * use the sampling median as methylation quantification result
     * @param samplingResult sampling result, 1 × (samplingTime-burnIn)
     * @return sampling median
     */
    private double distributionMedian(double[] samplingResult) {
        int keepNum = samplingResult.length;
        Arrays.sort(samplingResult);

        int medianIdx = keepNum / 2;
        if (keepNum % 2 == 0) {
            return (samplingResult[medianIdx] + samplingResult[medianIdx+1]) / 2;
        } else {
            return samplingResult[medianIdx];
        }
    }

    public void output() {

    }
}
