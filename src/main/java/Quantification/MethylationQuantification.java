package Quantification;

import ProbabilityCalculation.ProbabilityCalculator;
import SeqDataModel.BackgroundExpression;
import SeqDataModel.ReadsExpectation;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.*;
import java.util.stream.Collectors;

/**
 * Quantify the methylation level with component-wise MH sampling
 * reference https://theclevermachine.wordpress.com/2012/11/04/mcmc-multivariate-distributions-block-wise-component-wise-updates/
 */
public class MethylationQuantification {
    private MethylationLevelSampler methLevelSampler = null;
    private NonSpecificEnrichmentSampler nonSpecificEnrichmentSampler = null;
    private OverdispersionSampler ipOverdispersionSampler = null, inputOverdispersionSampler = null;
//    private BackgroundExpressionSampler[] geneBackgroundExpressionSamplers = null;
    private int individualNumber, geneNumber, samplingTime, burnIn, threadNumber;
    private int[][] ipReads, inputReads, ipNonPeak, inputNonPeak;
    private double a, b, c, d, w, k, r1, r2;
    private double[] backgroundExpressionMean, nonPeakBackgroundExpressionMean;
    private String outputFile;
    private ConcurrentHashMap<Integer, String> testRecord;

    public MethylationQuantification(int[][] ipReads, int[][] inputReads, int[][] ipNonPeak, int[][] inputNonPeak,
                                     double ipShape, double ipScale, double inputShape, double inputScale,
                                     double methLevelShape1, double methLevelShape2,
                                     double nonSpecificEnrich1, double nonSpecificEnrich2,
                                     int samplingTime, int burnIn, int threadNumber, String outputFile) {
        assert samplingTime > burnIn;
        // assert arrays contains same sample and gene number
        assert inputReads.length == ipReads.length && inputReads[0].length == ipReads[0].length;
        assert inputNonPeak.length == ipNonPeak.length && inputNonPeak[0].length == ipNonPeak[0].length;
        this.individualNumber = inputReads.length;
        this.geneNumber = inputReads[0].length;
        this.inputReads = inputReads;
        this.ipReads = ipReads;
        this.inputNonPeak = inputNonPeak;
        this.ipNonPeak = ipNonPeak;
        this.a = ipShape;
        this.b = ipScale;
        this.c = inputShape;
        this.d = inputScale;
        this.w = methLevelShape1;
        this.k = methLevelShape2;
        this.r1 = nonSpecificEnrich1;
        this.r2 = nonSpecificEnrich2;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.threadNumber = threadNumber;
        this.outputFile = outputFile;
    }

    /**
     * API function
     */
    public void runExecution() {
        this.initialize();
        this.generatePosteriorDistribution();
        this.output();
    }

    public ArrayList<String> getTestRecord() {
        ArrayList<String> result = new ArrayList<>();
        this.testRecord.entrySet().stream().sorted((o1, o2) -> o1.getKey() - o2.getKey()).forEach(x -> result.add(x.getValue()));
        return result;
    }

    /**
     * initialize sampling parameters with given parameters
     */
    private void initialize() {
        this.methLevelSampler = new MethylationLevelSampler(this.w, this.k);
        this.inputOverdispersionSampler = new OverdispersionSampler(this.c, this.d);
        this.ipOverdispersionSampler = new OverdispersionSampler(this.a, this.b);
        this.nonSpecificEnrichmentSampler = new NonSpecificEnrichmentSampler(this.r1, this.r2);
        BackgroundExpression be = new BackgroundExpression(this.ipReads, this.inputReads);
        // background expression expectation and standard deviation, shape 1 × geneNumber
        this.backgroundExpressionMean = be.geneBackgroundExp();
        double[] backgroundExpressionStd = be.geneExpressionStd();

        be = new BackgroundExpression(this.ipNonPeak, this.inputNonPeak);
        this.nonPeakBackgroundExpressionMean = be.geneBackgroundExp();

        // for each gene, initialize a background expression sampler. cancel background expression sampling
//        this.geneBackgroundExpressionSamplers = new BackgroundExpressionSampler[this.geneNumber];
//        BackgroundExpressionSampler sampler;
//        for (int i=0; i<this.geneNumber; i++) {
//            // shape and scale parameters of log-normal distribution
//            double scale = Math.log(backgroundExpressionMean[i]);
//            double shape = backgroundExpressionStd[i];
//            sampler = new BackgroundExpressionSampler(scale, shape);
//            // random init background expression from corresponding distribution
//            this.geneBackgroundExpressionSamplers[i] = sampler;
//        }
    }

    /**
     * MH sampling process
     */
    private void generatePosteriorDistribution() {
        this.testRecord = new ConcurrentHashMap<>();
        CountDownLatch countDownLatch = new CountDownLatch(this.geneNumber);
        ExecutorService threadPool = Executors.newFixedThreadPool(this.threadNumber);
        QuantifyTask task = (sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation,
                             nonPeakIPReads, nonPeakINPUTReads, nonPeakIPExpectation, nonPeakINPUTExpectation, geneIdx) -> {
            return () -> {
                try {
                    SamplingRecords sr = new SamplingRecords();
                    double[] nonSpecificEnrichmentRatio = new double[this.samplingTime];
                    nonSpecificEnrichmentRatio[0] = 0.1;
                    double[] methylationLevel = new double[this.samplingTime];
                    methylationLevel[0] = 0.1;
                    double[] ipOverdispersion = new double[this.samplingTime];
                    ipOverdispersion[0] = 1;
                    double[] inputOverdispersion = new double[this.samplingTime];
                    inputOverdispersion[0] = 1;
                    double[] backgroundExpression = new double[this.samplingTime];
                    backgroundExpression[0] = this.backgroundExpressionMean[geneIdx];

                    sr.setReads(sampleIPReads, sampleINPUTReads, nonPeakIPReads, nonPeakINPUTReads);
                    sr.setExpectations(sampleIPExpectation, sampleINPUTExpectation, nonPeakIPExpectation, nonPeakINPUTExpectation);
                    sr.setNonSpecificEnrichment(nonSpecificEnrichmentRatio);
                    sr.setMethylationLevels(methylationLevel);
                    sr.setIpOverdispersions(ipOverdispersion);
                    sr.setInputOverdispersions(inputOverdispersion);
                    sr.setBackgroundExpressions(backgroundExpression);
                    sr.setGeneBackgroundExp(this.backgroundExpressionMean[geneIdx]);
                    sr.setNonPeakExpression(this.nonPeakBackgroundExpressionMean[geneIdx]);

                    int time = 0;
                    while (time < this.samplingTime-1) {
                        time += 1;
                        // sampling non-specific enrichment ratio
                        this.nonSpecificEnrichmentRatioSampling(sr, time);
                        // sampling methylation level
                        this.methylationLevelSampling(sr, time);
                        // sampling INPUT data overdispersion
                        this.inputOverdispersionSampling(sr, time);
                        // sampling IP data overdispersion
                        this.ipOverdispersionSampling(sr, time);
                        // sampling background expression
//                        this.geneBackgroundExpressionSampling(sr, time);
                    }

                    double[] quantificationResult = this.quantify(sr);
                    String[] res = new String[quantificationResult.length];
                    for (int i=0; i<quantificationResult.length; i++) {
                        res[i] = String.valueOf(quantificationResult[i]);
                    }
                    String record = String.join("\t", res);
//                    System.out.println(geneIdx + "\t" + record);
                    this.testRecord.put(geneIdx, record);
                } catch (Exception e) {
                    e.printStackTrace();
                } finally {
                    System.out.println(geneIdx+1 + " / " + this.geneNumber);
                    countDownLatch.countDown();
                }
            };
        };

        double[] initMethLevel = new double[] {0.1};
        double[] initNonSpecificEnrich = new double[] {0.1};
        for (int geneIdx=0; geneIdx<this.geneNumber; geneIdx++) {
            int[][] sampleIPReads = new int[this.individualNumber][1];
            int[][] sampleINPUTReads = new int[this.individualNumber][1];
            int[][] nonPeakIPReads = new int[this.individualNumber][1];
            int[][] nonPeakINPUTReads = new int[this.individualNumber][1];
            for (int individualIdx=0; individualIdx<this.individualNumber; individualIdx++) {
                sampleIPReads[individualIdx][0] = this.ipReads[individualIdx][geneIdx];
                sampleINPUTReads[individualIdx][0] = this.inputReads[individualIdx][geneIdx];
                nonPeakIPReads[individualIdx][0] = this.ipNonPeak[individualIdx][geneIdx];
                nonPeakINPUTReads[individualIdx][0] = this.inputNonPeak[individualIdx][geneIdx];
            }
            double[][] sampleIPExpectation, sampleINPUTExpectation, nonPeakIPExpectation, nonPeakINPUTExpectation;

            ReadsExpectation re = new ReadsExpectation(sampleIPReads, sampleINPUTReads, nonPeakIPReads, nonPeakINPUTReads, initMethLevel, initNonSpecificEnrich);
            sampleIPExpectation = re.getIPReadsExepectation();
            sampleINPUTExpectation = re.getINPUTReadsExpectation();
            nonPeakIPExpectation = re.getIPNonPeakExepectation();
            nonPeakINPUTExpectation = re.getINPUTNonPeakExpectation();
            re = null;

            Runnable runnable = task.getTask(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation,
                                             nonPeakIPReads, nonPeakINPUTReads, nonPeakIPExpectation, nonPeakINPUTExpectation, geneIdx);
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
     * new sampling iteration of nonspecific enrichment ratio for each gene
     */
    private void nonSpecificEnrichmentRatioSampling(SamplingRecords samplingRecords, int time) {
        double prevNonSpecificEnrichment, curNonSpecificEnrichment, curPosteriorProba, prevPosteriorProba;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation, ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount,
                 prevIPExpectation, prevINPUTExpectation, prevNonPeakIPExpectation, prevNonPeakINPUTExpectation;
        double ipOverdispersion, inputOverdispersion, methLevel;
        boolean samplingRes;

        int[][] sampleIPReads = samplingRecords.getSampleIPReads();
        int[][] sampleINPUTReads = samplingRecords.getSampleINPUTReads();
        int[][] nonPeakIPReads = samplingRecords.getNonPeakIPReads();
        int[][] nonPeakINPUTReads = samplingRecords.getNonPeakINPUTReads();
        double[][] sampleIPExpectation = samplingRecords.getSampleIPExpectations();
        double[][] sampleINPUTExpectation = samplingRecords.getSampleINPUTExpectations();
        double[][] nonPeakIPExpectation = samplingRecords.getNonPeakIPExpectations();
        double[][] nonPeakINPUTExpectation = samplingRecords.getNonPeakINPUTExpectations();

        // gene's new non-specific enrichment ratio
        prevNonSpecificEnrichment = samplingRecords.getNonSpecificEnrichmentRatio(time-1);
        curNonSpecificEnrichment = this.nonSpecificEnrichmentSampler.randomSample(prevNonSpecificEnrichment);
        methLevel = samplingRecords.getMethylationLevelValue(time-1);

        // use the new methylation level to calculate IP and INPUT reads count expectations
        // shape individualNumber × geneNumber
        double[][] newIPReadsExpectation, newNonPeakIPExpectation;
        ReadsExpectation re = new ReadsExpectation(sampleIPReads, sampleINPUTReads, nonPeakIPReads, nonPeakINPUTReads, new double[] {methLevel}, new double[] {curNonSpecificEnrichment});
        newIPReadsExpectation = re.getIPReadsExepectation();
        newNonPeakIPExpectation = re.getIPNonPeakExepectation();
        re = null;

        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation,
                                             nonPeakIPReads, nonPeakINPUTReads, nonPeakIPExpectation, nonPeakINPUTExpectation);
        prevIPExpectation = readsData[2];
        prevINPUTExpectation = readsData[3];
        prevNonPeakIPExpectation = readsData[6];
        prevNonPeakINPUTExpectation = readsData[7];

        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, newIPReadsExpectation, sampleINPUTExpectation,
                                             nonPeakIPReads, nonPeakINPUTReads, newNonPeakIPExpectation, nonPeakINPUTExpectation);
        ipCount = readsData[0];
        inputCount = readsData[1];
        ipExpectation = readsData[2];
        inputExpectation = readsData[3];
        ipNonPeakCount = readsData[4];
        inputNonPeakCount = readsData[5];
        ipNonPeakExpCount = readsData[6];
        inputNonPeakExpCount = readsData[7];
        ipOverdispersion = samplingRecords.getIpOverdispersionValue(time-1);
        inputOverdispersion = samplingRecords.getInputOverdispersionValue(time-1);

        //        BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[geneIdx];
        BackgroundExpressionSampler sampler = null;
//        double curGeneBkgExp = backgroundExpressions[time-1];
        double curGeneBkgExp = samplingRecords.getGeneBackgroundExp();
        prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, prevIPExpectation, prevINPUTExpectation,
                                                          ipNonPeakCount, inputNonPeakCount, prevNonPeakIPExpectation, prevNonPeakINPUTExpectation,
                                                          ipOverdispersion, inputOverdispersion, methLevel, prevNonSpecificEnrichment, sampler, curGeneBkgExp);
        curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                                                         ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount,
                                                         ipOverdispersion, inputOverdispersion, methLevel, curNonSpecificEnrichment, sampler, curGeneBkgExp);

        samplingRes = this.nonSpecificEnrichmentSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
        // if samplingRes is true, means the new sampling value was accepted, also need to renew reads count expectation
        if (samplingRes) {
            samplingRecords.setNonSpecificEnrichmentRatio(curNonSpecificEnrichment, time);
//            this.renewReadsExpectation(sampleIPExpectation, ipExpectation);
            samplingRecords.setSampleIPExpectations(newIPReadsExpectation);
            samplingRecords.setNonPeakIPExpectations(newNonPeakIPExpectation);
        } else { // otherwise, the methylation level of new iteration is same as the previous iteration
            samplingRecords.setNonSpecificEnrichmentRatio(prevNonSpecificEnrichment, time);
        }

        readsData = null;
        newIPReadsExpectation = null;
        ipCount = null;
        inputCount = null;
        ipExpectation = null;
        inputExpectation = null;
    }

    /**
     * new sampling iteration of Methylation level for each gene
     */
    private void methylationLevelSampling(SamplingRecords samplingRecords, int time) {

        double curMethLevel, curPosteriorProba, prevMethLevel, prevPosteriorProba;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation,
                 ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount,prevIPExpectation, prevINPUTExpectation;
        double nonSpecificEnrichRatio, ipOverdispersion, inputOverdispersion;
        boolean samplingRes;

        int[][] sampleIPReads = samplingRecords.getSampleIPReads();
        int[][] sampleINPUTReads = samplingRecords.getSampleINPUTReads();
        int[][] nonPeakIPReads = samplingRecords.getNonPeakIPReads();
        int[][] nonPeakINPUTReads = samplingRecords.getNonPeakINPUTReads();
        double[][] sampleIPExpectation = samplingRecords.getSampleIPExpectations();
        double[][] sampleINPUTExpectation = samplingRecords.getSampleINPUTExpectations();
        double[][] nonPeakIPExpectation = samplingRecords.getNonPeakIPExpectations();
        double[][] nonPeakINPUTExpectation = samplingRecords.getNonPeakINPUTExpectations();

        // gene's new methylation level
        prevMethLevel = samplingRecords.getMethylationLevelValue(time-1);
        curMethLevel =  this.methLevelSampler.randomSample(prevMethLevel);

        nonSpecificEnrichRatio = samplingRecords.getNonSpecificEnrichmentRatio(time);
        // use the new methylation level to calculate IP and INPUT reads count expectations
        // shape individualNumber × geneNumber
        double[][] newIPReadsExpectation;
        ReadsExpectation re = new ReadsExpectation(sampleIPReads, sampleINPUTReads, nonPeakIPReads, nonPeakINPUTReads, new double[] {curMethLevel}, new double[] {nonSpecificEnrichRatio});
        newIPReadsExpectation = re.getIPReadsExepectation();
        re = null;

        // get IP and INPUT gene reads count expectation of each individual, shape 1 × individualNumber
        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation,
                                             nonPeakIPReads, nonPeakINPUTReads, nonPeakIPExpectation, nonPeakINPUTExpectation);
        prevIPExpectation = readsData[2];
        prevINPUTExpectation = readsData[3];

        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, newIPReadsExpectation, sampleINPUTExpectation,
                                             nonPeakIPReads, nonPeakINPUTReads, nonPeakIPExpectation, nonPeakINPUTExpectation);
        ipCount = readsData[0];
        inputCount = readsData[1];
        ipExpectation = readsData[2];
        inputExpectation = readsData[3];
        ipNonPeakCount = readsData[4];
        inputNonPeakCount = readsData[5];
        ipNonPeakExpCount = readsData[6];
        inputNonPeakExpCount = readsData[7];
        ipOverdispersion = samplingRecords.getIpOverdispersionValue(time-1);
        inputOverdispersion = samplingRecords.getInputOverdispersionValue(time-1);

//        BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[geneIdx];
        BackgroundExpressionSampler sampler = null;
//        double curGeneBkgExp = backgroundExpressions[time-1];
        double curGeneBkgExp = samplingRecords.getGeneBackgroundExp();
        prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, prevIPExpectation, prevINPUTExpectation,
                                                          ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount,
                                                          ipOverdispersion, inputOverdispersion, prevMethLevel, nonSpecificEnrichRatio, sampler, curGeneBkgExp);
        curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                                                         ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount,
                                                         ipOverdispersion, inputOverdispersion, curMethLevel, nonSpecificEnrichRatio, sampler, curGeneBkgExp);

        samplingRes = this.methLevelSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
        // if samplingRes is true, means the new sampling value was accepted, also need to renew reads count expectation
        if (samplingRes) {
            samplingRecords.setMethylationLevelValue(curMethLevel, time);
//            this.renewReadsExpectation(sampleIPExpectation, ipExpectation);
            samplingRecords.setSampleIPExpectations(newIPReadsExpectation);
        } else { // otherwise, the methylation level of new iteration is same as the previous iteration
            samplingRecords.setMethylationLevelValue(prevMethLevel, time);
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
    private void inputOverdispersionSampling(SamplingRecords samplingRecords, int time) {

        double curNonSpecificEnrich, curMethLevel, curINPUTOverdispersion, curPosteriorProba, prevINPUTOverdispersion, prevPosteriorProba, curIPOverdispersion;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation, ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount;
        boolean samplingRes;

        int[][] sampleIPReads = samplingRecords.getSampleIPReads();
        int[][] sampleINPUTReads = samplingRecords.getSampleINPUTReads();
        int[][] nonPeakIPReads = samplingRecords.getNonPeakIPReads();
        int[][] nonPeakINPUTReads = samplingRecords.getNonPeakINPUTReads();
        double[][] sampleIPExpectation = samplingRecords.getSampleIPExpectations();
        double[][] sampleINPUTExpectation = samplingRecords.getSampleINPUTExpectations();
        double[][] nonPeakIPExpectation = samplingRecords.getNonPeakIPExpectations();
        double[][] nonPeakINPUTExpectation = samplingRecords.getNonPeakINPUTExpectations();

        curNonSpecificEnrich = samplingRecords.getNonSpecificEnrichmentRatio(time);
        curMethLevel = samplingRecords.getMethylationLevelValue(time);
        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation,
                                             nonPeakIPReads, nonPeakINPUTReads, nonPeakIPExpectation, nonPeakINPUTExpectation);
        ipCount = readsData[0];
        inputCount = readsData[1];
        ipExpectation = readsData[2];
        inputExpectation = readsData[3];
        ipNonPeakCount = readsData[4];
        inputNonPeakCount = readsData[5];
        ipNonPeakExpCount = readsData[6];
        inputNonPeakExpCount = readsData[7];
        prevINPUTOverdispersion = samplingRecords.getInputOverdispersionValue(time-1);
        curINPUTOverdispersion = this.inputOverdispersionSampler.randomSample(prevINPUTOverdispersion);
        curIPOverdispersion = samplingRecords.getIpOverdispersionValue(time-1);

//        BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[geneIdx];
        BackgroundExpressionSampler sampler = null;
//        double curGeneBkgExp = backgroundExpressions[time-1];
        double curGeneBkgExp = samplingRecords.getGeneBackgroundExp();
        prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                                                          ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount,
                                                          curIPOverdispersion, prevINPUTOverdispersion, curMethLevel, curNonSpecificEnrich, sampler, curGeneBkgExp);
        curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                                                         ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount,
                                                         curIPOverdispersion, curINPUTOverdispersion, curMethLevel, curNonSpecificEnrich, sampler, curGeneBkgExp);

        samplingRes = this.inputOverdispersionSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
        // if samplingRes is true, means the new sampling value was accepted
        if (samplingRes) {
            samplingRecords.setInputOverdispersionValue(curINPUTOverdispersion, time);
        } else { // otherwise, rejected
            samplingRecords.setInputOverdispersionValue(prevINPUTOverdispersion, time);
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
    private void ipOverdispersionSampling(SamplingRecords samplingRecords, int time) {

        double curNonSpecificEnrich, curMethLevel, curIPOverdispersion, curPosteriorProba, prevIPOverdispersion, prevPosteriorProba, curINPUTOverdispersion;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation, ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount;
        boolean samplingRes;

        int[][] sampleIPReads = samplingRecords.getSampleIPReads();
        int[][] sampleINPUTReads = samplingRecords.getSampleINPUTReads();
        int[][] nonPeakIPReads = samplingRecords.getNonPeakIPReads();
        int[][] nonPeakINPUTReads = samplingRecords.getNonPeakINPUTReads();
        double[][] sampleIPExpectation = samplingRecords.getSampleIPExpectations();
        double[][] sampleINPUTExpectation = samplingRecords.getSampleINPUTExpectations();
        double[][] nonPeakIPExpectation = samplingRecords.getNonPeakIPExpectations();
        double[][] nonPeakINPUTExpectation = samplingRecords.getNonPeakINPUTExpectations();

        curNonSpecificEnrich = samplingRecords.getNonSpecificEnrichmentRatio(time);
        curMethLevel = samplingRecords.getMethylationLevelValue(time);
        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation,
                                             nonPeakIPReads, nonPeakINPUTReads, nonPeakIPExpectation, nonPeakINPUTExpectation);
        ipCount = readsData[0];
        inputCount = readsData[1];
        ipExpectation = readsData[2];
        inputExpectation = readsData[3];
        ipNonPeakCount = readsData[4];
        inputNonPeakCount = readsData[5];
        ipNonPeakExpCount = readsData[6];
        inputNonPeakExpCount = readsData[7];
        prevIPOverdispersion = samplingRecords.getIpOverdispersionValue(time-1);
        curIPOverdispersion = this.ipOverdispersionSampler.randomSample(prevIPOverdispersion);
        curINPUTOverdispersion = samplingRecords.getInputOverdispersionValue(time);

//        BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[geneIdx];
        BackgroundExpressionSampler sampler = null;
//        double curGeneBkgExp = backgroundExpressions[time-1];
        double curGeneBkgExp = samplingRecords.getGeneBackgroundExp();
        prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                                                          ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount,
                                                          prevIPOverdispersion, curINPUTOverdispersion, curMethLevel, curNonSpecificEnrich, sampler, curGeneBkgExp);
        curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                                                         ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount,
                                                         curIPOverdispersion, curINPUTOverdispersion, curMethLevel,curNonSpecificEnrich, sampler, curGeneBkgExp);

        samplingRes = this.ipOverdispersionSampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
        // if samplingRes is true, means the new sampling value was accepted
        if (samplingRes) {
            samplingRecords.setIpOverdispersionValue(curIPOverdispersion, time);
        } else { // otherwise, rejected
            samplingRecords.setIpOverdispersionValue(prevIPOverdispersion, time);
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
    private void geneBackgroundExpressionSampling(SamplingRecords samplingRecords, int time) {

        double curBackgroundExp, curPosteriorProba, prevBackgroundExp, prevPosteriorProba;
        double[][] readsData;
        double[] ipCount, inputCount, ipExpectation, inputExpectation, prevIPExpectation, prevINPUTExpectation,
                 ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount;
        double ipOverdispersion, inputOverdispersion;
        boolean samplingRes;

        int[][] sampleIPReads = samplingRecords.getSampleIPReads();
        int[][] sampleINPUTReads = samplingRecords.getSampleINPUTReads();
        int[][] nonPeakIPReads = samplingRecords.getNonPeakIPReads();
        int[][] nonPeakINPUTReads = samplingRecords.getNonPeakINPUTReads();
        double[][] sampleIPExpectation = samplingRecords.getSampleIPExpectations();
        double[][] sampleINPUTExpectation = samplingRecords.getSampleINPUTExpectations();
        double[][] nonPeakIPExpectation = samplingRecords.getNonPeakIPExpectations();
        double[][] nonPeakINPUTExpectation = samplingRecords.getNonPeakINPUTExpectations();

//        BackgroundExpressionSampler sampler = this.geneBackgroundExpressionSamplers[geneIdx];

        // gene's new background expression
//        curBackgroundExp = sampler.randomSample(backgroundExpressions[time-1]);
        prevBackgroundExp = samplingRecords.getBackgroundExpressionValue(time-1);
        curBackgroundExp = prevBackgroundExp;
        // use the new background expression value to calculate IP and INPUT reads count expectations
        // shape individualNumber × geneNumber
        double[] curMethLevel = new double[] {samplingRecords.getMethylationLevelValue(time)};
        double[] curNonSpecificEnrich = new double[] {samplingRecords.getNonSpecificEnrichmentRatio(time)};
        double[][] newIPReadsExpectation, newINPUTReadsExpectation;
        ReadsExpectation re = new ReadsExpectation(sampleIPReads, sampleINPUTReads, nonPeakIPReads, nonPeakINPUTReads, curMethLevel, curNonSpecificEnrich);
        newIPReadsExpectation = re.getIPReadsExepectation();
        newINPUTReadsExpectation = re.getINPUTReadsExpectation();
        re = null;

        // background expression sampler for corresponding gene
//        sampler = this.geneBackgroundExpressionSamplers[geneIdx];
        BackgroundExpressionSampler sampler = null;
//        prevBackgroundExp = backgroundExpressions[time-1];

        // get IP and INPUT gene reads count expectation of each individual, shape 1 × individualNumber
        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, sampleIPExpectation, sampleINPUTExpectation,
                                             nonPeakIPReads, nonPeakINPUTReads, nonPeakIPExpectation, nonPeakINPUTExpectation);
        prevIPExpectation = readsData[2];
        prevINPUTExpectation = readsData[3];

        readsData = this.getReadsDataForGene(sampleIPReads, sampleINPUTReads, newIPReadsExpectation, newINPUTReadsExpectation,
                                             nonPeakIPReads, nonPeakINPUTReads, nonPeakIPExpectation, nonPeakINPUTExpectation);
        ipCount = readsData[0];
        inputCount = readsData[1];
        ipExpectation = readsData[2];
        inputExpectation = readsData[3];
        ipNonPeakCount = readsData[4];
        inputNonPeakCount = readsData[5];
        ipNonPeakExpCount = readsData[6];
        inputNonPeakExpCount = readsData[7];
        ipOverdispersion = samplingRecords.getIpOverdispersionValue(time);
        inputOverdispersion = samplingRecords.getInputOverdispersionValue(time);

        double meth = samplingRecords.getMethylationLevelValue(time);
        double nonSpecificEnrich = samplingRecords.getNonSpecificEnrichmentRatio(time);
        prevPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, prevIPExpectation, prevINPUTExpectation,
                                                          ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount,
                                                          ipOverdispersion, inputOverdispersion, meth, nonSpecificEnrich, sampler, prevBackgroundExp);
        curPosteriorProba = this.logPosteriorProbability(ipCount, inputCount, ipExpectation, inputExpectation,
                                                         ipNonPeakCount, inputNonPeakCount, ipNonPeakExpCount, inputNonPeakExpCount,
                                                         ipOverdispersion, inputOverdispersion, meth, nonSpecificEnrich, sampler, curBackgroundExp);

//        samplingRes = sampler.getSamplingRes(curPosteriorProba, prevPosteriorProba, true);
//        // if samplingRes is true, means the new sampling value was accepted, also need to renew reads count expectation
//        if (samplingRes) {
//            backgroundExpressions[time] = curBackgroundExp;
//            this.renewReadsExpectation(sampleIPExpectation, ipExpectation);
//        } else { // otherwise, the methylation level of new iteration is same as the previous iteration
//            backgroundExpressions[time] = prevBackgroundExp;
//        }

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
                                           double[][] ipCountExpectation, double[][] inputCountExpectation,
                                           int[][] ipNonPeakReads, int[][] inputNonPeakReads,
                                           double[][] ipNonPeakExpectation, double[][] inputNonPeakExpectation) {
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

        double[] ipNonPeak = new double[this.individualNumber];
        double[] inputNonPeak = new double[this.individualNumber];
        double[] nonPeakIPExpectation = new double[this.individualNumber];
        double[] nonPeakINPUTExpectation = new double[this.individualNumber];
        for (int j=0; j<this.individualNumber; j++) {
            ipNonPeak[j] = ipNonPeakReads[j][0];
            inputNonPeak[j] = inputNonPeakReads[j][0];
            nonPeakIPExpectation[j] = ipNonPeakExpectation[j][0];
            nonPeakINPUTExpectation[j] = inputNonPeakExpectation[j][0];
        }

        return new double[][] {ipCount, inputCount, ipExpectation, inputExpectation,
                               ipNonPeak, inputNonPeak, nonPeakIPExpectation, nonPeakINPUTExpectation};
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
                                           double[] ipNonPeakCount, double[] inputNonPeakCount,
                                           double[] ipNonPeakExpCount, double[] inputNonPeakExpCount,
                                           double ipOverdispersion, double inputOverdispersion,
                                           double methLevel, double nonSpecifcEnrichment,
                                           BackgroundExpressionSampler sampler, double backgroundExp) {
        double nonSpecificEnrichProba = this.nonSpecificEnrichmentSampler.getLogDensity(nonSpecifcEnrichment);
        double methLevelProba = this.methLevelSampler.getLogDensity(methLevel);
        double ipOverdispersionProba = this.ipOverdispersionSampler.getLogDensity(ipOverdispersion);
        double inputOverdispersionProba = this.inputOverdispersionSampler.getLogDensity(inputOverdispersion);
//        double backgroundExpressionProba = sampler.getLogDensity(backgroundExp);

        double ipIndividualProba = ProbabilityCalculator.logNegativeProbability(ipReadsCount, ipExpectation, ipOverdispersion);
        double inputIndividualProba = ProbabilityCalculator.logNegativeProbability(inputReadsCount, inputExpectation, inputOverdispersion);
        double ipNonPeakProba = ProbabilityCalculator.logNegativeProbability(ipNonPeakCount, ipNonPeakExpCount, ipOverdispersion);
        double inputNonPeakProba = ProbabilityCalculator.logNegativeProbability(inputNonPeakCount, inputNonPeakExpCount, inputOverdispersion);
        return nonSpecificEnrichProba + methLevelProba + ipOverdispersionProba + inputOverdispersionProba + ipIndividualProba + inputIndividualProba + ipNonPeakProba + inputNonPeakProba;    // + backgroundExpressionProba;
    }

    /**
     * use the median of methylation level sampling result as methylation quantification result, shape 1 × geneNumber
     */
    private double[] quantify(SamplingRecords samplingRecords) {
        int keepNum = this.samplingTime - this.burnIn;
        double[] nonSpecificEnrichRemainValues = new double[keepNum];
        double[] methLevelRemainValues = new double[keepNum];
        double[] ipOverdispersionRemainValues = new double[keepNum];
        double[] inputOverdispersionRemainValues = new double[keepNum];

        double[] nonSpecificEnrichment = samplingRecords.getNonSpecificEnrichment();
        double[] methylationLevels = samplingRecords.getMethylationLevels();
        double[] ipOverdispersions = samplingRecords.getIpOverdispersions();
        double[] inputOverdispersions = samplingRecords.getInputOverdispersions();
        double[] backgroundExpressions = samplingRecords.getBackgroundExpressions();
        System.arraycopy(nonSpecificEnrichment, this.burnIn-1, nonSpecificEnrichRemainValues, 0, keepNum);
        System.arraycopy(methylationLevels, this.burnIn-1, methLevelRemainValues, 0, keepNum);
        System.arraycopy(ipOverdispersions, this.burnIn-1, ipOverdispersionRemainValues, 0, keepNum);
        System.arraycopy(inputOverdispersions, this.burnIn-1, inputOverdispersionRemainValues, 0, keepNum);

//        double[] geneBackgroundExpressionRemainValue = null;
//
//        geneBackgroundExpressionRemainValue = new double[keepNum];
//        System.arraycopy(backgroundExpressions, this.burnIn-1, geneBackgroundExpressionRemainValue, 0, keepNum);

        double nonSpecificEnrich, methylationValue, ipOverdispersionValue, inputOverdispersionValue, backgroundExpressionValue;
        nonSpecificEnrich = this.distributionMedian(nonSpecificEnrichRemainValues);
        methylationValue = this.distributionMedian(methLevelRemainValues);
        ipOverdispersionValue = this.distributionMedian(ipOverdispersionRemainValues);
        inputOverdispersionValue = this.distributionMedian(inputOverdispersionRemainValues);
//        backgroundExpressionValue = this.distributionMedian(geneBackgroundExpressionRemainValue);
        backgroundExpressionValue = samplingRecords.getGeneBackgroundExp();

        return new double[] {nonSpecificEnrich, methylationValue, ipOverdispersionValue, inputOverdispersionValue, backgroundExpressionValue};
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

    /**
     * output result
     */
    private void output() {
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(this.outputFile))));
            bfw.write("#index\tnon-specificEnrichment\tmethylationLevel\tipOverdispersion\tinputOverdispersion\tbackgroundExpression\n");
            ArrayList<Integer> sortedIdx = new ArrayList<>(this.testRecord.keySet().stream().sorted((o1, o2) -> o1 - o2).collect(Collectors.toList()));
            for (Integer idx: sortedIdx) {
                bfw.write(idx + "\t" + this.testRecord.get(idx));
                bfw.newLine();
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
